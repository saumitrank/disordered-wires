"""
This module contains functions to create various zigzag and armchair graphene ribbons with disorder
"""


import kwant
import numpy as np





def add_disorder(sys:kwant.builder.Builder, lat, Nimp:float=1.0, V:float=0.0, dis:int=0)->kwant.builder.Builder:
    """Add on-site or NN hopping disorder to a system

    Args:
        sys (kwant.builder.Builder): A KWANT Builder to which we add disorder
        lat (kwant.lattice.Monoatomic or kwant.lattice.Polyatomic): A lattice in KWANT
        Nimp (float, optional): Probability of having disorder on a site/bond. Defaults to 1.0.
        V (float, optional): Disorder strength. Disorder is chosen uniformly from [-1, 1]. Defaults to 0.0.
        dis (int, optional): Disorder type. 0 for on-site and 1 for hopping. Defaults to 0.

    Returns:
        kwant.builder.Builder: Builder with disorder added
    """
    #Get number of orbitals
    norbs = list(sys.sites())[0].family.norbs
    m = np.identity(norbs)
    
    #On-site disorder
    if (dis == 0):
        for site in sys.sites():
            if (np.random.random() <= Nimp):
                Vx = 2*V*np.random.random() - V
                sys[site] = sys[site] + Vx*m

    #NN Hopping disorder
    elif (dis == 1):
        for hopping_kind in lat.neighbors():
            pairs = hopping_kind(sys)
            for hop in pairs:
                if (np.random.random() <= Nimp):
                    Vx = 2*V*np.random.random() - V
                    sys[hop] = sys[hop] + Vx*m #Builder automatically takes care of hermitian conjugate
    
    return sys





def zigzag(params:dict)->kwant.system.FiniteSystem:
    """Create a zigzag ribbon with no spin and only NN hopping

    Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. 
        This contains:
        a (float): Lattice constant 
        t (float): NN hopping
        mu (float): Chemical potential in scattering region
        mul (float): Additional potential in leads
        L (int): Length
        W (int): Number of horizontal chains
        Nimp (float): Probability of bond/site having disorder (from 0 to 1)
        V (float): Strength of disorder chosen from uniform distribution in [-V, V]
        dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
        leads (bool): Whether or not there should be leads

    Returns:
        kwant.system.FiniteSystem: Zigzag graphene ribbon with/without leads attached
    """
    
    #Default values of parameters
    p = {'a' : 1.0, 't' : 1.0, 'mu' : 0.0, 'mul' : 0.0,
         'W' : 10, 'W' : 10, 'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0, 'leads':True}
    
    #Update values
    p.update(params)
    
    #Create lattice and builder
    lat = kwant.lattice.honeycomb(p['a'], norbs=1)
    B, A = lat.sublattices
    sys = kwant.Builder()
    
    #Special case: 1D chain
    if p['W']==1:
        for i in range(p['L']):
            sys[A(i, 0)] = p['mu']
            sys[B(i, 1)] = p['mu']

        sys[lat.neighbors()] = p['t']
    
    #For width larger than 1
    else:
        for i in range(p['L']):
            for j in range(p['W']+1):
                sys[A(i-j//2, j)] = p['mu']
                sys[B(i-j//2, j)] = p['mu']
                if j%2 == 0:
                    sys[A(p['L']-j//2, j)] = p['mu']
                    sys[B(p['L']-j//2, j)] = p['mu']
            

        sys[lat.neighbors()] = p['t']
        sys.eradicate_dangling()
        
    #Add disorder
    sys = add_disorder(sys, lat, p['Nimp'], p['V'], p['dis'])
    
    #Add Leads
    if (p['leads']):
        #Create left lead
        sym_left = kwant.TranslationalSymmetry(lat.vec((-1, 0)))
        lead_left = kwant.Builder(sym_left)
        
        for j in range(p['W']+1):
                lead_left[A(0, j)] = (p['mu'] + p['mul'])
                lead_left[B(0, j)] = (p['mu'] + p['mul'])
        
        lead_left[lat.neighbors()] = p['t']
        lead_left.eradicate_dangling()
        
        #Right lead is reversed version of left
        lead_right = lead_left.reversed()
        
        #Attach leads
        sys.attach_lead(lead_left)
        sys.attach_lead(lead_right)
    
    return sys.finalized()





def armchair(params:dict)->kwant.system.FiniteSystem:
    """Create arm-chair ribbon with no spin an NN hopping

    Args:
        params (dict): A dictionary of parameters for the arm-chair ribbon. 
        This contains:
        a (float): Lattice constant 
        t (float): NN hopping
        mu (float): Chemical potential in scattering region
        mul (float): Additional potential in leads
        L (int): Length
        W (int): Number of chains
        Nimp (float): Probability of bond/site having disorder (from 0 to 1)
        V (float): Strength of disorder chosen from uniform distribution in [-V, V]
        dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
        leads (bool): Whether or not there should be leads

    Returns:
        kwant.system.FiniteSystem: Arm-chair ribbon with / without leads attached
    """
    #Default values of params
    p = {'a' : 1.0, 't' : 1.0, 'mu' : 0.0, 'mul' : 0.0, 
         'L' : 10, 'W' : 10, 'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0, 'leads':True}
  
    #Update params
    p.update(params)

    #Define lattice and system
    lat = kwant.lattice.honeycomb(p['a'], norbs=1)
    B, A = lat.sublattices
    sys = kwant.Builder()
    
    #Create sites with chemical potential and add NN hopping
    for i in range(p['L']):
        for j in range(p['W']//2):
            sys[A(j-i//2, i)] = p['mu']
            sys[B(j-i//2, i)] = p['mu']
        if (p['W']%2 != 0) and (i%2 == 0):
            sys[A(j+1 - i//2, i)] = p['mu']
            sys[B(j+1 - i//2, i)] = p['mu']

    sys[lat.neighbors()] = p['t']
    sys.eradicate_dangling()
    
    #Add Disorder
    sys = add_disorder(sys, lat, p['Nimp'], p['V'], p['dis'])
    
    #Add Leads
    if(p['leads']):
        #Create top lead
        sym_top = kwant.TranslationalSymmetry(lat.vec((-1, 2)))
        lead_top = kwant.Builder(sym_top)
        
        for i in range(p['W']//2):
            lead_top[A(i, 0)] = (p['mu'] + p['mul'])
            lead_top[B(i, 0)] = (p['mu'] + p['mul'])
            lead_top[A(i, 1)] = (p['mu'] + p['mul'])
            lead_top[B(i, 1)] = (p['mu'] + p['mul'])
        if (p['W']%2 != 0):
            lead_top[A(i+1, 0)] = (p['mu'] + p['mul'])
            lead_top[B(i+1, 0)] = (p['mu'] + p['mul'])
                    
            lead_top[lat.neighbors()] = p['t']
            lead_top.eradicate_dangling()
            
            #Bottom lead is reveresed top lead
            lead_bot = lead_top.reversed()
            
            #Attach lead
            sys.attach_lead(lead_top)
            sys.attach_lead(lead_bot)
    
    return sys.finalized()





def kane_mele(params:dict)->kwant.system.FiniteSystem:
    """Create a zigzag Kane-Mele ribbon with NN and SOC

    Args:
        params (dict): A dictionary of parameters for the Kane-Mele ribbon. 
        This contains:
        a (float): Lattice constant 
        t (float): NN hopping
        t2 (float): NNN SOC
        mu (float): Chemical potential in scattering region
        mul (float): Additional potential in leads
        L (int): Length
        W (int): Number of chains
        Nimp (float): Probability of bond/site having disorder (from 0 to 1)
        V (float): Strength of disorder chosen from uniform distribution in [-V, V]
        dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
        leads (bool): Whether or not there should be leads

    Returns:
        kwant.system.FiniteSystem: Arm-chair ribbon with / without leads attached
    """
    
    #Default values
    p = {'a' : 1.0, 't' : 1.0, 't2' : 0.3, 't3' : 0.3, 'mu' : 0.0, 'mul' : 0.0, 
         'L' : 10, 'W' : 10, 'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0,
         'leads':True}
    
    #Update params
    p.update(params)  
    
    #Define Pauli matrices
    sigma = np.array((((1, 0), (0, 1)),
                      ((0, 1), (1, 0)),
                      ((0, -1j), (1j, 0)),
                      ((1, 0), (0, -1))))
    
    #Define lattice and sub-lattices with spin
    lat = kwant.lattice.honeycomb(p['a'], norbs=2)
    B, A = lat.sublattices
    sys = kwant.Builder()
    
    #Create sites and NN hopping
    #Special case: 1D chain
    if p['W']==1:
        for i in range(p['L']):
            sys[A(i, 0)] = p['mu']*sigma[0]
            sys[B(i, 1)] = p['mu']*sigma[0]

        sys[lat.neighbors()] = p['t']*sigma[0]

    #For W>1
    else:
        for i in range(p['L']):
            for j in range(p['W']+1):
                sys[A(i-j//2, j)] = p['mu']*sigma[0]
                sys[B(i-j//2, j)] = p['mu']*sigma[0]
                if j%2 == 0:
                    sys[A(p['L']-j//2, j)] = p['mu']*sigma[0]
                    sys[B(p['L']-j//2, j)] = p['mu']*sigma[0]
        
        sys[lat.neighbors()] = p['t']*sigma[0]
        sys.eradicate_dangling()    
    
    #Add disorder
    sys = add_disorder(sys, lat, p['Nimp'], p['V'], p['dis'])
    
    #NNN spin orbit term
    kinds = [(1, 0), (0, 1), (1, -1)]
    
    nuA = [1, 0, -1]  #Spin orbit sign for A to A
    
    for kind, nu in zip(kinds, nuA):
        sys[kwant.builder.HoppingKind(kind, A)] =  1j*p['t2']*nu*sigma[3]
        sys[kwant.builder.HoppingKind(kind, B)] = -1j*p['t2']*nu*sigma[3]
      
    #Add Leads
    if (p['leads']):
        #Create left leads
        sym_left = kwant.TranslationalSymmetry(lat.vec((-1, 0)))
        lead_left = kwant.Builder(sym_left)
        
        for j in range(p['W']+1):
                lead_left[A(0, j)] = (p['mu'] + p['mul'])*sigma[0]
                lead_left[B(0, j)] = (p['mu'] + p['mul'])*sigma[0]
        
        lead_left[lat.neighbors()] = p['t']*sigma[0]
        lead_left.eradicate_dangling()

        #NNN spin orbit term
        kinds = [(1, 0), (0, 1), (1, -1)]
        nuA = [1, 0, -1]  #Spin orbit sign for A to A
        
        for kind, nu in zip(kinds, nuA):
            lead_left[kwant.builder.HoppingKind(kind, A)] =  1j*p['t2']*nu*sigma[3]
            lead_left[kwant.builder.HoppingKind(kind, B)] = -1j*p['t2']*nu*sigma[3]
        
        lead_right = lead_left.reversed()
        
        #Attach leads
        sys.attach_lead(lead_left)
        sys.attach_lead(lead_right)

    return sys.finalized()




def zigzag_lead(params:dict)->kwant.system.InfiniteSystem:
    """Creates a right-facing zizgag nanoribbon lead

    Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. 
        This contains:
        a (float): Lattice constant 
        t (float): NN hopping
        W (int): Number of horizontal chains

    Returns:
        kwant.system.InfiniteSystem: Finalized zigzag lead
    """

    #Default values
    p = {'a' : 1.0, 't' : 1.0, 'mu' : 0.0, 'W' : 10}
    
    #Update params
    p.update(params)

    #Create H and V for clean system
    lat = kwant.lattice.honeycomb(p['a'], norbs=1)
    B, A = lat.sublattices

    sym = kwant.TranslationalSymmetry(lat.vec((1, 0))) # A right facing lead
    lead = kwant.Builder(sym)
    for j in range(p['W']+1):
        lead[A(0, j)] = (p['mu'])
        lead[B(0, j)] = (p['mu'])

    lead[lat.neighbors()] = p['t']
    lead.eradicate_dangling()
    
    return lead.finalized()




def armchair_lead(params:dict)->kwant.system.InfiniteSystem:
    """Creates a top-facing armchair nanoribbon lead

    Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. 
        This contains:
        a (float): Lattice constant 
        t (float): NN hopping
        W (int): Number of horizontal chains

    Returns:
        kwant.system.InfiniteSystem: Finalized armchair lead
    """

    #Default values
    p = {'a' : 1.0, 't' : 1.0, 'mu' : 0.0, 'W' : 10}
    
    #Update params
    p.update(params)

    #Create H and V for clean system
    lat = kwant.lattice.honeycomb(p['a'], norbs=1)
    B, A = lat.sublattices

    sym = kwant.TranslationalSymmetry(lat.vec((-1, 2))) # A top facing lead
    lead = kwant.Builder(sym)
    for i in range((p['W']+1)//2):
        lead[A(i, 0)] = (p['mu'])
        lead[B(i, 0)] = (p['mu'])
        if i != (p['W']//2):
            lead[A(i, 1)] = (p['mu'])
            lead[B(i, 1)] = (p['mu'])


    lead[lat.neighbors()] = p['t']
    lead.eradicate_dangling()
    
    return lead.finalized()
