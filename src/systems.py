import kwant
import numpy as np

def add_disorder(sys:kwant.builder.Builder, lat, Nimp:float=1.0, V:float=0.0, dis:int=0)->kwant.system.FiniteSystem:
    """Add on-site or NN hopping disorder to a system

    Args:
        sys (kwant.builder.Builder): _description_
        lat (_type_): _description_
        Nimp (float, optional): _description_. Defaults to 1.0.
        V (float, optional): _description_. Defaults to 0.0.
        dis (int, optional): _description_. Defaults to 0.

    Returns:
        kwant.system.FiniteSystem: _description_
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
        V (float): Strength of disorder chosen from uniform dsbn in [-V, V]
        dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
        leads (bool): Whether or not there should be leads

    Returns:
        kwant.system.FiniteSystem: Zigzag graphene ribbon with/without leads attached
    """
    
    #Default values of parameters
    p = {'a' : 1.0, 't' : 1.0, 'mu' : 0.0, 'mul' : 0.0,
         'L' : 10, 'W' : 10, 'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0, 'leads':True}
    
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