import kwant

def zigzag(params:dict)->kwant.system.FiniteSystem:
    """Create a zigzag ribbon with no spin and only NN hopping

    Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. This contains:
        a, t, mu, mul, Nx, Ny, Nimp, V, dis, leads, pbc
        a, t, mu, mul, L, W, Nimp, V, dis, leads

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
    # sys = add_disorder(sys, lat, p['Nimp'], p['V'], p['dis'])
    
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