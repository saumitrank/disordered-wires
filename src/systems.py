import kwant
import numpy

def zz(params:dict)->kwant.system.FiniteSystem:
    """Create a zigzag ribbon with no spin

    Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. This contains:
        a, t, mu, mul, Nx, Ny, Nimp, V, dis, leads, pbc

    Returns:
        kwant.system.FiniteSystem: Zigzag graphene ribbon with/without leads attached
    """

    #Default values of parameters
    p = {'a' : 1.0, 't' : 1.0, 'mu' : 0.0, 'mul' : 0.0, 
         'Nx' : 10, 'Ny' : 10, 'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0,
         'leads':True, 'pbc' : False}
    
    for key in params.keys():
        if key in p.keys():
            p[key] = params[key]
            
    lat = kwant.lattice.honeycomb(p['a'], norbs=1)
    B, A = lat.sublattices
    
    sys = kwant.Builder()
    
    #Usual stuff
    if p['Ny']==2:
        for i in range(p['Nx']):
            sys[A(i, 0)] = p['mu']
            sys[B(i, 1)] = p['mu']

        sys[lat.neighbors()] = p['t']
        
    else:
        for i in range(p['Nx']):
            for j in range(p['Ny']):
                sys[A(i-j//2, j)] = p['mu']
                sys[B(i-j//2, j)] = p['mu']
                if j%2 == 0:
                    sys[A(p['Nx']-j//2, j)] = p['mu']
                    sys[B(p['Nx']-j//2, j)] = p['mu']
            

        sys[lat.neighbors()] = p['t']
        sys.eradicate_dangling()
    
    #Add PBC (Along x)         
    if p['pbc']:
        p['leads'] = False
        
        sys[A(0, 0)] = p['mu']
        sys[B(p['Nx']-1, p['Ny'] - 1)] = p['mu']
        sys[A(0, 0), B(0, 1)] = p['t']
        sys[B(p['Nx']-1, p['Ny'] - 1), A(p['Nx'] - 1, p['Ny'] - 2)] = p['t']
        for j in range(p['Ny']-1):
            sys[A(0, j), B(p['Nx']-1, j+1)] = p['t']
        
    #Add disorder
    # sys = add_history(sys, lat, p['Nimp'], p['V'], p['dis'])
    
 
    #Add Leads
    if (p['leads']):
        sym_left = kwant.TranslationalSymmetry(lat.vec((-1, 0)))
        lead_left = kwant.Builder(sym_left)
        
        for j in range(p['Ny']):
                lead_left[A(0, j)] = (p['mu'] + p['mul'])
                lead_left[B(0, j)] = (p['mu'] + p['mul'])
        
        lead_left[lat.neighbors()] = p['t']
        lead_left.eradicate_dangling()
             
        lead_right = lead_left.reversed()
        
        sys.attach_lead(lead_left)
        sys.attach_lead(lead_right)

    
    return sys.finalized()


zz({})