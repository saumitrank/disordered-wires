import kwant
import numpy as np
import math

class rgf():
    """Recursive Green's function object that is based on a KWANT lead
    """




    def __init__(self, lead:kwant.system.InfiniteSystem):
        """Initialize the rgf object

        Args:
            lead (kwant.system.InfiniteSystem): Lead with which to caclulate properties
        """

        #Create inter-cell hopping with full unit cell
        self.lead = lead
        self.ham = lead.cell_hamiltonian()
        self.hop = np.zeros(self.ham.shape, dtype=np.complex128)
        self.hop[:, :lead.inter_cell_hopping().shape[1]] = lead.inter_cell_hopping()
        
        #Position of hoppings in ham and hop (so that we can add disorder here later)
        non_zero_pos = np.transpose(np.nonzero(self.ham)) 
        self.pairs_ham = np.array([(i[0], i[1]) for i in non_zero_pos if i[1]>i[0]]) #Only include upper right triangle
        
        non_zero_pos = np.transpose(np.nonzero(self.hop))
        self.pairs_hop = np.array([(i[0], i[1]) for i in non_zero_pos]) #Not hermitian

        


    def _add_disorder(self, params:dict):
        """Add disorder to a slice and return the inter and intra unit cell Hamiltonians

        Args:
            params (dict): Disorder parameters consisting of
            Nimp (float): Probability of bond/site having disorder (from 0 to 1)
            V (float): Strength of disorder chosen from uniform distribution in [-V, V]
            dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
            leads (bool): Whether or not there should be leads

        Returns:
            Inter and intra unit-cell Hamiltonians with disorder
        """
        #Default params
        p = {'Nimp':1.0, 'dis':0, 'V':0.0}

        #Update
        p.update(params)
        
        ham = np.copy(self.ham)
        hop = np.copy(self.hop)

        #On-site disorder
        if p['dis'] == 0:
            
            chances = np.array(np.random.random(ham.shape[0]) <= p['Nimp'])
            values = 2*p['V']*np.random.random(ham.shape[0]) - p['V']
            values = values*chances
            np.fill_diagonal(ham, values)

        #Hopping disorder
        if p['dis'] == 1 :
            
            #Add on the bonds within the unit cell
            chances = np.array(np.random.random(self.pairs_ham.shape[0]) <= p['Nimp'])
            values = 2*p['V']*np.random.random(self.pairs_ham.shape[0]) - p['V']
            values = values*chances
            
            ham[tuple(self.pairs_ham.T)] = ham[tuple(self.pairs_ham.T)] + values
            ham[tuple(self.pairs_ham[:, ::-1].T)] = ham[tuple(self.pairs_ham[:, ::-1].T)] + values.conjugate()
            
            #Add on the bonds connecting unit cells
            chances = np.array(np.random.random(self.pairs_hop.shape[0]) <= p['Nimp'])
            values = 2*p['V']*np.random.random(self.pairs_hop.shape[0]) - p['V']
            values = values*chances
            
            hop[tuple(self.pairs_hop.T)] = hop[tuple(self.pairs_hop.T)] + values
            
        return ham, hop
    
    
    

    def bands(self, momenta:np.array=np.array([0]))->np.array:
        """Calculate band energies at momenta

        Args:
            momenta (np.array, optional): Momenta at which to calculate energies. Defaults to np.array([0]).

        Returns:
            np.array: Array of energies for each momenta
        """
        momenta = np.array(momenta)
        
        #Define Hamiltonian
        H = self.ham[np.newaxis, :] \
            + self.hop[np.newaxis, :]*np.exp(-1j*momenta)[:, np.newaxis, np.newaxis] \
            + self.hop.T.conjugate()[np.newaxis, :]*np.exp(1j*momenta)[:, np.newaxis, np.newaxis]
    
        return np.linalg.eigvalsh(H)
    
    
    

    def density(self, params:dict)->np.array:
        """Calculate the DOS using RGF

        Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. 
        This contains:
        lengths (list): Lengths at which to calculate
        energies (list): Energies at which to calculate
        Nimp (float): Probability of bond/site having disorder (from 0 to 1)
        V (float): Strength of disorder chosen from uniform distribution in [-V, V]
        dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
        gamma_ratio (float): Ratio of imaginary part of energy to real part
        leads (bool): Whether or not the lead self-energy is included

        Returns:
            np.array: Array of densities of shape (len(lengths), len(energies))
        """
        #Default params
        p = {'lengths':[100], 'energies':[0.0],  'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0, 
             'gamma_ratio' : 1, 'leads':True}
        
        #Update
        p.update(params)

        energies = np.sort(p['energies'])
        lengths = np.sort(p['lengths']).astype(int)
        
        densities = []
        H, V = np.copy(self.ham), np.copy(self.hop)
        I = np.identity(H.shape[0], dtype=np.complex128)
        
        #Calculate self energy of leads from KWANT. Assuming they are symmetric
        sigma = np.zeros((len(energies), H.shape[0], H.shape[1]), dtype=np.complex128)
        for i, energy in enumerate(energies):
            S = self.lead.selfenergy(energy=energy)
            sigma[i, :S.shape[0], :S.shape[1]] = S
            
        #Start of Recursive procedure
        Z = (energies + 1j*np.abs(energies/p['gamma_ratio']))[:, np.newaxis, np.newaxis]*np.tile(I, (len(energies), 1, 1))
            
        #Initialize using lead self energy. F can be zero according to MacKinnon
        if p['leads']:
            R = np.linalg.inv(Z - H - sigma)
        else:
            R = np.linalg.inv(Z - H) #No Self energy
            
        F = np.zeros((len(energies), H.shape[0], H.shape[1]), dtype=np.complex128) 
        s = np.zeros(len(energies), dtype=np.complex128)
        H, V = self._add_disorder(p)
        
        for i in range(lengths[-1]):
            
            #Recursion relations
            R = np.linalg.inv(Z - H - V.T.conj() @ R @ V)
            s = s + np.trace(R @ (F + I), axis1=1, axis2=2)
            H, V = self._add_disorder(p)
            F = V.T.conj() @ R @ (F + I) @ R @ V
            
            #Add right side and compute density if we reach a given length
            if i+1 in lengths:
                
                if p['leads']:
                    R_final = np.linalg.inv(Z - H - sigma - V.T.conj() @ R @ V)
                    s_final = s + np.trace(R_final @ (F + I), axis1=1, axis2=2)
                    
                else:
                    s_final = np.copy(s) #No self energy
                
                densities.append(abs(-s_final.imag/(math.pi*H.shape[0]*(i+1))))
                
        return np.array(densities)
    



    def localization_length(self, params:dict)->tuple:
        """Calculate the localization length using RGF

        Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. 
        This contains:
        lengths (list): Lengths at which to calculate
        energies (list): Energies at which to calculate
        Nimp (float): Probability of bond/site having disorder (from 0 to 1)
        V (float): Strength of disorder chosen from uniform distribution in [-V, V]
        dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
        gamma_ratio (float): Ratio of imaginary part of energy to real part
        leads (bool): Whether or not the lead self-energy is included

        Returns:
            Two arrays for localization length and the error in it
        """
        #Default values
        p = {'lengths':[100], 'energies':[0.0],  'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0, 
             'gamma_ratio' : 1}
               
        #Update
        p.update(params)
        
        energies = np.sort(p['energies'])
        lengths = np.sort(p['lengths']).astype(int)
        
        xi = []
        xi_error = []
        H, V = np.copy(self.ham), np.copy(self.hop)
        I = np.identity(H.shape[0], dtype=np.complex128)
        
        #Calculate self energy of leads from KWANT. Assuming they are symmetric
        sigma = np.zeros((len(energies), H.shape[0], H.shape[1]), dtype=np.complex128)
        for i, energy in enumerate(energies):
            S = self.lead.selfenergy(energy=energy)
            sigma[i, :S.shape[0], :S.shape[1]] = S
            
        #Start of Recursive procedure
        Z = (energies + 1j*np.abs(energies/p['gamma_ratio']))[:, np.newaxis, np.newaxis]*np.tile(I, (len(energies), 1, 1))
            
        #Initialize using lead self energy
        R = np.linalg.inv(Z - H - sigma) # G_{1, 1}
        B = np.tile(I, (len(energies), 1, 1))
        c = np.zeros(len(energies), dtype=np.complex128)
        d = np.zeros(len(energies), dtype=np.complex128)
        H, V = self._add_disorder(p)
        
        for i in range(lengths[-1]):
            
            #Recursion relations
            R = np.linalg.inv(Z - H - V.T.conj() @ R @ V)
            b = np.linalg.norm(B @ R, ord='fro', axis=(1, 2))
            c = c + np.log(b)
            d = d + np.log(b)**2
            H, V = self._add_disorder(p)
            B = B @ R @ V / b[:, np.newaxis, np.newaxis]
            
            #Add right side and compute xi if we reach a given length
            if i+1 in lengths:
                R_final = np.linalg.inv(Z - H - sigma - V.T.conj() @ R @ V)       
                b_final = np.linalg.norm(B @ R_final, ord='fro', axis=(1, 2))
                c_final = c + np.log(b_final)
                d_final = d + np.log(b_final)**2
                
                xi.append(abs((i+1)/c_final))
                xi_error.append(np.sqrt(abs( (d_final/(i+1)) - (c_final/(i+1))**2 )))
                
        return np.array(xi), np.array(xi_error)
    
    


    def conductivity(self, params:dict)->np.array:
        """Calculate the conductivity using RGF (Not as reliable as S-matrix methods)

        Args:
        params (dict): A dictionary of parameters for the zigzag ribbon. 
        This contains:
        lengths (list): Lengths at which to calculate
        energies (list): Energies at which to calculate
        Nimp (float): Probability of bond/site having disorder (from 0 to 1)
        V (float): Strength of disorder chosen from uniform distribution in [-V, V]
        dis (0, 1): Disorder type. 0 for on-site, 1 for hopping
        gamma_ratio (float): Ratio of real part of energy to imaginary part
        leads (bool): Whether or not the lead self-energy is included

        Returns:
            np.array: numpy array of conductivities
        """
        
        #Default values
        p = {'lengths':[100], 'energies':[0.0],  'Nimp' : 1.0, 'V' : 0.0, 'dis' : 0, 
             'gamma_ratio' : 1}
               
        #Update
        p.update(params)
        
        energies = np.sort(p['energies'])
        lengths = np.sort(p['lengths']).astype(int)
        
        conductivities = []
        H, V = np.copy(self.ham), np.copy(self.hop)
        I = np.identity(H.shape[0], dtype=np.complex128)
        
        #Calculate self energy of leads from KWANT. Assuming they are symmetric
        sigma = np.zeros((len(energies), H.shape[0], H.shape[1]), dtype=np.complex128)
        for i, energy in enumerate(energies):
            S = self.lead.selfenergy(energy=energy)
            sigma[i, :S.shape[0], :S.shape[1]] = S
            
        #Start of Recursive procedure
        Z = (energies + 1j*np.abs(energies/p['gamma_ratio']))[:, np.newaxis, np.newaxis]*np.tile(I, (len(energies), 1, 1))
            
        #Initialize using lead self energy the rest are zero
        R = np.linalg.inv(Z - H - sigma)
        B = np.zeros((len(energies), H.shape[0], H.shape[1]), dtype=np.complex128)
        Cp = np.zeros((len(energies), H.shape[0], H.shape[1]), dtype=np.complex128)
        Cm = np.zeros((len(energies), H.shape[0], H.shape[1]), dtype=np.complex128)
        sxx = np.zeros(len(energies), dtype=np.complex128)
        H, V = self._add_disorder(p)
        
        for i in range(lengths[-1]):
            
            #B, C tranformations (to avoid using positions)
            B = B + 1j*(Cp + Cm) + (V.T.conj() @ (R - np.transpose(R, (0, 2, 1)).conj()) @ V)/2
            Cp = Cp - 1j*(V.T.conj() @ (R - np.transpose(R, (0, 2, 1)).conj()) @ V)/2
            Cm = Cm - 1j*(V.T.conj() @ (R - np.transpose(R, (0, 2, 1)).conj()) @ V)/2
                        
            #Recursion relations
            R = np.linalg.inv(Z - H - V.T.conj() @ R @ V)
            sxx = sxx + np.trace( (B @ R).real + Cp @ np.transpose(R, (0, 2, 1)).conj() @ Cm @ R , axis1=1, axis2=2)
            H, V = self._add_disorder(p)
            B = V.T.conj() @ R @ (B + 2*Cp @ np.transpose(R, (0, 2, 1)).conj() @ Cm) @ R @ V
            Cp = V.T.conj() @ R @ Cp @ np.transpose(R, (0, 2, 1)).conj() @ V
            Cm = V.T.conj() @ np.transpose(R, (0, 2, 1)).conj() @ Cm @ R @ V
            
            #Add right side and compute cond if we reach a given length
            if i+1 in lengths:
                B_final = B + 1j*(Cp + Cm) + (V.T.conj() @ (R - np.transpose(R, (0, 2, 1)).conj()) @ V)/2
                Cp_final = Cp - 1j*(V.T.conj() @ (R - np.transpose(R, (0, 2, 1)).conj()) @ V)/2
                Cm_final = Cm - 1j*(V.T.conj() @ (R - np.transpose(R, (0, 2, 1)).conj()) @ V)/2
                R_final = np.linalg.inv(Z - H - V.T.conj() @ R @ V)
                sxx_final = sxx + np.trace( (B_final @ R_final).real +  \
                                           Cp_final @ np.transpose(R_final, (0, 2, 1)).conj() @ Cm_final @ R_final , axis1=1, axis2=2)
                conductivities.append(4*sxx_final/((i+1)*H.shape[0]))

        return abs(np.array(conductivities).real)