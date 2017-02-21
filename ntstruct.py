import numpy as np
import os
import supertb as tb

def create_nt(n,m):
    """
    Creates a STRUCTURE object for an (n,m) Carbon Nanotube

    INPUT:
        n : integer
        m : integer
    OUTPUT:
        structure object
    """
    from pymatgen.io.cif import CifParser
    #home=os.getenv('HOME')
    os.system('gfortran CNT.F -o ./cnt')
    cntdir='./cnt'
    os.system('echo '+str(n)+'  '+str(m)+' > cntin')
    os.system('echo 1 >> cntin')
    os.system(cntdir+'< cntin > output')
    struct=CifParser('nanotube.cif').as_dict()
    #print struct
    a = np.float(struct['block_1']['_cell_length_a'])
    b = np.float(struct['block_1']['_cell_length_b'])
    c = np.float(struct['block_1']['_cell_length_c'])
    alpha = np.float(struct['block_1']['_cell_angle_alpha'])
    beta = np.float(struct['block_1']['_cell_angle_beta'])
    gamma = np.float(struct['block_1']['_cell_angle_gamma'])
    atoms = struct['block_1']['_atom_site_type_symbol']
    x = np.array(struct['block_1']['_atom_site_fract_x']).astype(np.float)
    y = np.array(struct['block_1']['_atom_site_fract_y']).astype(np.float)
    z = np.array(struct['block_1']['_atom_site_fract_z']).astype(np.float)
    
    fracs = np.vstack((x,y,z)).T
    lattice = tb.Lattice.from_parameters(a,b,c,alpha,beta,gamma)
    cnt = tb.Structure(lattice, ['C']*len(fracs), fracs, coords_are_cartesian=False)
    
    #clean
    os.system("rm nanotube.cif")
    os.system("rm output")
    os.system("rm cntin")
    return cnt
