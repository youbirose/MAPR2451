import numpy as np
import os
import supertb as tb


def read_cif(cifile):
    try:
       f = open(cifile,'r')
    except:
       error_msg = 'File '+cifile+' not found!'
       raise IOError(error_msg)
    lines = f.readlines()
    f.close()

    a = float(lines[1].strip().split()[1])
    b = float(lines[2].strip().split()[1])
    c = float(lines[3].strip().split()[1])
    alpha = float(lines[4].strip().split()[1])
    beta = float(lines[5].strip().split()[1])
    gamma = float(lines[6].strip().split()[1])
    lattice = tb.Lattice.from_parameters(a,b,c,alpha,beta,gamma)
    nat = len(lines)-12
    fracs = np.zeros((nat,3))
    specs = ['C']*nat
    for iat in range(nat):
        fracs[iat,:] = [float(x) for x in lines[12+iat].strip().split()[1:]]

    return tb.Structure(lattice, specs, fracs, coords_are_cartesian=False)

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
    cnt = read_cif('nanotube.cif')

    #struct=CifParser('nanotube.cif').as_dict()
    #a = np.float(struct['block_1']['_cell_length_a'])
    #b = np.float(struct['block_1']['_cell_length_b'])
    #c = np.float(struct['block_1']['_cell_length_c'])
    #alpha = np.float(struct['block_1']['_cell_angle_alpha'])
    #beta = np.float(struct['block_1']['_cell_angle_beta'])
    #gamma = np.float(struct['block_1']['_cell_angle_gamma'])
    #atoms = struct['block_1']['_atom_site_type_symbol']
    #x = np.array(struct['block_1']['_atom_site_fract_x']).astype(np.float)
    #y = np.array(struct['block_1']['_atom_site_fract_y']).astype(np.float)
    #z = np.array(struct['block_1']['_atom_site_fract_z']).astype(np.float)
    
    #fracs = np.vstack((x,y,z)).T
    #lattice = tb.Lattice.from_parameters(a,b,c,alpha,beta,gamma)
    #cnt = tb.Structure(lattice, ['C']*len(fracs), fracs, coords_are_cartesian=False)
    
    #clean
    os.system("rm nanotube.cif")
    os.system("rm output")
    os.system("rm cntin")
    return cnt
