import numpy as np
from itertools import permutations
from matscipy.neighbours import neighbour_list

supported_projection = [v+w for v,w in permutations('abc', 2)]

def projection_magnitude(v, dir ='c'):
    """
    Parameters
    ----------
    v: Numpy Array 
        Cartesian vector
    dir: String 
        Direction to check if v is on, either 'a', 'b', or 'c'
    thres: float (optional)
        Threshold for vector alignment along direction specified. 
        This is implemented as the dot product of v and unit vector 
        corresponding to dir.

    Given cartesian vector v and direction either a, b, or c, dir, calculate projected 
    vector magnitude on a plane.
    """
    dir = dir.lower()
    dir_idx = 'abc'.find(dir)
    ref_vec = [0,0,0]
    ref_vec[dir_idx] = 1
    proj_magnitude = v@ref_vec/np.linalg.norm(v)
    return proj_magnitude

def _projected_vector_angle(v, proj ='ab'):
    """
    Parameters
    ----------
    v: Numpy Array 
        cartesian vector
    proj: String
    Projection onto which v is projected, Must be a 2-product of 'a', 'b', or 'c'.

    Given cartesian vector v and projection, evaluate the angle made by opp and adj.
    """
    proj = proj.lower()
    ordered_idx = list(proj) +[axis for axis in 'abc' if axis not in proj]
    opp_idx, adj_idx, flatted_idx = ['abc'.find(i) for i in ordered_idx]
    v[flatted_idx] = 0
    #cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arctan(v[opp_idx]/v[adj_idx])
    return np.degrees(angle)

def vec_along_dir(v, dir ='c', thres = 0.9):
    """
    Parameters
    ----------
    v: Numpy Array 
        Cartesian vector
    dir: String 
        Direction to check if v is on, either 'a', 'b', or 'c'
    thres: float (optional)
        Threshold for vector alignment along direction specified. 
        This is implemented as the dot product of v and unit vector 
        corresponding to dir.

    Given cartesian vector v and direction either a, b, or c, dir, check if 
    v is aligned with dir through specified tolerance.
    """
    proj_magnitude = projection_magnitude(v, dir)
    return np.abs(proj_magnitude) > thres #, proj_magnitude >0

def vec_align_dir(v, dir ='c'):
    """
    Parameters
    ----------
    v: Numpy Array 
        Cartesian vector
    dir: String 
        Direction to check if v is on, either 'a', 'b', or 'c'
    thres: float (optional)
        Threshold for vector alignment along direction specified. 
        This is implemented as the dot product of v and unit vector 
        corresponding to dir.

    Given cartesian vector v and direction either a, b, or c, dir, check if 
    v is align or anti-aligned to direction.
    """
    proj_magnitude = projection_magnitude(v, dir)
    return proj_magnitude >0

def tilt_angle(atoms, proj, thres=0.9, absolute=True):
    '''
    Parameters
    ----------
    atoms: ASE Atoms 
        monochalcogenide atoms
    proj: String 
        Direction to check if v is on, either 'a', 'b', or 'c'
    thres: float (optional)
        Threshold for vector alignment along direction specified. 
        This is implemented as the dot product of v and unit vector 
        corresponding to dir.
    absolute: bool (optional)
        Whether to return positive-definite angles or not. 

    proj is the plane to project v. The proj string 
    has the following order: opposite, adjacent. If not absolute, 
    angle is positive when Se->Ge vector is aligned with the positive
    direction of the adjacent direction and negative when it's anti-aligned. 
    '''
    if proj not in supported_projection:
        raise ValueError(f"Projection {proj} not supported.")
    src, dst, dis_vec = neighbour_list(
        'ijD',
        atoms=atoms,
        cutoff={('Ge', 'Se'): 3.0}
    )
    Se_idx = [atom.index for atom in atoms if atom.symbol=='Se']
    tilt_angle = {}
    for idx in Se_idx:
        for D in dis_vec[src==idx]:
            along_dir = vec_along_dir(D, dir =proj[1], thres = thres)
            if along_dir:
                align_dir = vec_align_dir(D, dir =proj[1]) #align or anti align?
                proj_angle = _projected_vector_angle(D, proj)
                proj_angle = -proj_angle if not align_dir else proj_angle
                proj_angle = np.absolute(proj_angle) if absolute else proj_angle
                try: 
                    tilt_angle[idx]
                    print(f'Duplicate found at idx = {idx}.')
                    print(f'dis_vec = {dis_vec[src == idx]}.')
                    tilt_angle[idx] = tilt_angle[idx]+[proj_angle]
                except KeyError:
                    tilt_angle[idx] = [proj_angle]
    return tilt_angle
                
