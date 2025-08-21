from ase.spacegroup.spacegroup import Spacegroup
import numpy as np
from monochalcogenpy.utils import fill_abc, sort_axis_indices

class Spacegroup_MX(Spacegroup):
    """
    Custom Spacegroup class with rotatable symmetry operation.
    """
    def __init__(self, sg_no=31, orientation = 'ac'):
        """
        Parameters
        ----------
        sg_no: Int
            International Union of Crystallography number
            for desired monochalcogenide. In practice, only
            31 (monolayer) and 62 (bulk) are applicable. 
        orientation: Str
            orientation of the armchair direction followed by 
            the direction normal to the layer. For example 
            orientation ='ac' means armchair is along a and vdw along c. 
            orientation ='bc' means armchair is along b and vdw along c.

        Custom Spacegroup class with rotatable symmetry operation.
        Created because ASE's Spacegroup implementation forces orientation
        to be along 'cb' for spacegroup number=31 and 'ca' for 62. These default 
        orientations are used to rotate symmertry operation matrices given to 
        monochalcogenpy.crystal.crystal to create a crystal. 
        """
        super().__init__(sg_no)
        self.sg_no = sg_no
        if self.sg_no == 31:
            self.ref_orientation = 'cb'
        if self.sg_no == 62:
            self.ref_orientation = 'ca'
        self.orientation = orientation
        self.sort_index = sort_axis_indices(
            fill_abc(self.orientation),
            fill_abc(self.ref_orientation)
        )
    def get_symop(self):
        """
        Return rotated rotation and translation
        matrices
        """
        rots = self._sorted_rotations()
        trans = self._sorted_translations()
        return [(r,t) for r,t in zip(rots, trans)]
    def _sorted_translations(self):
        """
        Rotate default translation matrices to self.orientation. 
        """
        _, def_trans = zip(*super().get_symop())
        def_trans = np.array(def_trans)
        trans = def_trans[:,self.sort_index]
        return tuple(trans)
    def _sorted_rotations(self):
        """
        Rotate default rotation matrices to self.orientation. 
        """
        def_rots, _ = zip(*super().get_symop())
        rots = np.array(def_rots).copy()
        #trans = def_trans[:,self.sort_index]
        for i, arr in enumerate(rots):
            diag_elements = np.diag(arr)
            rots[i] = np.diag(diag_elements[self.sort_index])
        return tuple(rots)