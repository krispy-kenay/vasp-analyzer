import numpy as np

from .data.xml import xml_handler
from .data.incar import INCAR
from .data.kpoints import KPOINTS
from .data.dos import DOSC

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

Main VASP handler
- Most other classes are unified within this to give a simplified and unified interface

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class VASP:
    def __init__(self,
               file,
               total_dos=None,
               partial_dos=None,
               k_points=None,
               in_car=None,
               eferm=None
               ):
        """
        Usually this object should be set up empty (you are unlikely to want to set it up manually) and fetched from the 'vasprun.xml' file
        Args:
            total_dos: DOSC object of 'total' type
            partial_dos: DOSC object of 'partial' type
            k_points: KPOINTS object
            in_car: INCAR object
            eferm: fermi energy
        """
        if type(file) == str:
            self.files = [file]
        elif type(file) == list:
            self.files = file
        else:
            raise ValueError("Please provide either a string or list of files!")

        # Objects
        self.tdos = total_dos
        self.pdos = partial_dos
        self.kpoints = k_points
        self.incar = in_car

        # Properties
        self.efermi = eferm
        self.elements = None

        # Set up trigger
        self.trigger_projections = False
    
    def load_all(self):
        self.load(inc=True, dosc=True, kp=True, pro=True)
    
    def load(self,
             inc:bool=False,
             dosc:bool=False,
             kp:bool=False,
             pro:bool=False
             ):

        for file in self.files:
            xml = xml_handler(file + '/vasprun.xml')
            if len(self.files) != 1:
                modi = file + ' '
            else:
                modi = ''
            
            if inc == True and self.incar is None:
                self.incar = INCAR()
                self.incar.load(xml)
            if dosc == True and (self.tdos is None or self.pdos is None):
                self.tdos = DOSC(dostype='total')
                self.pdos = DOSC(dostype='partial')
                self.tdos.load(xml)
                self.efermi = self.tdos.efermi
                self.pdos.load(xml)
            if kp == True and self.kpoints is None:
                self.kpoints = KPOINTS()
                print('Loading KPOINTS')
                self.kpoints.load(xml, modifier=modi)
            if pro == True and self.trigger_projections == False:
                self.trigger_projections = True
                self.kpoints.load_site_projections(xml, modifier=modi)
    
    def get_dos(self, typ:str, subtyp:str=None):
        if self.tdos is None:
            self.load(dosc=True)
        if typ == 'total':
            return self.tdos.get_tdos()
        if typ == 'partial':
            if not subtyp:
                raise ValueError('Please provide a type for the partial DOS')
            else:
                return self.pdos.calc_dos(subtyp)
    




    # FINAL?
    # FINAL?
    # FINAL?
    def get_bs(self):
        self.load(kp=True)
        return self.kpoints.get_band_structure()
    
    def get_spin_delta_occ(self, cutoff=1e-3):
        self.load(kp=True)
        return self.kpoints.get_spin_difference_hl(mode='occupied', cutoff=cutoff)
    
    def get_spin_delta_tot(self, cutoff=1e-3):
        self.load(kp=True)
        return self.kpoints.get_spin_difference_hl(mode='total', cutoff=cutoff)

    def get_spin_delta_hob(self):
        self.load(kp=True)
        return self.kpoints.get_spin_difference_hl(mode='highest')
    
    def get_spin_delta_lub(self):
        self.load(kp=True)
        return self.kpoints.get_spin_difference_hl(mode='lowest')
    
    def get_spin_delta_hun(self):
        self.load(kp=True, inc=True)
        return self.kpoints.get_spin_difference_hl(mode='occupied_hungarian', moment=np.array(self.incar['MAGMOM']))






    def rotation_matrix(self, angle):
        theta = np.radians(angle)
        c, s = np.cos(theta), np.sin(theta)
        R = np.array(((c, -s, 0), (s, c, 0), (0,0,1)))
        return R
    
    def rotate_coords_weight(self, data):
        coords = data[:,0:3]
        weight = data[:,3]
        weight = weight.reshape(-1,1)
        angles = [60, 120, 180, 240, 300]
        out = data
        for angle in angles:
            rot = self.rotation_matrix(angle)
            new_coord = np.dot(coords, rot)
            new_data = np.concatenate([new_coord, weight], axis=1)
            out = np.vstack((out, new_data))
        return out