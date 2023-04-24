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

        self.tdos = total_dos
        self.pdos = partial_dos
        self.kpoints = k_points
        self.incar = in_car
        self.efermi = eferm
        self.elements = None

        # Set up triggers
        self.trigger_incar = False
        self.trigger_doscar = False
        self.trigger_kpoints = False
    
    def load_all(self):
        self.load(inc=True, dosc=True, kp=True, pro=True)
    
    def load(self,
             inc:bool=False,
             dosc:bool=False,
             kp:bool=False,
             pro:bool=False
             ):
        if not self.incar: self.incar = INCAR()
        if not self.tdos: self.tdos = DOSC(dostype='total')
        if not self.pdos: self.pdos = DOSC(dostype='partial')
        if not self.kpoints: self.kpoints = KPOINTS()

        for file in self.files:
            xml = xml_handler(file + '/vasprun.xml')
            if len(self.files) != 1:
                modi = file + ' '
            else:
                modi = ''
            
            if inc == True:
                self.trigger_incar = True
                self.incar.load(xml)
            if dosc == True:
                self.trigger_doscar = True
                self.tdos.load(xml)
                self.pdos.load(xml)
            if kp == True:
                self.trigger_kpoints = True
                self.kpoints.load(xml, modifier=modi)
            if pro == True:
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

    def from_vasprun(self, folder, file):
        if type(file) == str:
            xml = xml_handler(folder + file + '/vasprun.xml')
            self.load_INCAR(xml)
            self.load_DOSCAR(xml)
            self.load_KPOINTS(xml)
        elif type(file) == list:
            xml = xml_handler(folder + file[0] + '/vasprun.xml')
            self.load_INCAR(xml)
            self.load_DOSCAR(xml)
            kpointnum = xml.find(".//kpoints/generation/*[@name='divisions']")
            if not self.kpoints:
                try:
                    self.kpoints = KPOINTS(divisions=np.prod([int(i) for i in kpointnum.text.split()]))
                except:
                    self.kpoints = KPOINTS(divisions=1)
            for i in range(len(file)):
                xml = xml_handler(folder + file[i] + '/vasprun.xml')
                self.load_KPOINTS(xml, modifier=str(file[i])+' ')

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