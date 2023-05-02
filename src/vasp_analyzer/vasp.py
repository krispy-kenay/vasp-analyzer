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
               filepath=None,
               incar=INCAR(),
               eferm=None
               ):
        """
        Usually this object should be set up empty (you are unlikely to want to set it up manually) and fetched from the 'vasprun.xml' file
        Args:
            total_dos: DOSC object of 'total' type
            partial_dos: DOSC object of 'partial' type
            kpoints: KPOINTS object
            incar: INCAR object
            eferm: fermi energy
        """
        if filepath is not None: self.add_filepath(filepath) 
        else: self.files = None

        self.classes = {'incar': incar,
                        'kpoints': KPOINTS(),
                        'dos': DOSC(),
                        'projections': 0}
        
        self.triggers = {k: False for k in self.classes}

        # Properties
        self.efermi = eferm
        self.elements = None
    
    def __getattr__(self, key):
        if key in self.classes:
            return self.classes.get(key)
        else:
            raise KeyError("Nonexistent")

    ########################################
    # Adding Methods
    ########################################
    
    def add_filepath(self, filepath):
        if type(filepath) == str:
            self.files = [filepath]
        elif type(filepath) == list:
            self.files = filepath
        else:
            raise ValueError("Please provide either a string or list of filename strings!")


    ########################################
    # Getter Methods
    ########################################

    # For Kpoints
    def get_band_structure(self):
        self.load(kpoints=True)
        return self.kpoints.get_band_structure()
    
    def get_band_structure_matched(self):
        self.load(kpoints=True, projections=True)
        return self.kpoints.get_sorted_hun_bands()
    
    def get_spin_delta_occ(self, cutoff=1e-3, absolute=False):
        self.load(kpoints=True)
        return self.kpoints.get_spin_difference_hl(mode='occupied', cutoff=cutoff, absolute=absolute)
    
    def get_spin_delta_tot(self, cutoff=1e-3, absolute=False):
        self.load(kp=True)
        return self.kpoints.get_spin_difference_hl(mode='total', cutoff=cutoff, absolute=absolute)

    def get_spin_delta_hob(self):
        self.load(kp=True)
        return self.kpoints.get_spin_difference_hl(mode='highest')
    
    def get_spin_delta_lub(self):
        self.load(kp=True)
        return self.kpoints.get_spin_difference_hl(mode='lowest')
    
    def get_spin_delta_hun(self, absolute=False):
        self.load(kpoints=True, incar=True, projections=True)
        return self.kpoints.get_spin_difference_hl(mode='occupied_hungarian', moment=np.array(self.incar['MAGMOM']), absolute=absolute)

    # For DOS
    def get_dos_total(self):
        self.load(dos=True)
        return (self.dos.get_energy(), self.dos.get_tdos())
    
    def get_dos_element(self):
        self.load(dos=True)
        return (self.dos.get_energy(), self.dos.calc_dos('element'))
    
    def get_dos_orbital(self):
        self.load(dos=True)
        return (self.dos.get_energy(), self.dos.calc_dos('spd'))
    
    def get_dos_orbital_element(self):
        self.load(dos=True)
        return (self.dos.get_energy(), self.dos.calc_dos('element spd'))

    ########################################
    # Loading Methods
    ########################################

    def load_all(self):
        self.load(incar=True, dos=True, kpoints=True, projections=True)
    
    def load(self, **kwargs):
        
        keys = []
        for key, value in kwargs.items():
            if key in self.classes and value == True:
                if self.triggers[key] == False:
                    keys.append(key)
                    self.triggers[key] = True

        if keys != []:
            for i, file in enumerate(self.files):
                if 'vasprun.xml' not in file: xml = xml_handler(file + '/vasprun.xml')
                else: self.xml = xml_handler(file)

                if len(self.files) != 1: modi = 'Section ' + str(i) + ' '
                else: modi = ''

                for key in keys:
                    print("Loading section: " + str(key))
                    if key == 'projections':
                        self.classes['kpoints'].load_site_projections(xml, modifier=modi)
                    else:
                        self.classes[key].load(xml, modifier=modi)
