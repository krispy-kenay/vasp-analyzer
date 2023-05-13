import numpy as np

from .data.xml import xml_handler
from .data.incar import INCAR
from .data.kpoints import KPOINTS
from .data.dos import DOSC

from .graph.plotter import vplotter

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

Main VASP handler
- Most other classes are unified within this to give a simplified and unified interface

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class VASP:
    def __init__(self, filepath=None):
        """
        Main class which contains other classes and provides a simplified interface to access functions
        Args:
            filepath: can add filepath upon object creation (also possible later on)
        """
        if filepath is not None: self.add_filepath(filepath) 
        else: self.files = None

        self.reset()
        self.plotter = None

        # Properties
        self.efermi = None
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
    # Setter Methods
    ########################################

    def set_fermi_energy(self, energy:float):
        self.dos.efermi = energy
        self.kpoints.efermi = energy
        self.efermi = energy
    
    def set_elements(self, elements:list):
        self.dos.elements = elements
    
    # !CAREFUL! Resets all subclasses, which empties all previously stored content
    def reset(self):
        kp = KPOINTS()
        self.classes = {'incar': INCAR(),
                        'kpoints': kp,
                        'dos': DOSC(),
                        'projections': kp}
        
        self.triggers = {k: False for k in self.classes}

    ########################################
    # Getter Methods
    ########################################

    # Get properties
    def get_fermi_energy(self):
        self.load(dos=True)
        return self.dos.efermi
    
    # Density of states
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

    # Band structure
    def get_band_structure(self):
        self.load(kpoints=True)
        return self.kpoints.get_band_structure()
    
    def get_band_structure_matched(self):
        self.load(kpoints=True, projections=True)
        return self.kpoints.get_sorted_hun_bands()
    
    # Splitting
    def get_splitting_normal_occupied(self, cutoff=1e-3, absolute=False):
        self.load(kpoints=True)
        return self.kpoints.get_spin_difference_hl(mode='occupied', cutoff=cutoff, absolute=absolute)
    
    def get_splitting_normal_total(self, cutoff=1e-3, absolute=False):
        self.load(kpoints=True)
        return self.kpoints.get_spin_difference_hl(mode='total', cutoff=cutoff, absolute=absolute)

    def get_splitting_normal_hob(self):
        self.load(kpoints=True)
        return self.kpoints.get_spin_difference_hl(mode='highest')
    
    def get_splitting_normal_lub(self):
        self.load(kpoints=True)
        return self.kpoints.get_spin_difference_hl(mode='lowest')
    
    def get_splitting_hungarian_occupied(self, absolute=False):
        self.load(kpoints=True, incar=True, projections=True)
        return self.kpoints.get_spin_difference_hl(mode='occupied_hungarian', moment=np.array(self.incar['MAGMOM']), absolute=absolute)

    def get_splitting_hungarian_total(self, absolute=False):
        pass


    
    def get_plotter(self):
        pass

    ########################################
    # Plot Methods
    ########################################

    def get_plotter(self, size=None, width=None, height=None):
        return vplotter(size=size, width=width, height=height)
    
    def set_plotter(self, size=None, width=None, height=None):
        self.plotter = vplotter(size=size, width=width, height=height)

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
