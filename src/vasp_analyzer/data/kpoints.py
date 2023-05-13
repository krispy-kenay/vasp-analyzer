import numpy as np
from scipy.optimize import linear_sum_assignment

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

KPOINT object

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class KPOINT:
    def __init__(self, 
                 name:str, 
                 coords:list = None, 
                 weight:float = None, 
                 spin_up:np.ndarray = None, 
                 spin_down:np.ndarray = None,
                 spin_neut:np.ndarray = None,
                 spin_up_pro:np.ndarray = None,
                 spin_down_pro:np.ndarray = None,
                 spin_neut_pro:np.ndarray = None):
        """
        Args:
            name: k-point tag
            coords: Coordinates of the k-point
            weight: how heavily the k-point is weighed
            spin_up: energy vs. occupancy for spin up state
            spin_down: energy vs. occupancy for spin down state
        """
        self.name = name

        if coords is not None: self.add_coordinates(coords)
        else: self.coords = None
        if weight is not None: self.add_weight(weight)
        else: self.weight = None
        if spin_up is not None: self.add_spin(1, spin_up)
        else: self.spin_up = None
        if spin_down is not None: self.add_spin(-1, spin_down)
        else: self.spin_down = None
        if spin_neut is not None: self.add_spin(0, spin_neut)
        else: self.spin_neut = None
        
        if spin_up_pro is not None:
            self.add_spin_pro(1, spin_up)
        if spin_down_pro is not None:
            self.add_spin_pro(-1, spin_down)
            self.spin_down_site
        if spin_neut_pro is not None:
            self.add_spin_pro(0, spin_neut)
        
        self.spin_diff_abs = None

        self.spin_up_site = None
        self.spin_down_site = None
        self.spin_neut_site = None
    






    # FINAL?
    # FINAL?
    # FINAL?
    def __str__(self):
        name = "Name:"
        coords = "Coordinates:"
        wgt = "Weight:"
        return  f"{name:<15}{self.name:>40} \n{coords:<15}{str(self.coordinates):>40} \n{wgt:<15}{self.weight:>40}"
    
    ########################################
    # Adding Methods
    ########################################

    def add_coordinates(self, coords):
        self.coordinates = coords
    
    def add_spin(self, val, spin):
        if val == 0:
            self.spin_neut = spin
        if val == 1:
            self.spin_up = spin
        if val == -1 or val == 2:
            self.spin_down = spin
    
    def add_spin_pro(self, val, spin):
        if val == 0:
            self.spin_neut_pro = spin
        if val == 1:
            self.spin_up_pro = spin
        if val == -1 or val == 2:
            self.spin_down_pro = spin
    
    def add_weight(self, wgt):
        self.weight = wgt
    
    def add_site_projections(self, sign, band, data):
        dat = np.reshape(data, (-1))

        datlen = len(dat)
        index = int(band.replace('band ', '')) - 1
        length = self.get_len()

        if sign == 1:
            if self.spin_up_site is None:
                self.spin_up_site = np.zeros((length, datlen))
            self.spin_up_site[index] = dat
        
        if sign == -1:
            if self.spin_down_site is None:
                self.spin_down_site = np.zeros((length, datlen))
            self.spin_down_site[index] = dat
        
        if sign == 0:
            if self.spin_neut_site is None:
                self.spin_neut_site = np.zeros((length, datlen))
            self.spin_neut_site[index] = dat
    
    ########################################
    # Calculation Methods
    ########################################
    
    # Simple calculation of up spin energies - down spin energies (no weighing of orbital character)
    def calc_spin_diff(self, occupied=True, cutoff=1e-3, absolute=False):
        if self.spin_up is None and self.spin_down is None:
            raise ValueError("Either spin up or spin down information is missing! \n Make sure that you have set LORBIT = 2 and LSORBIT = True")
        
        if occupied == True: difference = self.spin_up[self.spin_up[:, 1] == 1, 0] - self.spin_down[self.spin_down[:, 1] == 1, 0]
        else: difference = self.spin_up - self.spin_down
        if absolute == True: diff = np.sum(np.abs(difference))
        else: diff = np.sum(difference)
        if abs(diff) < cutoff: diff = 0
        
        prefac = 1 / len(difference)
        diffg = diff * prefac
        return diffg
    
    # More complex calculation of up spin energies - down spin energies based on matching pairs of orbital characters (using the hungarian algorithm)
    def calc_hun_sort_diff(self, moment, occupied=True, cutoff=2, absolute=False):
        if self.spin_up is None or self.spin_down is None or self.spin_up_site is None or self.spin_down_site is None:
            raise ValueError("Either spin up, spin down or site projection information is missing! \n Make sure that you have set LORBIT = 2 and LSORBIT = True")

        if occupied == True:
            spin_up_filtered = self.spin_up[self.spin_up[:, 1] == 1]
            spin_down_filtered = self.spin_down[self.spin_down[:, 1] == 1]
        else:
            spin_up_filtered = self.spin_up
            spin_down_filtered = self.spin_down
        num_up = len(spin_up_filtered)
        num_down = len(spin_down_filtered)

        site_ups = self._reorder_sum_by_magmom(self.spin_up_site[:num_up], moments=moment, inv=False)
        site_downs= self._reorder_sum_by_magmom(self.spin_down_site[:num_down], moments=moment, inv=True)
        cost_matrix = - np.dot(site_ups, site_downs.T)

        abs_diff = np.abs(spin_up_filtered[:,0].reshape(-1,num_up) - spin_down_filtered[:,0].reshape(-1,num_down).T)

        cost_matrix += abs_diff/cutoff
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        

        difference = spin_up_filtered[row_ind] - spin_down_filtered[col_ind]
        if absolute == True: diff = np.sum(np.abs(difference[:,0]))
        else: diff = np.sum(difference[:,0])

        if abs(diff) < 1e-3: diff = 0

        prefac = 1 / num_up
        diffg = diff * prefac
        return diffg
    
    ########################################
    # Getter Methods
    ########################################

    # Return index of highest occupied band
    def get_hob(self):
        if self.spin_up is None and self.spin_down is None:
            raise ValueError("Either spin up or spin down information is missing! \n Make sure that you have set LORBIT = 2 and LSORBIT = True")
        
        ids = [np.argwhere(self.spin_up[:,1] == 1).argmax(), np.argwhere(self.spin_down[:,1] == 1).argmax()]
        return [self.spin_up[ids[0],0], self.spin_down[ids[1], 0]]
    
    # Return index of lowest unoccupied band
    def get_lub(self):
        if self.spin_up is None and self.spin_down is None:
            raise ValueError("Either spin up or spin down information is missing! \n Make sure that you have set LORBIT = 2 and LSORBIT = True")
        
        ids = [np.argwhere(self.spin_up[:,1] == 0).argmin(), np.argwhere(self.spin_down[:,1] == 0).argmin()]
        return [self.spin_up[ids[0],0], self.spin_down[ids[1], 0]]
    
    # Return length
    def get_len(self):
        if self.spin_up is not None:
            return len(self.spin_up)
        elif self.spin_down is not None:
            return len(self.spin_down)
        elif self.spin_neut is not None:
            return len(self.spin_neut)
        else:
            raise ValueError("Please create object with a valid eigenval before trying to access the length")

    ########################################
    # Static Methods
    ########################################

    @staticmethod
    def _reorder_sum_by_magmom(arrs, moments, inv=False):
        p_mask = moments > 0
        n_mask = moments < 0
        o_mask = moments == 0

        p_sums = [np.sum(arr.reshape(-1, 9)[p_mask], axis=0) for arr in arrs]
        n_sums = [np.sum(arr.reshape(-1, 9)[n_mask], axis=0) for arr in arrs]
        o_sums = [np.sum(arr.reshape(-1, 9)[o_mask], axis=0) for arr in arrs]

        if inv:
            outs = np.stack([np.concatenate([n, p, o]) for n, p, o in zip(n_sums, p_sums, o_sums)])
        else:
            outs = np.stack([np.concatenate([p, n, o]) for p, n, o in zip(p_sums, n_sums, o_sums)])
        return outs

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

KPOINTS file handler

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''


class KPOINTS:
    def __init__(self, divisions:int=1, name:str = None):
        """
        Args:
            divisions: Number of divisions between two points given in the KPOINTS file
        """
        self.name = name
        self.data = {}
        self.div = divisions
        self.basis = 0
        self.rec_basis = 0
        self.efermi = None
    
    def __iter__(self):
        return self
    
    def __getitem__(self, items):
        if type(items) == str:
            return self.data[items]
        if type(items) == int:
            keys = [*self.data]
            return self.data[keys[items]]
    
    def __str__(self):
        name = "NAME"
        weight = "WEIGHT"
        coordinates = "COORDINATES"
        string = f"{name:<15} {weight:<15} {coordinates:<40}\n"
        for key in self.data.keys():
            string += f"{key:<15} {str(self.data[key].weight):<15} {str(self.data[key].coordinates):<40}\n"
        return string

    ########################################
    # Adding Methods
    ########################################

    def add_point(self, kpoint):
        self.data[kpoint.name] = kpoint
    
    ########################################
    # Getter Methods
    ########################################
    
    def get_path(self):
        coord = np.stack([value.coordinates for value in self.data.values()])
        if self.rec_basis is not None:
            coord = np.dot(coord, self.rec_basis)
        X = coord[:,0]
        Y = coord[:,1]
        Z = coord[:,2]
        return (X, Y, Z)

    # Returns data nicely formatted for a Bandstructure plot
    def get_band_structure(self):
        bs_up = np.vstack([value.spin_up[:, 0] - self.efermi for value in self.data.values()])
        bs_down = np.vstack([value.spin_down[:, 0] - self.efermi for value in self.data.values()])
        kpoints = np.arange(0, len(self.data) + self.div, self.div) - 1
        kpoints[0] = 0
        return np.arange(len(self.data)), kpoints, bs_up.T, bs_down.T

    # Difference between highest occupied spin up and down energy
    def get_splitting_normal_highest(self):
        length = len(self.data)
        diff = [value.get_hob()[0] - value.get_hob()[1] for value in self.data.values()]
        coords = self.get_path()
        diff = np.array(diff).reshape(-1, 1)
        diff = diff / length
        return np.vstack([coords.T, diff.T]).T

    # Difference between lowest unoccupied spin up and down energy
    def get_splitting_normal_lowest(self):
        length = len(self.data)
        diff = [value.get_lub()[0] - value.get_lub()[1] for value in self.data.values()]
        coords = self.get_path()
        diff = np.array(diff).reshape(-1, 1)
        diff = diff / length
        return np.vstack([coords.T, diff.T]).T

    # Difference between spin up and down energy
    def get_splitting_normal(self, occupied:bool, absolute:bool, cutoff=1e-3):
        length = len(self.data)
        diff = [value.calc_spin_diff(occupied=occupied, cutoff=cutoff, absolute=absolute) for value in self.data.values()]
        X, Y, Z = self.get_path()
        diff = np.array(diff)
        diff = diff / length
        return (X, Y, Z, diff)
    
    # Difference between spin up and down energy with hungarian algorithm
    def get_splitting_hungarian(self, occupied:bool, absolute:bool, moment, cutoff=2):
        length = len(self.data)
        diff = [value.calc_hun_sort_diff(occupied=occupied, moment=moment, absolute=absolute, cutoff=cutoff) for value in self.data.values()]
        X, Y, Z = self.get_path()
        diff = np.array(diff)
        diff = diff / length
        return (X, Y, Z, diff)

    ########################################
    # Loading Methods
    ########################################

    
    def load(self, xml, modifier=''):
        spin_1 = xml.findall(".//calculation/eigenvalues/array/set/set[@comment='spin 1']/")
        spin_2 = xml.findall(".//calculation/eigenvalues/array/set/set[@comment='spin 2']/")
        kpointlist = xml.findall(".//kpoints/varray[@name='kpointlist']/")
        weights = xml.findall(".//kpoints/varray[@name='weights']/")
        kpointnum = xml.find(".//kpoints/generation/*[@name='divisions']")

        tree = xml.find('.//calculation/dos/i')
        self.efermi = float(tree.text.strip())
        
        try:
            self.div = np.prod([int(i) for i in kpointnum.text.split()])
        except:
            self.div = 1

        for i in range(len(weights)):
            spin1 = []
            spin2 = []
            for j in range(len(spin_1[i])):
                spin1.append([float(x) for x in spin_1[i][j].text.split()])
                spin2.append([float(x) for x in spin_2[i][j].text.split()])
            kpoint = KPOINT(name=modifier+spin_1[i].attrib['comment'], coords=[float(i) for i in kpointlist[i].text.split()], weight=float(weights[i].text.strip()), spin_up=np.array(spin1), spin_down=np.array(spin2))
            self.add_point(kpoint)
      
        basis = xml.iterfind('.//structure[@name="initialpos"]/crystal/varray[@name="rec_basis"]/')
        rec_vec = []
        for base in basis:
            rec_vec.append([float(i) for i in base.text.split()])
        self.rec_basis = np.array(rec_vec)
    
    def load_site_projections(self, xml, modifier=''):
        fields = xml.findall('.//projected/array/field')
        cols = [a.text.strip() for a in fields]
        set = xml.findall('.//projected/array/set/')
        momo = xml.findall('.//incar/v[@name="MAGMOM"]')
        for moms in momo:
            mom = [float(x) for x in moms.text.split()]

        for spins in set:
            spin = spins.attrib['comment']
            print(spin)
            if '1' in spin:
                spi = 1
            elif '2' in spin:
                spi = -1
            for point in spins:
                kpi = point.attrib['comment']
                for band in point:
                    bdi = band.attrib['comment']
                    out = np.array([list(map(lambda x: float(x), e.text.split())) for e in band])
                    #self.data[kpi].add_site_projections(spi, bdi, out, cols)
                    self.data[modifier+kpi].add_site_projections(spi, bdi, out, np.array(mom))