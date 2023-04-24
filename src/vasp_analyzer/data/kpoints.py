import numpy as np

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
    
    def add_coordinates(self, coords):
        self.coordinates = coords
    
    def add_spin(self, val, spin):
        if val == 0:
            self.spin_neut = spin
        if val == 1:
            self.spin_up = spin
        if val == -1:
            self.spin_down = spin
    
    def add_spin_pro(self, val, spin):
        if val == 0:
            self.spin_neut_pro = spin
        if val == 1:
            self.spin_up_pro = spin
        if val == -1:
            self.spin_down_pro = spin
    
    def add_weight(self, wgt):
        self.weight = wgt
    
    # Simple calculation of up spin energies - down spin energies (no weighing of orbital character)
    def calc_spin_diff(self, occupied=True, cutoff=1e-3):
        if self.spin_up is not None and self.spin_down is not None:
            if occupied == True:
                diff = np.sum(self.spin_up[self.spin_up[:, 1] == 1, 0] - self.spin_down[self.spin_down[:, 1] == 1, 0])
            else:
                diff = np.sum(self.spin_up - self.spin_down)
            if abs(diff) < cutoff:
                diff = 0
            return diff
        else:
            raise ValueError("Either spin up or spin down information is missing! \n Make sure that you have set LORBIT = 2 and LSORBIT = True")
    
    # Return index of highest occupied band
    def get_hob(self):
        ids = [np.argwhere(self.spin_up[:,1] == 1).argmax(), np.argwhere(self.spin_down[:,1] == 1).argmax()]
        return [self.spin_up[ids[0],0], self.spin_down[ids[1], 0]]
    
    # Return index of lowest unoccupied band
    def get_lub(self):
        ids = [np.argwhere(self.spin_up[:,1] == 0).argmin(), np.argwhere(self.spin_down[:,1] == 0).argmin()]
        return [self.spin_up[ids[0],0], self.spin_down[ids[1], 0]]

    










    
    def get_len(self):
        if self.spin_up is not None:
            return len(self.spin_up)
        elif self.spin_down is not None:
            return len(self.spin_down)
        elif self.spin_neut is not None:
            return len(self.spin_neut)
        else:
            raise ValueError("Please create object with a valid eigenval before trying to access the length")
    
    def calc_spin_del_hun(self):
        from scipy.optimize import linear_sum_assignment
        cost_matrix = np.zeros((len(self.spin_up[self.spin_up[:, 1] == 1]), len(self.spin_down[self.spin_down[:, 1] == 1])))

        for i in range(len(self.spin_up[self.spin_up[:, 1] == 1])):
            for j in range(len(self.spin_down[self.spin_down[:, 1] == 1])):
                if abs(self.spin_up[i, 0] - self.spin_down[j, 0]) > 1 or np.dot(self.spin_up_site[i], self.spin_down_site[j]) == 0:
                    cost_matrix[i, j] = 0
                else:
                    cost_matrix[i, j] = - np.dot(self.spin_up_site[i], self.spin_down_site[j])
        
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        diff = 0
        for q in range(len(row_ind)):
            y = row_ind[q]
            x = col_ind[q]
            #
            #
            # SHOULD DOT PRODUCT BE INCLUDED?????????
            # np.dot(self.spin_up_site[y], self.spin_down_site[x])
            #
            diff += (self.spin_up[y, 0] - self.spin_down[x, 0])
        
        prefac = 1 / (len(self.spin_up[self.spin_up[:, 1] == 1]))
        diffg = diff * prefac
        return diffg
    
    def calc_spin_del_hun_mod(self):
        from scipy.optimize import linear_sum_assignment
        cost_matrix = np.zeros((len(self.spin_up[self.spin_up[:, 1] == 1]), len(self.spin_down[self.spin_down[:, 1] == 1])))

        for i in range(len(self.spin_up[self.spin_up[:, 1] == 1])):
            for j in range(len(self.spin_down[self.spin_down[:, 1] == 1])):
                if abs(self.spin_up[i, 0] - self.spin_down[j, 0]) > 1 or np.dot(self.spin_up_site[i], self.spin_down_site[j]) == 0:
                    cost_matrix[i, j] = 0
                else:
                    cost_matrix[i, j] = - np.dot(self.spin_up_site[i], self.spin_down_site[j])
        
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        diff = 0
        for q in range(len(row_ind)):
            y = row_ind[q]
            x = col_ind[q]
            #
            #
            # SHOULD DOT PRODUCT BE INCLUDED?????????
            # np.dot(self.spin_up_site[y], self.spin_down_site[x])
            #
            diff += (self.spin_up[y, 0] - self.spin_down[x, 0])
        
        prefac = 1 / (len(self.spin_up[self.spin_up[:, 1] == 1]))
        diffg = diff * prefac
        return diffg
    
    def calc_hun_pos(self):
        from scipy.optimize import linear_sum_assignment
        cost_matrix = np.zeros((len(self.spin_up[self.spin_up[:, 1] == 1]), len(self.spin_down[self.spin_down[:, 1] == 1])))

        for i in range(len(self.spin_up[self.spin_up[:, 1] == 1])):
            for j in range(len(self.spin_down[self.spin_down[:, 1] == 1])):
                if abs(self.spin_up[i, 0] - self.spin_down[j, 0]) > 3 or np.dot(self.spin_up_site[i], self.spin_down_site[j]) == 0:
                    cost_matrix[i, j] = 0
                else:
                    cost_matrix[i, j] = - np.dot(self.spin_up_site[i], self.spin_down_site[j])
        
        row_ind, col_ind = linear_sum_assignment(cost_matrix)

        out = []
        for q in range(len(row_ind)):
            y = row_ind[q]
            x = col_ind[q]
            out.append([self.spin_up[y], self.spin_down[x]])
        return out
    



    
    @staticmethod
    def _arr_converter_sublattice(arr, moments, inv=False):
        #arr = arr.reshape((-1, 9))
        #plat = np.zeros((1,9))
        #nlat = np.zeros((1,9))
        #onat = np.zeros((1,9))
        p = np.sum(arr[np.argwhere(moments > 0).reshape(-1)], axis=0)
        n = np.sum(arr[np.argwhere(moments < 0).reshape(-1)], axis=0)
        o = np.sum(arr[np.argwhere(moments == 0).reshape(-1)], axis=0)
        
        if inv==False:
            out = np.concatenate([p, n, o])
        else:
            out = np.concatenate([n, p, o])
        return out.reshape((-1))
    
    @staticmethod
    def _arr_order_by_magmom(arr, moments, inv=False):
        p = arr[np.argwhere(moments > 0).reshape(-1)]
        n = arr[np.argwhere(moments < 0).reshape(-1)]
        o = arr[np.argwhere(moments == 0).reshape(-1)]
        
        if inv==False:
            out = np.concatenate([p, n, o],axis=0)
        else:
            out = np.concatenate([n, p, o],axis=0)
        return out.reshape((-1))
    
    def calc_gap_diff(self, efermi, cutoff=1e-3, mode_alt=False):
        id_up = [np.abs(self.spin_up[np.argwhere(self.spin_up[:,0] > efermi),0] - efermi).argmin(), np.abs(self.spin_up[np.argwhere(self.spin_up[:,0] < efermi),0] - efermi).argmin()]
        id_down = [np.abs(self.spin_down[np.argwhere(self.spin_down[:,0] > efermi),0] - efermi).argmin(), np.abs(self.spin_down[np.argwhere(self.spin_down[:,0] < efermi),0] - efermi).argmin()]
        
        data_up = [self.spin_up[np.argwhere(self.spin_up[:,0] > efermi),0][id_up[0]], self.spin_up[np.argwhere(self.spin_up[:,0] < efermi),0][id_up[1]]]
        data_down = [self.spin_down[np.argwhere(self.spin_down[:,0] > efermi),0][id_down[0]], self.spin_down[np.argwhere(self.spin_down[:,0] < efermi),0][id_down[1]]]
        
        if mode_alt == False:
            A = data_up[0] - data_up[1]
            B = data_down[0] - data_down[1]
        else:
            A = data_up[0] - data_down[1]
            B = data_up[1] - data_down[0]
        return float(A - B)
    
    def add_site_projections(self, sign, band, data, mom):
        # Sum over all ions
        #dat = np.sum(data, axis=0)

        # Sum over all orbitals
        #dat = np.sum(data, axis=1)

        # Full array
        #dat = np.reshape(data, (-1))

        # Temp sum over elements
        #dat1 = np.sum(data[:12], axis=0)
        #dat2 = np.sum(data[12:], axis=0)
        #dat = np.concatenate([dat1, dat2])

        # Sum by spin
        if sign == 1:
            dat = self._arr_converter_sublattice(data, mom)
        else:
            dat = self._arr_converter_sublattice(data, mom, inv=True)
        
        index = int(band.replace('band ', '')) - 1
        length = self.get_len()

        if sign == 1:
            if not self.spin_up_site:
                self.spin_up_site = [None] * length
            self.spin_up_site[index] = dat
        
        if sign == -1:
            if not self.spin_down_site:
                self.spin_down_site = [None] * length
            self.spin_down_site[index] = dat
        
        if sign == 0:
            if not self.spin_neut_site:
                self.spin_neut_site = [None] * length
            self.spin_neut_site[index] = dat

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

    def add_point(self, kpoint):
        self.data[kpoint.name] = kpoint

    def get_path(self):
        coord = []
        for key in self.data.keys():
            coord.append(self.data[key].coordinates)
        arr = np.array(coord)
        if type(self.rec_basis) != int:
            arr = np.dot(arr, self.rec_basis)
        return arr
    
    
    def get_tot_diff(self): #Using projections
        diff = []
        length = len(self.data)
        for key, value in self.data.items():
            val = value.calc_spin_del()
            diff.append(val / length)
        coords = self.get_path()
        diff = np.array(diff).reshape(-1,1)
        return np.concatenate([coords, diff], axis=1)
    
    def get_proj_diff(self): #Using projections
        diff = []
        length = len(self.data)
        for key, value in self.data.items():
            val = value.calc_spin_del_hun()
            diff.append(val / length)
        coords = self.get_path()
        diff = np.array(diff).reshape(-1,1)
        return np.concatenate([coords, diff], axis=1)

    
    def get_homb_lumb(self, efermi, cutoff=1e-3, mode_alt=False): # Trash I think
        diff = []
        for key, value in self.data.items():
            diff.append(value.calc_gap_diff(efermi, cutoff, mode_alt))
        coords = self.get_path()
        diff = np.array(diff).reshape(-1,1)
        return np.concatenate([coords, diff], axis=1)






    # FINAL?
    # FINAL?
    # FINAL?
    # Returns data nicely formatted for a Bandstructure plot
    def get_band_structure(self):
        bs_up = np.vstack([value.spin_up[:, 0] - self.efermi for value in self.data.values()])
        bs_down = np.vstack([value.spin_down[:, 0] - self.efermi for value in self.data.values()])
        kpoints = np.arange(0, len(self.data) + self.div, self.div) - 1
        kpoints[0] = 0
        return np.arange(len(self.data)), kpoints, bs_up.T, bs_down.T

    # Calculates difference between up spin energy and hdown spin energy (in different variations). 
    def get_spin_difference_hl(self, mode, cutoff=1e-3):
        length = len(self.data)
        if mode == 'highest':
            diff = [value.get_hob()[0] - value.get_hob()[1] for value in self.data.values()]
        elif mode == 'lowest':
            diff = [value.get_lub()[0] - value.get_lub()[1] for value in self.data.values()]
        elif mode == 'occupied':
            diff = [value.calc_spin_diff(occupied=True, cutoff=cutoff) for value in self.data.values()]
        elif mode == 'total':
            diff = [value.calc_spin_diff(occupied=False, cutoff=cutoff) for value in self.data.values()]
        elif mode == 'occupied_hungarian':
            diff = [value.calc_spin_del_hun() for value in self.data.values()]
        else:
            raise ValueError("provide a valid mode!")
        coords = self.get_path()
        diff = np.array(diff).reshape(-1, 1)
        diff = diff / length
        return np.vstack([coords.T, diff.T]).T







    
    def get_single_band(self, efermi, mode = 'highest'): # Useless??
        band_up = []
        band_down = []
        kpoints = []
        i = 0
        if mode == 'highest':
            for key, value in self.data.items():
                hobs = value.get_hob()
                band_up.append(hobs[0] - efermi)
                band_down.append(hobs[1] - efermi)
                if i == 0 or ((i+1) % self.div == 0):
                    kpoints.append(i)
                i += 1
            points = [*range(len(self.data))]
            return points, kpoints, np.array([band_up]), np.array([band_down])
    
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