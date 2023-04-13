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
        
        self.spin_diff = None
        self.spin_diff_abs = None

        self.spin_up_site = None
        self.spin_down_site = None
        self.spin_neut_site = None
    
    def __str__(self):
        name = "Name:"
        coords = "Coordinates:"
        wgt = "Weight:"
        return  f"{name:<15}{self.name:>40} \n{coords:<15}{str(self.coordinates):>40} \n{wgt:<15}{self.weight:>40}"
    
    def calc_spin_diff(self, efermi, cutoff=1e-3):
        if self.spin_diff == None:
            spin_diff = np.sum(self.spin_up[np.all(self.spin_up <= efermi, axis=1)] - self.spin_down[np.all(self.spin_down <= efermi, axis=1)])
            if abs(spin_diff) > cutoff:
                self.spin_diff = spin_diff
            else:
                self.spin_diff = 0
    
    def get_hobs(self):
        ids = [np.argwhere(self.spin_up[:,1] == 1).argmax(), np.argwhere(self.spin_down[:,1] == 1).argmax()]
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
    
    def calc_spin_del(self, absolute=False):
        diff = 0
        if self.spin_up is not None and self.spin_down is not None:
            for i in range(len(self.spin_up)):
                for j in range(len(self.spin_down)):
                    if self.spin_up[i, 1] == 1 and self.spin_down[j, 1] == 1:
                        diff += ((self.spin_up[i, 0]-self.spin_down[j, 0]) * np.dot(self.spin_up_site[i],self.spin_down_site[j]))
                        #* np.abs((self.spin_up[i, 0]-self.spin_down[j, 0]) * np.dot(self.spin_up_site[i],self.spin_down_site[j]))
        prefac = 1/(np.count_nonzero(self.spin_up[:,1]) * len(self.spin_up))
        diffg = diff * prefac
        return diffg
    
    def calc_spin_del2(self, absolute=False):
        ids = [np.argwhere(self.spin_up[:,1] == 1).argmax(), np.argwhere(self.spin_down[:,1] == 1).argmax()]
        diff = 0
        if self.spin_up is not None and self.spin_down is not None:
            for i in range(len(self.spin_up)):
                if self.spin_up[i, 1] == 1:
                    diff += ((self.spin_up[i, 0]-self.spin_down[ids[1], 0]) * np.dot(self.spin_up_site[i],self.spin_down_site[ids[1]]))

            for j in range(len(self.spin_down)):
                if self.spin_down[j, 1] == 1:
                    diff += ((self.spin_up[ids[0], 0]-self.spin_down[j, 0]) * np.dot(self.spin_up_site[ids[0]],self.spin_down_site[j]))

        prefac = 1/(np.count_nonzero(self.spin_up[:,1]) * len(self.spin_up))
        diffg = diff * prefac
        return diffg

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
    
    @staticmethod
    def _arr_converter_simplifier(data, cols):
        yx = np.vstack(np.nonzero(data)).T
        out = []
        for y, x in yx:
            string = 'ion '+ str(y+1) + ' : ' + str(data[y][x]) + ' ' + str(cols[x])
            out.append(string)
        return out
    
    def add_site_projections(self, sign, band, data):
        # Sum over all ions
        #dat = np.sum(data, axis=0)

        # Full array
        #dat = np.reshape(data, (-1))

        # Temp sum over elements
        dat1 = np.sum(data[:18], axis=0)
        dat2 = np.sum(data[18:], axis=0)
        dat = np.concatenate([dat1, dat2])
        
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
    
    def get_hob_diff(self):
        diff = []
        for key, value in self.data.items():
            hobs = value.get_hobs()
            diff.append(hobs[0] - hobs[1])
        coords = self.get_path()
        diff = np.array(diff).reshape(-1,1)
        return np.concatenate([coords, diff], axis=1)
    
    def get_tot_diff(self):
        diff = []
        for key, value in self.data.items():
            diff.append(value.calc_spin_del2())
        coords = self.get_path()
        diff = np.array(diff).reshape(-1,1)
        return np.concatenate([coords, diff], axis=1)

    
    def get_homb_lumb(self, efermi, cutoff=1e-3, mode_alt=False):
        diff = []
        for key, value in self.data.items():
            diff.append(value.calc_gap_diff(efermi, cutoff, mode_alt))
        coords = self.get_path()
        diff = np.array(diff).reshape(-1,1)
        return np.concatenate([coords, diff], axis=1)
    
    def get_band_structure(self, efermi):
        bs_up = []
        bs_down = []
        kpoints = []
        i = 0
        for key in self.data.keys():
            bs_up.append(self.data[key].spin_up[:,0] - efermi)
            bs_down.append(self.data[key].spin_down[:,0] - efermi)
            if i == 0 or ((i+1) % self.div == 0):
                kpoints.append(i)
            i += 1
        
        points = [*range(len(self.data))]

        return points, kpoints, np.array(bs_up).transpose(), np.array(bs_down).transpose()
    
    def get_single_band(self, efermi, mode = 'highest'):
        band_up = []
        band_down = []
        kpoints = []
        i = 0
        if mode == 'highest':
            for key, value in self.data.items():
                hobs = value.get_hobs()
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
        
        try:
            self.divisions = np.prod([int(i) for i in kpointnum.text.split()])
        except:
            self.divisions = 1

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
                    self.data[modifier+kpi].add_site_projections(spi, bdi, out)