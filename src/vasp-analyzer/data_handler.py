import numpy as np

'''
--- VARIABLES ---
'''

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

XML handler

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class xml_handler:
    def __init__(self, filename):
        from lxml import etree as ET
        #import xml.etree.ElementTree as ET
        self._tree = ET.parse(filename)
    
    def iterfind(self, finder):
        return self._tree.iterfind(finder)
    
    def iter(self, finder):
        return self._tree.iter(finder)
    
    def gener(self, generator, mode):
        data = []

    def find(self, finder):
        return self._tree.find(finder)
    
    def findall(self, finder):
        return self._tree.findall(finder)
    
    def close(self):
        self._tree = None

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

INCAR object

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class INCAR:
    def __init__(self,
                 ISTART=0,
                 ICHARG=2):
        self.container = {}
    
    def __setitem__(self, items, value):
        if type(value) == list:
            if len(value) == 1:
                value = value[0]
        if type(items) == str:
            self.container[items] = value
    
    def __getitem__(self, items):
        if type(items) == str:
            return self.container[items]
        else:
            raise ValueError("wrong input formatting")
    
    def __str__(self):
        string = ""
        for key in self.container.keys():
            out_key = str(key) + ":"
            string += f"{out_key:<15} {str(self.container[key]):>40}\n"
        return string
    
    def load(self):
        pass

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
    
    def add_hash_projections(self, sign, band, data):
        dat = np.sum(data, axis=0)
        #dat = np.reshape(data, (-1))
        
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
    
    def add_site_projections(self, sign, band, data, cols):
        dat = self._arr_converter_simplifier(data, cols)
        index = int(band.replace('band ', '')) - 1
        length = self.get_len()

        if sign == 1:
            if not self.spin_up_site:
                #self.spin_up_site = np.empty((self.get_len()))
                self.spin_up_site = [None] * length
            self.spin_up_site[index] = dat
        
        if sign == -1:
            if not self.spin_down_site:
                #self.spin_down_site = np.empty((self.get_len()))
                self.spin_down_site = [None] * length
            self.spin_down_site[index] = dat
        
        if sign == 0:
            if not self.spin_neut_site:
                #self.spin_neut_site = np.empty((self.get_len()))
                self.spin_neut_site = [None] * length
            self.spin_neut_site[index] = dat
            

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

KPOINTS object handler

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
    
    def load_site_projections(self, xml):
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
                    self.data[kpi].add_hash_projections(spi, bdi, out)

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
    
    def load(self, xml, verb=False, modifier=''):
        spin_1 = xml.findall(".//calculation/eigenvalues/array/set/set[@comment='spin 1']/")
        spin_2 = xml.findall(".//calculation/eigenvalues/array/set/set[@comment='spin 2']/")
        kpointlist = xml.findall(".//kpoints/varray[@name='kpointlist']/")
        weights = xml.findall(".//kpoints/varray[@name='weights']/")
        kpointnum = xml.find(".//kpoints/generation/*[@name='divisions']")
        
        try:
            self.divisions = np.prod([int(i) for i in kpointnum.text.split()])
        except:
            self.divisions = 1
        
        if len(kpointlist) >= 1 and verb == True:
            print("KPOINTS found!")

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
                

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

DOS object

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class DOS:
    def __init__(self, 
                 dostype:str,
                 ion:str,
                 spin:int,
                 columns:list,
                 data,
                 energies=None,
                 elements:str='undefined',
                 ):
        if dostype != 'total' and dostype != 'partial':
            raise TypeError("BandStructure needs to be either total or partial type!")
        
        if type(spin) == str:
            if spin == 'spin 1':
                self.spin = 1
            elif spin == 'spin 2':
                self.spin = -1
            else:
                self.spin = 0
        else:
            self.spin = spin
        
        self.type=dostype


        self.ion=int(ion.strip("ion "))
        self.elements=elements
        self.energy=energies

        self.container={}

        for i in range(len(columns)):
            self.container[columns[i]] = data[:,i]
    
    def get_spd(self):
        s = 0
        p = 0
        d = 0
        for key, value in self.container.items():
            if 's' in key:
                if type(s) == int:
                    s = np.copy(value)
                else:
                    s += value
            if 'p' in key:
                if type(p) == int:
                    p = np.copy(value)
                else:
                    p += value
            if 'd' in key or 'x2-y2' in key:
                if type(d) == int:
                    d = np.copy(value)
                else:
                    d += value
        return s, p, d

    def __getitem__(self, items):
        if items == 'type':
            return self.type
        elif items == 'ion':
            return self.ion
        elif items == 'spin':
            return self.spin
        elif items == 'elements':
            return self.elements
        elif items == 'energy':
            return self.energy
        elif type(items) == str:
            return self.container[items]
        else:
            data = self.container[items[0]]
            return data[items[1:]]
    
    def __str__(self):
      dtype = "Type:"
      ion = "Ion:"
      spin = "Spin:"
      element = "Element:"
      string = f"{dtype:<15} {self.type:>15}\n{ion:<15} {self.ion:>15}\n{spin:<15} {self.spin:>15}\n{element:<15} {self.elements:>15}"
      return string

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

DOS object handler

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class DOSC:
    def __init__(self, dostype:str):
        self.spin_up = {}
        self.spin_neut = {}
        self.spin_down = {}
        self.dostype = dostype
    
    def __getitem__(self, items):
        if items == 1:
          return self.spin_up
        if items == 0:
          return self.spin_neut
        if items == -1:
          return self.spin_down
    
    def add_DOS(self, dos):
        if dos['type'] != self.dostype:
            raise TypeError("types need to match between different DOS types!")
        if dos['spin'] == 1:
            self.spin_up[dos.ion] = dos
        if dos['spin'] == 0:
            self.spin_neut[dos.ion] = dos
        if dos['spin'] == -1:
            self.spin_down[dos.ion] = dos
    
    def get_spins(self):
      out = []
      if self.spin_up != {}:
        out.append(1)
      if self.spin_neut != {}:
        out.append(0)
      if self.spin_down != {}:
        out.append(-1)
      return out
    
    def calc_dos(self, typ):
        spins = self.get_spins()
        out = {}

        if 'spd' in typ:
            spd = ['s','p','d']
        else:
            spd = ['', '', '']
        
        if 'el' in typ:
            elements = self.get_elements()
            elements = [k + ' ' for k in elements]
        else:
            elements = ['']

        for i in spins:
            out[i] = {}
            for el in elements:
                for orb in spd:
                    out[i][el + orb] = 0
            for key, value in self[i].items():
                s, p, d = value.get_spd()
                dat_spd = [s,p,d]
                for j in range(3):
                    if 'el' in typ:
                        ele = value.elements + ' '
                    else:
                        ele = ''
                    out[i][ele + spd[j]] = self._combine_arr(out[i][ele + spd[j]], dat_spd[j])
        return out
    
    def get_energy(self, efermi):
       spins = self.get_spins()
       for key in self[spins[0]].keys():
          return self[spins[0]][key].energy - efermi
    
    def _combine_arr(self, to_add, add):
        if type(to_add) == int:
            to_add = np.copy(add)
        else:
            to_add += add
        return to_add

    def get_elements(self):
        spins = self.get_spins()
        out = []
        for i in spins:
            for key, value in self[i].items():
                if value.elements not in out:
                    out.append(value.elements)
        return out
    
    def get_tdos(self):
       spins = self.get_spins()

       out = {}
       for i in spins:
          out[i] = {}
          for key in self[i].keys():
             out[i]['total'] = self[i][key].container['total']
       return out

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

Main VASP handler
Most other classes are unified within this to give a simplified and unified interface

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class VASP:
    def __init__(self,
               total_dos=None,
               partial_dos=None,
               k_points=None,
               in_car=None,
               eferm=None,
               verbose=False
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
        self.tdos = total_dos
        self.pdos = partial_dos
        self.kpoints = k_points
        self.incar = in_car
        self.efermi = eferm
        self.elements = None
        self.verb = verbose

        # Set up triggers
        self.trigger_incar = False
        self.trigger_doscar = False
        self.trigger_kpoints = False
    
    def load(self,
             inc:bool=False,
             dosc:bool=False,
             kp:bool=False
             ):
        if inc == True:
            self.trigger_incar = True
            pass
        if dosc == True:
            self.trigger_doscar = True
            pass
        if kp == True:
            self.trigger_kpoints = True
            pass
    
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

    def load_INCAR(self, xml):
        self.incar = INCAR()

        data_big = xml.findall('.//incar/')
        if len(data_big) >= 1 and self.verb == True:
            print("INCAR found!")
        for data in data_big:
            field = data.attrib['name']
            try:
                types = data.attrib['type']
            except:
                types = 'float'
            if types == 'logical':
                if data.text == 'T':
                    self.incar[field] = True
                elif data.text == 'F':
                    self.incar[field] = False
            if types == 'string':
                self.incar[field] = str(data.text.strip())
            if types == 'int':
                self.incar[field] = [int(i) for i in data.text.split()]
            if types == 'float':
                self.incar[field] = [float(i) for i in data.text.split()]
            if types == 'list':
                self.incar[field] = [float(i) for i in data.text.split()]

    def load_DOSCAR(self, xml):
        tree = xml.iterfind('.//calculation/dos/')
        cat = [i.tag for i in tree]

        if 'i' in cat:
            tree = xml.find('.//calculation/dos/i')
            self.efermi = float(tree.text.strip())
            tree2 = xml.findall('.//atominfo/array[@name="atoms"]/set/rc/c')
            atoms = [q.text.split() for q in tree2]
            self.elements = np.array(atoms).reshape(-1,2)
        
        if 'total' in cat:
            if self.verb == True:
                print("Total DOS found!")
            self.tdos = DOSC(dostype='total')
            cols = [col.text.strip() for col in xml.iterfind('.//calculation/dos/total/array/field')]
            cols.remove('energy')
            spins = [j.attrib['comment'] for j in xml.find('.//calculation/dos/total/array/set')]
            for spin in spins:
                arr = np.array([k.text.split() for k in xml.findall('.//calculation/dos/total/array/set/set[@comment="' + spin + '"]/')], dtype=float)
                dos = DOS(dostype='total', ion='ion 0', spin=spin, columns=cols, data=arr[:,1:], energies=arr[:,0])
                self.tdos.add_DOS(dos=dos)
      
        if 'partial' in cat:
            if self.verb == True:
                print("Partial DOS found!")
            self.pdos = DOSC(dostype='partial')
            cols = [col.text.strip() for col in xml.iterfind('.//calculation/dos/partial/array/field')]
            cols.remove('energy')
            ions = [i.attrib['comment'] for i in xml.findall('.//calculation/dos/partial/array/set/')]
            spins = [j.attrib['comment'] for j in xml.find('.//calculation/dos/partial/array/set/set')]
            j = 0
            for ion in ions:
                for spin in spins:
                    arr = np.array([k.text.split() for k in xml.findall('.//calculation/dos/partial/array/set/set[@comment="' + ion + '"]/set[@comment="' + spin + '"]/')], dtype=float)
                    dos = DOS(dostype='partial', ion=ion, spin=spin, columns=cols, data=arr[:,1:], energies=arr[:,0], elements=atoms[j][0])
                    self.pdos.add_DOS(dos=dos)
                j += 2

    def load_KPOINTS(self, xml, modifier=''):
        spin_1 = xml.findall(".//calculation/eigenvalues/array/set/set[@comment='spin 1']/")
        spin_2 = xml.findall(".//calculation/eigenvalues/array/set/set[@comment='spin 2']/")
        kpointlist = xml.findall(".//kpoints/varray[@name='kpointlist']/")
        weights = xml.findall(".//kpoints/varray[@name='weights']/")
        kpointnum = xml.find(".//kpoints/generation/*[@name='divisions']")

        if not self.kpoints:
            try:
                self.kpoints = KPOINTS(divisions=np.prod([int(i) for i in kpointnum.text.split()]))
            except:
                self.kpoints = KPOINTS(divisions=1)
        
        if len(kpointlist) >= 1 and self.verb == True:
            print("KPOINTS found!")

        for i in range(len(weights)):
            spin1 = []
            spin2 = []
            for j in range(len(spin_1[i])):
                spin1.append([float(i) for i in spin_1[i][j].text.split()])
                spin2.append([float(i) for i in spin_2[i][j].text.split()])
            kpoint = KPOINT(name=modifier+spin_1[i].attrib['comment'], coords=[float(i) for i in kpointlist[i].text.split()], weight=float(weights[i].text.strip()), spin_up=np.array(spin1), spin_down=np.array(spin2))
            self.kpoints.add_point(kpoint)
      
        basis = xml.iterfind('.//structure[@name="initialpos"]/crystal/varray[@name="rec_basis"]/')
        rec_vec = []
        for base in basis:
            rec_vec.append([float(i) for i in base.text.split()])
      
        self.kpoints.rec_basis = np.array(rec_vec)

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