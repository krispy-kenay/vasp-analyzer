import numpy as np

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

DOS file handler

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class DOSC:
    def __init__(self, dostype:str):
        self.spin_up = {}
        self.spin_neut = {}
        self.spin_down = {}
        self.dostype = dostype
        self.efermi = None
        self.elements = None
    
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

    def load(self, xml):
        tree = xml.iterfind('.//calculation/dos/')
        cat = [i.tag for i in tree]

        if 'i' in cat:
            tree = xml.find('.//calculation/dos/i')
            self.efermi = float(tree.text.strip())
            tree2 = xml.findall('.//atominfo/array[@name="atoms"]/set/rc/c')
            atoms = [q.text.split() for q in tree2]
            self.elements = np.array(atoms).reshape(-1,2)
        
        if self.dostype == 'total':
            cols = [col.text.strip() for col in xml.iterfind('.//calculation/dos/total/array/field')]
            cols.remove('energy')
            spins = [j.attrib['comment'] for j in xml.find('.//calculation/dos/total/array/set')]
            for spin in spins:
                arr = np.array([k.text.split() for k in xml.findall('.//calculation/dos/total/array/set/set[@comment="' + spin + '"]/')], dtype=float)
                ats = np.array(atoms).reshape(-1,2)
                dos = DOS(dostype='total', ion='ion 0', spin=spin, columns=cols, data=arr[:,1:], energies=arr[:,0], elements=np.unique(ats[:,0]))
                self.add_DOS(dos=dos)

        if self.dostype == 'partial':
            cols = [col.text.strip() for col in xml.iterfind('.//calculation/dos/partial/array/field')]
            cols.remove('energy')
            ions = [i.attrib['comment'] for i in xml.findall('.//calculation/dos/partial/array/set/')]
            spins = [j.attrib['comment'] for j in xml.find('.//calculation/dos/partial/array/set/set')]
            j = 0
            for ion in ions:
                for spin in spins:
                    arr = np.array([k.text.split() for k in xml.findall('.//calculation/dos/partial/array/set/set[@comment="' + ion + '"]/set[@comment="' + spin + '"]/')], dtype=float)
                    dos = DOS(dostype='partial', ion=ion, spin=spin, columns=cols, data=arr[:,1:], energies=arr[:,0], elements=atoms[j][0])
                    self.add_DOS(dos=dos)
                j += 2