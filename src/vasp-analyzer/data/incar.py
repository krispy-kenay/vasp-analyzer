'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

INCAR object

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class INCAR:
    def __init__(self):
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

    def load(self, xml):
        data_big = xml.findall('.//incar/')

        for data in data_big:
            field = data.attrib['name']
            try:
                types = data.attrib['type']
            except:
                types = 'float'
            if types == 'logical':
                if data.text == 'T':
                    self[field] = True
                elif data.text == 'F':
                    self[field] = False
            if types == 'string':
                self[field] = str(data.text.strip())
            if types == 'int':
                self[field] = [int(i) for i in data.text.split()]
            if types == 'float':
                self[field] = [float(i) for i in data.text.split()]
            if types == 'list':
                self[field] = [float(i) for i in data.text.split()]