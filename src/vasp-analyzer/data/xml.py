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