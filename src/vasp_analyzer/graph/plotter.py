import plotly
import plotly.graph_objects as go
import numpy as np
from scipy import interpolate

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

VASP Plotter

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class vPlot:
    def __init__(self, size=None, width=None, height=None):
        """
        Plotly plotter object to facilitate plotting 3D interactive graphs. Includes helper functions to draw hexagonal brillouin zone, create surface and scatter plots and more. Relies on the offline module to prevent sending data to external cloud servers
        """
        plotly.offline.init_notebook_mode()
        self.fig = go.Figure()

        if size is not None:
            self.w = size
            self.h = size
        else:
            self.w = width
            self.h = height

        self.layout = go.Layout(
            width=width,
            height=height)

        self.traces = {}
        self.global_id = 3
        self._hidden_vars = []
    
    @staticmethod
    def interpolate_surface(xi, yi, wi, density:int=500, method:str='linear'):

        xxi = np.linspace(min(xi), max(xi), density)
        yyi = np.linspace(min(yi), max(yi), density)
        xx, yy = np.meshgrid(xxi, yyi)

        zz = interpolate.griddata((xi, yi), wi, (xx, yy), method=method)
        return xxi, yyi, zz
    
    @staticmethod
    def hex_to_RGB(hex_str):
        return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

    @staticmethod
    def get_color_gradient(c1, c2, n):
        assert n > 1
        c1_rgb = np.array(vPlot.hex_to_RGB(c1))/255
        c2_rgb = np.array(vPlot.hex_to_RGB(c2))/255
        mix_pcts = [x/(n-1) for x in range(n)]
        rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
        return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]