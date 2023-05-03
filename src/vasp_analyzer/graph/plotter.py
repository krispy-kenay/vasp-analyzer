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
    def __init__(self, width=None, height=None):
        """
        Plotly plotter object to facilitate plotting interactive graphs. Includes helper functions to draw hexagonal brillouin zone, create surface and scatter plots and more. Relies on the offline module to prevent sending data to external cloud servers
        """
        plotly.offline.init_notebook_mode()
        self.fig = go.Figure()

        self.w = width
        self.h = height
    
    @staticmethod
    def interpolate_surface(xi, yi, wi, density:int=500, method:str='linear'):

        xxi = np.linspace(min(xi), max(xi), density)
        yyi = np.linspace(min(yi), max(yi), density)
        xx, yy = np.meshgrid(xxi, yyi)

        zz = interpolate.griddata((xi, yi), wi, (xx, yy), method=method)
        return xxi, yyi, zz