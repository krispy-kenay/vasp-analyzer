import numpy as np
import plotly
import plotly.graph_objects as go
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
        Plotly plotter object to facilitate plotting 3D interactive graphs. Includes helper functions to draw hexagonal brillouin zone, create surface and scatter plots and more. Relies on the offline module to prevent sending data to external cloud servers
        """
        plotly.offline.init_notebook_mode()
        self.fig = go.Figure()

        self.w = width
        self.h = height
        self.layout = None