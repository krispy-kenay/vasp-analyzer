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
    
    ########################################
    # Adding Methods
    ########################################

    def add_bs(self, data:tuple,
               color_gradient:list = None,
               Type:str='bands',
               custom_labels:list=None,
               custom_labels_size:int=18,
               **kwargs):
        
        # Unpack data and check type
        if type(data) != tuple: raise ValueError("Provide data as a tuple of format (data, energies)")
        points, kpoints, bs_up, bs_down = data

        # Setup traces if not present
        if Type not in self.traces: self.traces[Type] = []
        else: raise Warning('Band structure is already present, use with caution!')

        # Add 2D tag
        self._hidden_vars.append('2D')

        # Handle color input
        if color_gradient is not None: colors = iter(self.get_color_gradient(color_gradient[0], color_gradient[1], len(bs_up)))
        else: colors = iter(self.get_color_gradient('#3A3736', '#3A3736', len(bs_up)))

        for i in range(len(bs_up)):
            c = next(colors)
            if i == 0: sl = True
            else: sl = False
            trace_up = go.Scatter(x=points, y=bs_up[i], line=dict(color=c, dash='solid'), mode='lines', name='spin up', legendgroup='1', showlegend=sl, legendgrouptitle_text='Bands', **kwargs)
            trace_down = go.Scatter(x=points, y=bs_down[i], line=dict(color=c, dash='5px,1px'), mode='lines', name='spin down', legendgroup='2', showlegend=sl, **kwargs)
            self.traces[Type].append(trace_up)
            self.traces[Type].append(trace_down)

        for kpoint in kpoints:
            trace_line = go.Scatter(x=[kpoint, kpoint], y=[np.min(bs_up) -10, np.max(bs_up)+10], mode='lines', line=dict(color='#3A3736', dash='dash'), showlegend=False)
            self.traces[Type].append(trace_line)
        
        # Handle labels for the kpoints
        if custom_labels is not None:
            self.fig.update_xaxes(tickmode='array', tickvals=kpoints, ticktext=custom_labels, tickfont=dict(size=custom_labels_size))
    
    def add_dos(self, data:tuple,
                color:str=None,
                color_dict:dict=None,
                fill:bool=False,
                Type:str='dos',
                **kwargs):
        
        # Unpack data and check type
        if type(data) != tuple: raise ValueError("Provide data as a tuple of format (data, energies)")
        energies, dos = data

        # Setup traces if not present
        if Type not in self.traces: self.traces[Type] = []

        # Add 2D tag
        self._hidden_vars.append('2D')
        
        # Reformat data
        spins = [*dos.keys()]
        tags = [list(dos[i].keys()) for i in spins]
        tags = np.unique(tags)

        # Check if should be filled
        if fill == True: filler = 'tozerox'
        else: filler = None

        for tag in tags:
            # Handle color input
            if color is not None: c = color
            elif color_dict is not None: c = color_dict[tag]
            else: c = 'lightgrey'
            # Assign Legend title only if it hasn't been before
            if 'dos' not in self._hidden_vars: 
                grptl = 'DOS'
                self._hidden_vars.append('dos')
            else: grptl = None
            # Handle the different spins
            for spin in spins:
                if spin == -1:
                    trace = go.Scatter(x=spin*dos[spin][tag], y=energies, line=dict(color=c, dash='5px,1px'), mode='lines', name=tag, fill=filler, legendgroup=str(self.global_id), showlegend=False, **kwargs)
                else:
                    trace = go.Scatter(x=spin*dos[spin][tag], y=energies, line=dict(color=c, dash='solid'), mode='lines', name=tag, fill=filler, legendgroup=str(self.global_id), legendgrouptitle_text=grptl, **kwargs)
                self.traces[Type].append(trace)
            # Increase global id (Hacky workaround for plotly)
            self.global_id += 1
    
    def add_3d_scatter_weight(self, data:tuple,
                           name:str='Data',
                           Type:str='scatter3D',
                           color_scale:str=None,
                           crange:tuple=None,
                           opacity:float=0.9,
                           size_marker:int=5):
        
        # Check data type and unpack into separate values
        if type(data) != tuple: raise ValueError("Provide data as a tuple of format (x, y, z, weight)")
        x, y, z, weight = data

        # Setup traces if not present
        if Type not in self.traces: self.traces[Type] = []

        # Add 3D tag
        self._hidden_vars.append('3D')

        # max, min and midpoints can be manually added (e.g. if scales should be the same between many different plots), if not then it inputs the default arguments
        if crange is None: crange = (min(weight), 0, max(weight))

        # Checks for custom colorscale, if none is provided it will automatically select a colorscale depending on the type of data
        if color_scale is None:
            if min(weight) >= 0: color_scale = 'inferno'
            else: color_scale = ['rgb(142, 181, 194)', 'rgb(69, 144, 185)', 'rgb(11, 102, 189)', 'rgb(41, 58, 143)', 'rgb(23, 28, 66)', 'rgb(23, 23, 23)', 'rgb(60, 9, 17)', 'rgb(120, 14, 40)', 'rgb(172, 43, 36)', 'rgb(196, 101, 72)', 'rgb(213, 157, 137)']
        
        # create a scatter trace with the data
        trace = go.Scatter3d(x=x, y=y, z=z, mode='markers', name=name,
            marker=dict(size=size_marker, color=weight, colorscale=color_scale, cmin=crange[0], cmid=crange[1], cmax=crange[2], opacity = opacity,
                colorbar = dict(title='Strength', xanchor='right', lenmode = 'fraction')))
        
        # add trace to the figure
        self.traces[Type].append(trace)

    ########################################
    # Plotting Methods
    ########################################

    # Automatic plotter that should get the correct plot type under default conditions (useful when not trying to modify much and just getting a plot)
    def plot(self, **kwargs):
        if self.plot_bands(sup_error=True, **kwargs) or self.plot_dos(sup_error=True, **kwargs) or self.plot_bsdos(sup_error=True, **kwargs) or self.plot_single(sup_error=True, **kwargs) or self.plot_multiple(sup_error=True, **kwargs):
            pass
        else:
            raise RuntimeError("Automatic plotting failed!")

    # Is called by all other plotting functions, checks for a file name and passes keyword arguments to the general layout.
    def _finish_plot(self, filename:str=None, render:str="vscode", show:bool=True, **kwargs):
        self.define_layout(**kwargs)

        if filename:
            if "html" in filename:
                self.fig.write_html(filename)
            else:
                self.fig.write_image(filename)
        
        if show == True:
            self.fig.show(renderer=render)

    # Plot density of states next to band structure with shared y-axis
    def plot_bsdos(self, sup_error=False, **kwargs):
        if 'bands' not in self.traces or 'dos' not in self.traces or len(self.traces) != 2:
            if sup_error == False: raise ValueError("Add only density of states and band structure!")
            else: return False

        self.fig.set_subplots(1,2, column_widths=[0.7, 0.3], subplot_titles=('Band Structure',  'Density of States'), shared_yaxes=True)
        self.fig.add_traces(self.traces['bands'], rows=1, cols=1)
        spacer1 = go.Scatter(x=[1,2,3], y=[1,2,3], visible='legendonly', mode='lines', name=' ', line=dict(color='rgba(5,5,5,0)'), legendgroup='spacer1')
        self.fig.add_trace(spacer1, row=1, col=1)
        self.fig.add_traces(self.traces['dos'], rows=1, cols=2)
        self.fig.update_xaxes(showgrid=False)

        self.fig['layout']['yaxis']['title']['text'] = 'E - Ef / eV'
        self.fig['layout']['xaxis']['title']['text'] = 'k Vector'
        self.fig['layout']['xaxis2']['title']['text'] = 'States / eV'

        self._finish_plot(**kwargs)

        return True

    # Plot density of states
    def plot_dos(self, sup_error=False, **kwargs):
        if 'dos' not in self.traces or len(self.traces) != 1:
            if sup_error == False: raise ValueError("Add only density of states!")
            else: return False
        
        self.fig.set_subplots(1,1, subplot_titles=('Density of States',))
        self.fig.add_traces(self.traces['dos'], rows=1, cols=1)
        self.fig.update_xaxes(showgrid=False)
        

        self.fig['layout']['yaxis']['title']['text'] = 'E - Ef / eV'
        self.fig['layout']['xaxis']['title'] = 'States / eV'

        self._finish_plot(**kwargs)

        return True
    
    # Plot band structure
    def plot_bands(self, sup_error=False, **kwargs):
        if 'bands' not in self.traces or len(self.traces) != 1:
            if sup_error == False: raise ValueError("Add only band structure!")
            else: return False

        self.fig.set_subplots(1,1, subplot_titles=('Band Structure',))
        self.fig.add_traces(self.traces['bands'], rows=1, cols=1)
        self.fig.update_xaxes(showgrid=False)

        self.fig['layout']['yaxis']['title']['text'] = 'E - Ef / eV'
        self.fig['layout']['xaxis']['title'] = 'k Vector'

        self._finish_plot(**kwargs)

        return True
    
    # Plot one type of data without specific adjustments
    def plot_single(self, sup_error=False, **kwargs):
        if len(self.traces) != 1:
            if sup_error == False: raise ValueError("Make sure only one type of data is added!")
            else: return False
        
        self.fig.add_traces(*self.traces.values())

        self._finish_plot(**kwargs)

        return True

    # Plot multiple sets of different data
    def plot_multiple(self, sup_error=False, custom_titles:tuple=None, **kwargs):
        if len(self.traces) == 0:
            if sup_error == False: raise ValueError("Add at least some data!")
            else: return False

        sums = len(self.traces)
        i = 0

        self.fig.set_subplots((sums//2) + (sums%2>0), 2, subplot_titles=custom_titles, vertical_spacing=0.18)

        for key, value in self.traces.items():
            row = (i // 2) + 1
            col = (i % 2) + 1 
            self.fig.add_traces(value, rows=row, cols=col)
            i += 1
        self._finish_plot(**kwargs)

        return True
    
    ########################################
    # Define Layout
    ########################################
    
    def define_layout(self, title:str='Title',
                              margins:dict={'l': 50, 'r': 20, 'b': 50, 't': 50, 'pad':5},
                              font_size:int=12,
                              font_family:str='Courier New, monospace',
                              font_color:str='rgba(60, 60, 60, 1)',
                              color_background:str='rgba(245, 245, 245, 1)',
                              color_gridlines:str='rgba(220, 220, 220, 1)',
                              color_grid:str='rgba(0, 0, 0, 0)',
                              width_line:int=1.5,
                              ylim:list=None,
                              xlim:list=None,
                              **kwargs):
        
        margins['t'] = 60 + font_size
        layout = go.Layout(
            autosize=True,
            width=self.w,
            height=self.h,
            margin=margins,
            paper_bgcolor=color_background,
            title=dict(text=title, yanchor='top', y=1, xanchor='center', x=0.5, pad={'t': 15}),
            titlefont=dict(family=font_family, size=font_size*2, color=font_color),
            font=dict(family=font_family, size=font_size, color=font_color),
            legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.01, tracegroupgap=5),
            plot_bgcolor=color_grid,**kwargs)
        self.fig.update_layout(layout)

        if '2D' in self._hidden_vars:
            self.fig.update_xaxes(dict(gridcolor=color_gridlines, zerolinecolor=color_gridlines, gridwidth=width_line, zerolinewidth=2*width_line),
                                  range=xlim)
            self.fig.update_yaxes(dict(gridcolor=color_gridlines, zerolinecolor=color_gridlines, gridwidth=width_line, zerolinewidth=2*width_line),
                                  range=ylim)
        if '3D' in self._hidden_vars:
            self.fig.update_scenes(dict(
                xaxis=dict(backgroundcolor=color_grid, gridcolor=color_gridlines, zerolinecolor=color_gridlines, gridwidth=width_line, zerolinewidth=2*width_line, showline=True, nticks=5),
                yaxis=dict(backgroundcolor=color_grid, gridcolor=color_gridlines, zerolinecolor=color_gridlines, gridwidth=width_line, zerolinewidth=2*width_line, showline=True, nticks=5),
                zaxis=dict(backgroundcolor=color_grid, gridcolor=color_gridlines, zerolinecolor=color_gridlines, gridwidth=width_line, zerolinewidth=2*width_line, showline=True, nticks=5)))


    ########################################
    # Static Methods
    ########################################

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