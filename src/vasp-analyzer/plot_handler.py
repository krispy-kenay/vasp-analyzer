import plotly
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

VARIABLES

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''


CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
CB91_Grad_BP = ['#2cbdfe', '#2fb9fc', '#33b4fa', '#36b0f8',
                '#3aacf6', '#3da8f4', '#41a3f2', '#449ff0',
                '#489bee', '#4b97ec', '#4f92ea', '#528ee8',
                '#568ae6', '#5986e4', '#5c81e2', '#607de0',
                '#6379de', '#6775dc', '#6a70da', '#6e6cd8',
                '#7168d7', '#7564d5', '#785fd3', '#7c5bd1',
                '#7f57cf', '#8353cd', '#864ecb', '#894ac9',
                '#8d46c7', '#9042c5', '#943dc3', '#9739c1',
                '#9b35bf', '#9e31bd', '#a22cbb', '#a528b9',
                '#a924b7', '#ac20b5', '#b01bb3', '#b317b1']

P_black = 'rgba(60, 60, 60, 1)'
P_black_hex = '#3A3736'
P_jasper = '#BF4E30'
P_coral = '#FF7F50'
P_r = '#f26419'

P_green = '#3A8B3C'

P_blue = '#4357AD'
P_violet = '#9F5975'
P_tomato = '#FA5B3D'

P_burgundy = '#D7263D'
P_silver = '#5D89BA'

col = {'Fe ':P_silver, 
       'O ':P_burgundy,
       'Cr ':P_green,
       's':P_blue,
       'p':P_violet,
       'd':P_tomato,
       'Fe s':'#5070B4',
       'Fe p':'#7E7198',
       'Fe d':'#AC727C',
       'O s':'#8D3F75',
       'O p':'#BB4059',
       'O d':'#E9413D',
       'Cr s':'#3F7175',
       'Cr p':'#6D7259',
       'Cr d':'#9A733D'
       }

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

3D Plotter

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class hPlot:
    def __init__(self, width=None, height=None):
        """
        Plotly plotter object to facilitate plotting 3D interactive graphs. Includes helper functions to draw hexagonal brillouin zone, create surface and scatter plots and more. Relies on the offline module to prevent sending data to external cloud servers
        """
        plotly.offline.init_notebook_mode()
        self.fig = go.Figure()

        self.w = width
        self.h = height
        self.layout = None

    def draw_hexagon(self, basis_vector, verbose:bool=False):
        converter = []
        fac = np.sqrt(3)/3
        rotation = np.pi/3
        scaling_x = np.linalg.norm(basis_vector[0])*fac
        scaling_y = np.linalg.norm(basis_vector[1])*fac
        scaling_z = np.linalg.norm(basis_vector[2])*0.5

        upper_hex = []
        lower_hex = []
        for n in range(7):
            x = scaling_x*np.cos(n*np.pi/3 + rotation)
            y = scaling_y*np.sin(n*np.pi/3 + rotation)
            z = scaling_z
            upper_hex.append([x,y,z])
            lower_hex.append([x,y,-z])
        
        converter.append(np.array(upper_hex))
        converter.append(np.array(lower_hex))

        for n in range(6):
            x = scaling_x*np.cos(n*np.pi/3 + rotation)
            y = scaling_y*np.sin(n*np.pi/3 + rotation)
            z_up = scaling_z
            z_down = -scaling_z
            coord = [[x,y,z_up],[x,y,z_down]]
            converter.append(np.array(coord))
        
        traces = []
        for conv in converter:
            trace = go.Scatter3d(x=conv[:,0], y=conv[:,1], z=conv[:,2], mode='lines', showlegend=False, line=dict(
                color='black',
                dash='dash'
                ))
            traces.append(trace)
        
        if verbose == False:
            for trc in traces:
                self.fig.add_trace(trc)
        else:
            return traces
    
    def add_3d_scatter_weight(self,
                           data,
                           name:str='Data',
                           color_scale:list=[(0,"red"), (0.5,"rgba(0,0,0,0.1)"), (1,"blue")],
                           opacity:float=0.9,
                           size_marker:int=5,
                           midpoint:float=0,
                           verbose:bool=False):
        trace = go.Scatter3d(
            x=data[:,0],
            y=data[:,1],
            z=data[:,2],
            mode='markers',
            name=name,
            marker=dict(
                size=size_marker,
                color=data[:,3],
                colorscale=color_scale,
                cmid = midpoint,
                opacity = opacity,
                colorbar = dict(
                    title='Strength',
                    xanchor='right',
                    lenmode = 'fraction'
                    )
                )
                )
        
        if verbose == False:
            self.fig.add_trace(trace)
        else:
            return trace
    
    def add_surface_weight(self,
                           data,
                           z_explicit:float=None,
                           density:int=500,
                           method:str='linear',
                           color_scale:list=[(0,"red"), (0.5,"rgba(200,200,200,1)"), (1,"blue")],
                           midpoint:float=0,
                           contour:bool=True,
                           contour_step:int=8):
        z_list = np.unique(data[:,2])
        for z_expl in z_list:
            xx, yy, zz = self.interpolate_surface(data, z_exp=z_expl, density=density, meth=method)
            trace = go.Surface(
            visible=False,
            x=xx,
            y=yy,
            z=zz,
            colorscale=color_scale,
            cmid=midpoint,
            contours = {
                "z": {"show": contour, "start": float(np.nanmin(zz)), "end": float(np.nanmax(zz)), "size": float((np.nanmax(zz)-np.nanmin(zz))/contour_step), 'usecolormap':True}})
            self.fig.add_trace(trace)
        self.fig.data[0].visible = True

        steps = []
        for i in range(len(self.fig.data)):
            step = dict(method='update', args=[{'visible': [False]*len(self.fig.data)}])
            step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
            steps.append(step)
        
        sliders = [dict(
            active=10,
            currentvalue={"prefix": "height-step: " + str(z_list[i])},
            pad={"t": 50},
            steps=steps)]
        
        self.fig.update_layout(
            sliders=sliders
            )


    def define_layout(self,
                      title_text:str='unknown',
                      margins:dict={'l': 0, 'r': 0, 'b': 0, 't': 0},
                      color_background:str='rgba(245,245,245,1)',
                      color_gridlines:str=P_black,
                      color_grid:str='rgba(235, 235, 235, 1)',
                      font_family:str='Courier New, monospace',
                      font_size:int=12,
                      font_color:str=P_black
                      ):
        layout = go.Layout(
            width=self.w,
            height=self.h,
            margin=margins,
            paper_bgcolor=color_background,
            title=dict(text=title_text, yanchor='top', y=0.95, xanchor='center', x=0.5),
            titlefont=dict(family=font_family, size=font_size*2, color=font_color),
            font=dict(family=font_family, size=font_size, color=font_color),
            legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.01),
            scene=dict(
                xaxis=dict(backgroundcolor=color_grid, gridcolor=color_gridlines, showline=True, nticks=5),
                yaxis=dict(backgroundcolor=color_grid, gridcolor=color_gridlines, showline=True, nticks=5),
                zaxis=dict(backgroundcolor=color_grid, gridcolor=color_gridlines, showline=True, nticks=5))
        )
        self.fig.update_layout(layout)
        self.fig.update_xaxes(showline=True)
        self.fig.update_layout(yaxis_range=[-1.2,1.2])
        self.layout = 1
    
    def interpolate_surface(self, data, z_ind=0, z_exp=None, density=600, meth='linear'):
        if z_exp:
            z = z_exp
        else:
            z_list = np.unique(data[:,2])
            z = z_list[z_ind]
        loc = np.argwhere(data[:,2] == z)
        xy_data = data[loc, np.array([0,1,3])]

        xi = np.linspace(min(xy_data[:,0]), max(xy_data[:,0]), density)
        yi = np.linspace(min(xy_data[:,1]), max(xy_data[:,1]), density)
        xi, yi = np.meshgrid(xi, yi)

        zi = interpolate.griddata((xy_data[:,0] ,xy_data[:,1]), xy_data[:,2], (xi, yi), method=meth)
        return xi, yi, zi
    
    def plot(self,
             filename:str=None,
             render:str="vscode"):
        if not self.layout:
            self.define_layout()
        if filename:
            if "html" in filename:
                self.fig.write_html("pictures/" + filename)
            else:
                self.fig.write_image("pictures/" + filename)
        self.fig.show(renderer=render)
        #plotly.offline.iplot(self.fig)

'''
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

3D Plotter

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
'''

class vPlot:
    def __init__(self,size=10):
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec

        golden = (1+ np.sqrt(5))/2
        self.fig = plt.figure(figsize=(size, size))
        #heights = [1]
        #widths = [golden, 1]
        #self.gs = GridSpec(1, 2, figure=self.fig, width_ratios=widths, height_ratios=heights)
        self.gs = GridSpec(1, 3, figure=self.fig)
        self.bs_trigger = False
        self.dos_trigger = False
        self.bz_trigger = False
        

        plt.rcParams['figure.facecolor'] = '#F5F5F5'
        plt.rcParams['axes.linewidth'] = 0.5
        plt.rcParams['axes.spines.left'] = False
        plt.rcParams['axes.spines.right'] = False
        plt.rcParams['axes.spines.top'] = False
    
    def init_dos(self):
        self.ax2 = self.fig.add_subplot(self.gs[0, -1], label='DOS')
        self.dos_trigger = True
    
    def init_bs(self):
        self.ax1 = self.fig.add_subplot(self.gs[0, :-1], label='BS')
        self.bs_trigger = True
    
    def init_bz(self):
        self.ax3 = self.fig.add_subplot(self.gs[0, :], label='BZ', projection='3d')
        self.bz_trigger = True
    
    def add_bs(self,
               points,
               kpoints,
               bs_up,
               bs_down,
               c_grad:list = None,
               c_single:str = None
               ):
        
        if self.bs_trigger == False:
            self.init_bs()
        
        if c_grad:
            colors = iter(self.get_color_gradient(c_grad[0], c_grad[1], len(bs_up)))
        elif c_single:
            colors = iter([c_single, c_single])

        for i in range(len(bs_up)):
            c = next(colors)
            self.ax1.plot(points, bs_up[i], c=c, linestyle='-')
            self.ax1.plot(points, bs_down[i], c=c, linestyle=':')
        
        self.ax1.vlines(kpoints, np.min(bs_up), np.max(bs_up), colors=P_black_hex, linestyle='--')
    
    def add_dos(self,
                dos:dict,
                energies:list,
                c_list:list=None,
                c_grad:list=None,
                c_dict:dict=None,
                fill:bool=False):
        
        if self.dos_trigger == False:
            self.init_dos()
        
        spins = [*dos.keys()]
        tags = [list(dos[i].keys()) for i in spins]
        tags = np.unique(tags)
        
        if c_grad != None:
            colors = iter(self.get_color_gradient(c_grad[0], c_grad[1], len(tags)))
        elif c_list != None:
            colors = iter(c_list)
        else:
            colors = iter(self.get_color_gradient('#3A3736', '#3A3736', len(tags)))
        
        for tag in tags:
            c = next(colors)
            for spin in spins:
                if spin == -1:
                    style = (0, (5, 1))
                else:
                    style = '-'
                if c_dict:
                    self.ax2.plot(spin*dos[spin][tag], energies, c=c_dict[tag], linestyle=style, label=tag, linewidth=2)
                else:
                    self.ax2.plot(spin*dos[spin][tag], energies, c=c, linestyle=style, label=tag, linewidth=2)

                
                if fill == True:
                    self.ax2.fill_between(spin*dos[spin][tag], energies, color=c, alpha=0.8)
    
    def add_kpoints(self, 
                    coords, 
                    weights:list=0,
                    colors=None):

        if self.bz_trigger == False:
            self.init_bz()
        
        if type(weights) != int:
            self.kpoint_weight_add(coords, weights, colors)
        else:
            #self.ax3.scatter3D(coords[:,0],coords[:,1],coords[:,2],c = weights, cmap = 'coolwarm')
            self.ax3.scatter3D(coords[:,0],coords[:,1],coords[:,2])
    
    def kpoint_weight_add(self, coords, weights, colors):
        mm = max(np.abs(weights))
        for i in range(len(weights)):
            if weights[i] > 0:
                self.ax3.scatter3D(coords[i,0], coords[i,1], coords[i,2], c = colors[0], alpha=weights[i]/mm)
            elif weights[i] < 0:
                self.ax3.scatter3D(coords[i,0], coords[i,1], coords[i,2], c = colors[1], alpha=abs(weights[i]/mm))
            else:
                self.ax3.scatter3D(coords[i,0], coords[i,1], coords[i,2], c = 'black', alpha=colors[2])

    def prop_bs(self, x_lim, y_lim, title, x_label,y_label, fntsize):
        self.ax1.set_xlabel(x_label, fontsize=fntsize)
        self.ax1.set_ylabel(y_label, fontsize=fntsize)
        self.ax1.set_xlim(x_lim)
        self.ax1.set_ylim(y_lim)

        self.ax1.set_xticklabels([])
        self.ax1.grid(visible=True, which='major', axis='y')

        self.ax1.set_title(title, fontsize=fntsize)

        self.ax1.axhline(linewidth=2, color=P_black_hex)

        self.ax1.set_facecolor('#F5F5F5')

    def prop_dos(self, x_lim, y_lim, title, x_label,y_label, fntsize):
        self.ax2.set_xlabel(x_label, fontsize=fntsize)
        self.ax2.set_ylabel(y_label, fontsize=fntsize)

        self.ax2.set_ylim(y_lim)
        if x_lim == None:
            lin = self.ax2.get_lines()
            xd = lin[0].get_xdata()
            yd = lin[0].get_ydata()
            lo, hi = self.ax2.get_ylim()
            x_displayed = xd[((yd>lo) & (yd<hi))]
            h = np.max(x_displayed) - np.min(x_displayed)
            bot = np.min(x_displayed)-0.1*h
            top = np.max(x_displayed)+0.1*h
            self.ax2.set_xlim(-top, top)
        else:
            self.ax2.set_xlim(x_lim)
        
        self.ax2.grid(visible=True, which='major', axis='y')
        if y_label == "":
            self.ax2.set_yticklabels([])
            #self.ax2.set_yticks([])
        

        self.ax2.axvline(linewidth=2, color=P_black_hex) #adds thick red line @ y=0
        self.ax2.axhline(linewidth=2, color=P_black_hex)

        self.ax2.set_xticklabels([])
        

        self.ax2.set_title(title, fontsize=fntsize)

        handles, labels = self.ax2.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        self.ax2.legend(by_label.values(), by_label.keys(), fontsize=fntsize/1.5)
        
        self.ax2.set_facecolor('#F5F5F5')
    
    def prop_bz(self, azimuth, elevation):
        self.ax3.set_facecolor('#F5F5F5')

        self.ax3.xaxis.pane.fill = False
        self.ax3.yaxis.pane.fill = False
        self.ax3.zaxis.pane.fill = False

        self.ax3.xaxis.pane.set_edgecolor('w')
        self.ax3.yaxis.pane.set_edgecolor('w')
        self.ax3.zaxis.pane.set_edgecolor('w')

        self.ax3.view_init(azim=azimuth, elev=elevation)
    
    def hex_to_RGB(self, hex_str):
        return [int(hex_str[i:i+2], 16) for i in range(1,6,2)]

    def get_color_gradient(self, c1, c2, n):
        assert n > 1
        c1_rgb = np.array(self.hex_to_RGB(c1))/255
        c2_rgb = np.array(self.hex_to_RGB(c2))/255
        mix_pcts = [x/(n-1) for x in range(n)]
        rgb_colors = [((1-mix)*c1_rgb + (mix*c2_rgb)) for mix in mix_pcts]
        return ["#" + "".join([format(int(round(val*255)), "02x") for val in item]) for item in rgb_colors]
    
    def draw_hex_brillouin(self, rotation:float=np.pi/3, scaling_x=1, scaling_y=1, scaling_z=1, basis_vector=0):
        fac = np.sqrt(3)/3
        if type(basis_vector) != int:
            scaling_x = np.linalg.norm(basis_vector[0])*fac
            scaling_y = np.linalg.norm(basis_vector[1])*fac
            scaling_z = np.linalg.norm(basis_vector[2])*0.5
        else:
            scaling_x = scaling_x * fac
            scaling_y = scaling_y * fac
            scaling_z = scaling_z * 0.5
        
        for n in range(6):
            ls = ':'
            #scaling = np.sqrt(3)/3
            #rotation = np.pi/3
            # Lower Hexagon
            self.ax3.plot3D([scaling_x*np.cos(n*np.pi/3 + rotation), scaling_x*np.cos((n+1)*np.pi/3 + rotation)], [scaling_y*np.sin(n*np.pi/3 + rotation), scaling_y*np.sin((n+1)*np.pi/3 + rotation)], [-scaling_z,-scaling_z], color='k', linestyle=ls)
            # Upper Hexagon
            self.ax3.plot3D([scaling_x*np.cos(n*np.pi/3 + rotation), scaling_x*np.cos((n+1)*np.pi/3 + rotation)], [scaling_y*np.sin(n*np.pi/3 + rotation), scaling_y*np.sin((n+1)*np.pi/3 + rotation)], [scaling_z,scaling_z], color='k', linestyle=ls)
            # Vertical connections
            self.ax3.plot3D([scaling_x*np.cos(n*np.pi/3 + rotation), scaling_x*np.cos(n*np.pi/3 + rotation)], [scaling_y*np.sin(n*np.pi/3 + rotation), scaling_y*np.sin(n*np.pi/3 + rotation)], [-scaling_z,scaling_z], color='k', linestyle=ls)
    
    def plot_all(self,
                 filename=None,
                 transparent=False,
                 fontsize=20,
                 y_lim:list=None,
                 y_label:str="$E-E_F$ / eV",
                 bs_xlim:list=None,
                 bs_xlabel:str="$k$ Vector",
                 bs_title:str='Band Structure',
                 dos_xlim:list=None,
                 dos_xlabel:str="states / eV",
                 dos_title:str='Density of States',
                 azimuth:float=0,
                 elevation:int=90
                 ):        

        if (self.bs_trigger == True and self.dos_trigger == True):
            self.prop_bs(bs_xlim,y_lim, bs_title, bs_xlabel, y_label, fontsize)
            self.prop_dos(dos_xlim,y_lim, dos_title, dos_xlabel, "", fontsize)
        
        elif self.bs_trigger == True:
            self.prop_bs(bs_xlim,y_lim, bs_title, bs_xlabel, y_label, fontsize)
        elif self.dos_trigger == True:
            self.prop_dos(dos_xlim,y_lim, dos_title, dos_xlabel, y_label, fontsize)
        elif self.bz_trigger == True:
            self.prop_bz(azimuth, elevation)
        if filename:
            plt.savefig('pictures/' + filename, transparent=transparent, bbox_inches = "tight", dpi=1000)
            print('Saved as ' + filename + '.png')
        plt.show()