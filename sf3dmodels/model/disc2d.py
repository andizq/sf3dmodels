"""
2D Disc models
==============
Classes: Rosenfeld2d, General2d, Velocity, Intensity, Cube, Tools
"""
#TODO in show(): Scroll over channels, 2nd slidder for several cubes
#TODO in Cube: Add convolve() function  
from ..utils.constants import G, kb
from ..utils import units as u
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import numbers
import warnings
import os

#warnings.filterwarnings("error")
__all__ = ['Cube', 'Tools', 'Intensity', 'Velocity', 'General2d', 'Rosenfeld2d']

matplotlib.rcParams['font.family'] = 'monospace'
matplotlib.rcParams['font.weight'] = 'normal'
matplotlib.rcParams['lines.linewidth'] = 1.5
matplotlib.rcParams['axes.linewidth'] = 3.0
matplotlib.rcParams['xtick.major.width']=1.6
matplotlib.rcParams['ytick.major.width']=1.6

SMALL_SIZE = 10
MEDIUM_SIZE = 15
BIGGER_SIZE = 22

matplotlib.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
matplotlib.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
matplotlib.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
matplotlib.rc('xtick', labelsize=MEDIUM_SIZE-2)    # fontsize of the tick labels
matplotlib.rc('ytick', labelsize=MEDIUM_SIZE-2)    # fontsize of the tick labels
matplotlib.rc('legend', fontsize=SMALL_SIZE-1)    # legend fontsize
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

params = {'xtick.major.size': 6.5,
          'ytick.major.size': 6.5
          }

matplotlib.rcParams.update(params)


class Cube(object):
    def __init__(self, nchan, channels, data):
        self.nchan = nchan
        self.channels = channels
        self.data = data
        self._interactive = self.cursor

    @property
    def interactive(self): 
        return self._interactive
          
    @interactive.setter 
    def interactive(self, func): 
        print('Setting interactive function to', func) 
        self._interactive = func

    @interactive.deleter 
    def interactive(self): 
        print('Deleting interactive function') 
        del self._interactive
    
    def _plot_region(self, x0, x1, y0, y1, ax, extent=None, stat_func=np.mean):
        print (x0,x1,y0,y1)
        def get_ji(x,y):
            pass
        if extent is None:
            j0, i0 = int(x0), int(y0)
            j1, i1 = int(x1), int(y1)
        else: 
            nz, ny, nx = np.shape(self.data)
            dx = extent[1] - extent[0]
            dy = extent[3] - extent[2]
            j0 = int(nx*(x0-extent[0])/dx)
            i0 = int(ny*(y0-extent[2])/dy)
            j1 = int(nx*(x1-extent[0])/dx)
            i1 = int(ny*(y1-extent[2])/dy)

        slice_cube = self.data[:,i0:i1,j0:j1]
        spectrum = np.array([stat_func(chan) for chan in slice_cube])

        v0, v1 = self.channels[0], self.channels[-1]
        if np.logical_or(np.isinf(spectrum), np.isnan(spectrum)).all(): return False
        else:
            ax.fill_between(self.channels, spectrum, alpha=0.1)
            return ax.step(self.channels, spectrum, where='mid', linewidth=2.5, label=r'%d,%d'%(x0,x1))
   
    def box(self, fig, ax, extent=None, stat_func=np.mean):
        from matplotlib.widgets import RectangleSelector
        import matplotlib.patches as patches
        
        def onselect(eclick, erelease):
            if eclick.inaxes is ax[0]:
                plot_spec = self._plot_region(eclick.xdata, erelease.xdata, eclick.ydata, erelease.ydata, 
                                              ax[1], extent=extent, stat_func=stat_func) 
                if plot_spec:
                    print('startposition: (%f, %f)' % (eclick.xdata, eclick.ydata))
                    print('endposition  : (%f, %f)' % (erelease.xdata, erelease.ydata))
                    print('used button  : ', eclick.button)
                    xc, yc = eclick.xdata, eclick.ydata #Left, bottom corner
                    dx, dy = erelease.xdata-eclick.xdata, erelease.ydata-eclick.ydata
                    rect = patches.Rectangle((xc,yc), dx, dy, lw=2, edgecolor=plot_spec[0].get_color(), facecolor='none')
                    ax[0].add_patch(rect)
                    ax[1].legend()
                    fig.canvas.draw()
                    fig.canvas.flush_events()

        def toggle_selector(event):
            print('Key pressed.')
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print('RectangleSelector deactivated.')
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print('RectangleSelector activated.')
                toggle_selector.RS.set_active(True)

        rectprops = dict(facecolor='none', edgecolor = 'white',
                         alpha=0.8, fill=False)

        lineprops = dict(color='white', linestyle='-',
                         linewidth=3, alpha=0.8)

        toggle_selector.RS = RectangleSelector(ax[0], onselect, drawtype='box', rectprops=rectprops, lineprops=lineprops)
        cid = fig.canvas.mpl_connect('key_press_event', toggle_selector)

    def ellipse(self):
        pass

    def convolve(self, kernel=None):
        pass

    def _plot_spectrum(self, x, y, ax, extent=None):
        def get_ji(x,y):
            pass
        if extent is None:
            j, i = int(x), int(y)
        else: 
            nz, ny, nx = np.shape(self.data)
            dx = extent[1] - extent[0]
            dy = extent[3] - extent[2]
            j = int(nx*(x-extent[0])/dx)
            i = int(ny*(y-extent[2])/dy)
            
        spectrum = self.data[:,i,j]
        v0, v1 = self.channels[0], self.channels[-1]
        if np.logical_or(np.isinf(spectrum), np.isnan(spectrum)).all(): return False
        else:
            ax.fill_between(self.channels, spectrum, alpha=0.1)
            return ax.step(self.channels, spectrum, where='mid', linewidth=2.5, label=r'%d,%d'%(x,y))
        
    def cursor(self, fig, ax, extent=None):

        def onclick(event):
            if event.button==3: 
                print ('Right click. Disconnecting click event...')
                fig.canvas.mpl_disconnect(cid)
            elif event.inaxes is ax[0]:
                plot_spec = self._plot_spectrum(event.xdata, event.ydata, ax[1], extent=extent) 
                if plot_spec:
                    print('%s click: button=%d, xdata=%f, ydata=%f' %
                          ('double' if event.dblclick else 'single', event.button,
                           event.xdata, event.ydata))
                    ax[0].scatter(event.xdata, event.ydata)#, color=plot_spec[0].get_color())
                    ax[1].legend()
                    fig.canvas.draw()
                    fig.canvas.flush_events()

        cid = fig.canvas.mpl_connect('button_press_event', onclick)

    def show(self, extent=None, chan_init=20, cursor_grid=True, 
             int_unit=r'Brightness Temperature [K]', pos_unit='au', vel_unit=r'km s$^{-1}$',
             **kwargs):
        from matplotlib.widgets import Slider, Cursor
        v0, v1 = self.channels[0], self.channels[-1]
        dv = v1-v0
        fig, ax = plt.subplots(ncols=2, figsize=(12,5))
        axchan = plt.axes([0.2, 0.9, 0.24, 0.05], facecolor='0.7')
        
        y0, y1 = ax[1].get_position().y0, ax[1].get_position().y1
        axcbar = plt.axes([0.48, y0, 0.03, y1-y0])
        max_data = np.max(self.data)
        ax[0].set_xlabel(pos_unit)
        ax[0].set_ylabel(pos_unit)
        ax[1].set_xlabel('l.o.s velocity [%s]'%vel_unit)
        ax[1].tick_params(direction='in', right=True, labelright=False, labelleft=False)
        axcbar.tick_params(direction='out')
        ax[1].set_ylabel(int_unit, labelpad=15)
        ax[1].yaxis.set_label_position('right')
        ax[1].set_xlim(v0-0.1, v1+0.1)
        ax[1].set_ylim(-1, max_data)
        ax[1].grid(lw=1.5, ls=':')
        cmap = plt.get_cmap('hot')
        cmap.set_bad(color=(0.9,0.9,0.9))

        img = ax[0].imshow(self.data[chan_init], cmap=cmap, extent=extent, origin='lower left', vmin=-1, vmax=max_data)
        cbar = plt.colorbar(img, cax=axcbar)
        current_chan = ax[1].axvline(self.channels[chan_init], color='black', lw=2, ls='--')
        text_chan = ax[1].text((self.channels[chan_init]-v0)/dv, 1.02, #Converting xdata coords to Axes coords 
                               '%4.1f km/s'%self.channels[chan_init], ha='center', 
                               color='black', transform=ax[1].transAxes)

        self.interactive(fig, ax, extent=extent, **kwargs)

        if cursor_grid: cg=Cursor(ax[0], useblit=True, color='lime', linewidth=1.5)
        
        

        slider_chan = Slider(axchan, 'Channel', 0, self.nchan-1, 
                             valstep=1, valinit=chan_init, valfmt='%2d', color='dodgerblue')        

        def update(val):
            chan = int(val)
            vchan = self.channels[chan]
            img.set_data(self.data[chan])
            current_chan.set_xdata(vchan)
            text_chan.set_x((vchan-v0)/dv)
            text_chan.set_text('%4.1f km/s'%vchan)
            fig.canvas.draw_idle()
            
        slider_chan.on_changed(update)
        plt.show()

    def make_fits(self, output, **kw_header):
        from astropy.io import fits
        hdr = fits.Header()
        hdr.update(**kw_header)
        data = np.where(np.isfinite(self.data), self.data, 0)
        fits.writeto(output, data, hdr, overwrite=True)
    
    def make_gif(self, folder='./movie/', extent=None, velocity2d=None, 
                 unit=r'Brightness Temperature [K]',
                 gif_command='convert -delay 10 *int2d* cube_channels.gif'):
        cwd = os.getcwd()
        if folder[-1] != '/': folder+='/'
        os.system('mkdir %s'%folder)
        max_data = np.max(self.data)

        clear_list, coll_list = [], []
        fig, ax = plt.subplots()
        contour_color = 'red'
        cmap = plt.get_cmap('binary')
        cmap.set_bad(color=(0.9,0.9,0.9))
        ax.plot([None],[None], color=contour_color, linestyle='--', linewidth=2, label='Near side') 
        ax.plot([None],[None], color=contour_color, linestyle=':', linewidth=2, label='Far side') 
        ax.set_xlabel('au')
        ax.set_ylabel('au')
        for i in range(self.nchan):
            vchan = self.channels[i]
            int2d = ax.imshow(self.data[i], cmap=cmap, extent=extent, origin='lower left', vmax=max_data)
            cbar = plt.colorbar(int2d)
            cbar.set_label(unit)
            if velocity2d is not None:
                vel_near=ax.contour(velocity2d['near'], levels=[vchan], colors=contour_color, linestyles='--', linewidths=1.3, extent = extent)
                vel_far=ax.contour(velocity2d['far'], levels=[vchan], colors=contour_color, linestyles=':', linewidths=1.3, extent = extent)
                coll_list = [vel_near, vel_far]
            text_chan = ax.text(0.7, 1.02, '%4.1f km/s'%vchan, color='black', transform=ax.transAxes)
            ax.legend(loc='upper left')
            plt.savefig(folder+'int2d_chan%04d'%i)
            #print ('Saved channel %d'%i)
            #plt.cla()
            clear_list = [cbar, int2d, text_chan]
            for obj in clear_list: obj.remove()
            for obj in coll_list: 
                for coll in obj.collections:
                    coll.remove()
        plt.close()
        os.chdir(folder)
        print ('Making movie...')
        os.system(gif_command)
        os.chdir(cwd)


class Tools:
    @staticmethod
    def _rotate_sky_plane(x, y, ang):
        xy = np.array([x,y])
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        rot = np.array([[cos_ang, -sin_ang],
                        [sin_ang, cos_ang]])
        return np.dot(rot, xy)

    @staticmethod
    def _rotate_sky_plane3d(x, y, z, ang, axis='z'):
        xyz = np.array([x,y,z])
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        if axis == 'x':
            rot = np.array([[1, 0, 0],
                            [0, cos_ang, -sin_ang],
                            [0, sin_ang, cos_ang]])
        if axis == 'y':
            rot = np.array([[cos_ang, 0, -sin_ang],
                            [0, 1, 0],
                            [sin_ang, 0, cos_ang]])
            
        if axis == 'z':
            rot = np.array([[cos_ang, -sin_ang , 0],
                            [sin_ang, cos_ang, 0], 
                            [0, 0, 1]])
        return np.dot(rot, xyz)

    @staticmethod
    def _project_on_skyplane(x, y, z, cos_incl, sin_incl):
        x_pro = x
        y_pro = y * cos_incl - z * sin_incl
        z_pro = y * sin_incl + z * cos_incl
        return x_pro, y_pro, z_pro

    @staticmethod
    def _compute_prop(grid, prop_funcs, prop_kwargs):
        n_funcs = len(prop_funcs)
        props = [{} for i in range(n_funcs)]
        for side in ['near', 'far']:
            x, y, z, R, phi = grid[side]
            coord = {'x': x, 'y': y, 'z': z, 'phi': phi, 'R': R}
            for i in range(n_funcs): props[i][side] = prop_funcs[i](coord, **prop_kwargs[i])
        return props


class Linewidth:
    @property
    def linewidth_func(self): 
        return self._linewidth_func
          
    @linewidth_func.setter 
    def linewidth_func(self, linewidth): 
        print('Setting linewidth function to', linewidth) 
        self._linewidth_func = linewidth

    @linewidth_func.deleter 
    def linewidth_func(self): 
        print('Deleting linewidth function') 
        del self._linewidth_func

    @staticmethod
    def linewidth_powerlaw(coord, L0=0.1, p=-0.4, q=0.3, R0=100*u.au, z0=100*u.au): #A=600.0, p=-0.4, q=0.3): #
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        z = coord['z']        
        A = L0*R0**-p*z0**-q
        return A*R**p*abs(z)**q


class Velocity:
    @property
    def velocity_func(self): 
        return self._velocity_func
          
    @velocity_func.setter 
    def velocity_func(self, vel): 
        print('Setting velocity function to', vel) 
        self._velocity_func = vel

    @velocity_func.deleter 
    def velocity_func(self): 
        print('Deleting velocity function') 
        del self._velocity_func

    @staticmethod
    def keplerian(coord, Mstar=u.Msun):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        return np.sqrt(G*Mstar/R) 
    
    @staticmethod
    def keplerian_vertical(coord, Mstar=u.Msun):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        if 'r' not in coord.keys(): r = np.hypot(R, coord['z'])
        else: r = coord['r']
        return np.sqrt(G*Mstar/r**3)*R


class Intensity:        
    @property
    def use_temperature(self):
        return self._use_temperature

    @use_temperature.setter 
    def use_temperature(self, use): 
        use = bool(use)
        print('Setting use_temperature var to', use)
        if use: self.line_profile = self.line_profile_temp
        else: self.line_profile = self.line_profile_v_sigma
        self._use_temperature = use

    @use_temperature.deleter 
    def use_temperature(self): 
        print('Deleting use_temperature var') 
        del self._use_temperature

    @property
    def line_profile(self): 
        return self._line_profile
          
    @line_profile.setter 
    def line_profile(self, profile): 
        print('Setting line profile function to', profile) 
        self._line_profile = profile

    @line_profile.deleter 
    def line_profile(self): 
        print('Deleting intensity function') 
        del self._line_profile
    
    @property
    def intensity_func(self): 
        return self._intensity_func
          
    @intensity_func.setter 
    def intensity_func(self, intensity): 
        print('Setting intensity function to', intensity) 
        self._intensity_func = intensity

    @intensity_func.deleter 
    def intensity_func(self): 
        print('Deleting intensity function') 
        del self._intensity_func

    @staticmethod
    def intensity_powerlaw(coord, I0=30.0, R0=100*u.au, p=-0.4, z0=100*u.au, q=0.3): #A=600.0, p=-0.4, q=0.3): #
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        z = coord['z']        
        A = I0*R0**-p*z0**-q
        return A*R**p*abs(z)**q
        
    @staticmethod
    def nuker(coord, I0=30.0, Rt=100*u.au, alpha=-0.5, gamma=0.1, beta=0.2):
        if 'R' not in coord.keys(): R = np.hypot(coord['x'], coord['y'])
        else: R = coord['R'] 
        A = I0*Rt**gamma
        return A*(R**-gamma) * (1+(R/Rt)**alpha)**((gamma-beta)/alpha)

    @staticmethod
    def line_profile_temp(v_chan, v, T, v_turb=0.0, mmol=2*u.amu):
        v_sigma = np.sqrt(2*kb*T/mmol + v_turb**2) * 1e-3 #in km/s
        return 1/(np.sqrt(np.pi)*v_sigma) * np.exp(-((v-v_chan)/v_sigma)**2)

    @staticmethod
    def line_profile_v_sigma(v_chan, v, v_sigma, mmol=2*u.amu):
        return 1/(np.sqrt(np.pi)*v_sigma) * np.exp(-((v-v_chan)/v_sigma)**2)

    def get_channel(self, velocity2d, intensity2d, temperature2d, v_chan, **kwargs):                    
        vel2d, temp2d, int2d = velocity2d, {}, {}
        if isinstance(temperature2d, numbers.Number): temp2d['near'] = temp2d['far'] = temperature2d
        else: temp2d = temperature2d
        if isinstance(intensity2d, numbers.Number): int2d['near'] = int2d['far'] = intensity2d
        else: int2d = intensity2d
    
        v_near = self.line_profile(v_chan, vel2d['near'], temp2d['near'], **kwargs)
        v_far = self.line_profile(v_chan, vel2d['far'], temp2d['far'], **kwargs)

        v_near_clean = np.where(np.isnan(v_near), -np.inf, v_near)
        v_far_clean = np.where(np.isnan(v_far), -np.inf, v_far)
        
        #int2d_near = int2d['near'] * v_near_clean / v_near_clean.max()
        int2d_near = np.where(np.isnan(int2d['near']), -np.inf, int2d['near'] * v_near_clean / v_near_clean.max())
        int2d_far = np.where(np.isnan(int2d['far']), -np.inf, int2d['far'] * v_far_clean / v_far_clean.max())
        
        vmap_full = np.array([v_near_clean, v_far_clean]).max(axis=0)
        int2d_full = np.array([int2d_near, int2d_far]).max(axis=0)

        return int2d_full

    def get_cube(self, vchan0, vchan1, velocity2d, intensity2d, temperature2d, nchan=30, folder='./movie_channels/', **kwargs):
        vel2d, temp2d, int2d = velocity2d, {}, {}
        line_profile = self.line_profile
        channels = np.linspace(vchan0, vchan1, num=nchan)
        cube = []

        if isinstance(temperature2d, numbers.Number): temp2d['near'] = temp2d['far'] = temperature2d
        else: temp2d = temperature2d
        if isinstance(intensity2d, numbers.Number): int2d['near'] = int2d['far'] = intensity2d
        else: int2d = intensity2d

        int2d_near_nan = np.isnan(int2d['near'])
        int2d_far_nan = np.isnan(int2d['far'])
        vel2d_near_nan = np.isnan(vel2d['near']) 
        vel2d_far_nan = np.isnan(vel2d['far'])

        for i, v_chan in enumerate(channels):    
            v_near = line_profile(v_chan, vel2d['near'], temp2d['near'], **kwargs)
            v_far = line_profile(v_chan, vel2d['far'], temp2d['far'], **kwargs)
            v_near_clean = np.where(vel2d_near_nan, -np.inf, v_near)
            v_far_clean = np.where(vel2d_far_nan, -np.inf, v_far)
            
            int2d_near = np.where(int2d_near_nan, -np.inf, int2d['near'] * v_near_clean / v_near_clean.max())
            int2d_far = np.where(int2d_far_nan, -np.inf, int2d['far'] * v_far_clean / v_far_clean.max())        
            #vmap_full = np.array([v_near_clean, v_far_clean]).max(axis=0)
            int2d_full = np.array([int2d_near, int2d_far]).max(axis=0)

            cube.append(int2d_full)

        return Cube(nchan, channels, np.array(cube))

    @staticmethod
    def make_channels_movie(vchan0, vchan1, velocity2d, intensity2d, temperature2d, nchans=30, folder='./movie_channels/', **kwargs):
        import matplotlib.pyplot as plt
        channels = np.linspace(vchan0, vchan1, num=nchans)
        int2d_cube = []
        for i, vchan in enumerate(channels):
            int2d = Intensity.get_channel(velocity2d, intensity2d, temperature2d, vchan, **kwargs)
            int2d_cube.append(int2d)
            extent = [-600, 600, -600, 600]
            plt.imshow(int2d, cmap='binary', extent=extent, origin='lower left', vmax=np.max(temperature2d['near']))
            plt.xlabel('au')
            plt.ylabel('au')
            plt.text(200, 500, '%.1f km/s'%vchan)
            cbar = plt.colorbar()
            cbar.set_label(r'Brightness Temperature [K]')
            plt.contour(velocity2d['near'], levels=[vchan], colors='red', linestyles='--', linewidths=1.3, extent = extent)
            plt.contour(velocity2d['far'], levels=[vchan], colors='red', linestyles=':', linewidths=1.3, extent = extent)
            plt.plot([None],[None], color='red', linestyle='--', linewidth=2, label='Near side') 
            plt.plot([None],[None], color='red', linestyle=':', linewidth=2, label='Far side') 
            plt.legend(loc='upper left')
            plt.savefig(folder+'int2d_chan%04d'%i)
            print ('Saved channel %d'%i)
            plt.close()

        os.chdir(folder)
        print ('Making channels movie...')
        os.system('convert -delay 10 *int2d* cube_channels.gif')
        return np.array(int2d_cube)

   
class General2d(Velocity, Intensity, Linewidth, Tools):
    def __init__(self, grid):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self._velocity_func = General2d.keplerian
        self._intensity_func = General2d.intensity_powerlaw
        self._linewidth_func = General2d.linewidth_powerlaw
        self._line_profile = General2d.line_profile_temp
        self._use_temperature = True

    def make_model(self, incl, z_func, PA=0.0, get_2d=True, z_far=None, int_kwargs={}, vel_kwargs={}, lw_kwargs=None):
        from scipy.interpolate import griddata
        #*************************************
        #MAKE TRUE GRID FOR NEAR AND FAR SIDES
        cos_incl, sin_incl = np.cos(incl), np.sin(incl)

        x_init, y_init = self.grid.XYZ[:2]
        x_true, y_true = x_init, y_init
        phi_true = np.arctan2(y_true, x_true)        
        R_true = np.hypot(x_init, y_init)
        z_true = z_func(R_true)

        if z_far is not None: z_true_far = z_far(R_true) 
        else: z_true_far = -z_true
            
        grid_true = {'near': [x_true, y_true, z_true, R_true, phi_true], 
                     'far': [x_true, y_true, z_true_far, R_true, phi_true]}

        #*******************************
        #COMPUTE PROPERTIES ON TRUE GRID
        avai_kwargs = [vel_kwargs, int_kwargs, lw_kwargs]
        avai_funcs = [self.velocity_func, self.intensity_func, self.linewidth_func]
        true_kwargs = [isinstance(kwarg, dict) for kwarg in avai_kwargs]
        prop_kwargs = [kwarg for i, kwarg in enumerate(avai_kwargs) if true_kwargs[i]]
        prop_funcs = [func for i, func in enumerate(avai_funcs) if true_kwargs[i]]
        props = self._compute_prop(grid_true, prop_funcs, prop_kwargs)
        #Positive vel is positive along z, i.e. pointing to the observer, for that reason imposed a (-) factor to convert to the standard convention: (+) receding  
        if true_kwargs[0]:
            ang_fac = sin_incl * np.cos(phi_true) #np.sin(incl+0*np.pi/2)*np.cos(phi_true) 
            props[0]['near'] *= ang_fac 
            props[0]['far'] *= ang_fac

        #***********************************
        #PROJECT PROPERTIES ON THE SKY PLANE
        #grid_axes = np.meshgrid(*self.grid.XYZgrid[:2]) #Avoiding this to keep backward compatibility with python2
        grid_axes = np.meshgrid(self.grid.XYZgrid[0], self.grid.XYZgrid[1]) 
        
        for side in ['near', 'far']:
            xt, yt, zt = grid_true[side][:3]
            x_pro, y_pro, z_pro = self._project_on_skyplane(xt, yt, zt, cos_incl, sin_incl)
            if PA: x_pro, y_pro = self._rotate_sky_plane(x_pro, y_pro, PA)             
            for prop in props: 
                prop[side] = griddata((x_pro, y_pro), prop[side], (grid_axes[0], grid_axes[1]))
        #*************************************

        return props
    
class Rosenfeld2d(Velocity, Intensity, Linewidth, Tools):
    """
    Host class for the Rosenfeld+2013 model which describes the velocity field of a flared disc in 2D. 
    This model assumes a (Keplerian) double cone to account for the near and far sides of the disc 
    and solves analytical equations to find the line-of-sight velocity v_obs projected on the sky-plane from both sides. 
    
    Parameters
    ----------
    grid : array_like, shape (nrows, ncols)
       (x', y') map of the sky-plane onto which the disc velocity field will be projected.

    Attributes
    ----------
    velocity_func : function(coord, **kwargs) 
       Velocity function describing the kinematics of the disc. The argument coord is a dictionary
       of coordinates (e.g. 'x', 'y', 'z', 'r', 'R', 'theta', 'phi') where the function will be evaluated. 
       Additional arguments are optional and depend upon the function definition, e.g. Mstar=1.0*u.Msun
    """

    def __init__(self, grid):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self._velocity_func = Rosenfeld2d.keplerian
        self._intensity_func = Rosenfeld2d.intensity_powerlaw
        self._linewidth_func = Rosenfeld2d.linewidth_powerlaw
        self._line_profile = General2d.line_profile_temp
        self._use_temperature = True

    def _get_t(self, A, B, C):
        t = []
        for i in range(self.grid.NPoints):
            p = [A, B[i], C[i]]
            t.append(np.sort(np.roots(p)))
        return np.array(t)

    def make_model(self, incl, psi, PA=0.0, get_2d=True, int_kwargs={}, vel_kwargs={}, lw_kwargs=None):
        """
        Executes the Rosenfeld+2013 model.
        The sum of incl+psi must be < 90, otherwise the quadratic equation will have imaginary roots as some portions of the cone (which has finite extent)
        do not intersect with the sky plane.  

        Parameters
        ----------
        incl : scalar
           Inclination of the disc midplane with respect to the x'y' plane; pi/2 radians is edge-on.
    
        psi : scalar
           Opening angle of the cone describing the velocity field of the gas emitting layer in the disc; 
           0 radians returns the projected velocity field of the disc midplane (i.e no conic emission). 

        PA : scalar, optional
           Position angle in radians. Measured from North (+y) to East (-x).

        get_2d : bool, optional
           If True regrids the resulting velocity field into a 2D map and stores it in the attribute 'velocity2d'. 

        Attributes
        ----------
        velocity : array_like, size (n,)
           Velocity field computed using the Rosenfeld+2013 model.

        velocity2d : array_like, size (nx, ny)
           If set get_2d=True: Velocity field computed using the Rosenfeld+2013 model, reshaped to 2D to facilitate plotting.
        """
        if PA: x_plane, y_plane = Rosenfeld2d._rotate_sky_plane(self.grid.XYZ[0], self.grid.XYZ[1], -PA)
        else: x_plane, y_plane = self.grid.XYZ[:2]

        cos_incl = np.cos(incl)
        sin_incl = np.sin(incl)
        y_plane_cos_incl = y_plane/cos_incl

        #**********************
        #ROSENFELD COEFFICIENTS
        fac = -2*np.sin(psi)**2
        A = np.cos(2*incl) + np.cos(2*psi)
        B = fac * 2*(sin_incl/cos_incl) * y_plane
        C = fac * (x_plane**2 + (y_plane_cos_incl)**2)
        t = self._get_t(A,B,C).T

        #****************************
        #ROSENFELD CONVERSION X<-->X'
        x_true_near = x_plane
        y_true_near = y_plane_cos_incl + t[1]*sin_incl
            
        x_true_far = x_plane
        y_true_far = y_plane_cos_incl + t[0]*sin_incl
        
        #np.hypot 2x faster than np.linalg.norm([x,y], axis=0)
        R_true_near = np.hypot(x_true_near, y_true_near) 
        R_true_far = np.hypot(x_true_far, y_true_far)

        z_true_near = t[1] * cos_incl
        z_true_far = t[0] * cos_incl 

        phi_true_near = np.arctan2(y_true_near, x_true_near)        
        phi_true_far = np.arctan2(y_true_far, x_true_far)        

        #****************************
            
        grid_true =  {'near': [x_true_near, y_true_near, z_true_near, R_true_near, phi_true_near], 
                      'far': [x_true_far, y_true_far, z_true_far, R_true_far, phi_true_far]}

        #*******************************
        #COMPUTE PROPERTIES ON TRUE GRID
        avai_kwargs = [vel_kwargs, int_kwargs, lw_kwargs]
        avai_funcs = [self.velocity_func, self.intensity_func, self.linewidth_func]
        true_kwargs = [isinstance(kwarg, dict) for kwarg in avai_kwargs]
        prop_kwargs = [kwarg for i, kwarg in enumerate(avai_kwargs) if true_kwargs[i]]
        prop_funcs = [func for i, func in enumerate(avai_funcs) if true_kwargs[i]]
        props = self._compute_prop(grid_true, prop_funcs, prop_kwargs)
        #Positive vel is positive along z, i.e. pointing to the observer, for that reason imposed a (-) factor to convert to the standard convention: (+) receding  
        if true_kwargs[0]:
            ang_fac_near = -sin_incl * np.cos(phi_true_near)
            ang_fac_far = -sin_incl * np.cos(phi_true_far)
            props[0]['near'] *= ang_fac_near 
            props[0]['far'] *= ang_fac_far
                
        #*************************************

        return [{side: prop[side].reshape(self.grid.Nodes[:2]) for side in ['near', 'far']} for prop in props]

