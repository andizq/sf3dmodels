"""
2D Disc models
==============
Classes: Rosenfeld2d, General2d, Velocity, Intensity, Cube, Tools
"""
#TODO in show(): Perhaps use text labels on line profiles to distinguish prof from more than 2 cubes  
from ..utils.constants import G, kb
from ..utils import units as sfu
from astropy.convolution import Gaussian2DKernel, convolve
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import ticker
import numpy as np
import matplotlib
import itertools
import warnings
import numbers
import emcee
import copy
import time
import os

from multiprocessing import Pool
os.environ["OMP_NUM_THREADS"] = "1"

#warnings.filterwarnings("error")
__all__ = ['Cube', 'Tools', 'Intensity', 'Velocity', 'General2d', 'Rosenfeld2d']
path_file = os.path.dirname(os.path.realpath(__file__))+'/'

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
matplotlib.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of axes title
matplotlib.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of x and y labels
matplotlib.rc('xtick', labelsize=MEDIUM_SIZE-2)    # fontsize of y tick labels
matplotlib.rc('ytick', labelsize=MEDIUM_SIZE-2)    # fontsize of x tick labels
matplotlib.rc('legend', fontsize=SMALL_SIZE-1)    # legend fontsize
matplotlib.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of figure title

params = {'xtick.major.size': 6.5,
          'ytick.major.size': 6.5
          }

matplotlib.rcParams.update(params)
sigma2fwhm = np.sqrt(8*np.log(2))

hypot_func = lambda x,y: np.sqrt(x**2 + y**2) #Slightly faster than np.hypot<np.linalg.norm<scipydistance. Checked precision up to au**2 orders and it all seems ok.

class InputError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        
    def __str__(self):
        return '%s --> %s'%(self.expression, self.message)

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
  
    @staticmethod
    def get_beam_from(beam, grid, distance=None, frac_pixels=1.0):
        """
        beam must be str pointing to fits file to extract beam from header or radio_beam Beam object.
        If Beam object is provided, distance must be set.
        #frac_pixels: number of averaged pixels on the data to reduce size and time
        """
        from radio_beam import Beam
        from astropy.io import fits
        from astropy import units as u
        if isinstance(beam, str):
            header = fits.getheader(beam)
            beam = Beam.from_fits_header(header)
            pix_scale = header['CDELT2'] * u.Unit(header['CUNIT2'])
        elif isinstance(beam, Beam):
            if distance is None: InputError(distance, 'Wrong input distance. Please provide a value for the distance (in pc) to transform grid pix to arcsec')
            pix_arcsec = np.arctan(grid.step[0] / (distance*sfu.pc)) #dist*ang=projdist
            pix_scale = (pix_arcsec*u.radian).to(u.arcsec)  
        else: InputError(beam, 'beam object must either be str or Beam instance')

        x_stddev = ((beam.major/pix_scale) / sigma2fwhm).value / frac_pixels 
        y_stddev = ((beam.minor/pix_scale) / sigma2fwhm).value / frac_pixels
        print (x_stddev, beam.major)
        angle = (90*u.deg+beam.pa).to(u.radian).value
        gauss_kern = Gaussian2DKernel(x_stddev, y_stddev, angle) 

        #gauss_kern = beam.as_kernel(pix_scale) #as_kernel() is slowing down the run when used in astropy.convolve
        return beam, gauss_kern
    
    @staticmethod
    def _get_tb(I, nu, beam):
        """
        nu in GHz
        Intensity in mJy/beam
        beam object from radio_beam
        """
        from astropy import units as u
        return (1222.0*I/(nu**2*(beam.minor/1.0).to(u.arcsecond)*(beam.major/1.0).to(u.arcsecond))).value


class Residuals:
    pass


class PlotTools:
    @staticmethod
    def mod_nticks_cbars(cbars, nbins=5):
        for cb in cbars:
            cb.locator = ticker.MaxNLocator(nbins=nbins)
            cb.update_ticks()

    @staticmethod
    def mod_major_ticks(ax, axis='both', nbins=6):
        ax.locator_params(axis=axis, nbins=nbins)

    @staticmethod
    def mod_minor_ticks(ax):
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2)) #1 minor tick per major interval
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        
    @staticmethod
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap

    @staticmethod
    def get_cmap_from_color(color, lev=3):
        cmap = matplotlib.colors.to_rgba(color)
        newcolors = np.tile(cmap, lev).reshape(lev,4) #Repeats the colour lev times
        newcolors[:,-1] = np.linspace(0.25, 0.95, lev) #Modifies alpha only
        new_cmap = ListedColormap(newcolors)
        return new_cmap

    @classmethod
    def make_up_ax(cls, ax, xlims=(None, None), ylims=(None, None), 
                   mod_minor=True, mod_major=True, **kwargs_tick_params):
        kwargs_t = dict(labeltop=True, labelbottom=False, top=True, right=True, which='both', direction='in')
        kwargs_t.update(kwargs_tick_params)
        if mod_major: cls.mod_major_ticks(ax)
        if mod_minor: cls.mod_minor_ticks(ax)
        ax.set_xlim(*xlims)
        ax.set_ylim(*ylims)
        ax.tick_params(**kwargs_t)
    
    @staticmethod
    def append_stddev_panel(ax, prop):
        gauss = lambda x, A, mu, sigma: A*np.exp(-(x-mu)**2/(2.*sigma**2))
        ax[-1].set_xlim(-0.2, 1.2)
        ax1_ylims = ax[-2].get_ylim()
        for axi in ax[:-1]: axi.tick_params(which='both', right=False, labelright=False)
        ax[-1].tick_params(which='both', top=False, bottom=False, labelbottom=False, 
                           left=False, labelleft=False, right=True, labelright=True)
        ax[-1].yaxis.set_label_position('right')
        ax[-1].spines['left'].set_color('0.6')
        ax[-1].spines['left'].set_linewidth(3.5)
    
        prop_mean = np.mean(prop)
        prop_std = np.std(prop)
        prop_x = np.linspace(prop_mean-4*prop_std, prop_mean+4*prop_std, 100)
        prop_pars =  [1.0, prop_mean, prop_std]
        prop_y = gauss(prop_x, *prop_pars)
        ax[-1].plot(prop_y, prop_x, color='limegreen', lw=3.5)
        #ax[-1].plot([-0.2, 1.0], [prop_mean]*2, color='0.6', lw=2.5)
        #for axi in ax[:-1]: axi.axhline(prop_mean, color='0.6', lw=2.5)
    
        for i in range(0,4): 
            prop_stdi = prop_mean+i*prop_std
            gauss_prop_stdi = gauss(prop_stdi, *prop_pars)
            ax[-1].plot([-0.2, gauss_prop_stdi], [prop_stdi]*2, color='0.6', ls=':', lw=2.)
            for axi in ax[:-1]: axi.axhline(prop_stdi, color='0.6', ls=':', lw=2.)
            if prop_stdi < ax1_ylims[-1] and i>0:
                ax[-1].text(gauss_prop_stdi+0.2, prop_stdi, r'%d$\sigma$'%i, 
                            fontsize=14, ha='center', va='center', rotation=-90)

        for axi in ax: axi.set_ylim(*ax1_ylims)

       
class Canvas3d:
    pass


class Contours(PlotTools):
    @staticmethod
    def emission_surface(ax, R, phi, R_lev=None, phi_lev=None, extent=None, proj_offset=None, X=None, Y=None, kwargs_R={}, kwargs_phi={}):
        kwargs_phif = dict(linestyles=':', linewidths=1.0, colors='k')
        kwargs_Rf = dict(linewidths=1.4, colors='k')
        kwargs_phif.update(kwargs_phi)        
        kwargs_Rf.update(kwargs_R)

        near_nonan = ~np.isnan(R['near'])

        Rmax = np.max(R['near'][near_nonan])
        if extent is None:
            extent = np.array([-Rmax, Rmax, -Rmax, Rmax])/sfu.au
        kwargs_phif.update({'extent': extent})
        kwargs_Rf.update({'extent': extent})

        if R_lev is None: R_lev = np.linspace(0.06, 0.97, 4)*Rmax
        else: R_lev = np.sort(R_lev)
        if phi_lev is None: phi_lev = np.linspace(-np.pi*0.95, np.pi, 11, endpoint=False)

        #Splitting phi into pos and neg to try and avoid ugly contours close to -pi and pi
        phi_lev_neg = phi_lev[phi_lev<0] 
        phi_lev_pos = phi_lev[phi_lev>0]
        phi_neg_near = np.where((phi['near']<0) & (R['near']>R_lev[0]) & (R['near']<R_lev[-1]), phi['near'], np.nan)
        phi_pos_near = np.where((phi['near']>0) & (R['near']>R_lev[0]) & (R['near']<R_lev[-1]), phi['near'], np.nan)
        phi_neg_far = np.where((phi['far']<0) & (R['far']>R_lev[0]) & (R['far']<R_lev[-1]), phi['far'], np.nan)
        phi_pos_far = np.where((phi['far']>0) & (R['far']>R_lev[0]) & (R['far']<R_lev[-1]), phi['far'], np.nan)

        if proj_offset is not None: #For 3d projections
            ax.contour(X, Y, R['near'], offset=proj_offset, levels=R_lev, **kwargs_Rf)
            ax.contour(X, Y, np.where(near_nonan, np.nan, R['far']), offset=proj_offset, levels=R_lev, **kwargs_Rf)
            ax.contour(X, Y, phi_pos_near, offset=proj_offset, levels=phi_lev_pos, **kwargs_phif)
            ax.contour(X, Y, phi_neg_near, offset=proj_offset, levels=phi_lev_neg, **kwargs_phif)
            ax.contour(X, Y, np.where(near_nonan, np.nan, phi_pos_far), offset=proj_offset, levels=phi_lev_pos, **kwargs_phif)
            ax.contour(X, Y, np.where(near_nonan, np.nan, phi_neg_far), offset=proj_offset, levels=phi_lev_neg, **kwargs_phif)
            
        else:
            ax.contour(R['near'], levels=R_lev, **kwargs_Rf)
            ax.contour(np.where(near_nonan, np.nan, R['far']), levels=R_lev, **kwargs_Rf)
            ax.contour(phi_pos_near, levels=phi_lev_pos, **kwargs_phif)
            ax.contour(phi_neg_near, levels=phi_lev_neg, **kwargs_phif)
            ax.contour(np.where(near_nonan, np.nan, phi_pos_far), levels=phi_lev_pos, **kwargs_phif)
            ax.contour(np.where(near_nonan, np.nan, phi_neg_far), levels=phi_lev_neg, **kwargs_phif)

    #The following method can be optimised if the contour finding process is separated from the plotting
    # by returning coords_list and inds_cont first, which will allow the user use the same set of contours to plot different props.
    @staticmethod
    def prop_along_coords(ax, prop, coords, coord_ref, coord_levels, 
                          ax2=None, X=None, Y=None, 
                          PA=0,
                          acc_threshold=0.05, 
                          color_bounds=[np.pi/5, np.pi/2],
                          colors=['k', 'dodgerblue', (0,1,0), (1,0,0)],
                          lws=[2, 0.5, 0.2, 0.2], lw_ax2_factor=1,
                          subtract_quadrants=False):
        """
        Compute radial/azimuthal contours according to the model disc geometry 
        to get and plot information from the input 2D property ``prop``.    

        Parameters
        ----------
        ax : `matplotlib.axes` instance, optional
           ax instance to make the plot. 

        prop : array_like, shape (nx, ny)
           Input 2D field to extract information along the computed contours.
        
        coords : list, shape (2,)
           coords[0] [array_like, shape (nx, ny)], is the coordinate 2D map onto which contours will be computed using the input ``coord_levels``;
           coords[1] [array_like, shape (nx, ny)], is the coordinate 2D map against which the ``prop`` values are plotted. The output plot is prop vs coords[1]       
           
        coord_ref : scalar
           Reference coordinate (referred to ``coords[0]``) to highlight among the other contours.
           
        coord_levels : array_like, shape (nlevels,)
           Contour levels to be extracted from ``coords[0]``.
        
        ax2 : `matplotlib.axes` instance (or list of instances), optional
           Additional ax(s) instance(s) to plot the location of contours in the disc. 
           If provided, ``X`` and ``Y`` must also be passed.
           
        X : array_like, shape (nx, ny), optional
           Meshgrid of the model x coordinate (see `numpy.meshgrid`). Required if ax2 instance(s) is provided.

        Y : array_like, shape (nx, ny), optional
           Meshgrid of the model y coordinate (see `numpy.meshgrid`). Required if ax2 instance(s) is provided.
        
        PA : scalar, optional
           Reference position angle.
           
        acc_threshold : float, optional 
           Threshold to accept contours at constant coords[0].

        color_bounds : array_like, shape (nbounds,), optional
           Colour bounds with respect to the reference contour coord_ref.
           
        colors : array_like, shape (nbounds+2,), optional
           Contour colors. (i=0) is reserved for the reference contour coord_ref, 
           (i>0) for contour colors according to the bounds in color_bounds. 
           
        lws : array_like, shape (nbounds+2), optional
           Contour linewidths. Similarly, (i=0) is reserved for coord_ref and 
           (i>0) for subsequent bounds.

        subtract_quadrants : bool, optional
           If True, subtract residuals by folding along the projected minor axis of the disc.
        """
        from skimage import measure 

        coord_list, lev_list, resid_list, color_list = [], [], [], []
        if np.sum(coord_levels==coord_ref)==0: coord_levels = np.append(coord_levels, coord_ref)
        for lev in coord_levels:
            contour = measure.find_contours(coords[0], lev) #, fully_connected='high', positive_orientation='high')
            if len(contour)==0:
                print ('no contours found for phi =', lev)
                continue
            inds_cont = np.round(contour[-1]).astype(np.int)
            inds_cont = [tuple(f) for f in inds_cont]
            first_cont = np.array([coords[0][i] for i in inds_cont])
            second_cont = np.array([coords[1][i] for i in inds_cont])
            prop_cont = np.array([prop[i] for i in inds_cont])
            corr_inds = np.abs(first_cont-lev) < acc_threshold
            if lev == coord_ref: zorder=10
            else: zorder=np.random.randint(0,10)

            lw = lws[-1]
            color = colors[-1]
            for i,bound in enumerate(color_bounds):
                if lev == coord_ref: 
                    lw = lws[0]
                    color = colors[0]
                    zorder = 10
                    break
                if np.abs(coord_ref - lev) < bound: 
                    lw = lws[i+1]
                    color = colors[i+1]
                    break

            if subtract_quadrants:
                if lev < color_bounds[0]: continue
                ref_pos = PA+90
                ref_neg = PA-90
                angles = second_cont[corr_inds]
                prop_ = prop_cont[corr_inds]
                angles_pos = angles[angles>=0]
                angles_neg = angles[angles<0]
                relative_diff_pos = ref_pos - angles_pos
                relative_diff_neg = ref_neg - angles_neg
                angle_diff_pos, prop_diff_pos = [], []
                angle_diff_neg, prop_diff_neg = [], []

                for i,diff in enumerate(relative_diff_pos):
                    #Finding where the difference matches that of the current analysis angle
                    #The -1 flips the sign so that the number on the other side of the symmetry axis is found                
                    ind = np.argmin(np.abs(-1*relative_diff_pos - diff))  
                    mirror_ind = angles==angles_pos[ind]
                    current_ind = angles==angles_pos[i]
                    prop_diff = prop_[current_ind][0] - prop_[mirror_ind][0]
                    angle_diff_pos.append(angles_pos[i])
                    prop_diff_pos.append(prop_diff)
                angle_diff_pos = np.asarray(angle_diff_pos)
                prop_diff_pos = np.asarray(prop_diff_pos)

                if len(angle_diff_pos)>1:
                    ind_sort_pos = np.argsort(angle_diff_pos)
                    plot_ang_diff_pos = angle_diff_pos[ind_sort_pos]
                    plot_prop_diff_pos = prop_diff_pos[ind_sort_pos]
                    ax.plot(plot_ang_diff_pos, plot_prop_diff_pos, color=color, lw=lw, zorder=zorder)
                    coord_list.append(plot_ang_diff_pos)
                    resid_list.append(plot_prop_diff_pos)
                    color_list.append(color)
                    lev_list.append(lev)
                else: 
                    plot_ang_diff_pos = []
                    plot_prop_diff_pos = []

                for i,diff in enumerate(relative_diff_neg):
                    ind = np.argmin(np.abs(-1*relative_diff_neg - diff))
                    mirror_ind = angles==angles_neg[ind]
                    current_ind = angles==angles_neg[i]
                    prop_diff = prop_[current_ind][0] - prop_[mirror_ind][0]
                    angle_diff_neg.append(angles_neg[i])
                    prop_diff_neg.append(prop_diff)
                angle_diff_neg = np.asarray(angle_diff_neg)
                prop_diff_neg = np.asarray(prop_diff_neg)

                if len(angle_diff_neg)>1:
                    ind_sort_neg = np.argsort(angle_diff_neg)    
                    plot_ang_diff_neg = angle_diff_neg[ind_sort_neg]
                    plot_prop_diff_neg = prop_diff_neg[ind_sort_neg]
                    ax.plot(plot_ang_diff_neg, plot_prop_diff_neg, color=color, lw=lw, zorder=zorder)
                    coord_list.append(plot_ang_diff_neg)
                    resid_list.append(plot_prop_diff_neg)
                    color_list.append(color)
                    lev_list.append(lev)
                else: 
                    plot_ang_diff_neg = []
                    plot_prop_diff_neg = []

                """
                if len(angle_diff_pos)>1 or len(angle_diff_neg)>1:
                    coord_list.append(np.append(plot_ang_diff_pos, plot_ang_diff_neg))
                    resid_list.append(np.append(plot_prop_diff_pos, plot_prop_diff_neg))
                    color_list.append(color)
                    lev_list.append(lev)
                """
            else: ax.plot(second_cont[corr_inds], prop_cont[corr_inds], color=color, lw=lw, zorder=zorder)

            if ax2 is not None:
                x_cont = np.array([X[i] for i in inds_cont])
                y_cont = np.array([Y[i] for i in inds_cont])
            if isinstance(ax2, matplotlib.axes._subplots.Axes): 
                ax2.plot(x_cont[corr_inds], y_cont[corr_inds], color=color, lw=lw*lw_ax2_factor)
            elif isinstance(ax2, list):
                for axi in ax2: 
                    if isinstance(axi, matplotlib.axes._subplots.Axes): axi.plot(x_cont[corr_inds], y_cont[corr_inds], color=color, lw=lw*lw_ax2_factor)

        return [np.asarray(tmp) for tmp in [coord_list, resid_list, color_list, lev_list]]


class Cube(object):
    def __init__(self, nchan, channels, data, beam=False, beam_kernel=False, tb={'nu': False, 'beam': False}):
        self.nchan = nchan
        self.channels = channels
        self.data = data
        self.point = self.cursor
        self._interactive = self.cursor
        self._interactive_path = self.curve
        if beam: self.beam = beam
        if beam_kernel: self.beam_kernel = beam_kernel
        if isinstance(tb, dict):
            if tb['nu'] and tb['beam']: self.data = Tools._get_tb(self.data, tb['nu'], tb['beam'])

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

    @property
    def interactive_path(self): 
        return self._interactive_path
          
    @interactive_path.setter 
    def interactive_path(self, func): 
        print('Setting interactive_path function to', func) 
        self._interactive_path = func

    @interactive_path.deleter 
    def interactive_path(self): 
        print('Deleting interactive_path function') 
        del self._interactive_path

    def ellipse(self):
        pass
    
    def _plot_spectrum_region(self, x0, x1, y0, y1, ax, extent=None, compare_cubes=[], stat_func=np.mean, **kwargs):
        kwargs_spec = dict(where='mid', linewidth=2.5, label=r'x0:%d,x1:%d'%(x0,x1))
        kwargs_spec.update(kwargs)
        v0, v1 = self.channels[0], self.channels[-1]
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
        ncubes = len(compare_cubes)
        if ncubes > 0: 
            slice_comp = [compare_cubes[i].data[:,i0:i1,j0:j1] for i in range(ncubes)]
            cubes_spec = [np.array([stat_func(chan) for chan in slice_comp[i]]) for i in range(ncubes)]

        if np.logical_or(np.isinf(spectrum), np.isnan(spectrum)).all(): return False
        else:
            plot_spec = ax.step(self.channels, spectrum, **kwargs_spec)
            if ncubes > 0:
                alpha = 0.2
                dalpha = -alpha/ncubes
                for i in range(ncubes):
                    ax.fill_between(self.channels, cubes_spec[i], color=plot_spec[0].get_color(), step='mid', alpha=alpha)
                    alpha+=dalpha
            else: ax.fill_between(self.channels, spectrum, color=plot_spec[0].get_color(), step='mid', alpha=0.2)
            return plot_spec
   
    def box(self, fig, ax, extent=None, compare_cubes=[], stat_func=np.mean, **kwargs):
        from matplotlib.widgets import RectangleSelector
        
        def onselect(eclick, erelease):
            if eclick.inaxes is ax[0]:
                plot_spec = self._plot_spectrum_region(eclick.xdata, erelease.xdata, eclick.ydata, erelease.ydata, 
                                                       ax[1], extent=extent, compare_cubes=compare_cubes, 
                                                       stat_func=stat_func, **kwargs) 
                                                       
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
        return toggle_selector.RS

    def _plot_spectrum_cursor(self, x, y, ax, extent=None, compare_cubes=[], **kwargs):
        kwargs_spec = dict(where='mid', linewidth=2.5, label=r'%d,%d'%(x,y))
        kwargs_spec.update(kwargs)
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
            #plot_fill = ax.fill_between(self.channels, spectrum, alpha=0.1)
            plot_spec = ax.step(self.channels, spectrum, **kwargs_spec)
            ncubes = len(compare_cubes)
            if ncubes > 0:
                alpha = 0.2
                dalpha = -alpha/ncubes
                for cube in compare_cubes: 
                    ax.fill_between(self.channels, cube.data[:,i,j], color=plot_spec[0].get_color(), step='mid', alpha=alpha)
                    alpha+=dalpha
            else: ax.fill_between(self.channels, spectrum, color=plot_spec[0].get_color(), step='mid', alpha=0.2)
            return plot_spec
        
    #def point(self, *args, **kwargs):
     #   return self.cursor(*args, **kwargs)

    def cursor(self, fig, ax, extent=None, compare_cubes=[], **kwargs):
        def onclick(event):
            if event.button==3: 
                print ('Right click. Disconnecting click event...')
                fig.canvas.mpl_disconnect(cid)
            elif event.inaxes is ax[0]:
                plot_spec = self._plot_spectrum_cursor(event.xdata, event.ydata, ax[1], extent=extent, 
                                                       compare_cubes=compare_cubes, **kwargs) 
                if plot_spec:
                    print('%s click: button=%d, xdata=%f, ydata=%f' %
                          ('double' if event.dblclick else 'single', event.button,
                           event.xdata, event.ydata))
                    ax[0].scatter(event.xdata, event.ydata, color=plot_spec[0].get_color())
                    ax[1].legend()
                    fig.canvas.draw()
                    fig.canvas.flush_events()

        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        return cid

    def _plot_beam(self, ax):
        x_fwhm = self.beam_kernel.model.x_fwhm
        y_fwhm = self.beam_kernel.model.y_fwhm
        ny_pix, nx_pix = np.shape(self.data[0])
        ellipse = patches.Ellipse(xy = (0.05,0.05), angle = 90+self.beam.pa.value,
                                  width=x_fwhm/nx_pix, height=y_fwhm/ny_pix, lw=1, fill=True, 
                                  fc='gray', ec='k', transform=ax.transAxes)
        ax.add_artist(ellipse)
        
    def surface(self, ax, *args, **kwargs): return Contours.emission_surface(ax, *args, **kwargs)

    def show(self, extent=None, chan_init=20, compare_cubes=[], cursor_grid=True, cmap='gnuplot2_r',
             int_unit=r'Intensity [mJy beam$^{-1}$]', pos_unit='au', vel_unit=r'km s$^{-1}$',
             show_beam=False, surface={'args': (), 'kwargs': {}}, **kwargs):
        from matplotlib.widgets import Slider, Cursor, Button
        v0, v1 = self.channels[0], self.channels[-1]
        dv = v1-v0
        fig, ax = plt.subplots(ncols=2, figsize=(12,5))
        plt.subplots_adjust(wspace=0.25)

        y0, y1 = ax[1].get_position().y0, ax[1].get_position().y1
        axcbar = plt.axes([0.47, y0, 0.03, y1-y0])
        max_data = np.max([self.data]+[comp.data for comp in compare_cubes])
        ax[0].set_xlabel(pos_unit)
        ax[0].set_ylabel(pos_unit)
        ax[1].set_xlabel('l.o.s velocity [%s]'%vel_unit)
        ax[1].tick_params(direction='in', right=True, labelright=False, labelleft=False)
        axcbar.tick_params(direction='out')
        ax[1].set_ylabel(int_unit, labelpad=15)
        ax[1].yaxis.set_label_position('right')
        ax[1].set_xlim(v0-0.1, v1+0.1)
        vmin, vmax = -1*max_data/100, 0.98*max_data#0.8*max_data#
        ax[1].set_ylim(vmin, vmax)
        ax[1].grid(lw=1.5, ls=':')
        cmap = plt.get_cmap(cmap)
        cmap.set_bad(color=(0.9,0.9,0.9))

        if show_beam and self.beam_kernel: self._plot_beam(ax[0])

        img = ax[0].imshow(self.data[chan_init], cmap=cmap, extent=extent, origin='lower left', vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(img, cax=axcbar)
        img.cmap.set_under('w')
        current_chan = ax[1].axvline(self.channels[chan_init], color='black', lw=2, ls='--')
        text_chan = ax[1].text((self.channels[chan_init]-v0)/dv, 1.02, #Converting xdata coords to Axes coords 
                               '%4.1f %s'%(self.channels[chan_init], vel_unit), ha='center', 
                               color='black', transform=ax[1].transAxes)

        if cursor_grid: cg = Cursor(ax[0], useblit=True, color='lime', linewidth=1.5)

        def get_interactive(func):
            return func(fig, ax, extent=extent, compare_cubes=compare_cubes, **kwargs)
        
        interactive_obj = [get_interactive(self.interactive)]
        #***************
        #SLIDERS
        #***************
        def update_chan(val):
            chan = int(val)
            vchan = self.channels[chan]
            img.set_data(self.data[chan])
            current_chan.set_xdata(vchan)
            text_chan.set_x((vchan-v0)/dv)
            text_chan.set_text('%4.1f %s'%(vchan, vel_unit))
            fig.canvas.draw_idle()

        def update_cubes(val):
            i = int(slider_cubes.val)
            chan = int(slider_chan.val)
            vchan = self.channels[chan]
            if i==0: img.set_data(self.data[chan])
            else: img.set_data(compare_cubes[i-1].data[chan])
            current_chan.set_xdata(vchan)
            text_chan.set_x((vchan-v0)/dv)
            text_chan.set_text('%4.1f km/s'%vchan)
            fig.canvas.draw_idle()

        ncubes = len(compare_cubes)
        if ncubes>0:
            axcubes = plt.axes([0.2, 0.90, 0.24, 0.025], facecolor='0.7')
            axchan = plt.axes([0.2, 0.95, 0.24, 0.025], facecolor='0.7')
            slider_cubes = Slider(axcubes, 'Cube id', 0, ncubes, 
                                  valstep=1, valinit=0, valfmt='%1d', color='dodgerblue')                                  
            slider_chan = Slider(axchan, 'Channel', 0, self.nchan-1, 
                                 valstep=1, valinit=chan_init, valfmt='%2d', color='dodgerblue')        
            slider_cubes.on_changed(update_cubes)
            slider_chan.on_changed(update_cubes)
        else: 
            axchan = plt.axes([0.2, 0.9, 0.24, 0.05], facecolor='0.7')
            slider_chan = Slider(axchan, 'Channel', 0, self.nchan-1, 
                                 valstep=1, valinit=chan_init, valfmt='%2d', color='dodgerblue')        
            slider_chan.on_changed(update_chan)
    
        #*************
        #BUTTONS
        #*************
        def go2cursor(event):
            if self.interactive == self.cursor or self.interactive == self.point: return 0
            interactive_obj[0].set_active(False)
            self.interactive = self.cursor
            interactive_obj[0] = get_interactive(self.interactive)
        def go2box(event):
            if self.interactive == self.box: return 0
            fig.canvas.mpl_disconnect(interactive_obj[0])
            self.interactive = self.box
            interactive_obj[0] = get_interactive(self.interactive)
        def go2trash(event):
            print ('Cleaning interactive figure...')
            plt.close()
            chan = int(slider_chan.val)
            self.show(extent=extent, chan_init=chan, compare_cubes=compare_cubes, 
                      cursor_grid=cursor_grid, int_unit=int_unit, pos_unit=pos_unit, 
                      vel_unit=vel_unit, surface=surface, **kwargs)
        def go2surface(event):
            self.surface(ax[0], *surface['args'], **surface['kwargs'])
            fig.canvas.draw()
            fig.canvas.flush_events()
            
        box_img = plt.imread(path_file+'button_box.png')
        cursor_img = plt.imread(path_file+'button_cursor.jpeg')
        trash_img = plt.imread(path_file+'button_trash.jpg') 
        surface_img = plt.imread(path_file+'button_surface.png') 
        axbcursor = plt.axes([0.05, 0.709, 0.05, 0.05])
        axbbox = plt.axes([0.05, 0.65, 0.05, 0.05])
        axbtrash = plt.axes([0.05, 0.591, 0.05, 0.05], frameon=True, aspect='equal')
        bcursor = Button(axbcursor, '', image=cursor_img)
        bcursor.on_clicked(go2cursor)
        bbox = Button(axbbox, '', image=box_img)
        bbox.on_clicked(go2box)
        btrash = Button(axbtrash, '', image=trash_img, color='white', hovercolor='lime')
        btrash.on_clicked(go2trash)
        if len(surface['args'])>0:
            axbsurf = plt.axes([0.005, 0.689, 0.07, 0.07], frameon=True, aspect='equal')
            bsurf = Button(axbsurf, '', image=surface_img)
            bsurf.on_clicked(go2surface)
        plt.show()

        
    """
    #Lasso functions under development
    def _plot_lasso(self, ax, x, y, chan, color=False, show_path=True, extent=None, compare_cubes=[], **kwargs): 
        if len(self._lasso_path) == 0: return
        #for i in range(len(self.lasso_path))
        if extent is None:
            j = x.astype(np.int)
            i = y.astype(np.int)
        else: 
            nz, ny, nx = np.shape(self.data)
            dx = extent[1] - extent[0]
            dy = extent[3] - extent[2]
            j = (nx*(x-extent[0])/dx).astype(np.int)
            i = (ny*(y-extent[2])/dy).astype(np.int)
        
        if color: self._plot_path = ax[1].step(np.arange(len(i)), self.data[chan,i,j], color=color)
        else: self._plot_path = ax[1].step(np.arange(len(i)), self.data[chan,i,j])
        self._plot_color = self._plot_path[0].get_color()
        if show_path: self._path_on_cube = ax[0].plot(x,y, color=self._plot_color)
        else: self._path_on_cube = None

    def lasso(self, fig, ax, chan, color=False, show_path=True, extent=None, compare_cubes=[], **kwargs): 
        from matplotlib.widgets import LassoSelector
        canvas = ax[0].figure.canvas        
        def onselect(verts):
            #path = Path(verts)
            canvas.draw_idle()
            self._lasso_path.append(np.array(verts).T)
            self._plot_lasso(ax, *np.array(verts).T, chan, color, show_path, extent, compare_cubes, **kwargs)
            print (verts)
        def disconnect():
            self._lasso_obj.disconnect_events()
            canvas.draw_idle()
        self._lasso_obj = LassoSelector(ax[0], onselect, lineprops={'color': 'lime'})
        def onclick(event):
            if event.button == 3:
                print ('Right click. Disconnecting click event...')
                disconnect()
                fig.canvas.draw()
        cid = fig.canvas.mpl_connect('button_press_event', onclick) 
    """

    def curve(self, ax, x, y, chan, color=False, show_path=True, extent=None, compare_cubes=[], **kwargs): 
        kwargs_curve = dict(linewidth=2.5)#, label=r'x0:%d,x1:%d'%(x0,x1))
        kwargs_curve.update(kwargs)

        if extent is None:
            j = x.astype(np.int)
            i = y.astype(np.int)
        else: 
            nz, ny, nx = np.shape(self.data)
            dx = extent[1] - extent[0]
            dy = extent[3] - extent[2]
            j = (nx*(x-extent[0])/dx).astype(np.int)
            i = (ny*(y-extent[2])/dy).astype(np.int)

        pix_ids = np.arange(len(i))
        path_val = self.data[chan,i,j]
        if color: plot_path = ax[1].step(pix_ids, path_val, where='mid', color=color, **kwargs_curve)
        else: plot_path = ax[1].step(pix_ids, path_val, where='mid', **kwargs_curve)
        plot_color = plot_path[0].get_color()
        if show_path: path_on_cube = ax[0].plot(x, y, color=plot_color, **kwargs_curve)
        else: path_on_cube = None

        cube_fill = []
        plot_fill = None
        ncubes = len(compare_cubes)        
        if ncubes > 0:
            alpha = 0.2
            dalpha = -alpha/ncubes
            for cube in compare_cubes: 
                cube_fill.append(ax[1].fill_between(pix_ids, cube.data[chan,i,j], color=plot_color, step='mid', alpha=alpha))
                alpha+=dalpha
        else: plot_fill = ax[1].fill_between(pix_ids, path_val, color=plot_color, step='mid', alpha=0.2)

        return path_on_cube, plot_path, plot_color, plot_fill, cube_fill

    def show_path(self, x, y, extent=None, chan_init=20, compare_cubes=[], cursor_grid=True,
                  int_unit=r'Intensity [mJy beam$^{-1}$]', pos_unit='au', vel_unit=r'km s$^{-1}$',
                  show_beam=False, **kwargs):
        from matplotlib.widgets import Slider, Cursor, Button
        v0, v1 = self.channels[0], self.channels[-1]
        dv = v1-v0
        fig, ax = plt.subplots(ncols=2, figsize=(12,5))
        plt.subplots_adjust(wspace=0.25)

        y0, y1 = ax[1].get_position().y0, ax[1].get_position().y1
        axcbar = plt.axes([0.47, y0, 0.03, y1-y0])
        max_data = np.max(self.data)
        ax[0].set_xlabel(pos_unit)
        ax[0].set_ylabel(pos_unit)
        ax[1].set_xlabel('Pixel id along path')
        ax[1].tick_params(direction='in', right=True, labelright=False, labelleft=False)
        axcbar.tick_params(direction='out')
        ax[1].set_ylabel(int_unit, labelpad=15)
        ax[1].yaxis.set_label_position('right')
        #ax[1].set_xlim(v0-0.1, v1+0.1)
        #ax[1].set_ylim(-1, max_data)
        vmin, vmax = -max_data/30, max_data
        ax[1].set_ylim(vmin, vmax)
        ax[1].grid(lw=1.5, ls=':')
        cmap = plt.get_cmap('brg')
        cmap.set_bad(color=(0.9,0.9,0.9))

        if show_beam and self.beam_kernel: self._plot_beam(ax[0])

        img = ax[0].imshow(self.data[chan_init], cmap=cmap, extent=extent, origin='lower left', vmin=vmin, vmax=vmax)
        cbar = plt.colorbar(img, cax=axcbar)
        text_chan = ax[1].text(0.15, 1.04, #Converting xdata coords to Axes coords 
                               r'v$_{\rmchan}$=%4.1f %s'%(self.channels[chan_init], vel_unit), ha='center', 
                               color='black', transform=ax[1].transAxes)

        if cursor_grid: cg = Cursor(ax[0], useblit=True, color='lime', linewidth=1.5)
        box_img = plt.imread(path_file+'button_box.png')
        cursor_img = plt.imread(path_file+'button_cursor.jpeg')

        def get_interactive(func, chan=chan_init, color=False, show_path=True):
            return func(ax, x, y, chan, color=color, show_path=show_path, extent=extent, compare_cubes=compare_cubes, **kwargs)

        interactive_obj = [get_interactive(self.interactive_path)]
        #***************
        #SLIDERS
        #***************
        def update_chan(val):
            chan = int(val)
            vchan = self.channels[chan]
            img.set_data(self.data[chan])
            text_chan.set_text(r'v$_{\rmchan}$=%4.1f %s'%(vchan, vel_unit))
            path_on_cube, plot_path, plot_color, plot_fill, cube_fill = interactive_obj[0]
            plot_path[0].remove()
            if plot_fill is not None: plot_fill.remove()
            for cbfill in cube_fill: cbfill.remove()
            interactive_obj[0] = get_interactive(self.interactive_path, chan, color=plot_color, show_path=False)
            fig.canvas.draw_idle()

        def update_cubes(val):
            i = int(slider_cubes.val)
            chan = int(slider_chan.val)
            vchan = self.channels[chan]
            if i==0: img.set_data(self.data[chan])
            else: img.set_data(compare_cubes[i-1].data[chan])
            text_chan.set_text(r'v$_{\rmchan}$=%4.1f %s'%(vchan, vel_unit))
            path_on_cube, plot_path, plot_color, plot_fill, cube_fill = interactive_obj[0]
            plot_path[0].remove()
            if plot_fill is not None: plot_fill.remove()
            for cbfill in cube_fill: cbfill.remove()
            interactive_obj[0] = get_interactive(self.interactive_path, chan, color=plot_color, show_path=False)
            fig.canvas.draw_idle()

        ncubes = len(compare_cubes)
        if ncubes>0:
            axcubes = plt.axes([0.2, 0.90, 0.24, 0.025], facecolor='0.7')
            axchan = plt.axes([0.2, 0.95, 0.24, 0.025], facecolor='0.7')
            slider_cubes = Slider(axcubes, 'Cube id', 0, ncubes, 
                                  valstep=1, valinit=0, valfmt='%1d', color='dodgerblue')                                  
            slider_chan = Slider(axchan, 'Channel', 0, self.nchan-1, 
                                 valstep=1, valinit=chan_init, valfmt='%2d', color='dodgerblue')        
            slider_cubes.on_changed(update_cubes)
            slider_chan.on_changed(update_cubes)
        else: 
            axchan = plt.axes([0.2, 0.9, 0.24, 0.05], facecolor='0.7')
            slider_chan = Slider(axchan, 'Channel', 0, self.nchan-1, 
                                 valstep=1, valinit=chan_init, valfmt='%2d', color='dodgerblue')        
            slider_chan.on_changed(update_chan)

        plt.show()

        """
        self._path_on_cube, self._plot_path, self._plot_color = None, None, None
        self._lasso_path = []
        self.interactive_path(fig, ax, chan_init, color=False, show_path=True, extent=extent, compare_cubes=compare_cubes, **kwargs)

        def get_interactive(func, chan=chan_init, color=False, show_path=True):
            #func(fig, ax, chan, color=color, show_path=show_path, extent=extent, compare_cubes=compare_cubes, **kwargs)
            if func == self.lasso:
                return self._plot_lasso(ax, True, True, chan, color=color, show_path=show_path, extent=extent, compare_cubes=compare_cubes, **kwargs)
        
        #interactive_obj = [get_interactive(self.interactive_path)]
        #print (interactive_obj)
        #***************
        #SLIDERS
        #***************
        def update_chan(val):
            chan = int(val)
            vchan = self.channels[chan]
            img.set_data(self.data[chan])
            current_chan.set_xdata(vchan)
            text_chan.set_x((vchan-v0)/dv)
            text_chan.set_text('%4.1f km/s'%vchan)
            #path_on_cube, plot_path, plot_color = interactive_obj[0]
            if self._path_on_cube is not None: 
                self._plot_path[0].remove()
                get_interactive(self.interactive_path, chan, color=self._plot_color, show_path=False)
            fig.canvas.draw_idle()
        """

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


class Height:
    @property
    def z_near_func(self): 
        return self._z_near_func
          
    @z_near_func.setter 
    def z_near_func(self, near): 
        print('Setting near-side height function to', near) 
        self._z_near_func = near

    @z_near_func.deleter 
    def z_near_func(self): 
        print('Deleting near-side height function') 
        del self._z_near_func

    @property
    def z_far_func(self): 
        return self._z_far_func
          
    @z_far_func.setter 
    def z_far_func(self, far): 
        print('Setting far-side height function to', far) 
        self._z_far_func = far

    @z_far_func.deleter 
    def z_far_func(self): 
        print('Deleting far-side height function') 
        del self._z_far_func

    psi0 = 15*np.pi/180
    @staticmethod
    def z_cone(coord, psi=psi0):
        R = coord['R'] 
        z = np.tan(psi) * R
        return z

    @staticmethod
    def z_cone_neg(coord, psi=psi0):
        return -Height.z_cone(coord, psi)


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
    def linewidth_powerlaw(coord, L0=0.2, p=-0.4, q=0.3, R0=100*sfu.au, z0=100*sfu.au): #A=600.0, p=-0.4, q=0.3): #
        if 'R' not in coord.keys(): R = hypot_func(coord['x'], coord['y'])
        else: R = coord['R'] 
        z = coord['z']        
        A = L0*R0**-p*z0**-q
        return A*R**p*np.abs(z)**q


class Lineslope:
    @property
    def lineslope_func(self): 
        return self._lineslope_func
          
    @lineslope_func.setter 
    def lineslope_func(self, lineslope): 
        print('Setting lineslope function to', lineslope) 
        self._lineslope_func = lineslope

    @lineslope_func.deleter 
    def lineslope_func(self): 
        print('Deleting lineslope function') 
        del self._lineslope_func

    @staticmethod
    def lineslope_powerlaw(coord, Ls=5.0, p=0.0, q=0.0, R0=100*sfu.au, z0=100*sfu.au): #A=600.0, p=-0.4, q=0.3): #
        if p==0.0 and q==0.0:
            return Ls
        else:
            if 'R' not in coord.keys(): R = hypot_func(coord['x'], coord['y'])
            else: R = coord['R'] 
            z = coord['z']        
            A = Ls*R0**-p*z0**-q
            return A*R**p*np.abs(z)**q


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
    def keplerian(coord, Mstar=1.0, vel_sign=1, vsys=0):
        Mstar *= sfu.MSun
        if 'R' not in coord.keys(): R = hypot_func(coord['x'], coord['y'])
        else: R = coord['R'] 
        return vel_sign*np.sqrt(G*Mstar/R) * 1e-3
    
    @staticmethod
    def keplerian_vertical(coord, Mstar=1.0, vel_sign=1, vsys=0):
        Mstar *= sfu.MSun
        if 'R' not in coord.keys(): R = hypot_func(coord['x'], coord['y'])
        else: R = coord['R'] 
        if 'r' not in coord.keys(): r = hypot_func(R, coord['z'])
        else: r = coord['r']
        return vel_sign*np.sqrt(G*Mstar/r**3)*R * 1e-3 #to km/s


class Intensity:   
    @property
    def beam_info(self):
        return self._beam_info

    @beam_info.setter 
    def beam_info(self, beam_info): 
        print('Setting beam_info var to', beam_info)
        self._beam_info = beam_info

    @beam_info.deleter 
    def beam_info(self): 
        print('Deleting beam_info var') 
        del self._beam_info     

    @property
    def beam_kernel(self):
        return self._beam_kernel

    @beam_kernel.setter 
    def beam_kernel(self, beam_kernel): 
        print('Setting beam_kernel var to', beam_kernel)
        x_stddev = beam_kernel.model.x_stddev.value
        y_stddev = beam_kernel.model.y_stddev.value
        self._beam_area = 2*np.pi*x_stddev*y_stddev
        self._beam_kernel = beam_kernel

    @beam_kernel.deleter 
    def beam_kernel(self): 
        print('Deleting beam_kernel var') 
        del self._beam_kernel     

    @property
    def beam_from(self):
        return self._beam_from

    @beam_from.setter 
    def beam_from(self, file): 
        #Rework this, missing beam kwargs info
        print('Setting beam_from var to', file)
        if file: self.beam_info, self.beam_kernel = Tools.get_beam_from(file) #Calls beam_kernel setter
        self._beam_from = file

    @beam_from.deleter 
    def beam_from(self): 
        print('Deleting beam_from var') 
        del self._beam_from     

    @property
    def use_temperature(self):
        return self._use_temperature

    @use_temperature.setter 
    def use_temperature(self, use): 
        use = bool(use)
        print('Setting use_temperature var to', use)
        if use: self.line_profile = self.line_profile_temp
        #else: self.line_profile = self.line_profile_v_sigma
        self._use_temperature = use

    @use_temperature.deleter 
    def use_temperature(self): 
        print('Deleting use_temperature var') 
        del self._use_temperature

    @property
    def use_full_channel(self):
        return self._use_full_channel

    @use_full_channel.setter 
    def use_full_channel(self, use): #Needs remake, there is now a new kernel (line_profile_bell) 
        use = bool(use)
        print('Setting use_full_channel var to', use)
        if use: 
            if self.use_temperature: self.line_profile = self.line_profile_temp_full
            else: self.line_profile = self.line_profile_v_sigma_full
        else: 
            if self.use_temperature: self.line_profile = self.line_profile_temp
            else: self.line_profile = self.line_profile_v_sigma
        self._use_full_channel = use

    @use_full_channel.deleter 
    def use_full_channel(self): 
        print('Deleting use_full_channel var') 
        del self._use_full_channel

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
    def intensity_powerlaw(coord, I0=30.0, R0=100*sfu.au, p=-0.4, z0=100*sfu.au, q=0.3): #A=600.0, p=-0.4, q=0.3): #
        if 'R' not in coord.keys(): R = hypot_func(coord['x'], coord['y'])
        else: R = coord['R'] 
        z = coord['z']        
        A = I0*R0**-p*z0**-q
        return A*R**p*np.abs(z)**q
        
    @staticmethod
    def nuker(coord, I0=30.0, Rt=100*sfu.au, alpha=-0.5, gamma=0.1, beta=0.2):
        if 'R' not in coord.keys(): R = hypot_func(coord['x'], coord['y'])
        else: R = coord['R'] 
        A = I0*Rt**gamma
        return A*(R**-gamma) * (1+(R/Rt)**alpha)**((gamma-beta)/alpha)

    @staticmethod
    def line_profile_subchannel(line_profile_func, v_chan, v, v_sigma, b_slope, channel_width=0.1, **kwargs): #Not currently used
        half_chan = 0.5*channel_width
        v0 = v_chan - half_chan
        v1 = v_chan + half_chan
        nsub = 10
        vsub = np.linspace(v0, v1, nsub)
        dvsub = vsub[1]-vsub[0]
        J = 0
        for vs in vsub:
            J += line_profile_func(vs, v, v_sigma, b_slope, **kwargs) 
        J = J * dvsub/channel_width
        return J
        
    @staticmethod
    def line_profile_temp(v_chan, v, T, dum, v_turb=0.0, mmol=2*sfu.amu):
        v_sigma = np.sqrt(kb*T/mmol + v_turb**2) * 1e-3 #in km/s
        #return 1/(np.sqrt(np.pi)*v_sigma) * np.exp(-((v-v_chan)/v_sigma)**2)
        return np.exp(-0.5*((v-v_chan)/v_sigma)**2)

    @staticmethod
    def line_profile_temp_full(v_chan, v, T, dum, v_turb=0, mmol=2*sfu.amu, channel_width=0.1):
        v_sigma = np.sqrt(kb*T/mmol + v_turb**2) * 1e-3 #in km/s
        half_chan = 0.5*channel_width
        v0 = v_chan - half_chan
        v1 = v_chan + half_chan
        nsub = 10
        vsub = np.linspace(v0, v1, nsub)
        dvsub = vsub[1]-vsub[0]
        J = 0
        for vs in vsub:
            J += np.exp(-0.5*((v-vs)/v_sigma)**2)
        J = J * dvsub/channel_width
        return J

    @staticmethod
    def line_profile_v_sigma(v_chan, v, v_sigma, dum):
        #return 1/(np.sqrt(2*np.pi)*v_sigma) * np.exp(-0.5*((v-v_chan)/v_sigma)**2)
        #return np.where(np.abs(v-v_chan) < 0.5*v_sigma, 1, 0)
        return np.exp(-0.5*((v-v_chan)/v_sigma)**2)
    
    @staticmethod
    def line_profile_v_sigma_full(v_chan, v, v_sigma, dum, channel_width=0.1):
        half_chan = 0.5*channel_width
        v0 = v_chan - half_chan
        v1 = v_chan + half_chan
        nsub = 10
        vsub = np.linspace(v0, v1, nsub)
        dvsub = vsub[1]-vsub[0]
        J = 0
        for vs in vsub:
            J += np.exp(-0.5*((v-vs)/v_sigma)**2)
        J = J * dvsub/channel_width
        return J

    @staticmethod
    def line_profile_bell(v_chan, v, v_sigma, b_slope):
        return 1/(1+np.abs((v-v_chan)/v_sigma)**(2*b_slope))        

    @staticmethod
    def line_profile_bell_full(v_chan, v, v_sigma, b_slope, channel_width=0.1):
        half_chan = 0.5*channel_width
        v0 = v_chan - half_chan
        v1 = v_chan + half_chan
        nsub = 10
        vsub = np.linspace(v0, v1, nsub)
        dvsub = vsub[1]-vsub[0]
        J = 0
        for vs in vsub:
            J += 1/(1+np.abs((v-vs)/v_sigma)**(2*b_slope))
        J = J * dvsub/channel_width
        return J

    def get_line_profile(self, v_chan, vel2d, linew2d, lineb2d, **kwargs):
        if self.subpixels:
            v_near, v_far = [], []
            for i in range(self.subpixels_sq):
                v_near.append(self.line_profile(v_chan, vel2d[i]['near'], linew2d['near'], lineb2d['near'], **kwargs))
                v_far.append(self.line_profile(v_chan, vel2d[i]['far'], linew2d['far'], lineb2d['far'], **kwargs))

            integ_v_near = np.sum(np.array(v_near), axis=0) * self.sub_dA / self.pix_dA
            integ_v_far = np.sum(np.array(v_far), axis=0) * self.sub_dA / self.pix_dA
            return integ_v_near, integ_v_far
        
        else: 
            v_near = self.line_profile(v_chan, vel2d['near'], linew2d['near'], lineb2d['near'], **kwargs)
            v_far = self.line_profile(v_chan, vel2d['far'], linew2d['far'], lineb2d['far'], **kwargs)
            return v_near, v_far 

    def get_channel(self, velocity2d, intensity2d, linewidth2d, lineslope2d, v_chan, **kwargs):                    
        vel2d, int2d, linew2d, lineb2d = velocity2d, {}, {}, {}

        if isinstance(intensity2d, numbers.Number): int2d['near'] = int2d['far'] = intensity2d
        else: int2d = intensity2d
        if isinstance(linewidth2d, numbers.Number): linew2d['near'] = linew2d['far'] = linewidth2d
        else: linew2d = linewidth2d
        if isinstance(lineslope2d, numbers.Number): lineb2d['near'] = lineb2d['far'] = lineslope2d
        else: lineb2d = lineslope2d
    
        v_near, v_far = self.get_line_profile(v_chan, vel2d, linew2d, lineb2d, **kwargs)

        v_near_clean = np.where(np.isnan(v_near), -np.inf, v_near)
        v_far_clean = np.where(np.isnan(v_far), -np.inf, v_far)
        
        #int2d_near = int2d['near'] * v_near_clean / v_near_clean.max()
        int2d_near = np.where(np.isnan(int2d['near']), -np.inf, int2d['near'] * v_near_clean)# / v_near_clean.max())
        int2d_far = np.where(np.isnan(int2d['far']), -np.inf, int2d['far'] * v_far_clean)# / v_far_clean.max())
        
        #vmap_full = np.array([v_near_clean, v_far_clean]).max(axis=0)
        int2d_full = np.array([int2d_near, int2d_far]).max(axis=0)
        
        if self.beam_kernel:
            inf_mask = np.isinf(int2d_full)
            int2d_full = np.where(inf_mask, 0.0, int2d_full) 
            int2d_full = self._beam_area*convolve(int2d_full, self.beam_kernel, preserve_nan=False)

        return int2d_full

    #def get_cube(self, vchan0, vchan1, velocity2d, intensity2d, linewidth2d, lineslope2d, nchan=30, tb={'nu': False, 'beam': False}, **kwargs):
    def get_cube(self, vchannels, velocity2d, intensity2d, linewidth2d, lineslope2d, 
                 nchan=None, tb={'nu': False, 'beam': False}, return_data_only=False, **kwargs):
        vel2d, int2d, linew2d, lineb2d = velocity2d, {}, {}, {}
        line_profile = self.line_profile
        if nchan is None: nchan=len(vchannels)
        #channels = np.linspace(vchan0, vchan1, num=nchan)
        cube = []

        if isinstance(intensity2d, numbers.Number): int2d['near'] = int2d['far'] = intensity2d
        else: int2d = intensity2d
        if isinstance(linewidth2d, numbers.Number): linew2d['near'] = linew2d['far'] = linewidth2d
        else: linew2d = linewidth2d
        if isinstance(lineslope2d, numbers.Number): lineb2d['near'] = lineb2d['far'] = lineslope2d
        else: lineb2d = lineslope2d

        int2d_near_nan = np.isnan(int2d['near'])#~int2d['near'].mask#
        int2d_far_nan = np.isnan(int2d['far'])#~int2d['far'].mask#
        if self.subpixels:
            vel2d_near_nan = np.isnan(vel2d[self.sub_centre_id]['near'])
            vel2d_far_nan = np.isnan(vel2d[self.sub_centre_id]['far'])
        else:
            vel2d_near_nan = np.isnan(vel2d['near'])#~vel2d['near'].mask#
            vel2d_far_nan = np.isnan(vel2d['far'])#~vel2d['far'].mask#

        #for i, v_chan in enumerate(vchannels):
        #viter = iter(vchannels)
        #for _ in itertools.repeat(None, nchan):
        for i in range(nchan):
            v_near, v_far = self.get_line_profile(vchannels[i], vel2d, linew2d, lineb2d, **kwargs)
            v_near_clean = np.where(vel2d_near_nan, -np.inf, v_near)
            v_far_clean = np.where(vel2d_far_nan, -np.inf, v_far)
            
            int2d_near = np.where(int2d_near_nan, -np.inf, int2d['near'] * v_near_clean)# / v_near_clean.max())
            int2d_far = np.where(int2d_far_nan, -np.inf, int2d['far'] * v_far_clean)# / v_far_clean.max())        
            #vmap_full = np.array([v_near_clean, v_far_clean]).max(axis=0)
            int2d_full = np.array([int2d_near, int2d_far]).max(axis=0)

            if self.beam_kernel:
                inf_mask = np.isinf(int2d_full)
                int2d_full = np.where(inf_mask, 0.0, int2d_full)
                int2d_full = self._beam_area*convolve(int2d_full, self.beam_kernel, preserve_nan=False)

            cube.append(int2d_full)
            
        if return_data_only: return np.asarray(cube)
        else: return Cube(nchan, vchannels, np.asarray(cube), beam=self.beam_info, beam_kernel=self.beam_kernel, tb=tb)

    @staticmethod
    def make_channels_movie(vchan0, vchan1, velocity2d, intensity2d, linewidth2d, nchans=30, folder='./movie_channels/', **kwargs):
        import matplotlib.pyplot as plt
        channels = np.linspace(vchan0, vchan1, num=nchans)
        int2d_cube = []
        for i, vchan in enumerate(channels):
            int2d = Intensity.get_channel(velocity2d, intensity2d, linewidth2d, lineslope2d, vchan, **kwargs)
            int2d_cube.append(int2d)
            extent = [-600, 600, -600, 600]
            plt.imshow(int2d, cmap='binary', extent=extent, origin='lower left', vmax=np.max(linewidth2d['near']))
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


class Mcmc:
    @staticmethod
    def _get_params2fit(mc_params, boundaries):
        header = []
        kind = []
        params_indices = {}
        boundaries_list = []
        check_param2fit = lambda val: val and isinstance(val, bool)
        i = 0
        for key in mc_params:
            if isinstance(mc_params[key], dict):
                params_indices[key] = {}
                for key2 in mc_params[key]: 
                    if check_param2fit(mc_params[key][key2]):
                        header.append(key2)
                        kind.append(key)
                        boundaries_list.append(boundaries[key][key2])
                        params_indices[key][key2] = i
                        i+=1
            else: raise InputError(mc_params, 'Wrong input parameters. Base keys in mc_params must be categories; parameters of a category must be within a dictionary as well.')

        return header, kind, len(header), boundaries_list, params_indices
    
    @staticmethod
    def plot_walkers(samples, best_params, nstats=None, header=None, kind=None, tag=''):
        npars, nsteps, nwalkers = samples.shape
        if kind is not None:
            ukind, neach = np.unique(kind, return_counts=True)
            ncols = len(ukind)
            nrows = np.max(neach)
        else:
            ukind = [''] 
            ncols = 1
            nrows = npars
            kind = ['' for i in range(nrows)] 
        
        if header is not None:
            if len(header) != npars: raise InputError(header, 'Number of headers must be equal to number of parameters')
            
        kind_col = {ukind[i]: i for i in range(ncols)}
        col_count = np.zeros(ncols).astype('int')

        fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(3.*ncols, 3*nrows))

        print (header, kind, samples.shape)
        x0_hline = 0
        for k, key in enumerate(kind):
            j = kind_col[key]
            i = col_count[j] 
            for walker in samples[k].T:
                if ncols == 1: axij = ax[i]
                else: axij = ax[i][j]
                axij.plot(walker, alpha=0.1, lw=1.0, color='k')
                if header is not None: 
                    #axij.set_ylabel(header[k])
                    axij.text(0.06, 0.84, header[k], va='baseline', fontsize=MEDIUM_SIZE-2, transform=axij.transAxes, rotation=90)
            if i==0: axij.set_title(key, pad=10)
            if nstats is not None: 
                axij.axvline(nsteps-nstats, ls=':', lw=2, color='r')
                x0_hline = nsteps-nstats

            axij.plot([x0_hline, nsteps], [best_params[k]]*2, ls='-', lw=3, color='dodgerblue')
            axij.text((nsteps-1)+0.024*nsteps, best_params[k], '%.3f'%best_params[k], va='center', color='dodgerblue', fontsize=SMALL_SIZE+2, rotation=90) 
            axij.tick_params(axis='y', which='major', labelsize=SMALL_SIZE, rotation=45)
            axij.set_xlim(None, nsteps-1 + 0.01*nsteps)
            col_count[j]+=1

        for j in range(ncols):
            i_last = col_count[j]-1
            if ncols==1: ax[i_last].set_xlabel('Steps')
            else: ax[i_last][j].set_xlabel('Steps')
            if i_last < nrows-1: #Remove empty axes
                for k in range((nrows-1)-i_last): ax[nrows-1-k][j].axis('off')

        plt.tight_layout()
        plt.savefig('mc_walkers_%s_%dwalkers_%dsteps.png'%(tag, nwalkers, nsteps))
        plt.close()


    @staticmethod
    def plot_corner(samples, labels=None, quantiles=None):
        """Plot the corner plot to check for covariances."""
        import corner
        quantiles = [0.16, 0.5, 0.84] if quantiles is None else quantiles
        corner.corner(samples, labels=labels, title_fmt='.4f', bins=30,
                      quantiles=quantiles, show_titles=True)
        
    def ln_likelihood(self, new_params, **kwargs):
        for i in range(self.mc_nparams):
            if not (self.mc_boundaries_list[i][0] < new_params[i] < self.mc_boundaries_list[i][1]): return -np.inf
            else: self.params[self.mc_kind[i]][self.mc_header[i]] = new_params[i]

        vel2d, int2d, linew2d, lineb2d = self.make_model(**kwargs)

        lnx2=0    
        model_cube = self.get_cube(self.channels, vel2d, int2d, linew2d, lineb2d, nchan=self.nchan, return_data_only=True)#, tb = {'nu': 230, 'beam': self.beam_info})
        for i in range(self.nchan):
            model_chan = model_cube[i] #model_cube.data[i] #self.get_channel(vel2d, int2d, linew2d, lineb2d, self.channels[i])
            mask_data = np.isfinite(self.data[i])
            mask_model = np.isfinite(model_chan)
            data = np.where(np.logical_and(mask_model, ~mask_data), 0, self.data[i])
            model = np.where(np.logical_and(mask_data, ~mask_model), 0, model_chan)
            mask = np.logical_and(mask_data, mask_model)
            lnx =  np.where(mask, np.power((data - model)/self.noise_stddev, 2), 0) 
            #lnx = -0.5 * np.sum(lnx2[~np.isnan(lnx2)] * 0.00001)# * self.ivar)
            lnx2 += -0.5 * np.sum(lnx)
            
        #print (new_params, "\nLOG LIKELIHOOD %.4e"%lnx2)
        return lnx2 if np.isfinite(lnx2) else -np.inf
    
     
class General2d(Height, Velocity, Intensity, Linewidth, Lineslope, Tools, Mcmc): #Inheritance should only be from Intensity and Mcmc, the others contain just staticmethods...
    def __init__(self, grid, prototype=False, subpixels=False, beam=None, kwargs_beam={}):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self.prototype = prototype

        self._beam_info = False
        self._beam_from = False #Should be deprecated
        self._beam_kernel = False
        self._beam_area = False 
        if beam is not None: 
            self.beam_info, self.beam_kernel = Tools.get_beam_from(beam, grid, **kwargs_beam)

        self._z_near_func = General2d.z_cone
        self._z_far_func = General2d.z_cone_neg
        self._velocity_func = General2d.keplerian
        self._intensity_func = General2d.intensity_powerlaw
        self._linewidth_func = General2d.linewidth_powerlaw
        self._lineslope_func = General2d.lineslope_powerlaw
        self._line_profile = General2d.line_profile_bell
        self._use_temperature = False
        self._use_full_channel = False
 
        x_true, y_true = grid.XYZ[:2] 
        self.phi_true = np.arctan2(y_true, x_true) #grid.rRTP[3] 
        self.R_true = hypot_func(x_true, y_true) #grid.rRTP[1] #Slightly different as in the grid object the pixels R=0 actually take the closest-neighbour value. Current approach masks r,R=0
        self.x_true, self.y_true = x_true, y_true
        self.mesh = np.meshgrid(grid.XYZgrid[0], grid.XYZgrid[1])

        if subpixels and isinstance(subpixels, int):
            if subpixels%2 == 0: subpixels+=1 #If input even becomes odd to contain pxl centre
            pix_size = grid.step[0]
            dx = dy = pix_size / subpixels
            centre = int(round((subpixels-1)/2.))
            centre_sq = int(round((subpixels**2-1)/2.))
            x_shift = np.arange(0, subpixels*dx, dx) - dx*centre
            y_shift = np.arange(0, subpixels*dy, dy) - dy*centre
            sub_x_true = [x_true + x0 for x0 in x_shift]
            sub_y_true = [y_true + y0 for y0 in y_shift]
            self.sub_R_true = [[hypot_func(sub_x_true[j], sub_y_true[i]) for j in range(subpixels)] for i in range(subpixels)]
            self.sub_phi_true = [[np.arctan2(sub_y_true[i], sub_x_true[j]) for j in range(subpixels)] for i in range(subpixels)]
            self.sub_x_true = sub_x_true
            self.sub_y_true = sub_x_true
            self.sub_dA = dx*dy
            self.pix_dA = pix_size**2
            self.sub_centre_id = centre_sq
            self.subpixels = subpixels
            self.subpixels_sq = subpixels**2
        else: self.subpixels=False

        #Get and print default parameters for default functions
        self.categories = ['velocity', 'orientation', 'intensity', 'linewidth', 'lineslope', 'height_near', 'height_far']

        self.mc_params = {'velocity': {'Mstar': True, 
                                       'vel_sign': 1,
                                       'vsys': 0},
                          'orientation': {'incl': True, 
                                          'PA': True},
                          'intensity': {'I0': True, 
                                        'p': True, 
                                        'q': False},
                          'linewidth': {'L0': True, 
                                        'p': True, 
                                        'q': 0.1}, 
                          'lineslope': {'Ls': False, 
                                        'p': False, 
                                        'q': False},
                          'height_near': {'psi': True},
                          'height_far': {'psi': True},
                          }
        
        self.mc_boundaries = {'velocity': {'Mstar': [0.05, 5.0],
                                           'vsys': [-10, 10]},
                              'orientation': {'incl': [-np.pi/3, np.pi/3], 
                                              'PA': [-np.pi, np.pi]},
                              'intensity': {'I0': [0, 100], 
                                            'p': [-5.0, 5.0], 
                                            'q': [0, 5.0]},
                              'linewidth': {'L0': [0.005, 5.0], 
                                            'p': [-5.0, 5.0], 
                                            'q': [-5.0, 5.0]},
                              'lineslope': {'Ls': [0.005, 100], 
                                            'p': [-5.0, 5.0], 
                                            'q': [-5.0, 5.0]},
                              'height_near': {'psi': [0, np.pi/2]},
                              'height_far': {'psi': [0, np.pi/2]}
                              }

        if prototype:
            self.params = {}
            for key in self.categories: self.params[key] = {}
            print ('Available categories for prototyping:', self.params)

        else: 
            self.mc_header, self.mc_kind, self.mc_nparams, self.mc_boundaries_list, self.mc_params_indices = General2d._get_params2fit(self.mc_params, self.mc_boundaries)
            print ('Default parameter header for mcmc model fitting:', self.mc_header)
            print ('Default parameters to fit and fixed parameters:', self.mc_params)

    def run_mcmc(self, data, channels, p0_mean='optimize', p0_stddev=1e-3, noise_stddev=1.0,
                 nwalkers=30, nsteps=100, frac_stats=0.5, frac_stddev=1e-3, mc_layers=1, z_mirror=False, 
                 custom_header={}, custom_kind={}, tag='',
                 plot_walkers=True, plot_corner=True, **kwargs_model): #p0 from 'optimize', 'min', 'max', list of values.
        self.data = data
        self.channels = channels
        self.nchan = len(channels)
        self.noise_stddev = noise_stddev

        kwargs_model.update({'z_mirror': z_mirror})
        if z_mirror: 
            for key in self.mc_params['height_far']: self.mc_params['height_far'][key] = 'height_near'
        self.mc_header, self.mc_kind, self.mc_nparams, self.mc_boundaries_list, self.mc_params_indices = General2d._get_params2fit(self.mc_params, self.mc_boundaries)
        self.params = copy.deepcopy(self.mc_params)

        print ('Parameter header set for mcmc model fitting:', self.mc_header)
        print ('Parameters to fit and fixed parameters:', self.mc_params)            
        print ('Number of mc parameters:', self.mc_nparams)
        print ('Kind of parameters:', self.mc_kind)
        print ('Parameter boundaries:', self.mc_boundaries_list)
        
        if p0_mean == 'optimize': p0_mean = optimize_p0()
        if isinstance(p0_mean, (list, tuple, np.ndarray)): 
            if len(p0_mean) != self.mc_nparams: raise InputError(p0_mean, 'Length of input p0_mean must be equal to number of parameters to fit: %d'%self.mc_nparams)
            else: pass
        print ('Mean for initial guess p0:', p0_mean)

        nstats = int(round(frac_stats*(nsteps-1))) #python2 round returns float, python3 returns int
        ndim = self.mc_nparams

        p0_stddev = [frac_stddev*(self.mc_boundaries_list[i][1] - self.mc_boundaries_list[i][0]) for i in range(self.mc_nparams)]
        print ('p0 pars stddev:', p0_stddev)
        p0 = np.random.normal(loc=p0_mean,
                              scale=p0_stddev,
                              size=(nwalkers, ndim)
                              )

        with Pool() as pool:
            sampler = emcee.EnsembleSampler(nwalkers, ndim, self.ln_likelihood, pool=pool, kwargs=kwargs_model)                                                        
            start = time.time()
            sampler.run_mcmc(p0, nsteps, progress=True)
            end = time.time()
            multi_time = end - start
            print("Multiprocessing took {0:.1f} seconds".format(multi_time))

        samples = sampler.chain[:, -nstats:] #3d matrix, shape (nwalkers, nstats, npars)
        samples = samples.reshape(-1, samples.shape[-1]) #2d matrix, shape (nwalkers*nstats, npars). With the -1 np guesses the dimensionality
        best_params = np.median(samples, axis=0)
        self.best_params = best_params
        print ('Median from parameter walkers for the last %d steps:'%nstats, list(zip(self.mc_header, best_params)))

        #************
        #PLOTTING
        #************
        for key in custom_header: self.mc_header[key] = custom_header[key]
        if plot_walkers: 
            Mcmc.plot_walkers(sampler.chain.T, best_params, header=self.mc_header, kind=self.mc_kind, nstats=nstats, tag=tag)
        if plot_corner: 
            Mcmc.plot_corner(samples, labels=self.mc_header)
            plt.savefig('mc_corner_%s_%dwalkers_%dsteps.png'%(tag, nwalkers, nsteps))
            plt.close()

    @staticmethod
    def orientation(incl=np.pi/4, PA=0.0):
        return incl, PA

    def get_projected_coords(self, z_mirror=False, R_inner=0, R_disc=None, 
                             R_nan_val=0, phi_nan_val=10*np.pi, z_nan_val=0):

        from scipy.interpolate import griddata
        #*************************************
        #MAKE TRUE GRID FOR NEAR AND FAR SIDES
        if self.prototype: print ('Getting projected coords for prototype model:', self.params)
        
        incl, PA = General2d.orientation(**self.params['orientation'])
        cos_incl, sin_incl = np.cos(incl), np.sin(incl)

        z_true = {}
        z_true['near'] = self.z_near_func({'R': self.R_true}, **self.params['height_near'])

        if z_mirror: z_true['far'] = -z_true['near']
        else: z_true['far'] = self.z_far_func({'R': self.R_true}, **self.params['height_far']) 
            
        grid_true = {'near': [self.x_true, self.y_true, z_true['near'], self.R_true, self.phi_true], 
                     'far': [self.x_true, self.y_true, z_true['far'], self.R_true, self.phi_true]}
        
        #***********************************
        #PROJECT PROPERTIES ON THE SKY PLANE        
        R, phi, z = {}, {}, {}
        for side in ['near', 'far']:
            xt, yt, zt = grid_true[side][:3]
            x_pro, y_pro, z_pro = self._project_on_skyplane(xt, yt, zt, cos_incl, sin_incl)
            if PA: x_pro, y_pro = self._rotate_sky_plane(x_pro, y_pro, PA)             
            R[side] = griddata((x_pro, y_pro), self.R_true, (self.mesh[0], self.mesh[1]), method='linear')
            x_grid = griddata((x_pro, y_pro), xt, (self.mesh[0], self.mesh[1]), method='linear')
            y_grid = griddata((x_pro, y_pro), yt, (self.mesh[0], self.mesh[1]), method='linear')
            phi[side] = np.arctan2(y_grid, x_grid) 
            #Since this one is periodic it has to be recalculated, otherwise the interpolation will screw up things at the boundary -np.pi->np.pi
            # When plotting contours there seems to be in any case some sort of interpolation, so there is still problems at the boundary
            #phi[side] = griddata((x_pro, y_pro), self.phi_true, (self.mesh[0], self.mesh[1]), method='linear')
            z[side] = griddata((x_pro, y_pro), z_true[side], (self.mesh[0], self.mesh[1]), method='linear')
            #r[side] = hypot_func(R[side], z[side])
            if R_disc is not None: 
                for prop in [R, phi, z]: prop[side] = np.where(np.logical_and(R[side]<R_disc, R[side]>R_inner), prop[side], np.nan)
            
        R_nonan, phi_nonan, z_nonan = None, None, None
        if R_nan_val is not None: R_nonan = {side: np.where(np.isnan(R[side]), R_nan_val, R[side]) for side in ['near', 'far']}
        if phi_nan_val is not None: phi_nonan = {side: np.where(np.isnan(phi[side]), phi_nan_val, phi[side]) for side in ['near', 'far']}
        if z_nan_val is not None: z_nonan = {side: np.where(np.isnan(z[side]), z_nan_val, z[side]) for side in ['near', 'far']}

        return R, phi, z, R_nonan, phi_nonan, z_nonan
        
    def make_model(self, z_mirror=False, R_inner=0, R_disc=None):
                   
        from scipy.interpolate import griddata
        #*************************************
        #MAKE TRUE GRID FOR NEAR AND FAR SIDES
        if self.prototype: print ('Prototype model:', self.params)
        
        incl, PA = General2d.orientation(**self.params['orientation'])
        int_kwargs = self.params['intensity']
        vel_kwargs = self.params['velocity']
        lw_kwargs = self.params['linewidth']
        ls_kwargs = self.params['lineslope']

        cos_incl, sin_incl = np.cos(incl), np.sin(incl)

        z_true = self.z_near_func({'R': self.R_true}, **self.params['height_near'])

        if z_mirror: z_true_far = -z_true
        else: z_true_far = self.z_far_func({'R': self.R_true}, **self.params['height_far']) 
            
        grid_true = {'near': [self.x_true, self.y_true, z_true, self.R_true, self.phi_true], 
                     'far': [self.x_true, self.y_true, z_true_far, self.R_true, self.phi_true]}

        #*******************************
        #COMPUTE PROPERTIES ON SKY GRID #This will no longer be necessary as all four functions will always be called
        avai_kwargs = [vel_kwargs, int_kwargs, lw_kwargs, ls_kwargs]
        avai_funcs = [self.velocity_func, self.intensity_func, self.linewidth_func, self.lineslope_func]
        true_kwargs = [isinstance(kwarg, dict) for kwarg in avai_kwargs]
        prop_kwargs = [kwarg for i, kwarg in enumerate(avai_kwargs) if true_kwargs[i]]
        prop_funcs = [func for i, func in enumerate(avai_funcs) if true_kwargs[i]]
       
        if self.subpixels:
            subpix_vel = []
            for i in range(self.subpixels):
                for j in range(self.subpixels):
                    z_true = self.z_near_func({'R': self.sub_R_true[i][j]}, **self.params['height_near'])
                    
                    if z_mirror: z_true_far = -z_true
                    else: z_true_far = self.z_far_func({'R': self.sub_R_true[i][j]}, **self.params['height_far']) 

                    subpix_grid_true = {'near': [self.sub_x_true[j], self.sub_y_true[i], z_true, self.sub_R_true[i][j], self.sub_phi_true[i][j]], 
                                        'far': [self.sub_x_true[j], self.sub_y_true[i], z_true_far, self.sub_R_true[i][j], self.sub_phi_true[i][j]]}
                    subpix_vel.append(self._compute_prop(subpix_grid_true, [self.velocity_func], [vel_kwargs])[0])

            ang_fac = sin_incl * np.cos(self.phi_true) 
            for i in range(self.subpixels_sq):
                for side in ['near', 'far']: subpix_vel[i][side] *= ang_fac

            props = self._compute_prop(grid_true, prop_funcs[1:], prop_kwargs[1:])
            props.insert(0, subpix_vel)
            
        else: 
            props = self._compute_prop(grid_true, prop_funcs, prop_kwargs)
            if true_kwargs[0]: #Positive vel is positive along z, i.e. pointing to the observer, for that reason imposed a (-) factor to convert to the standard convention: (+) receding  
                ang_fac = sin_incl * np.cos(self.phi_true) 
                props[0]['near'] *= ang_fac 
                props[0]['far'] *= ang_fac
                props[0]['near'] += vel_kwargs['vsys']
                props[0]['far'] += vel_kwargs['vsys']

        #***********************************
        #PROJECT PROPERTIES ON THE SKY PLANE        
        x_pro_dict = {}
        y_pro_dict = {}
        z_pro_dict = {}
        for side in ['near', 'far']:
            xt, yt, zt = grid_true[side][:3]
            x_pro, y_pro, z_pro = self._project_on_skyplane(xt, yt, zt, cos_incl, sin_incl)
            if PA: x_pro, y_pro = self._rotate_sky_plane(x_pro, y_pro, PA)             
            if R_disc is not None: R_grid = griddata((x_pro, y_pro), self.R_true, (self.mesh[0], self.mesh[1]), method='linear')
            x_pro_dict[side] = x_pro
            y_pro_dict[side] = y_pro
            z_pro_dict[side] = z_pro

            if self.subpixels:
                for i in range(self.subpixels_sq): #Subpixels are projected on the same plane where true grid is projected
                    props[0][i][side] = griddata((x_pro, y_pro), props[0][i][side], (self.mesh[0], self.mesh[1]), method='linear') #subpixels velocity
                for prop in props[1:]:
                    prop[side] = griddata((x_pro, y_pro), prop[side], (self.mesh[0], self.mesh[1]), method='linear')
                    if R_disc is not None: prop[side] = np.where(np.logical_and(R_grid<R_disc, R_grid>R_inner), prop[side], np.nan) #Todo: allow for R_in as well
            else:
                for prop in props:
                    if not isinstance(prop[side], numbers.Number): prop[side] = griddata((x_pro, y_pro), prop[side], (self.mesh[0], self.mesh[1]), method='linear')
                    if R_disc is not None: prop[side] = np.where(np.logical_and(R_grid<R_disc, R_grid>R_inner), prop[side], np.nan)
            
        """
        grid_axes_3d = np.meshgrid(self.grid.XYZgrid[0], self.grid.XYZgrid[1], self.grid.XYZgrid[0])
        from functools import reduce  
        x_full = np.array(reduce(np.append, list(x_pro_dict.values())))
        y_full = np.array(reduce(np.append, list(y_pro_dict.values())))
        z_full = np.array(reduce(np.append, list(z_pro_dict.values())))
        z_grid = griddata((x_full, y_full, z_full), z_full, (grid_axes_3d[0], grid_axes_3d[1], grid_axes_3d[2]), method='linear') 
        """
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
       Additional arguments are optional and depend upon the function definition, e.g. Mstar=1.0*sfu.Msun
    """

    def __init__(self, grid):
        self.flags = {'disc': True, 'env': False}
        self.grid = grid
        self._velocity_func = Rosenfeld2d.keplerian
        self._intensity_func = Rosenfeld2d.intensity_powerlaw
        self._linewidth_func = Rosenfeld2d.linewidth_powerlaw
        self._line_profile = General2d.line_profile_v_sigma
        self._use_temperature = False
        self._use_full_channel = False


    def _get_t(self, A, B, C):
        t = []
        for i in range(self.grid.NPoints):
            p = [A, B[i], C[i]]
            t.append(np.sort(np.roots(p)))
        return np.array(t)

    def make_model(self, incl, psi, PA=0.0, int_kwargs={}, vel_kwargs={}, lw_kwargs=None, ls_kwargs=None):
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
        R_true_near = hypot_func(x_true_near, y_true_near) 
        R_true_far = hypot_func(x_true_far, y_true_far)

        z_true_near = t[1] * cos_incl
        z_true_far = t[0] * cos_incl 

        phi_true_near = np.arctan2(y_true_near, x_true_near)        
        phi_true_far = np.arctan2(y_true_far, x_true_far)        

        #****************************
            
        grid_true =  {'near': [x_true_near, y_true_near, z_true_near, R_true_near, phi_true_near], 
                      'far': [x_true_far, y_true_far, z_true_far, R_true_far, phi_true_far]}

        #*******************************
        #COMPUTE PROPERTIES ON TRUE GRID
        avai_kwargs = [vel_kwargs, int_kwargs, lw_kwargs, ls_kwargs]
        avai_funcs = [self.velocity_func, self.intensity_func, self.linewidth_func, self.lineslope_func]
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
