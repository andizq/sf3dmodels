from __future__ import print_function
import os
import time
import copy
import inspect
import itertools
import numpy as np
from . import GridInit
from . import GridSet
from ..utils.units import au, pc, amu
from ..utils.constants import temp_cmb
from ..utils.prop import propTags
from ..tools import formatter
from .. import Model

#****************
#MAKE GRID (Included: Random)
#****************
class Grid(object):
    def __init__(self):
        self.kind = {'random_weighted': True}
        
    def _accept_point(self, val, norm, power):
        flag = np.random.random()
        val = (val / norm)**power
        if val >= flag: return True
        else: return False

    def _sph_to_cart(self, r,th,phi):
        x = r * np.sin(th) * np.cos(phi)
        y = r * np.sin(th) * np.sin(phi)
        z = r * np.cos(th)
        return x,y,z

    def random(self, function=None, r_size=100*au, normalization=1e16, power=0.5, npoints=50000, kwargs_func={}):
        #Include an option to define a certain r to compute the normalization.
        x,y,z = np.zeros((3,npoints))
        rh,Rh,th,ph = np.zeros((4,npoints))
        n = 0
        twopi = 2*np.pi
        while (n < npoints):
            #i = np.random.randint(r_size)
            r = np.random.uniform(0,r_size)
            phi = np.random.uniform(0,twopi)
            theta = np.random.uniform(0,np.pi)
            R = r * np.sin(theta)
            kwargs_func.update({'loc': {'r': r, 'R': R, 'theta': theta, 'phi': phi, 'z': r*np.cos(theta)}})
            val = function(**kwargs_func)
            if self._accept_point(val,normalization,power): 
                x[n],y[n],z[n] = self._sph_to_cart(r, theta, phi)
                rh[n],Rh[n],th[n],ph[n] = r, R, theta, phi
                n+=1
            else: continue
        GRID = Model.Struct( XYZ = np.array([x,y,z]), NPoints = npoints)
        GRID.rRTP = [rh,Rh,th,ph]
        return GRID

#*************
#OVERLAP GRIDS
#*************
class NeighbourRegularGrid(object):
    """
    Class for getting nearest neighbours in regular grids. 
    """
    def _get_nearest_id_lime(self,i_list,j_list,k_list,ns):
        nx,ny,nz = ns
        return np.array(i_list)*ny*nz + np.array(j_list)*nz + np.array(k_list)

    def _get_nearest_id_radmc3d(self,i_list,j_list,k_list,ns):
        nx,ny,nz = ns
        return np.array(k_list)*ny*nx + np.array(j_list)*nx + np.array(i_list) 

    def _neighbour1d(self,x,xma,Nx):
        """
        Returns the index of the nearest value from the array ``xma`` to the scalar ``x``. 
        """
        j = np.where(xma < x)[0]
        if len(j) == 0: j = 0 #left-border
        else: j = j[-1] #closest from the left
        if j+1 == Nx: return j #right-border
        else: return j if (x-xma[j]) < (xma[j+1]-x) else j+1 #compare to the closest from the right

class Overlap(NeighbourRegularGrid):
    """
    Host class with functions to overlap submodels either from files or from prop objects into a single regular grid.
    
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure where submodels will be merged in.

    Attributes
    ----------
    min_values : dict, optional
       Dictionary containing the base minimum values for the final-overlaped physical properties.
    """
    def __init__(self, GRID):
        self.GRID = GRID
        self.min_values = {'dens_H': 1e3,
                           'dens_H2': 1e3,
                           'dens_Hplus': 1e3,
                           'dens_ion': 1e3,
                           'dens_e': 1e3,
                           'temp_gas': temp_cmb,
                           'temp_dust': temp_cmb,
                           'gtdratio': 1e2,
                           }

    def _get_files_in_folder(self,folder):
        """
        Returns the list of .dat files in folder.
        """
        num=int(np.loadtxt(os.popen("ls -1 %s*.dat| wc -l"%folder)))
        allfiles=os.popen("ls -1 %s*.dat"%folder).read().split('\n',num)[:-1]
        return allfiles

    def fromfiles(self, columns, 
                  submodels = 'all',
                  weighting_dens = 'all', 
                  rt_code = 'lime',
                  folder = './Subgrids'): 
                  #weighting_dens = {'Lime': ['dens_H', 'dens_H2'],
                  #                  'Radmc3d': ['dens_ion']}):
        """
        Overlaps submodels from files. 
        The input files structure is the same as that of the files written by `~sf3dmodels.rt.MakeDatatab.submodel`.
        
        Parameters
        ----------
        columns : array_like, shape (ncols,)
           List object containing the column names of the input submodel files. 

        submodels : 'all' or array_like, optional 
           If 'all': reads all the '.dat' files within ``folder``.\n
           If array_like: File names to be considered by the algorithm.\n
           Defaults to 'all'.
        
        weighting_dens : str, optional
           Density column name for weighting the non-density properties. See equations in the **Notes** section.\n
           If 'all': The algorithm takes the sum of all the density columns multiplied by their respective atomic mass.\n
           Defaults to 'all'.
                   
        rt_code : 'lime' or 'radmc3d', optional
           Radiative transfer code that is going to be used later with the output prop.

        folder : str, optional
           Folder name were the submodel files are located. Defaults to './Subgrids'.

        Returns
        -------
        final_dict : dict
           Dictionary containing the overlaped properties. The dictionary keys are the physical properties from the input ``columns``.

        Notes
        -----
        The overlaping is computed via a direct addition for the densities ( :math:`\\rho`) and a weighted-by-density average for 
        the remaining properties (:math:`X`) as follows:
 
        .. math:: \\rho &= \\sum^{N}_{i=0} \\rho_i,\n\n
                  X &= \\frac{\\sum^{N}_{i=0} \\rho_{wi} X_i } {\\rho_w},\n
        
        where :math:`N` is the number of submodels involved in the process and :math:`\\rho_w` the weighting density specified via ``weighting_dens``.
        """

        #***************************
        #PREPARING AND READING FILES
        #***************************
        func_name = inspect.stack()[0][3]
        print ("Running function '%s'..."%func_name)

        if folder[-1] != '/': folder += '/'
        allfiles = self._get_files_in_folder(folder)
        if submodels == 'all': files = allfiles
        elif isinstance(submodels, list) or isinstance(submodels, np.ndarray): files = [folder + sub for sub in submodels]
        else: raise TypeError("Invalid type: %s for 'submodels'. Please provide a valid 'submodels' object: list, np.ndarray or str 'all'"%type(submodels))
        
        #If wish to preserve data types and formats use np.genfromtxt, but it's slower and the output is a 1-D array of tuples, where each tuple is a row from the table. 
        data = [np.loadtxt(file, dtype=None) for file in files] 
        detected = [file.split(folder)[1] for file in allfiles]
        read = [file.split(folder)[1] for file in files]
        
        nfiles = len(files)
        print ('Files detected (%d):'%len(allfiles), detected, 
               '\nFiles to merge in grid (%d):'%nfiles, read)

        data_dicts = [{columns[i]: data[j][:,i] for i in range(len(columns))} for j in range(nfiles)]

        if weighting_dens == 'all': 
            weighting_dens = 'dens_mass'
            columns = np.append(columns,weighting_dens)
        elif weighting_dens not in columns: raise ValueError("The weighting column '%s' is not amongst the written columns"%weighting_dens[i], columns)  
                
        #***************************
        #DEFINING DICTS 
        #***************************        
        GRID = self.GRID
        ntotal = GRID.NPoints
        nrows = [len(data[nf]) for nf in range(nfiles)]
        nx, ny, nz = GRID.Nodes
        xgrid, ygrid, zgrid = GRID.XYZcentres 
        cm3_to_m3 = 1e6
        
        coords, densities, velocities, others = [], [], [], []
        for col in columns: 
            kind = propTags.get_prop_kind(col)
            if kind == 'grid': coords.append(col) 
            elif kind == 'density': densities.append(col)
            elif kind == 'velocity': velocities.append(col)
            else: others.append(col)

        tmp_dict = {}
        val0_dict = {}
        
        for col in densities+others:
            tmp_dict[col] = np.zeros(ntotal)
            val0_dict[col] = 0            
        for col in velocities:
            tmp_dict[col] = np.zeros(ntotal) #np.random.normal(scale=10*amu,size=ntotal)
            val0_dict[col] = 0

        if weighting_dens == 'dens_mass': 
            for nf in range(nfiles):
                data_dicts[nf][weighting_dens] = np.zeros(nrows[nf])
                for col in densities:
                    if col != weighting_dens: data_dicts[nf][weighting_dens] += data_dicts[nf][col] * propTags.get_dens_mass(col)
    
            tmp_dict[weighting_dens] += -1*amu
            val0_dict[weighting_dens] += -1*amu
        
        else:
            tmp_dict[weighting_dens] += -1
            val0_dict[weighting_dens] += -1
        
        partial_dicts = [copy.deepcopy(tmp_dict) for _ in range(nfiles)]

        #***************************
        #FILLING EACH FILE's DICT 
        #***************************
        others_tmp_dicts = [{col: data_dicts[j][col] * data_dicts[j][weighting_dens] for col in velocities+others} for j in range(nfiles)]
        num_list = [] 
        num_uniques = []
        if rt_code == 'lime': get_id = self._get_nearest_id_lime
        elif rt_code == 'radmc3d': get_id = self._get_nearest_id_radmc3d
        else: raise ValueError("The value '%s' in rt_code is invalid. Please choose amongst the following: 'lime', 'radmc3d'"%rt_code)
        for nf in range(nfiles): 
            i_list, j_list, k_list = [], [], []
            xiter = iter(data_dicts[nf]['x'])
            yiter = iter(data_dicts[nf]['y'])
            ziter = iter(data_dicts[nf]['z'])
            for _ in itertools.repeat(None, nrows[nf]): 
                i_list.append(self._neighbour1d(next(xiter),xgrid,nx))
                j_list.append(self._neighbour1d(next(yiter),ygrid,ny))
                k_list.append(self._neighbour1d(next(ziter),zgrid,nz))
            num_list.append(get_id(i_list,j_list,k_list,GRID.Nodes))
            num_uniques.append(np.unique(num_list[nf], return_counts=True))

        for nf in range(nfiles):
            rowiter = iter(np.arange(nrows[nf]))
            num = num_list[nf]
            for _ in itertools.repeat(None, nrows[nf]): 
                row = next(rowiter)
                for col in densities: partial_dicts[nf][col][num[row]] += data_dicts[nf][col][row]
                for col in velocities+others: partial_dicts[nf][col][num[row]] += others_tmp_dicts[nf][col][row]
            
            for col in velocities+others:
                partial_dicts[nf][col][num_uniques[nf][0]] /= partial_dicts[nf][weighting_dens][num_uniques[nf][0]]         
            for col in densities: 
                #partial_dicts[nf][col][num_uniques[nf][0]] -= val0_dict[col] #Commented to avoid nans in the final division on final_dict
                partial_dicts[nf][col][num_uniques[nf][0]] /= num_uniques[nf][1] 

            print ('Finished merging for: %s'%files[nf])
        
        #***************************
        #FILLING FINAL GLOBAL DICT
        #***************************
        print ('Computing combined physical properties...')
        
        final_dict = {}        
        for col in densities: final_dict[col] = np.sum([partial_dicts[nf][col] for nf in range(nfiles)], axis=0)
        for col in velocities+others: final_dict[col] = np.sum([partial_dicts[nf][weighting_dens] * partial_dicts[nf][col] for nf in range(nfiles)], axis=0) / final_dict[weighting_dens]

        #******************************************
        #FILLING DICT WITH min_values and 0's
        #******************************************
        #final_unique = np.unique(reduce(np.append, [num_uniques[nf][0] for nf in range(nfiles)]))
        #mask_empty = np.ones(ntotal, dtype=bool)
        #mask_empty[final_unique] = False

        for col in densities+others:
            if col in self.min_values: 
                #final_dict[col][mask_empty] = empty_cells[col] 
                final_dict[col] = np.where(final_dict[col] < self.min_values[col], self.min_values[col], final_dict[col])
                print ('Using constant minimum value %.3e'%self.min_values[col], "for column '%s'."%col)
            else: 
                final_dict[col] = np.where(final_dict[col] < 0.0, 0.0, final_dict[col])
                #final_dict[col] = np.where(final_dict[col] < 0., 0., final_dict[col])
                print ("Using constant minimum value 0.0 for column '%s'."%col)

        if weighting_dens == 'dens_mass': _ = final_dict.pop(weighting_dens)
 
        return final_dict
    
    def fromprops(self, props):
        """
        In preparation.
        """
        pass

#****************
#GRID AROUND AXIS
#****************
class RandomGridAroundAxis(object):
    """
    Base class for Random grids around a given axis.
    """
    def __init__(self, pos_c, axis, z_min, z_max, dx, mirror=False):
        self.pos_c = np.asarray(pos_c)
        self.axis = np.asarray(axis)
        self.pos_f = self.pos_c + self.axis
        self.z_min = z_min
        self.z_max = z_max
        self.dx = dx        
        self.mirror = mirror
        self._axis()

    def _axis(self):
        self.z_dir = self.axis/np.linalg.norm(self.axis) #Axis direction
        z_seg = (self.z_max - self.z_min) * self.z_dir #The segment length is zmax-zmin  
        self.dr = self.dx * 2.**-1 #3.**-1#2.5**-1

        self._cart_th = np.zeros(3) #Initializing an empty cartesian vector.
        self._cart_th[np.argmin(self.axis)] = 1. #Make 1 the component where the axis vector is shorter. 
        self._axis_th = np.cross(self.z_dir, self._cart_th) #The reference axis for theta is the cross product of axis and cart.

    def _grid(self, func_width, width_pars, R_min=None, dummy_frac=0.0):
        
        #Guess the number of random grid points to generate: (approximation to a rectangular region)
        # Number of divisions along the main segment, 
        #  times the number of divisions along an axis perpendicular to the segment,
        #   times the number of divisions along an axis perpendicular to both of the segments above.
        #    times an exageration factor (in dr) to increase the number of grid points.
        #     Twice to account for the other branch of the jet.

        mean_w = np.mean(func_width(np.linspace(self.z_min,self.z_max,num=100), *width_pars))
        z_seg_mag =  self.z_max - self.z_min #Long-axis length
        
        self.NPoints = int(z_seg_mag/self.dr * (mean_w/self.dr)**2 ) 

        mirror_int = 1
        if self.mirror: mirror_int = 2
        
        npoints = self.NPoints
        print ('Number of grid points: %d'%(mirror_int*npoints))
        
        z = np.random.uniform(self.z_min, self.z_max, size=npoints) #Random z's along long axis
        z_vec = z[:,None]*self.z_dir #Random vectors along the long axis
        width = func_width(z, *width_pars)
        
        rand_vec = np.random.uniform(-1,1,size=(npoints,3)) #Random vector to generate the perpendicular vector to the long axis
        
        if R_min is None: R_min=0 
        R = np.random.uniform(R_min, width) #Random R from the long axis
        R_dir = np.cross(self.z_dir, rand_vec)        
        R_dir /= np.linalg.norm(R_dir, axis=1, keepdims=True) #Perpendicular (random) unit vector to the long axis
        R_vec = R[:,None]*R_dir

        theta = np.sign(np.dot(np.cross(self._axis_th,R_dir), self.z_dir))*np.arccos(np.dot(self._axis_th,R_dir.T)) 
        theta_dir = np.cross(self.z_dir,R_dir)
        theta_dir /= np.linalg.norm(theta_dir, axis=1, keepdims=True)

        r_vec = z_vec + R_vec #Vectors from the object centre to the generated points
        r_dir = r_vec / np.linalg.norm(r_vec, axis = 1)[np.newaxis].T

        r_real = self.pos_c + r_vec #Positions from the origin of coordinates

        self.ndummies = int(round(npoints*dummy_frac))        
        if self.ndummies > 0:
            rand_ind = np.random.choice(np.arange(npoints), size=self.ndummies, replace=False)
            R_dummy = np.random.uniform(width[rand_ind], np.max(width))
            R_dummy_vec = R_dummy[:,None]*R_dir[rand_ind]
            r_dummy_vec = z_vec[rand_ind] + R_dummy_vec
            r_dummy_real = self.pos_c + r_dummy_vec
            print ('Number of dummy points: %d\nNew number of grid points: %d'
                   %(self.ndummies, self.ndummies+mirror_int*npoints))
            
        if self.mirror: 
            r_real_n = self.pos_c - r_vec #Mirror point to real position from the origin of coordinates
            self.grid = np.append(r_real, r_real_n, axis=0)
            if self.ndummies > 0:
                r_dummy_real_n = self.pos_c - r_dummy_vec 
                self.grid_dummy = np.append(r_dummy_real, r_dummy_real_n, axis=0)
            self.r_dir = np.append(r_dir, -1*r_dir, axis=0) 
            self.width = np.append(width, width)
            self.z = np.append(z, -1*z)
            self.z_dir = np.repeat([self.z_dir], 2*npoints, axis=0)
            self.R = np.append(R, R)
            self.R_dir = np.append(R_dir, -1*R_dir, axis=0)
            self.theta = np.append(theta, theta-np.sign(theta)*np.pi)
            self.theta_dir = np.append(theta_dir, -1*theta_dir, axis=0)
        else: 
            self.grid = r_real
            if self.ndummies > 0: self.grid_dummy = r_dummy_real
            self.r_dir = r_dir
            self.width = width
            self.z = z
            self.z_dir = np.repeat([self.z_dir], npoints, axis=0)
            self.R = R
            self.R_dir = R_dir
            self.theta = theta
            self.theta_dir = theta_dir
        
        self.NPoints = mirror_int*(npoints+self.ndummies)

#************
#GRID BUILDER 
#************
class Build_r(GridSet):
    _base = """
    Computes the spherical coordinate :math:`r` on each cartesian grid point.

       .. math:: r = \\sqrt{x^2 + y^2 + z^2}
    """
    _pars = GridSet._pars
    _returns = """
    Attributes
    ----------
    r : `numpy.ndarray`, shape ``(GRID.NPoints,)``
       Array of :math:`r` coordinate values in `input units`.
    """
    __doc__ = _base + _pars + _returns

    def __init__(self, GRID):
        super(Build_r, self).__init__(GRID)
        self._r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self._r_max = np.max(self._r_grid)
        self.GRID.r = self._r_grid
        self._set_flag('r')

class Build_theta(GridSet):
    _base = """
    Computes the polar angle :math:`\\theta` on each cartesian grid point.
       .. math:: 
          \\theta = \\cos^{-1}\\left(\\frac{z}{r}\\right) 
    """
    _pars = GridSet._pars
    _returns = """
    Attributes
    ----------    
    theta : `numpy.ndarray`, shape ``(GRID.NPoints,)`` 
       Array of angles in `radians`, in the range `[0, pi/2]` if :math:`z>=0`, 
       and `[pi/2, 0]` if :math:`z<0`; 
       where :math:`\\theta=` `pi/2` on the plane :math:`z=0`.
    """
    __doc__ = _base + _pars + _returns

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)
        
class Build_phi(GridSet):
    _base = """
    Computes the azimuthal angle :math:`\\phi` on each cartesian grid point. Uses `numpy.arctan2`.

       .. math:: \\phi = \\tan^{-1}\\left(\\frac{y}{x}\\right) + 2\\pi
    """
    _pars = GridSet._pars
    _returns = """
    Attributes
    ----------
    phi : `numpy.ndarray`, shape ``(GRID.NPoints,)``
       Array of angles in `radians`, in the range `[0, pi/2]`.
    """
    __doc__ = _base + _pars + _returns

    def __init__(self, GRID):
        self.GRID = GRID
        self.r_grid = np.linalg.norm(self.GRID.XYZ, axis = 0)
        self.r_max = np.max(self.r_grid)

