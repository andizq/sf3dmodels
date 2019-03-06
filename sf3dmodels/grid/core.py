from __future__ import print_function
import os
import time
import copy
import inspect
import itertools
import numpy as np
from . import GridInit
from . import GridSet
from ..utils.units import pc
from ..utils.constants import temp_cmb
from ..rt import propTags
from ..tools import formatter

#*************
#OVERLAP GRIDS
#*************

class Overlap(object):
    """
    Host class with functions to overlap submodels either from files or from a list of prop objects, into a unique uniform grid.
    
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure where submodels will be merged in. 
    """
    def __init__(self, GRID):
        self.GRID = GRID
        
    def _get_id_lime(self,i_list,j_list,k_list,ns):
        nx,ny,nz = ns
        return np.array(i_list)*ny*nz + np.array(j_list)*nz + np.array(k_list)

    def _get_id_radmc3d(self,i_list,j_list,k_list,ns):
        nx,ny,nz = ns
        return np.array(k_list)*ny*nx + np.array(j_list)*nx + np.array(i_list) 

    def _mindistance(self,x,xma,Nx):
        """
        Returns index of the nearest value in array xma to scalar x. 
        """
        j = np.where(xma < x)[0][-1]
        if j+1 == Nx: return j
        else: return j if (x-xma[j]) < (xma[j+1]-x) else j+1

    def _get_files_in_folder(self,folder):
        """
        Returns the list of .dat files in folder.
        """
        num=int(np.loadtxt(os.popen("ls -1 %s*.dat| wc -l"%folder)))
        allfiles=os.popen("ls -1 %s*.dat"%folder).read().split('\n',num)[:-1]
        return allfiles

    def fromfiles(self, columns, 
                  submodels = 'all',
                  weighting_dens = None, 
                  min_values = {'dens_H': 1e6,
                                'dens_H2': 1e6,
                                'dens_ion': 1e3,
                                'dens_e': 1e3,
                                'temp_gas': temp_cmb,
                                'temp_dust:': temp_cmb,
                                'gtdratio': 1e2},
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
           Array object containing the column names from input file. 

        submodels : 'all' or array_like, optional 
           If 'all': reads all the '.dat' files within ``folder``.\n
           If array_like: File names to be considered by the algorithm.\n
           Defaults to 'all'.
        
        weighting_dens : str, optional
           Density for weighting the non-density properties. See equations in the **Notes** section.\n
           If None: The algorithm takes the 4th column to weight the non-density properties.
           Defaults to None.
        
        min_values : dict, optional
           Dictionary containing the base minimum values for the final-overlaped physical properties.
                   
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
        The overlaping is computed via a standard average for the densities ( :math:`\\rho`) and a weighted-by-density average for 
        the remaining properties (:math:`X`) as follows:
 
        .. math:: \\rho &= \\frac{1}{N}\\sum^{N}_{i=0} \\rho_i,\n\n
                  X &= \\frac{\\sum^{N}_{i=0} \\rho_{wi} X_i } {\\rho_w},\n
        
        where :math:`N` is the number of submodels involved in the process and :math:`\\rho_w` the weighting density specified via ``weighting_dens``.
        """

        #***************************
        #PREPARING AND READING FILES
        #***************************

        func_name = inspect.stack()[0][3]
        print ("Running function '%s'..."%func_name)

        if weighting_dens is None: weighting_dens = columns[4] #First column after the grid columns
        elif weighting_dens not in columns: raise ValueError("The weighting column '%s' is not amongst the written columns"%weighting_dens, columns)  

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
    
        #***************************
        #DEFINING DICTS 
        #***************************

        GRID = self.GRID
        ntotal = GRID.NPoints
        nrows = [len(data[nf]) for nf in range(nfiles)]
        nx, ny, nz = GRID.Nodes
        xgrid, ygrid, zgrid = GRID.XYZcentres 
        cm3_to_m3 = 1e6
        
        data_dicts = [{columns[i]: data[j][:,i] for i in range(len(columns))} for j in range(nfiles)]
        coords, densities, others = [], [], []
        for col in columns: 
            kind = propTags.get_prop_kind(col)
            if kind == 'grid': coords.append(col) 
            elif kind == 'density': densities.append(col)
            else: others.append(col)
            
        tmp_dict = {}
        val0_dict = {}
        
        for col in densities+others:
            tmp_dict[col] = np.zeros(ntotal)
            val0_dict[col] = 0
        
        tmp_dict[weighting_dens] += -1
        val0_dict[weighting_dens] += -1
        
        partial_dicts = [copy.deepcopy(tmp_dict) for _ in range(nfiles)]

        #***************************
        #FILLING EACH FILE's DICT 
        #***************************
        others_tmp_dicts = [{col: data_dicts[j][col] * data_dicts[j][weighting_dens] for col in others} for j in range(nfiles)]
        num_list = [] 
        num_uniques = []
        if rt_code == 'lime': get_id = self._get_id_lime
        elif rt_code == 'radmc3d': get_id = self._get_id_radmc3d
        else: raise ValueError("The value '%s' in rt_code is invalid. Please choose amongst the following: 'lime', 'radmc3d'"%rt_code)
        for nf in range(nfiles): 
            i_list, j_list, k_list = [], [], []
            xiter = iter(data_dicts[nf]['x'])
            yiter = iter(data_dicts[nf]['y'])
            ziter = iter(data_dicts[nf]['z'])
            for _ in itertools.repeat(None, nrows[nf]): 
                i_list.append(self._mindistance(next(xiter),xgrid,nx))
                j_list.append(self._mindistance(next(yiter),ygrid,ny))
                k_list.append(self._mindistance(next(ziter),zgrid,nz))
            num_list.append(get_id(i_list,j_list,k_list,GRID.Nodes))
            num_uniques.append(np.unique(num_list[nf], return_counts=True))

        for nf in range(nfiles):
            rowiter = iter(np.arange(nrows[nf]))
            num = num_list[nf]
            for _ in itertools.repeat(None, nrows[nf]): 
                row = next(rowiter)
                for col in densities: partial_dicts[nf][col][num[row]] += data_dicts[nf][col][row]
                for col in others: partial_dicts[nf][col][num[row]] += others_tmp_dicts[nf][col][row]
            
            for col in others:
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
        for col in others: final_dict[col] = np.sum([partial_dicts[nf][weighting_dens] * partial_dicts[nf][col] for nf in range(nfiles)], axis=0) / final_dict[weighting_dens]

        #******************************************
        #FILLING DICT WITH min_values and 0's
        #******************************************
        #final_unique = np.unique(reduce(np.append, [num_uniques[nf][0] for nf in range(nfiles)]))
        #mask_empty = np.ones(ntotal, dtype=bool)
        #mask_empty[final_unique] = False
        for col in densities+others:
            if col in min_values: 
                #final_dict[col][mask_empty] = empty_cells[col] 
                final_dict[col] = np.where(final_dict[col] < min_values[col], min_values[col], final_dict[col])
                print ('Using constant minimum value %.3e'%min_values[col], "for column '%s'."%col)
            else: 
                final_dict[col] = np.where(final_dict[col] < 0.0, 0.0, final_dict[col])
                #final_dict[col] = np.where(final_dict[col] < 0., 0., final_dict[col])
                print ("Using constant minimum value 0.0 for column '%s'."%col)

        return final_dict

    def fromprops(self, props):
        """
        In preparation.
        """
        pass

#*************
#GRID BUILDER
#*************

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
       Array of :math:`r` in `input units`.
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

