from __future__ import print_function
import os
import inspect
import itertools
import numpy as np
from ..utils.units import cm, amu
from ..utils.prop import propTags
from ..tools import formatter

"""
class Emissivity(object):
    def __init__(self, GRID, prop):
"""     

#******************************************
#WRITING DATA (LIME v1.9.5, RADMC-3D v0.41)
#******************************************
class MakeDatatab(object):
    """
    Base class for rt (radiative transfer) classes. 
    """
    def __init__(self, GRID):
        self.GRID = GRID
        super(MakeDatatab, self).__init__()

    @classmethod
    def _get_class_name(cls): 
        return cls.__name__
    @classmethod
    def _get_available_props(cls, self):
        """
        Returns the sorted available physical properties which can be read by the RT codes.
        """
        parents = np.array(cls.__bases__)
        if cls == Lime or (parents == Lime).any(): 
            av_props = ['dens_H2', 'dens_p_H2', 'dens_o_H2', 'dens_e', 'dens_H', 'dens_He', 'dens_Hplus', 
                        'temp_gas', 'temp_dust',
                        'vel_x', 'vel_y', 'vel_z',
                        'abundance',
                        'gtdratio',
                        'doppler'
                        ]   
            #cls._col_ids(self)
            #self._propnames = np.array(list(self.sf3d_header))[np.argsort(list(self.sf3d_header.values()))]
            #print (self._propnames, self.sf3d_header)
 
        if cls == Radmc3d or (parents == Radmc3d).any(): 
            av_props = ['dens_ion', 'dens_e', 'dens_dust',
                        'temp_gas', 'temp_dust', 
                        'vel_x', 'vel_y', 'vel_z',
                        'microturbulence'
                        ]
        #if cls == Polaris or (parents == Polaris).any(): pass
        return av_props
        
    def _write_npoints_header(self, folder='./', lime_npoints = True, lime_header = True):
        """
        Writes the npoints.dat and header.dat files for Lime.
        """
        if folder[-1] != '/': folder += '/'

        if lime_npoints:
            try: Ns = self.GRID.Nodes
            except AttributeError: Ns = [0,0,0]
            size_file = folder+'npoints.dat'
            sfile = open(size_file,'w') 
            print ('Writing grid size into %s'%size_file)
            sfile.write("%d %d %d %d %d"%(len(self.columns), Ns[0], Ns[1], Ns[2], self.GRID.NPoints))
            sfile.close()
    
        if lime_header:
            header_file = folder+'header.dat'
            print ('Writing columns header into %s'%header_file)
            colswritten = np.insert(np.array([0,4242]), 1, self.prop_id).T
            np.savetxt(header_file, colswritten, fmt = '%d')

        
    def _prepare_prop(self, prop):
        """
        Prepare the prop object to write its content in ordered columns according to _col_ids().
        """
        prop_id, prop_list, prop_keys, abund_id, abund_list, abund_keys = [], [], [], [], [], []

        for key in prop: 
            if key in self.sf3d_header: 
                prop_id.append(self.sf3d_header[key])            
                prop_list.append(prop[key])
                prop_keys.append(key)
            elif ('abundance' in key) and ('abundance' in self.sf3d_header):
                abund_id.append( int(key.split('abundance')[-1][-1]) ) 
                abund_list.append(prop[key])
                abund_keys.append(key)            
            else: raise KeyError("The property '%s' is invalid for the class %s. Please make sure your property is amongst the following keys:"
                                 %(key,self._get_class_name()), self._get_available_props(self))
        
        arg_sorted = np.argsort(prop_id)        
        prop_id_sorted = np.array(prop_id)[arg_sorted]
        prop_list_sorted = np.array(prop_list)[arg_sorted]
        prop_keys_sorted = np.array(prop_keys)[arg_sorted].tolist()

        if len(abund_id) != 0:
            abund_sorted = np.argsort(abund_id)
            abund_list_sorted = np.array(abund_list)[abund_sorted]
            abund_keys_sorted = np.array(abund_keys)[abund_sorted]
            abund_id2insert = len(prop_id_sorted[prop_id_sorted <= self.sf3d_header['abundance']])
            prop_id_sorted = np.insert(prop_id_sorted, abund_id2insert, [self.sf3d_header['abundance']] * len(abund_id))
            prop_list_sorted = np.insert(prop_list_sorted, abund_id2insert, abund_list_sorted, axis=0)
            prop_keys_sorted = np.array(prop_keys_sorted[0:abund_id2insert] + abund_keys_sorted.tolist() + prop_keys_sorted[abund_id2insert:]) #For memory reasons, the np.insert was cutting off longer strings (than the longest property before abundances, e.g. temp_gas [8 characters]) such like abundance0 and abundance1 to be abundanc.

        prop_list_sorted = prop_list_sorted.tolist()
        if len(np.shape(prop_list_sorted)) == 1: self.prop_list = [prop_list_sorted]
        else: self.prop_list = prop_list_sorted
        self.n = len(prop_keys_sorted)
        self.prop_id = np.array(prop_id_sorted)
        self.prop_keys = np.array(prop_keys_sorted)
        self.prop_header = {prop_keys_sorted[i]: prop_id_sorted[i] for i in range(self.n)}
        self.prop = prop
        
    def submodel(self, prop, output = '0.dat', fmt = '%.6e', folder = './Subgrids', lime_npoints=False, lime_header=False):        
        """
        Writes a preliminary model. 
        This method **must** be used either from the radiative transfer class `~sf3dmodels.rt.Lime` or `~sf3dmodels.rt.Radmc3d`. 
        
        Motivation: the written submodel will tipically suffer further modifications from an external code or will be read back
        by sf3dmodels and merged in a global grid with other submodels. See `~sf3dmodels.grid.Overlap`. 
        
        Parameters
        ----------
        prop : dict
           Dictionary containing the model physical properties.\n
           See the **Notes** section of `~sf3dmodels.rt.Lime` or `~sf3dmodels.rt.Radmc3d` for available physical properties.
           
        output : str, optional
           Submodel output file. Defaults to '0.dat'
        
        fmt : str or dict, optional
            If str: Format string for all the columns. Defaults to '%.6e'.\n
            If dict: the input keys of the dictionary must be present in the ``prop`` keys. The keys in ``prop`` which
            were not specified in the ``fmt`` dict, will be written with the default format '%.6e'.

        folder : str, optional
           Sets the folder to write the files in. Defaults to './Subgrids'.
           
        Attributes
        ----------
        prop_header : dict
           Dictionary relating the prop keys to their characteristic id according to the 
           invoked radiative transfer class.

        prop_keys : array_like
           Array object containing the column names of the physical properties written into file.

        columns : array_like
           Array object containing all the column names written into file. 
           The prop columns are sorted based on their characteristic id from the radiative transfer class.
        
        Returns
        -------
        Data file with sorted columns according to the id's table at the **Notes** section 
        of `~sf3dmodels.rt.Lime` or `~sf3dmodels.rt.Radmc3d`. 
        """

        
        os.system('mkdir %s'%folder)
        if folder[-1] != '/': folder += '/'
        file_path = '%s%s'%(folder,output)
        x,y,z = self.GRID.XYZ
        self.id = np.arange(self.GRID.NPoints)

        xyz_dict = {'x': x, 'y': y, 'z': z}
        prop.update(xyz_dict)
        
        self._prepare_prop(prop)
        #self.columns = np.append(['id','x','y','z'], self.prop_keys)
        self.columns = np.append(['id'], self.prop_keys)

        #fmt_string = formatter(self.prop_keys, fmt, base = '%d %.8e %.8e %.8e')
        fmt_string = formatter(self.prop_keys, fmt, base = '%d')
        tmp_write = []
        #list2write = iter(np.array([self.id,x,y,z] + self.prop_list).T)
        list2write = iter(np.array([self.id] + self.prop_list).T)
        print ('Writing Submodel data in %s'%file_path)
        for _ in itertools.repeat(None, self.GRID.NPoints): tmp_write.append( fmt_string % tuple(next(list2write)) )
        file_data = open(file_path, 'w')        
        file_data.writelines(tmp_write)
        
        #self.prop_id = np.insert(self.prop_id, 0, [self.sf3d_header[coord] for coord in ['x','y','z']]) 

        self._write_npoints_header(folder=folder, lime_npoints=lime_npoints, lime_header=lime_header)
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')
                

class Lime(MakeDatatab):
    """
    Base class for lime related objects.
    Contains base functions to write formatted files to perform Radiative Transfer calculations with `LIME`_. 
    
    - Input and output units: SI (metres, kilograms, seconds).

    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure in which the model was computed.

    Notes
    -----
    Available grid [0,1,2,3] and physical properties for Lime, their id and description:

    +-----------+----------+---------------------------------+
    | prop key  |    id    |Description                      |
    +===========+==========+=================================+
    | id        | 0        | cell id                         |
    +-----------+----------+---------------------------------+
    | x         | 1        | x coordinate                    |
    +-----------+----------+---------------------------------+
    | y         | 2        | y coordinate                    |
    +-----------+----------+---------------------------------+
    | z         | 3        | z coordinate                    |
    +-----------+----------+---------------------------------+
    | dens_H2   | 4        |Molecular Hydrogen density       |
    +-----------+----------+---------------------------------+
    | dens_p_H2 | 5        |ParaHydrogen density             |
    +-----------+----------+---------------------------------+
    | dens_o_H2 | 6        |OrthoHydrogen density            |
    +-----------+----------+---------------------------------+
    | dens_e    | 7        |Electronic density               |
    +-----------+----------+---------------------------------+
    | dens_H    | 8        |Atomic Hydrogen density          |
    +-----------+----------+---------------------------------+
    | dens_He   | 9        |Helium density                   |
    +-----------+----------+---------------------------------+
    |dens_Hplus | 10       |Ionized Hydrogen density         |
    +-----------+----------+---------------------------------+
    |temp_gas   | 11       |Gas temperature                  |
    +-----------+----------+---------------------------------+
    |temp_dust  | 12       |Dust temperature                 |
    +-----------+----------+---------------------------------+
    |vel_x      | 13       |Velocity along x-axis            |
    +-----------+----------+---------------------------------+
    |vel_y      | 14       |Velocity along y-axis            |
    +-----------+----------+---------------------------------+
    |vel_z      | 15       |Velocity along z-axis            |
    +-----------+----------+---------------------------------+
    |abundance_0|   16     |                                 |
    +-----------+          |                                 |
    |abundance_1|          |                                 |
    +-----------+          |                                 |
    |    ...    |          |Molecular abundance(s)           |
    +-----------+          |                                 |
    |abundance_n|          |                                 |
    +-----------+----------+---------------------------------+
    |gtdratio   | 17       |Gas-to-dust (mass) ratio         |
    +-----------+----------+---------------------------------+
    |doppler    | 18       |Doppler broadening via turbulence|
    +-----------+----------+---------------------------------+
    |           | 4242     |End of file in header.dat        |
    +-----------+----------+---------------------------------+
    """

    def __init__(self, GRID):
        print ('Set LIME format...')
        self._col_ids()
        super(Lime, self).__init__(GRID)
        
    def _col_ids(self):
        
        CP_H2 = 1
        CP_p_H2 = 2
        CP_o_H2 = 3
        CP_e = 4
        CP_H = 5
        CP_He = 6
        CP_Hplus = 7

        base = dict(SF3D_id =            0,
                    SF3D_x =             1,
                    SF3D_y =             2,
                    SF3D_z =             3,
                    SF3D_dens_H2 =       CP_H2 + 3,     
                    SF3D_dens_p_H2 =     CP_p_H2 + 3,
                    SF3D_dens_o_H2 =     CP_o_H2 + 3,
                    SF3D_dens_e =        CP_e + 3,
                    SF3D_dens_H =        CP_H + 3,
                    SF3D_dens_He =       CP_He + 3,
                    SF3D_dens_Hplus =    CP_Hplus + 3,
                    SF3D_temp_gas =      11,
                    SF3D_temp_dust =     12,
                    SF3D_vel_x =         13,
                    SF3D_vel_y =         14,
                    SF3D_vel_z =         15,
                    SF3D_abundance =     16,
                    SF3D_gtdratio =      17,
                    SF3D_doppler =       18,
                    SF3D_max_cols =      19)

        self.sf3d_header = {key.split('_',1)[1]: base[key] for key in base}
        
    def finalmodel(self, prop, fmt = '%.6e', folder = './'):
        """
        Writes the final model into file. The output is ready to be read with `LIME`_ 
        
        Parameters
        ----------
        prop : dict
           Dictionary containing the model physical properties.\n
           See the **Notes** section above for a list with the available input properties.
           
        fmt : str or dict, optional
            If str: Format string for all the columns. Defaults to '%.6e'.\n
            If dict: the input keys of the dictionary must be present in the ``prop`` keys. The keys in ``prop`` which
            were not specified in the ``fmt`` dict, will be written with the default format '%.6e'.

        folder : str, optional
           Sets the folder to write the files in. Defaults to './'.
           
        Attributes
        ----------
        prop_header : dict
           Dictionary relating the input prop keys to their characteristic id according to the Lime requirements.

        prop_keys : array_like
           Array object containing the column names of the physical properties written into file.

        columns : array_like
           Array object containing all the column names written into file. 
           The prop columns are sorted based on their characteristic id from the Lime class.


        Examples
        --------
        Have a look at the following folders on the GitHub repository: 
           `examples/single_main/ <https://github.com/andizq/star-forming-regions/tree/master/examples/single_main>`_ \n
           `examples/hamburger_tapering/ <https://github.com/andizq/star-forming-regions/tree/master/examples/hamburger_tapering>`_ \n
           `examples/two_sources/ <https://github.com/andizq/star-forming-regions/tree/master/examples/two_sources>`_ \n
                      
        .. image:: ../../images/pIntensitypoints.png
           :width: 180px
           :height: 170px
           :align: left

        .. image:: ../../images/global_grid_temp.png
           :width: 260px
           :height: 210px

        .. image:: ../../images/lime_single_source.png
           :width: 180px
           :height: 180px

        
        Returns
        -------
        'datatab.dat' : file
           Data file made from the ``prop``'s content. The columns are sorted according to the id's table in the **Notes** section of the `Lime` class. 
        'npoints.dat' : file
           File containing the number of columns in 'datatab.dat', the number of cells along x,y,z, and the total number of cells.
        'header.dat' : file
           File specifying the column ids according to the table in the **Notes** section of the `Lime` class.
        'x.dat' : file
        'y.dat' : file
        'z.dat' : file     
        """

        self._prepare_prop(prop)
        self.id = np.arange(self.GRID.NPoints)
        self.columns = np.append(['id'], self.prop_keys)

        if folder[-1] != '/': folder += '/'

        fmt_string = formatter(self.prop_keys, fmt, base = '%d')
        tmp_write = []
        #if type(self.prop_list) == np.ndarray: self.prop_list = self.prop_list.tolist()
        list2write = iter(np.array([self.id] + self.prop_list).T)
        files = [folder + tag for tag in ['datatab.dat','x.dat','y.dat','z.dat']]
        
        print ('Writing Global grid data into %s'%files[0])        
        for _ in itertools.repeat(None, self.GRID.NPoints): tmp_write.append( fmt_string % tuple(next(list2write)) )
        file_data = open(files[0],'w')
        file_data.writelines(tmp_write)

        Ns = self.GRID.Nodes
        size_file = folder+'npoints.dat'
        sfile = open(size_file,'w') 
        print ('Writing grid size into %s'%size_file)
        sfile.write("%d %d %d %d %d"%(len(self.columns), Ns[0], Ns[1], Ns[2], self.GRID.NPoints))
        sfile.close()

        header_file = folder+'header.dat'
        print ('Writing columns header into %s'%header_file)
        colswritten = np.insert(np.array([0,4242]), 1, self.prop_id).T
        np.savetxt(header_file, colswritten, fmt = '%d')

        print ('Writing grid coordinates into:')
        for i in range(1,4):
            print ('%s'%files[i])
            np.savetxt(files[i], self.GRID.XYZcentres[i-1], fmt = '%.8e')

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')


#*****************************
#WRITING DATA (RADMC-3D v0.41)
#*****************************

class Radmc3d(MakeDatatab): #RADMC-3D uses the cgs units system
    """
    Base class for radmc3d related objects. Contains base functions to write formatted files to 
    perform Radiative Transfer calculations with `RADMC-3D`_.

    - Input units: SI (metres, kilograms, seconds).\n
    - Output units: cgs (centimetres, grams, seconds).

    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure in which the model was computed.
       
    nphot : int, optional
       Number of photon packages when computing dust temperature via Monte Carlo.

    Notes
    -----
    Available grid [0,1,2,3] and physical properties for Radmc3d, their ids and description:

    +---------------+----------+--------------------------+
    | prop key      |    id    |Description               |
    +===============+==========+==========================+
    | id            | 0        | cell id                  |
    +---------------+----------+--------------------------+
    | x             | 1        | x coordinate             |
    +---------------+----------+--------------------------+
    | y             | 2        | y coordinate             |
    +---------------+----------+--------------------------+
    | z             | 3        | z coordinate             |
    +---------------+----------+--------------------------+
    | dens_ion      | 4        |Ionized Hydrogen density  |
    +---------------+----------+--------------------------+
    | dens_e        | 5        |Electronic density        |
    +---------------+----------+--------------------------+
    | dens_dust     | 6        |Dust density              |
    +---------------+----------+--------------------------+
    |temp_gas       | 7        |Gas temperature           |
    +---------------+----------+--------------------------+
    |temp_dust      | 8        |Dust temperature          |
    +---------------+----------+--------------------------+
    |vel_x          | 9        |Velocity along x-axis     |
    +---------------+----------+--------------------------+
    |vel_y          | 10       |Velocity along y-axis     |
    +---------------+----------+--------------------------+
    |vel_z          | 11       |Velocity along z-axis     |
    +---------------+----------+--------------------------+
    |microturbulence| 12       |Microturbulent velocity   |
    +---------------+----------+--------------------------+
    """    

    def __init__(self, GRID, nphot = 1000000):
        print ('Set RADMC-3D format...')
        self.GRID = GRID
        self.nphot = nphot
        self._col_ids()
        super(Radmc3d, self).__init__(GRID)

    def _col_ids(self):
        
        base = dict(id =                   0,
                    x =                    1,
                    y =                    2,
                    z =                    3,
                    dens_ion =             4,
                    dens_e =               5,
                    dens_dust =            6,
                    temp_gas =             7,
                    temp_dust =            8,
                    vel_x =                9,
                    vel_y =                10,
                    vel_z =                11,
                    microturbulence =      12,
                    max_cols =             13)

        self.sf3d_header = base
        
    def write_amr_grid(self, iformat = 1, 
                       grid_style = 0, 
                       coord_system = 0,
                       grid_info = 0, 
                       include_dim = [1,1,1]):
        """
        Writes the file 'amr_grid.inp' for radmc3d. 
        
        Parameters
        ----------
        iformat : int, optional
           Defaults to 1.
        grid_style : int, optional
           Defaults to 0.
        coord_system : int, optional
           Defaults to 0.
        grid_info : int, optional
           Defaults to 0.
        include_dim : array_like, optional
           Defaults to [1,1,1]
        """
        nx,ny,nz = self.GRID.Nodes
        self.nn = nx * ny * nz
        xi, yi, zi = np.array(self.GRID.XYZgrid) * cm #from m to cm

        #if not self.amr_grid:
        with open('amr_grid.inp','w+') as f:
            f.write('%d\n'%iformat)                   # iformat
            f.write('%d\n'%grid_style)                # AMR grid style  (0=regular grid, no AMR)
            f.write('%d\n'%coord_system)              # Coordinate system
            f.write('%d\n'%grid_info)                 # grid_info
            f.write('%d %d %d\n'%tuple(include_dim))  # Include x,y,z coordinate
            f.write('%d %d %d\n'%(nx,ny,nz))          # Size of grid

            tmp = ['%13.6e '*(n+1) for n in [nx,ny,nz]]
            f.write((tmp[0]+'\n')%tuple(xi)) # X values (cell walls)
            f.write((tmp[1]+'\n')%tuple(yi)) # Y values (cell walls)
            f.write((tmp[2]+'\n')%tuple(zi)) # Z values (cell walls)
            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')
        
    def write_electron_numdens(self, dens_e, fmt = '%13.6e'):
        """
        Writes the file 'electron_numdens.inp' for radmc3d. 
        
        Parameters
        ----------
        dens_ion : list or array_like, shape(n,)  
           The model electronic number density.
        
        fmt : str
           Format string for numbers in the output file.
        """
        dens_e = np.asarray(dens_e)*cm**-3
        with open('electron_numdens.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = dens_e.ravel(order='F') # Create a 1-D view, fortran-style indexing
            dens_e.tofile(f, sep='\n', format=fmt)
            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_ion_numdens(self, dens_ion, fmt = '%13.6e'):
        """
        Writes the file 'ion_numdens.inp' for radmc3d. 
        
        Parameters
        ----------
        dens_ion : list or array_like, shape(n,)  
           The model ionized hydrogen number density.
        
        fmt : str
           Format string for numbers in the output file.
        """
        dens_ion = np.asarray(dens_ion)*cm**-3
        with open('ion_numdens.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = dens_ion.ravel(order='F') # Create a 1-D view, fortran-style indexing
            dens_ion.tofile(f, sep='\n', format=fmt)
            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_dust_density(self, dens_dust, nrspec = 1, fmt = '%13.6e'):
        """
        Writes the file 'dust_density.inp' for radmc3d. 
        
        Parameters
        ----------
        dens_dust : list or array_like, shape(n,)  
           The model dust density (in :math:`kg/m^3`).
        
        fmt : str
           Format string for numbers in the output file.
        """
        dens_dust = np.asarray(dens_dust)*1e3*cm**-3
        with open('dust_density.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            f.write('%d\n'%nrspec)                                          # Number of species
            #data = dens_e.ravel(order='F') # Create a 1-D view, fortran-style indexing
            dens_dust.tofile(f, sep='\n', format=fmt)
            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_gas_temperature(self, temp_gas, fmt = '%13.6e'):
        """
        Writes the file 'gas_temperature.inp' for radmc3d. 
        
        Parameters
        ----------
        temp_gas : list or array_like, shape(n,)  
           The model gas temperature.
        
        fmt : str
           Format string for numbers in the output file.
        """
        with open('gas_temperature.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = temp_gas.ravel(order='F') # Create a 1-D view, fortran-style indexing
            temp_gas.tofile(f, sep='\n', format=fmt)
            f.close()
            
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_dust_temperature(self, temp_dust, fmt = '%13.6e'):
        """
        Writes the file 'dust_temperature.inp' for radmc3d. 
        
        Parameters
        ----------
        temp_dust : list or array_like, shape(n,)  
           The model dust temperature.
        
        fmt : str
           Format string for numbers in the output file.
        """
        with open('dust_temperature.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            temp_dust.tofile(f, sep='\n', format=fmt)
            f.close()
            
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')


    def write_microturbulence(self, microturbulence, fmt = '%13.6e'):
        """
        Writes the file 'microturbulence.inp' for radmc3d. 
        See the section 7.4.8 at the `RADMC-3D`_ manual for further information. 
        
        Parameters
        ----------
        microturbulence : list or array_like, shape(n,)  
           The model microturbulent broadening.
        
        fmt : str
           Format string for numbers in the output file.
        """
        microturb = np.asarray(microturbulence) * cm
        with open('microturbulence.inp','w+') as f:
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = temp_gas.ravel(order='F') # Create a 1-D view, fortran-style indexing
            microturb.tofile(f, sep='\n', format=fmt)
            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_gas_velocity(self, vel, fmt = '%13.6e'):
        """
        Writes the file 'gas_velocity.inp' for radmc3d.
        
        Parameters
        ----------
        vel : list or array_like, shape(3,n)  
           The model velocity in x,y,z: [vel_x, vel_y, vel_z]
        
        fmt : str
           Format string for numbers in the output file.
        """
        vel2wrt = np.asarray(vel).T * cm
        tmp_write = []
        tmp = ((fmt+'\t')*3)[:-1] + '\n'
        
        with open('gas_velocity.inp','w+') as f: 
            f.write('1\n')                                          # Format number
            f.write('%d\n'%self.nn)                                 # Nr of cells
            #data = temp_gas.ravel(order='F') # Create a 1-D view, fortran-style indexing
            for i in range(self.nn): tmp_write.append(tmp % tuple(vel2wrt[i]))         
            f.writelines(tmp_write)
            #vel2wrt.tofile(f, sep='\n', format=fmt)
            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_radmc3d_control(self, **kwargs):
        """
        Writes the control file 'radmc3d.inp'.
        
        Parameters
        ----------
        **kwargs : control parameters for radmc3d. 
           For example: write_radmc3d_control(incl_freefree = 1, incl_dust = 0, setthreads = 4)
        
        Notes
        -----
        Have a look at the `RADMC-3D`_ manual, section A.1, for a comprehensive list of the available control parameters and their default values.
        """
        with open('radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(self.nphot))
            for key in kwargs: f.write('{} = {}\n'.format(key, kwargs[key]))
            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def _write_lam(self, file, lam, nxx):
        len_lam = len(lam)
        if len_lam - 1 == len(nxx):
            lam_list = [np.logspace(np.log10(lam[i]),
                                    np.log10(lam[i+1]),
                                    nxx[i], endpoint=False) 
                        for i in range(len_lam-2)]
            lam_list.append(np.logspace(np.log10(lam[-2]),
                                        np.log10(lam[-1]),
                                        nxx[-1], endpoint=True))
            lam_list = np.concatenate(lam_list)
            tmp = '%.6e'+'\n'
            for value in lam_list: file.write(tmp%(value))

        else: raise ValueError("Wrong length(s) for input list(s): len(lam)-1 must be equal to len(nxx)")
        
    def write_wavelength_micron(self, lam = [1e-1,5e2,2e4,4e4,3e5], nxx = [50,50,50,50], fmt = '%13.6e'):
        """
        Writes the file 'wavelength_micron.inp' for radmc3d.
        
        Parameters
        ----------
        lam : list or array_like,  
           Boundaries of the wavelength grid, in microns.
        
        nxx : list or array_like, length: len(lam)-1 
           Number of wavelengths within each pair of boundaries. 

        fmt : str
           Format string for numbers in the output file.
        """
        nlam = np.sum(nxx)
        with open('wavelength_micron.inp','w+') as f:
            f.write('%d\n'%(nlam))
            self._write_lam(f, lam, nxx)
            f.close()
        """
        len_lam = len(lam)
        if len_lam - 1 == len(nxx):
            lam_list = [np.logspace(np.log10(lam[i]),
                                    np.log10(lam[i+1]),
                                    nxx[i], endpoint=False) 
                        for i in range(len_lam-2)]
            lam_list.append(np.logspace(np.log10(lam[-2]),
                                        np.log10(lam[-1]),
                                        nxx[-1], endpoint=True))
            lam_list = np.concatenate(lam_list)
            nlam = lam_list.size
            with open('wavelength_micron.inp','w+') as f:
                f.write('%d\n'%(nlam))
                tmp = fmt+'\n'
                for value in lam_list: f.write(tmp%(value))
                f.close()

        else: sys.exit("ERROR: Wrong length(s) for input list(s): len(lam)-1 must be equal to len(nxx)")
        """
        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')

    def write_stars(self, nstars = 1, pos = [[0.,0.,0.]], 
                    rstars = [6.96e8], mstars = [1.99e30],
                    lam = [1e-1,5e2,2e4,4e4,3e5], nxx = [50,50,50,50], 
                    flux = [[-5780]], fmt = '%.6e'):
        """
        Writes the file 'stars.inp' for radmc3d. Defaults to Sun's properties.
        
        Parameters
        ----------
        nstars : scalar
           Number of stars invoked.
        
        pos : list of lists or array_like, shape (`nstars`,3)
           [x,y,z] position of each star in meters.
        
        rstars : list or array_like, length `nstars`
           Radii of the stars in meters.
        
        mstars : list or array_like, length `nstars`
           Mass of the stars in kilograms.

        lam : list or array_like,  
           Wavelength intervals, in microns.
        
        nxx : list or array_like, length: len(lam)-1 
           Number of wavelenghts in each interval 

        flux : list of lists or array_like, shape (`nstars`, number of wavelenghts)
           Flux from the stars at each wavelength in cgs. If the first number of the list(s) is negative, it is
           assumed to be the blackbody temperature of the star and the flux is computed accordingly. 

        fmt : str
           Format string for the flux values.

        Notes
        -----
        Have a look at the `RADMC-3D`_ manual, section A.7, for further information about the input parameters and scnearios.
        """
        nlam = np.sum(nxx)
        rstars = np.array(rstars)*cm
        mstars = np.array(mstars)*1e3
        with open('stars.inp','w+') as f:
            f.write('2\n')
            f.write('%d %d\n'%(nstars, nlam))
            for i in range(nstars): f.write('%.2e %.2e %.2e %.2e %.2e\n'
                                            %(rstars[i], mstars[i], 
                                              pos[i][0]*cm, pos[i][1]*cm, pos[i][2]*cm))
            self._write_lam(f, lam, nxx)

            for i in range(nstars): 
                if flux[i][0] < 0:
                    f.write((fmt%flux[i][0])+'\n')
                else: 
                    for j in range(nlam):
                        f.write((fmt%flux[i][j])+'\n')

            f.close()

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------')
        

class Radmc3dDefaults(Radmc3d):
    """
    Hosts predefined radiative transfer modes for `RADMC-3D`_: freefree, recombination lines.
    
    Parameters
    ----------
    GRID : `~sf3dmodels.Model.Struct`
       Grid structure in which the model was computed.
    """

    def __init__(self, GRID, nphot = 1000000):
        super(Radmc3dDefaults, self).__init__(GRID, nphot=nphot)
        

    def _freefree_keys(rstr = False):
        """
        Mandatory and optional prop object keys for the method `freefree`.

        Parameters
        ----------
        rstr : bool, optional
           If True returns the string of keys separated by commas. If False returns the list of keys.
        """
        mand = ['dens_e', 'dens_ion', 'temp_gas']
        mandstr = str(mand)[1:-1]
        opt = ['temp_dust']
        optstr = str(opt)[1:-1]
        if rstr: return {'mand': mandstr, 'opt': optstr} 
        else: return {'mand': mand, 'opt': opt} 

    def freefree(self, prop, fmt = '%13.6e', folder = './', 
                 kwargs_control = {'scattering_mode_max': 0,
                                   'incl_dust': 0,
                                   'camera_incl_stars': 0},
                 kwargs_wavelength = {'nxx': [20,20,20],
                                      'lam': [5e2,2e4,4e4,3e5]}
                 ): 
        """
        Writes the files required by `RADMC-3D`_ to compute free-free emission from the input model.

        Parameters
        ----------
        prop : dict
           Dictionary containing the model physical properties.\n
           Mandatory keys: 'dens_e', 'dens_ion', 'temp_gas'.\n
           Optional keys: 'temp_dust'.\n
           Additional keys will not be taken into account.
        
        fmt : str, optional
           Format string for numbers.

        folder : str, optional
           Sets the folder to write the files in.
        
        kwargs_control : dict, optional 
           Optional dictionary containing additional keys to be written into the radmc3d control file 'radmc3d.inp'.\n
           - If you set 'incl_dust': 1, then prop['dens_dust'] and prop['temp_dust'] must be provided. \n
           - Fixed keys for this method: 'incl_freefree': 1 \n
           Have a look at the `RADMC-3D`_ docs for all the available control commands.  

        kwargs_wavelength : {'nxx': array_like, 'lam': array_like}, optional
           Optional dictionary containing the number of points and boundaries for the grid of wavelengths to consider.
           For further information and default values take a look at the method `~sf3dmodels.Model.Radmc3d.write_wavelength_micron`.
        
        Examples
        --------
        Have a look at the following folders on the GitHub repository: 
           `examples/ionized_constant/ <https://github.com/andizq/star-forming-regions/tree/master/examples/ionized_constant>`_ \n
           `examples/ionized_powerlaw/ <https://github.com/andizq/star-forming-regions/tree/master/examples/ionized_powerlaw>`_ \n
           `examples/ionized_env+disc/ <https://github.com/andizq/star-forming-regions/tree/master/examples/ionized_env+disc>`_ \n
           `examples/outflows/ <https://github.com/andizq/star-forming-regions/tree/master/examples/outflows>`_ \n
                      
        .. image:: ../../examples/ionized_env+disc/img_HII_env+disc_.png
           :width: 220px
           :height: 170px

        .. image:: ../../examples/ionized_env+disc/sed_env+disc.png
           :width: 200px
           :height: 170px

        .. image:: ../../examples/outflows/img_outflows.png
           :width: 220px
           :height: 170px
           
        Returns
        -------
        data files : files
           amr_grid.inp, electron_numdens.inp, ion_numdens.inp, gas_temperature.inp [, dust_temperature.inp]
        """ 
        _keys = list(self._freefree_keys())
        
        func_name = inspect.stack()[0][3]
        print ('Setting %s mode for RADMC-3D...'%func_name)
        
        self.write_amr_grid()
        self.write_electron_numdens(prop['dens_e'], fmt=fmt)
        self.write_ion_numdens(prop['dens_ion'], fmt=fmt)
        self.write_gas_temperature(prop['temp_gas'], fmt=fmt)
        if 'temp_dust' in prop: self.write_dust_temperature(prop['temp_dust'], fmt=fmt)

        kwargs_control['incl_freefree'] = 1
        self.write_radmc3d_control(**kwargs_control)
        self.write_wavelength_micron(**kwargs_wavelength)

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')

    def _recomblines_keys(rstr = False):
        """
        Mandatory and optional prop object keys for the method `recomblines`.
        
        Parameters
        ----------
        rstr : bool, optional
           If True returns the string of keys separated by commas. If False returns the list of keys.  
        """
        mand = ['dens_e', 'dens_ion', 'temp_gas', 'vel']
        mandstr = str(mand)[1:-1]
        opt = ['temp_dust', 'microturbulence']
        optstr = str(opt)[1:-1]
        if rstr: return {'mand': mandstr, 'opt': optstr} 
        else: return {'mand': mand, 'opt': opt} 

    def recomblines(self, prop, transition, fmt = '%13.6e', folder = './', 
                    kwargs_control = {},
                    kwargs_wavelength = {'nxx': [20,20,20,20],
                                         'lam': [1e-1,5e2,2e4,4e4,3e5]}
                    ):
        """
        Writes the files required by `RADMC-3D`_ to compute recombination lines emission from the input model.

        Parameters
        ----------
        prop : dict
           Dictionary containing the model physical properties.\n
           Mandatory keys: 'dens_e', 'dens_ion', 'temp_gas', 'vel'.\n
           Optional keys: 'temp_dust', 'microturbulence'.\n
           Additional keys will not be taken into account.
        
        transition : [level_up, level_down]
           array_like object of scalar integers.\n
           Upper and lower levels of the transition to be computed.

        fmt : str, optional
           Format string for numbers.

        folder : str, optional
           Sets the folder to write the files in.
        
        kwargs_control : dict, optional 
           Dictionary containing additional keys to be written into the radmc3d control file 'radmc3d.inp'.\n
           - If you set 'incl_dust': 1, then prop['temp_dust'] and prop['dens_dust'] must be provided, and you need to have a dustopac.inp file as well. \n
           - Default keys for this method: 'scattering_mode_max': 0, 'camera_incl_stars': 0, 'incl_dust': 0, 'incl_freefree': 1, 'incl_lines': 1, 'lines_profile': 1,'userdef_nonlte': 1 \n
           Have a look at the `RADMC-3D`_ docs for all the available control commands.  

        kwargs_wavelength : {'nxx': array_like, 'lam': array_like}, optional
           Dictionary containing the desired range of wavelengths to be written into file.
           For further information and default values take a look at the method `~Radmc3d.write_wavelength_micron`.

        Examples
        --------
        Have a look at the `examples/recomblines/ <https://github.com/andizq/star-forming-regions/tree/master/examples/recomblines>`_ folder on the GitHub repository. 

        .. image:: ../../examples/recomblines/img_channel.png
           :width: 335px
           :height: 260px

        .. image:: ../../examples/recomblines/img_spectrum.png
           :width: 335px
           :height: 260px
        
        Returns
        -------
        data files : files
           amr_grid.inp, radmc3d.inp, electron_numdens.inp, ion_numdens.inp, gas_temperature.inp, gas_velocity.inp [, dust_temperature.inp, microturbulence.inp]

        Warnings
        --------
        The free version of radmc3d does not include the recombination lines module. 
        If you want to get this extension you will need to request it directly to the developer: `Peters+2012`_.   
        """ 
        _keys = list(self._recomblines_keys())

        func_name = inspect.stack()[0][3]
        print ('Setting %s mode for RADMC-3D...'%func_name)

        self.write_amr_grid()
        self.write_electron_numdens(prop['dens_e'], fmt=fmt)
        self.write_ion_numdens(prop['dens_ion'], fmt=fmt)
        self.write_gas_temperature(prop['temp_gas'])
        self.write_gas_velocity([prop['vel_x'],prop['vel_y'],prop['vel_z']], fmt=fmt) #--> Create this function (DONE)
        if 'temp_dust' in prop: self.write_dust_temperature(prop['temp_dust'], fmt=fmt)
        if 'dens_dust' in prop: self.write_dust_density(prop['dens_dust'], fmt=fmt)
        if 'microturbulence' in prop: self.write_microturbulence(prop['microturbulence'], fmt=fmt)

        kwargs_tmp = {}
        kwargs_tmp['scattering_mode_max'] = 0
        kwargs_tmp['incl_dust'] = 0
        kwargs_tmp['camera_incl_stars'] = 0
        kwargs_tmp['incl_freefree'] = 1
        kwargs_tmp['incl_lines'] = 1
        kwargs_tmp['lines_profile'] = 1
        kwargs_tmp['userdef_nonlte'] = 1 #User-defined non-lte mode.
        kwargs_tmp['userdef_nup'] = transition[0]
        kwargs_tmp['userdef_ndown'] = transition[1]        
        kwargs_tmp.update(kwargs_control)
        kwargs_tmp['lines_mode'] = -10 #User-defined populations: see chapter 7 on the radmc3d manual.

        self.write_radmc3d_control(**kwargs_tmp)
        self.write_wavelength_micron(**kwargs_wavelength)

        print ('%s is done!'%inspect.stack()[0][3])
        print ('-------------------------------------------------\n-------------------------------------------------')

Radmc3dRT = Radmc3dDefaults #Backwards compatibility


