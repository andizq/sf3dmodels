from __future__ import print_function
import numpy as np
from copy import copy, deepcopy

from ..Model import Struct


class UniqueCells(object):
    """
    Finds non-repeated cells on the input AREPO dataset.
    
    Parameters
    ----------
    data : dict
       Dictionary containing the physical data of the AREPO snapshot.
    
    header : dict
       Dictionary containing the header information of the AREPO snapshot.


    Notes
    -----
    This tool computes the search over the gas particles.

    """

    def __init__(self, data, header):
        self.origdata = deepcopy(data)
        self.data = data
        self.header = header
        
    def mergemass(self):
        """
        Merges the mass of replicated cells into the survivor twin cell.
        
        Returns
        -------
        indices: 1-D `numpy.ndarray`, length: *number of non-repeated cells*
           The indices of the computed non-repeated cells.
 
        Notes
        -----
        The initial data dictionary is modified on its sorting and on its ``data['mass']`` array. 
        If you want to keep stored your original data you can do so via the copy library:
        
        UPDATE THIS doc: I have defined the attribute self.origdata to save time to the user...

        .. code-block:: python

           ...
           >>> from copy import deepcopy
           >>> data_orig = deepcopy(data)
           >>> uni = UniqueCells(data, header) 
           >>> inds = uni.mergemass() 


        """

        data_pos = self.data['pos']
        ngas = self.header['num_particles'][0]
   
        li_x = np.array(data_pos[:,0])[0:ngas]
        li_y = np.array(data_pos[:,1])[0:ngas]
        li_z = np.array(data_pos[:,2])[0:ngas]
        
        r = np.linalg.norm([li_x, li_y, li_z], axis = 0)
        id_sort = np.argsort(r)

        r = r[id_sort]
        li_x = li_x[id_sort]
        li_y = li_y[id_sort]
        li_z = li_z[id_sort]

        for prop in self.data: self.data[prop][0:ngas] = self.data[prop][id_sort]  
        data_mass = self.data['mass']

        _,id_pos = np.unique(r, return_index = True) #First guess of unique points. The returned indices are sorted by less to greater r

        ind_diff = id_pos[1:] - id_pos[:-1] #Difference between item i+1 - item i to look for indices holes in the array
        ind_diff = np.append(ind_diff, ngas - id_pos[-1]) #Appending difference of last item
    
        ind_holes, = np.where(ind_diff > 1) #Holes indicate repeated cells
        #ind_noholes, = np.where(ind_diff == 1)

        id_new = []; new_list = []; merged_mass = []
        for ind in ind_holes: #iterates over holes
            id0 = id_pos[ind]
            for id in np.arange(id0+1, id0+ind_diff[ind]): #iterate over neighbor holes to check which of them are identical to id0
                if (li_x[id0],li_y[id0],li_z[id0]) == (li_x[id],li_y[id],li_z[id]): #if identical merge mass into id0
                    data_mass[id0] += data_mass[id]            
                else: #if not...
                    new_atall = True
                    for ch in np.arange(id0+1, id): #check the hole 'id' against the previous holes (!= id0) to avoid accepting duplicates
                        new = (li_x[ch],li_y[ch],li_z[ch]) != (li_x[id],li_y[id],li_z[id]) #True if different to ch
                        new_list.append(new)
                        new_atall &= new  
                    if new_atall: id_new.append(id) #True if new=True when compared to all the previous holes. Then the hole 'id' is a new point that must not be rejected
                    else: #Means id is repeated against one (or more than one) of its previous neighboring holes 
                          # (holes that in fact should have been added one step before to the list 'id_new' of new ids). 
                          #  Then, merge the mass of the repeated holes.
                        for i_n in range(len(new_list)): #Iterates over the list of neighboring holes that stores whether the point id is identical to them or not
                            new_point = new_list[i_n]
                            if not new_point: #True if the hole is twin of id 
                                data_mass[id0+1 + i_n] += data_mass[id]
                                break #Only merge with the first repeated hole which is the only one alive.
                    new_list = []

        id_new = np.array(id_new)
        id_pos = np.append(id_pos, id_new)
        np.random.shuffle(id_pos) #Shuffle ids to destroy the initial sorting by r

        print ("initial number of points", len(li_x))
        print ("final (non-repeated) points", len(id_pos))
        print ("merged mass", np.sum(data_mass[id_pos]) - np.sum(self.origdata['mass'][id_pos]))

        #NOTE: the final total mass d_obj.mass[id_pos].sum() may slightly differ to d_obj_orig.mass[:ngas].sum() due to precision approximations
        # made during the merging process. Converting the objects into a same bigger data type (say float64) which avoids rounding problem makes the 
        #  two quantities to coincide.

        return id_pos


class ReadSnapshot(UniqueCells):
    """
    Reads a binary snapshot from AREPO and stores its data and header information 
    into python dictionaries.

    Parameters
    ----------
    snap_path: str
       Path to the snapshot file to read in.

    Returns
    -------
    data: dict
       Dictionary containing the physical data of the snapshot.

    header: dict
       Dictionary containing the header information of the snapshot.    
       
    Notes
    -----
    This tool uses the main functionalities of the AREPO standard python library **arepy**.


    """

    def __init__(self, 
                 path, 
                 flags = {'mc_tracer': True,
                          'time_steps': True,
                          'sgchem': True,
                          'variable_metallicity': False}
                 ):
        try:
            from arepy.read_write import binary_read as rsnap
            from arepy.read_write import binary_write as wsnap
            from arepy.utility import cgs_constants as cgs
        except ImportError: 
            print ("The arepy package was not found, please double check whether it is \
                    in your machine and its path is defined in the $PYTHONPATH env variable")

        for key in flags.keys(): rsnap.io_flags[key] = flags[key]
        data, header = rsnap.read_snapshot(path)        

        super(ReadSnapshot, self).__init__(data, header)
        return data, header

    """
    def mergemass(self):
        return self.mergemass()
    """
