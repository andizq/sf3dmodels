from .units import amu

class PropTags(object):
    """
    Tag the physical properties. Contains functions for naming and tagging the grid and physical properties.
    """
    def __init__(self):
        self.prop_kind  = {'id': 'grid',
                           'x': 'grid',
                           'y': 'grid',
                           'z': 'grid',
                           'dens_H2': 'density', 
                           'dens_p_H2': 'density', 
                           'dens_o_H2': 'density', 
                           'dens_e': 'density', 
                           'dens_H': 'density', 
                           'dens_He': 'density', 
                           'dens_Hplus': 'density', 
                           'dens_dust': 'density',
                           'temp_gas': 'temperature', 
                           'temp_dust': 'temperature',
                           'vel_x': 'velocity', 
                           'vel_y': 'velocity', 
                           'vel_z': 'velocity',
                           'abundance': 'abundance',
                           'gtdratio': 'gtdratio',
                           'dens_ion': 'density', 
                           'microturbulence': 'velocity',
                           'doppler': 'velocity',
                           'dens_mass': 'density'
                           }

        self.dens_mass = {'dens_H2': 2.0159*amu,
                          'dens_p_H2': 2.0159*amu,
                          'dens_o_H2': 2.0159*amu,
                          'dens_e': 5.486e-4*amu,
                          'dens_H': 1.00794*amu,
                          'dens_He': 4.0026*amu,
                          'dens_Hplus': 1.00739*amu,
                          'dens_ion': 1.00739*amu,
                          'dens_dust': 1 # 1 as the dust density should already be in kg/m3
                          }

    def get_prop_kind(self,prop_name):
        """
        Get what kind of physical property 'prop_name' is.
        """
        if 'abundance' in prop_name and prop_name[-1].isdigit(): prop_kind = self.prop_kind['abundance']
        else: prop_kind = self.prop_kind[prop_name]
        
        return prop_kind

    def get_dens_mass(self,prop_name):
        """
        Returns the 1-particle mass associated to the input number density. Equivalent to dens_mass[prop_name].
        """
        return self.dens_mass[prop_name]

propTags = PropTags()
