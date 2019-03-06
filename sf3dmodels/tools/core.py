import numpy as np

__all__ = ['formatter']

def formatter(prop_keys, fmt = '%.6e', base = '', default = '%.6e'):
    """
    Formatter function for a given set of properties.
    
    Parameters
    ----------
    fmt : str or dict, optional
       If str: Format string for all the columns. Defaults to '%.6e'.\n
       If dict: the input keys of the dictionary must be present in the ``prop_keys`` list. The keys from ``prop_keys`` which
       were not specified in the ``fmt`` dict, will return the default format ``default``.
    
    base : str, optional
       Preceding format(s) to the final-computed format string. Defaults to None.
    
    default : str, optional
       Default format for non-specified keys if fmt is dict.
    
    Returns
    -------
    fmt_string : str
       String with the computed set of formats, to be used on output prints or files.
    """
    prop_keys = np.asarray(prop_keys)
    n = len(prop_keys)
    nvec = range(n)
    fmt_string = base

    if isinstance(fmt, str): #If a single format is provided
        print ("Using format '%s' for all the properties"%fmt) 
        for i in nvec: fmt_string += ' '+fmt #Same format for all properties

    elif isinstance(fmt, dict): #If a dict of formats
        fmt_list = np.array([default] * n)
        for key in fmt:
            if key in prop_keys: fmt_list = np.where(np.array(key) == prop_keys, fmt[key], fmt_list)
            else: raise KeyError("The property '%s' provided in 'fmt' was not defined in the 'prop' object of physical properties."%key)
        print ('Using formats {} for properties {}'.format(fmt_list, prop_keys)) 
        for f in fmt_list: fmt_string += ' '+f

    else: raise TypeError("Invalid type %s for 'fmt'. Please provide a valid 'fmt' object: str, list or np.ndarray"%type(fmt))
    fmt_string += '\n'
    return fmt_string
