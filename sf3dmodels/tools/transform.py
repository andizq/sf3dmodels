import numpy as np

#******************
#Useful TOOLS
#******************

def spherical2cartesian(r = None, theta = None, phi = None):
    """
    Converts spherical coordinates to rectangular cartesian coordinates.
    
    Parameters
    ----------
    r : scalar or array_like, shape (n,) 
       Radial distance from the origin of coordinates to the point in question.
    
    theta : scalar or array_like, shape (n,) 
       Polar angle referred to the positive `z` - axis.
    
    phi : scalar or array_like, shape (n,) 
       Azimuthal angle referred to the positive `x` - axis.
       
    Returns
    -------
    out : `numpy.ndarray`, shape (3,n) 
       Array-like object with the resulting x,y,z transformation.

    See Also
    --------
    cartesian2spherical
    
    Notes
    -----
    Transformation equations:
    
    .. math:: x &= r\\sin(\\theta)\\cos(\\phi) \n 
              y &= r\\sin(\\theta)\\sin(\\phi) \n
              z &= r\\cos(\\theta)

    """

    r = np.asarray(r)
    theta = np.asarray(theta)
    phi = np.asarray(phi)
    
    return np.array([
            r * np.sin(theta) * np.cos(phi),
            r * np.sin(theta) * np.sin(phi),
            r * np.cos(theta)
            ])

def cartesian2spherical(x=None, y=None, z=None):
    """
    Under development
    """
    return 0
