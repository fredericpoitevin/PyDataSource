from __future__ import absolute_import
from . import PyDataSource

class Camera(PyDataSource.Detector):
    """Camera Detector.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)

        self.add.property(center_xpos, doc='Camera X center_of_mass')
        self.add.property(center_ypos, doc='Camera Y center_of_mass')
        self.add.projection('corr', axis='x', name='img_x', axis_name='ximg')
        self.add.projection('corr', axis='y', name='img_y', axis_name='yimg')
        self.add.count('corr', name='img_count')

def center_xpos(self): 
    """
    X center of mass of camera image
    """
    from scipy import ndimage
    return ndimage.measurements.center_of_mass(self.calib)[1]

def center_ypos(self):
    """
    Y center of mass of camera image
    """
    from scipy import ndimage
    return ndimage.measurements.center_of_mass(self.calib)[0]

