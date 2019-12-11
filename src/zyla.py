from __future__ import absolute_import
from . import PyDataSource
import os

from scipy import ndimage

class Zyla(PyDataSource.Detector):
    """Zyla Camera Detector.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)

        self.add.property(center_xpos, doc='Camera X center_of_mass')
        self.add.property(center_ypos, doc='Camera Y center_of_mass')
        next(self)
        self.add.projection('calib', axis='x', name='calib_x')
        self.add.projection('calib', axis='y', name='calib_y')
        self.add.count('calib', name='calib_count')


def center_xpos(self): 
    return ndimage.measurements.center_of_mass(self.calib)[1]

def center_ypos(self):
    return ndimage.measurements.center_of_mass(self.calib)[0]

