from __future__ import absolute_import
from . import PyDataSource
import os

from pylab import *

class Gasdet(PyDataSource.Detector):
    """FEEGasDetEnergy Functions.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)
        self.add.property(Energy, 'energy') 


def Energy(self):
    """Calculated Mean Energy of f_11 and f_12, in mJ.
    """
    return (self.f_11_ENRC+self.f_12_ENRC)/2.



