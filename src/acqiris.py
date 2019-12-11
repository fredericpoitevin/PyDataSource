from __future__ import absolute_import
from . import PyDataSource
import os

from pylab import *

class Acqiris(PyDataSource.Detector):
    """Acqiris Functions.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)

        if hasattr(self.configData, 'nbrChannels'):
            self.add.parameter(nchannels=self.configData.nbrChannels)
        else: 
            self.add.parameter(nchannels=4)

    def _on_init(self):

        nchannels = self.nchannels
        xattrs = {'doc': 'Acqiris max value for 4 diode channel waveforms',
                  'unit': 'ADU'}
        self.add.property(peak_height, **xattrs)
        self._xarray_info['dims'].update({'peak_height': (['ch'], nchannels, xattrs)}) 
        
        xattrs = {'doc': 'Acqiris time of max value for 4 diode channel waveforms',
                  'unit': 's'}
        self.add.property(peak_time, **xattrs)
        self._xarray_info['dims'].update({'peak_time': (['ch'], nchannels, xattrs)}) 
        
def peak_height(self):
    """Max value of each waveform.
    """
    return np.array([max(self.waveform[ch]) for ch in range(self.nchannels)])

def peak_time(self):
    """Time of max of each waveform.
    """
    return np.array([self.wftime[ch][self.waveform[ch].argmax()] for ch in range(self.nchannels)])


