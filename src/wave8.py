from __future__ import absolute_import
from . import PyDataSource
import os

from pylab import *

class Wave8(PyDataSource.Detector):
    """Wave8 Functions.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)

        if False and hasattr(self.configData, 'NChannels'):
            self.add.parameter(nchannels=self.configData.NChannels)
        else: 
            self.add.parameter(nchannels=8)

        self.add.property(peaks)
        self.add.property(waveforms)
        self.add.property(count)
        bkrange = kwargs.get('bkrange', [500,600])
        self.add.parameter(bkrange=bkrange)

    def _update_xarray_info(self):

        nchannels = self.nchannels
        # Length used to be scalar -- now is array
        if isinstance(self.configData.Length, int):
            length = self.configData.Length
        else:
            length = self.configData.Length[0]

        xattrs = {'doc': 'BeamMonitor summed intensity for first 4 diode channel waveforms',
                  'unit': 'ADU'}
        self._xarray_info['dims'].update({'count': ([], (), xattrs)}) 
        
        xattrs = {'doc': 'BeamMonitor peak intensity for 8 diode channel waveforms',
                  'unit': 'ADU'}
        self._xarray_info['dims'].update({'peaks': (['ch'], nchannels, xattrs)}) 
        
        xattrs = {'doc': 'BeamMonitor waveforms for 8 diode channels',
                  'unit': 'ADU'}
        self._xarray_info['dims'].update({'waveforms': (['ch', 't'], (nchannels, length), xattrs)}) 
        self._xarray_info['coords'].update({'t': np.arange(length)}) 

def peaks(self):
    """Max value of each waveform.
    """
    return np.array([max(self.waveforms[ch]) for ch in range(self.nchannels)])

def count(self):
    """Sum of first four waveform peaks [ADU]
    """
    return sum(self.peaks[0:3])

def waveforms(self):
    """Array of 8 waveforms with background subtration from mean within self.bkrange.
    """
    wfs = []
    for ch in range(self.nchannels):
        wf = self.evtData.data_u32[ch]
        if len(wf) == 0:
            wf = self.evtData.data_u16[ch]
        back = wf[self.bkrange[0]:self.bkrange[1]].mean()
        wfs.append(-1.*wf+back)

    return np.array(wfs)


