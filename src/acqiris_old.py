from __future__ import print_function
from __future__ import absolute_import
from . import PyDataSource
import os

from pylab import *

class Acqiris(PyDataSource.Detector):
    """Acqiris Functions.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)
        
        for ch in range(self.nchannels):
            name = 'Ch'+str(ch+1) 
            setattr(self, name, Channel(self,ch,name))
    #        self.add.peak(ichannel=ch)

    @property
    def nchannels(self):
        """Number of Acqiris channels.
        """
        return self.configData.nbrChannels 

    @property
    def ch_sum(self):
        """Dictionary of max value from waveforms for acqiris channels.
           (no doc in show_info)
        """
        return [getattr(self,ch).waveform.sum() 
                for ch in self._channel_dict]

    @property
    def ch_max(self):
        """Dictionary of max value from waveforms for acqiris channels.
           (no doc in show_info)
        """
        return [getattr(self,ch).waveform.max() 
                for ch in self._channel_dict]

    @property
    def ch_min(self):
        """Dictionary of max value from waveforms for acqiris channels.
           (no doc in show_info)
        """
        return [getattr(self,ch).waveform.min() 
                for ch in self._channel_dict]

    @property
    def ch_std(self):
        """Dictionary of max value from waveforms for acqiris channels.
           (no doc in show_info)
        """
        return [getattr(self,ch).waveform.std() 
                for ch in self._channel_dict]

    @property
    def _channel_dict(self):
        """Dictionary of acqiris psana data attribute index values 
           for acqiris channels.
           (no doc in show_info)
        """
        return {'Ch'+str(num+1):num for num in range(self.nchannels)}


class Channel(object):
    """Channel class for Acqiris Detector.
       Warning -- potential problem with order loading functions.  
       Currently loads according to .__dir__ order (i.e., alphabetically)
       instead of seqentially.
    """
    def __init__(self,acqiris,channel,name):
        self._acqiris = acqiris
        self._channel = channel
        self._name = name
        self.baseline = 0.

    @property
    def waveform(self):
        """Waveform of channel.
        """
        return self._acqiris.waveform[self._channel]
    
    @property
    def min(self):
        """Min of Waveform of channel.
        """
        return self.waveform.min()

    @property
    def max(self):
        """Max of Waveform of channel.
        """
        return self.waveform.max()
   
    @property
    def std(self):
        """Standard deviation of Waveform of channel.
        """
        return self.waveform.std()


    @property
    def sum(self):
        """Sum of Waveform of channel.
        """
        return self.waveform.sum()-self.baseline*self.waveform.shape[0]

    def show_info(self):
        doc = 'Acqiris waveform [min,max,std]'
        name = self._name
        attrs = ['min', 'max', 'std']
        value = list(map(int,[getattr(self.waveform,attr)() for attr in attrs]))
        print('{:8s} {:26} {:}'.format(name, value, doc))

