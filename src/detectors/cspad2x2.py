from __future__ import division
import PyDataSource
import os

from pylab import *

class Cspad2x2(PyDataSource.Detector):
    """Cspad2x2 Detector Class.
    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)
       
        if 'roi' in kwargs:
            self.make_roi(kwargs['roi'])
        
        gain = kwargs.get('gain', 30)
        self.add.parameter(gain=gain)

    def make_roi(self, roi, xaxis=None, yaxis=None, **kwargs):       
        self.add.parameter(roi=roi)
        if not xaxis or len(xaxis) != roi[1][1]-roi[1][0]:
            xaxis = np.arange(roi[1][0],roi[1][1])    
        
        if not yaxis or len(yaxis) != roi[2][1]-roi[2][0]:
            yaxis = np.arange(roi[2][0],roi[2][1])    
        
        self.add.parameter(ximg=xaxis)
        self.add.parameter(yimg=yaxis)

        self._xarray_info['coords'].update({'ximg': self.ximg, 'yimg': self.yimg})
        self._xarray_info['dims'].update({'img': (['ximg', 'yimg'], (self.ximg.size, self.yimg.size))})

        self._xarray_info['dims'].update({'xrays': ([], ())})

    @property
    def img(self):
        if not self._xarray_info: 
            self._update_xarray_info()

        try:
            if self.detector.calib is not None:
                img = self.detector.calib
            else:
                return None 

            if hasattr(self, 'roi'):
                roi = self.roi
                img = img[roi[0],roi[1][0]:roi[1][1],roi[2][0]:roi[2][1]]
                #img = img[roi[0][0]:roi[0][1],roi[1][0]:roi[1][1]]
            
            return img

        except:
            return None

    @property
    def xrays(self):
        img = self.img
        if img is not None:
            return self.img.sum()/self.gain
        else:
            return None

