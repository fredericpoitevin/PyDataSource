from __future__ import absolute_import
from . import PyDataSource
import os

class Timetool(PyDataSource.Detector):
    """
    Timetool Functions.

    Paramters
    ---------
    roi : tuple
        region of interest 
        default:  ((sig_roi_lo.row, sig_roi_hi.row+1), (sig_roi_lo.column, sig_roi_hi.column+1))

    """

    def __init__(self,*args,**kwargs):

        PyDataSource.Detector.__init__(self,*args,**kwargs)
        roi = kwargs.get('roi')
        if not roi:
            conf = self.configData
            roi = ((conf.sig_roi_lo.row, conf.sig_roi_hi.row+1),
                   (conf.sig_roi_lo.column, conf.sig_roi_hi.column+1))
        self.add.roi('raw', roi=roi, name='sig', doc='Signal ROI of raw data', projection=True, )

    def add_evtData(self, attrs=None):
        import numpy as np
        #if 'projected_signal' not in self._xarray_info['dims']:
        if True:
            while not hasattr(self, 'evtData'):
                next(self._ds.events)

            if attrs is None:
                attrs = ['amplitude', 'nxt_amplitude','position_fwhm','position_pixel',
                         'position_time', 'ref_amplitude', 'projected_signal']
            for attr in attrs:
                ainfo = self.evtData._attr_info.get(attr,{})
                self.add.property('.'.join(['evtData',attr]), attr, doc=ainfo.get('doc'), unit=ainfo.get('unit'))

            if 'projected_signal' in attrs:
                config_attrs = {attr: item for attr, item in self.configData._all_values.items()}
                for attr in self.epicsData._attrs:
                    config_attrs.update({attr: getattr(self.epicsData, attr)})

                self._xarray_info['dims'].update(
                            {'projected_signal': (['X'], self.configData.signal_projection_size, config_attrs)} )

                self._xarray_info['coords'].update({'X': np.arange(self.configData.signal_projection_size)})


