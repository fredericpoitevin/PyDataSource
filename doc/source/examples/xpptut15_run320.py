from __future__ import division
#%pylab
import PyDataSource
ds  = PyDataSource.DataSource(exp='xpptut15',run=320)
ds.configData
evt = next(ds.events)
evt.cspad.add.psplot('image')

evt.cspad.add.parameter(gain=1/35.)
evt.cspad.add.count('calib', limits=(200,10000), gain=evt.cspad.gain, unit='photons', doc='Counts in rings')
bins = np.arange(-20,200,1)/evt.cspad.gain
evt.cspad.add.histogram('calib', bins=bins, gain=evt.cspad.gain, unit='photons', doc='Gain corrected histogram', publish=True)

bins = np.arange(-20,200,1)/evt.cspad.gain
evt.cspad.add.histogram('calib', bins=bins, gain=evt.cspad.gain, unit='photons', doc='Gain corrected histogram', publish=True)
#evt.cspad.add.peak('calib_hist2', doc='Peak in cspad photon histogram')
evt.cspad.show_info()

evt.opal_1.add.psplot('raw')
evt.opal_1.add.projection('raw', axis='x', roi=((0,300),(1024,400)), publish=True)
evt.opal_1.show_info()

import numpy as np
import PyDataSource
ds = PyDataSource.DataSource(exp='cxitut13',run=10)
ds.configData
evt = next(ds.events)
evt.DscCsPad.add.psplot('image')
#Jump to event with higher Background  
evt.next(14)
evt.DscCsPad.add.histogram('calib', bins=np.arange(-20,150,1), publish=True)

evt.DscCsPad.add.parameter(gain=1/23.)
bins = np.arange(-20,150,1)*evt.DscCsPad.gain
evt.DscCsPad.add.histogram('calib', bins=bins, gain=evt.DscCsPad.gain, name='photon_hist', unit='photons', doc='Gain corrected histogram', publish=True)
evt.DscCsPad.add.count('calib', limits=(12,10000), gain=evt.DscCsPad.gain, name='photon_count', unit='photons', doc='Photon Count')
evt.next(22)
evt.DscCsPad.show_info()
#Jump to event with Diffraction Rings
evt.DscCsPad.add.projection('calib', axis='r', publish=True)

