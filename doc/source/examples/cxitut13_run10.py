from pylab import *
#import numpy as np
import PyDataSource
ds = PyDataSource.DataSource(exp='cxitut13',run=10)
ds.configData
evt = next(ds.events)
evt.DscCsPad.add.parameter(gain=1/23.)
bins = np.arange(-20,150,1)*evt.DscCsPad.gain
evt.DscCsPad.add.histogram('calib', bins=bins, gain=evt.DscCsPad.gain, name='photon_hist', unit='photons', doc='Gain corrected histogram')
evt.DscCsPad.add.count('calib', limits=(12,10000), gain=evt.DscCsPad.gain, name='photon_count', unit='photons', doc='Photon Count')
evt.DscCsPad.add.projection('calib', axis='r')
evt.DscCsPad.show_info()

x = ds.to_xarray(nevents=10, max_size=1e9)

x.DscCsPad_image.mean(dim='time').plot(vmin=0,vmax=50)
plt.show()

x.DscCsPad_photon_hist.mean(dim='time').plot()

attrs = [a for a in x.keys() if a.startswith('FEEGasDetEnergy')]
df = x.reset_coords()[attrs].to_dataframe()
df.plot()
df.plot(subplots=True, sharex=True)
df.describe().T


