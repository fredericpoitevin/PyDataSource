import PyDataSource

ds = PyDataSource.DataSource(exp='xpptut15',run=280)

evt  = next(ds.events)

evt.ACQ1.show_info()

import matplotlib.pyplot as plt
plt.plot(evt.ACQ1.wftime[0], evt.ACQ1.waveform[0])
plt.show()

