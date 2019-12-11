from __future__ import absolute_import
# banner used in loading with idatasource
import argparse
import sys
import time
from .psmessage import Message

def banner(ds, base='ds', time0=None):
    message = Message()
    ds = getattr(sys.modules['__main__'], base)
    message.add("*"*80)
    message.add("")
    message.add('Data loaded for the psana data_source = {:}'.format(str(ds)))
    if ds.events:
        message.add('Total Events = {:}'.format(ds.nevents))
    
    if time0:
        message.add('Load time = {:5.1f} sec'.format(time.time()-time0))
    
    message.add("")
    ds.configData.show_info()
    evt = next(ds.events)
    message.add("")
    message.add("*"*80)
    message.add('{:} is a python friendly (i.e, tab accessible) form of the psana.DataSource.'.format(base))
    message.add('Event data is accessible through aliases (psana get and keys methods are also preserved),  e.g.,')
    message.add("")
    message.add('In [1]: evt = {:}.events.next()'.format(base))
    message.add('In [2]: evt.Evr.eventCodes')
    message.add('Out[2]: {:}'.format(str(evt.Evr.eventCodes)))
    message.add("")
    message.add("The current event is also available in the {:}.events object, e.g.,".format(base)) 
    message.add("")
    message.add('In [3]: ds.events.current.Evr.eventCodes')
    message.add('Out[3]: {:}'.format(str(ds.events.current.Evr.eventCodes)))
    message.add("")
    if ds.data_source.smd:
        message.add('For offline data (using smd small data access), one can also iterate through "calib cycle" steps.  e.g.,')
        message.add("In [4]: for events in {:}.steps:".format(base))
        message.add("            for evt in events:")
        message.add("                # Do something with events...")
        message.add("")
    
    message.add("*"*80)

    return message

