from __future__ import print_function
from __future__ import absolute_import
# standard python modules
import os
import operator
import time
import traceback
import numpy as np
from .xarray_utils import *
#from pylab import *

def runlist_to_str(runs, portable=False):
    """Convert list of runs to string representation.
    """
    if portable:
        strmulti='-'
        strsingle='_'
    else:
        strmulti=':'
        strsingle=','

    runs = np.array(sorted(set(runs)))
    runsteps = runs[1:] - runs[0:-1]
    runstr = '{:}'.format(runs[0])
    for i in range(len(runsteps)):    
        if i == 0:
            if runsteps[i] == 1:
                grouped_runs = True
            else:
                grouped_runs = False
            
            if i == len(runsteps)-1:
                runstr += '{:}{:}'.format(strsingle,runs[i+1])

        elif i == len(runsteps)-1:
            if grouped_runs:
                if runsteps[i] == 1:
                    runstr += '{:}{:}'.format(strmulti,runs[i+1])
                else:
                    runstr += '{:}{:}{:}{:}'.format(strmulti,runs[i],strsingle,runs[i+1])
            else:
                runstr += '{:}{:}:{:}'.format(strsingle,runs[i],runs[i+1])

        elif i > 0:
            if runsteps[i] == 1:
                if not grouped_runs:
                    runstr += '{:}{:}'.format(strsingle,runs[i])

                grouped_runs = True
            else:
                if grouped_runs:
                    runstr += '{:}{:}'.format(strmulti,runs[i])
                else:
                    runstr += '{:}{:}'.format(strsingle,runs[i])

                grouped_runs = False

    return runstr

def runstr_to_array(runstr, portable=False):
    """Convert run string to list.
    """
    if portable:
        strmulti='-'
        strsingle='_'
    else:
        strmulti=':'
        strsingle=','
    
    runs = []
    for item in runstr.split(strsingle):
        runrange = item.split(strmulti)
        if len(runrange) == 2:
            for run in range(int(runrange[0]),int(runrange[1])+1):
                runs.append(run)
        else:
            runs.append(int(runrange[0]))

    return np.array(sorted(runs))

def read_netcdfs(files, dim='time', transform_func=None, engine='h5netcdf'):
    """
    Read netcdf files and return concatenated xarray Dataset object

    Parameters
    ----------
    files : str or list
        File name(s) which may have '*' and '?' to be used for matching
    dim : str
        Diminsion along which to concatenate Dataset objects
        Default = 'time'
    transform_func : object
        Method to transform each Dataset before concatenating
    engine : str
        Engine for loading files.  default = 'h5netcdf'

    """
    import glob
    import xarray as xr
    def process_one_path(path):
        # use a context manager, to ensure the file gets closed after use
        with xr.open_dataset(path, engine=engine) as ds:
            # transform_func should do some sort of selection or
            # aggregation
            if transform_func is not None:
                ds = transform_func(ds)
            # load all data from the transformed dataset, to ensure we can
            # use it after closing each original file
            ds.load()
            return ds

    paths = sorted(glob.glob(files))
    datasets = []
    for p in paths:
        try:
            xo = process_one_path(p)
            datasets.append(xo)
            xattrs = xo.attrs
        except:
            print('Cannot open {:}'.format(p))
    
    x = resort(xr.concat(datasets, dim))
    x.attrs.update(**xattrs)
    x.attrs.pop('ichunk')
    x.attrs['files'] = paths
    return x 


#ds = PyDataSource.DataSource('exp=cxij4915:run=49:smd')
def open_h5netcdf(file_name=None, path='', file_base=None, exp=None, run=None, 
        h5folder='scratch', subfolder='nc', chunk=False, summary=False):
    """
    Open hdf5 file with netcdf4 convention using builtin xarray engine h5netcdf.
    """
    import xarray as xr
    if exp:
        instrument = exp[0:3]
    
    if not file_name and not path:
        path = '/reg/d/psdm/{:}/{:}/{:}/{:}'.format(instrument, exp, h5folder, subfolder)

    if chunk and run:
        file_names = '{:}/run{:04}_c*.nc'.format(path, int(run))
        if True:
            x = read_netcdfs(file_names)
        else:
            import glob
            files = glob.glob(file_names)
            xo = xr.open_dataset(files[0], engine='h5netcdf')
            x = resort(xr.open_mfdataset(file_names, engine='h5netcdf'))
            x.attrs.update(**xo.attrs)
            x.attrs.pop('ichunk')
            x.attrs['files'] = files

        return x

    else:
        if not file_name:
            if not file_base:
                file_base = 'run{:04}'.format(int(run))
                if summary:
                    file_base += '_sum'

            file_name = os.path.join(path,file_base+'.nc')

        return xr.open_dataset(file_name, engine='h5netcdf')

def to_h5netcdf(xdat=None, ds=None, file_name=None, path=None, 
        h5folder='scratch', subfolder='nc', **kwargs):
    """Write hdf5 file with netcdf4 convention using builtin xarray engine h5netcdf.
    """
    if xdat:
        if xdat.__class__.__name__ == 'DataSource':
            # 1st arg is actually PyDataSource.DataSource
            ds = xdat
            xdat = None

    if not xdat:
        if not ds:
            from . import PyDataSource
            ds = PyDataSource.DataSource(**kwargs)

        xdat = to_xarray(ds, **kwargs)
    
    if not path:
        path = '/reg/d/psdm/{:}/{:}/{:}/{:}/'.format(xdat.instrument,xdat.experiment,h5folder,subfolder)

    if not os.path.isdir(path):
        os.mkdir(path)

    if not file_name:
        if 'ichunk' in xdat.attrs:
            file_name = '{:}/run{:04}_c{:03}.nc'.format(path, int(xdat.run[0]), xdat.attrs['ichunk'])

        else:
            # add in sets for mulit run
            file_base = xdat.attrs.get('file_base')
            if not file_base:
                if 'stat' in xdat.dims:
                    run = xdat.run
                    file_base = 'run{:04}_sum'.format(int(run))
                else:
                    run = sorted(set(xdat.run.values))[0]
                    file_base = 'run{:04}'.format(int(run))
                
            file_name = os.path.join(path,file_base+'.nc')
    
    xdat.to_netcdf(file_name, engine='h5netcdf')

    return xdat

#def to_summary(x, dim='time', groupby='step', 
#        save_summary=False,
#        normby=None,
#        omit_list=['run', 'sec', 'nsec', 'fiducials', 'ticks'],
#        stats=['mean', 'std', 'min', 'max', 'count']):
#    """
#    Summarize a run.
#    
#    Parameters
#    ----------
#    x : xarray.Dataset
#        input xarray Dataset
#    dim : str
#        dimension to summarize over [default = 'time']
#    groupby : str
#        coordinate to groupby [default = 'step']
#    save_summary : bool
#        Save resulting xarray Dataset object
#    normby : list or dict
#        Normalize data by attributes
#    omit_list : list
#        List of Dataset attributes to omit
#    stats : list
#        List of statistical operations to be performed for summary.
#        Default = ['mean', 'std', 'min', 'max', 'count']
#
#    """
#    import xarray as xr
#    xattrs = x.attrs
#    data_attrs = {attr: x[attr].attrs for attr in x}
#
#    if 'Damage_cut' in x:
#        x = x.where(x.Damage_cut).dropna(dim)
#   
#    coords = [c for c in x.coords if c != dim and c not in omit_list and dim in x.coords[c].dims] 
#    x = x.reset_coords(coords)
#
#    if isinstance(normby, dict):
#        for norm_attr, attrs in normby.items():
#            x = normalize_data(x, attrs, norm_attr=norm_attr)
#    elif isinstance(normby, list):
#        x = normalize_data(x, normby)
#
#    # With new xarray 0.9.1 need to make sure loaded otherwise h5py error
#    x.load()
#
#    if groupby:
#        x = x.groupby(groupby)
#
#    dsets = [getattr(x, func)(dim=dim) for func in stats]
#    x = xr.concat(dsets, stats).rename({'concat_dim': 'stat'})
#    for attr,val in xattrs.items():
#        x.attrs[attr] = val
#    for attr,item in data_attrs.items():
#        if attr in x:
#            x[attr].attrs.update(item)
#
#    for c in coords:                                                       
#        x = x.set_coords(c)
#    
#    x = resort(x)
#    
#    if save_summary:
#        to_h5netcdf(x)
#
#    return x

def add_steps(x, attr, name=None):
    vals = getattr(x, attr).values
    steps = np.sort(list(set(vals)))
    asteps = np.digitize(vals, steps)
    if not name:
        name = attr+'_step'
 
    x.coords[name] = (['time'], asteps)


def add_index(x, attr, name=None, nbins=8, bins=None, percentiles=None):

    if not bins:
        if not percentiles:
            percentiles = (arange(nbins+1))/float(nbins)*100.

        bins = np.percentile(x[attr].to_pandas().dropna(), percentiles)

    if not name:
        name = attr+'_index'

    #per = [percentiles[i-1] for i in np.digitize(x[attr].values, bins)]
    x[name] = (['time'], np.digitize(x[attr].values, bins))

# Need to add in 'chunking based on steps'
def to_xarray(ds=None, nevents=None, max_size=10001, 
        xbase=None, 
        publish=False,
        store_data=[],
        chunk_steps=True,
        ichunk=None,
        nchunks=24,
        #code_flags={'XrayOff': [162], 'XrayOn': [-162], 'LaserOn': [183, -162], 'LaserOff': [184, -162]},
        code_flags={'XrayOff': [162], 'XrayOn': [-162]},
        drop_unused_codes=True,
        pvs=[], epics_attrs=[], 
        eventCodes=None,  
        save=None, **kwargs):
    """
    Build xarray object from PyDataSource.DataSource object.
       
    Parameters
    ----------
    max_size : uint
        Maximum array size of data objects to build into xarray.
    ichunk: int
        chunk index (skip ahead nevents*ichunk)
    pvs: list
        List of pvs to be loaded vs time
    epics_attrs: list
        List of epics pvs to be saved as run attributes based on inital value 
        of first event.
    code_flags : dict
        Dictionary of event code flags. 
        Default = {'XrayOff': [162], 'XrayOn': [-162]}
    eventCodes : list
        List of event codes 
        Default is all event codes in DataSource
    drop_unused_codes : bool
        If true drop unused eventCodes [default=True]

    Example
    -------
    import PyDataSource
    ds = PyDataSource.DataSource(exp='xpptut15',run=200)
    evt = ds.events.next()
    evt.opal_1.add.projection('raw', axis='x', roi=((0,300),(1024,400)))
    evt.cs140_rob.add.roi('calib',sensor=1,roi=((104,184),(255,335)))
    evt.cs140_rob.add.count('roi')

    """
    try:
        import xarray as xr
    except:
        raise Exception('xarray package not available. Use for example conda environment with "source conda_setup"')

    if not ds:
        PyDataSource
        ds = PyDataSource.DataSource(**kwargs)
  
    adat = {}
    ds.reload()
    evt = ds.events.next(publish=publish, init=publish)
    dtime = evt.EventId
    if not eventCodes:
        eventCodes = sorted(ds.configData._eventcodes.keys())
    
    if hasattr(ds.configData, 'ScanData') and ds.configData.ScanData:
        nsteps = ds.configData.ScanData.nsteps
    else:
        nsteps = 1
    
    ievent0 = 0
    if not nevents:
        if ichunk is not None:
            istep = ichunk-1
            if chunk_steps and nsteps > 1:
                ievent_start = ds.configData.ScanData._scanData['ievent_start'][istep]
                ievent_end = ds.configData.ScanData._scanData['ievent_end'][istep]
                nevents = ievent_end-ievent_start+1
                ievent0 = ievent_start
            else:
                nevents = int(np.ceil(ds.nevents/float(nchunks)))
                ievent0 = (ichunk-1)*nevents
            
            print('Do {:} of {:} events for {:} chunk'.format(nevents, ds.nevents, ichunk))
        else:
            nevents = ds.nevents
 
    neventCodes = len(eventCodes)
    det_funcs = {}
    epics_pvs = {}
    for pv in pvs:
        epics_pvs[pv] = {} 

    for srcstr, src_info in ds.configData._sources.items():
        det = src_info['alias']
        nmaxevents = 100
        ievt = 0
        try:
            while det not in evt._attrs and ievt < nmaxevents:
                evt.next(publish=publish, init=publish)
                ievt += 1
            
            detector = getattr(evt,det)
            if hasattr(detector, '_update_xarray_info'):
                print('updating', srcstr, det)
                detector._update_xarray_info()

            # Note that the a and b objects here link the det_funcs dictionary 
            # to the adat dictionary of xarray.Dataset objects.
            # Thus, updating the det_funcs dictionary updates the xarray.Dataset objects.
            adat[det] = {}
            det_funcs[det] = {}
            xarray_dims = detector._xarray_info.get('dims')
            if xarray_dims is not None: 
                for attr,item in sorted(xarray_dims.items(), key=operator.itemgetter(0)):
                    # Only save data with less than max_size total elements
                    alias = det+'_'+attr
                    if len(item) == 3:
                        attr_info = item[2]
                    else:
                        attr_info = {}

                    det_funcs[det][attr] = {'alias': alias, 'det': det, 'attr': attr, 'attr_info': attr_info}
                    if np.product(item[1]) <= max_size or alias in store_data:
                        a = [det+'_'+name for name in item[0]]
                        a.insert(0, 'time')
                        try:
                            b = list(item[1])
                        except:
                            b = [item[1]]
                        b.insert(0, nevents)

                        adat[det][alias] = (a, tuple(b))
                        det_funcs[det][attr]['event'] = {'dims': a, 'shape': b}

        except:
            print('ERROR loading', srcstr, det)
            traceback.print_exc()

    axdat = {}
    atimes = {}
    btimes = []
    if not xbase:
        xbase = xr.Dataset()
    
    # Experiment Attributes
    xbase.attrs['data_source'] = str(ds.data_source)
    xbase.attrs['run'] = ds.data_source.run
    for attr in ['instrument', 'experiment', 'expNum', 'calibDir']:
        xbase.attrs[attr] = getattr(ds, attr)

    xbase.coords['time'] = np.zeros(nevents, dtype=dtime.datetime64.dtype)
    
    ttypes = {'sec': 'int32', 
              'nsec': 'int32', 
              'fiducials': 'int32', 
              'ticks': 'int32', 
              'run': 'int32'}
    
    # explicitly order EventId coords in desired order 
    print('Begin processing {:} events'.format(nevents))
    for attr in ['sec', 'nsec', 'fiducials', 'ticks', 'run']:
        #dtyp = ttypes[attr]
        #xbase.coords[attr] = (['time'], np.zeros(nevents,dtype=dtyp))
        dtyp = int
        xbase.coords[attr] = (['time'], np.zeros(nevents,dtype=int))

    xbase.coords['step'] = (['time'], np.empty(nevents,dtype=int))
    
    # Event Codes -- earlier bool was not supported but now is. 
    for code in eventCodes:
        xbase.coords['ec{:}'.format(code)] = ('time', np.zeros(nevents, dtype=bool))

    for attr, ec in code_flags.items():
        xbase.coords[attr] = ('time', np.zeros(nevents, dtype=bool))
        xbase.coords[attr].attrs['doc'] = 'Event code flag: True if all positive and no negative "codes" are in eventCodes'
        xbase.coords[attr].attrs['codes'] = ec

    xbase.attrs['event_flags'] = code_flags.keys()

    xbase.coords['steps'] = range(nsteps)
    xbase.coords['codes'] = eventCodes
 
    # Scan Attributes -- cannot put None or dicts as attrs in netcdf4
    # e.g., pvAliases is a dict
    if hasattr(ds.configData, 'ScanData') and ds.configData.ScanData:
        if ds.configData.ScanData.nsteps == 1:
            attrs = ['nsteps']
        else:
            attrs = ['nsteps', 'pvControls', 'pvMonitors', 'pvLabels']
            #xbase.coords['pvControls'] = ds.configData.ScanData.pvControls
            for attr, vals in ds.configData.ScanData.control_values.items():
                alias = ds.configData.ScanData.pvAliases[attr]
                xbase.coords[alias+'_steps'] = (['steps'], vals) 
                xbase.coords[alias] = ('time', np.zeros(nevents, dtype=bool))

        for attr in attrs:
            val = getattr(ds.configData.ScanData, attr)
            if val:
                xbase.attrs[attr] = val 

   
    for srcstr, item in sorted(ds.configData._sources.items(), key=operator.itemgetter(0)):
        det = item['alias']
        if det in adat:
            atimes[det] = []
            axdat[det] = xr.Dataset()
            for attr, item in sorted(adat[det].items(), key=operator.itemgetter(0)):
                axdat[det][attr] = (item[0], np.zeros(item[1]))
            
            axdat[det].coords['steps'] = range(nsteps)
            axdat[det].coords['codes'] = eventCodes

    for srcstr, srcitem in sorted(ds.configData._sources.items(), key=operator.itemgetter(0)):
        src_info = ds.configData._sources.get(srcstr).copy()
        src_info['src'] = str(src_info.get('src'))
        det = srcitem.get('alias')
        readoutCode = srcitem.get('eventCode')
        detector = getattr(evt, det)
        det_func = det_funcs.get(det)
        if det_func and det in adat:
            while det not in evt._attrs:
                evt.next(publish=publish, init=publish)

            detector._update_xarray_info()
            # Make a config object for each detector
            # In future check if a detector config changes during a run and 
            # create coords for each attr that changes in steps during the run.
            #axdat[det].coords['steps'] = range(nsteps)

            config_info = detector._xarray_info.get('attrs')
#            if attrs:
#                axdat[det].coords[det+'_config'] = ([det+'_steps'], range(nsteps))
#                axdat[det].coords[det+'_config'].attrs.update(attrs)

            for attr, attr_func in det_func.items():
                alias = attr_func.get('alias')
                attrs_info = attr_func.get('attr_info', {})
                if detector._tabclass == 'evtData':
                    if detector.evtData is not None:
                        item = detector.evtData._attr_info.get(attr)
                        if item:
                            attrs_info.update({a: item[a] for a in ['doc', 'unit']})
                    else:
                        print('No data for {:} in {:}'.format(str(detector), attr))
                
                else:
                    if detector._calib_class is not None:
                        item = detector.calibData._attr_info.get(attr)
                        if item is not None:
                            attrs_info.update(item)
                    
                    elif detector.detector is not None:
                        item = detector.detector._attr_info.get(attr)
                        if item is not None:
                            attrs_info.update(item)

                if src_info:
                    attrs_info.update(src_info)

                # Make sure no None attrs
                for a, aitm in attr_info.items():
                    if aitm is None:
                        attr_info.update({a, ''})

                if 'event' in attr_func:
                    axdat[det][alias].attrs.update(attrs_info)
                
            coords = detector._xarray_info.get('coords')
            if coords:
                for coord, item in sorted(coords.items(), key=operator.itemgetter(0)):
                    try:
                        if isinstance(item, tuple):
                            dims = [det+'_'+dim for dim in item[0]]
                            vals = item[1]
                            axdat[det].coords[det+'_'+coord] = (dims,vals)
                        else:
                            axdat[det].coords[det+'_'+coord] = item
                    except:
                        print(det, coord, item)

    ds.reload()
    print('xarray Dataset configured')

    time0 = time.time()
    igood = -1
    aievt = {}
    aievents = {}

    # keep track of events for each det
    for srcstr, srcitem in ds.configData._sources.items():
        det = srcitem.get('alias')
        aievt[det] = -1
        aievents[det] = []
  
    if ichunk is not None:
        print('Making chunk {:}'.format(ichunk))
        print('Starting with event {:} of {:}'.format(ievent0,ds.nevents))
        print('Analyzing {:} events'.format(nevents))
        xbase.attrs['ichunk'] = ichunk
        # Need to update to jump to event.
        if ichunk > 1 and not chunk_steps:
            for i in range(ievent0):
                evt = next(ds.events)
        
            print('Previous event before current chunk:', evt)

    if ichunk is not None:
        evtformat = '{:10.1f} sec, Event {:} of {:} in chunk with {:} accepted'
    else:
        evtformat = '{:10.1f} sec, Event {:} of {:} with {:} accepted'
    
    #for ievent in range(ds.nevents+1):
    for ievt in range(nevents):
        ievent = ievent0+ievt
        if ievt > 0 and (ievt % 100) == 0:
            print(evtformat.format(time.time()-time0, ievt, nevents, igood+1))
        
        if ichunk > 0 and chunk_steps:
            if ievt == 0:
                # on first event skip to the desired step
                for i in range(ichunk):
                    step_events = next(ds.steps)
            
            try:
                evt = next(step_events)
            except:
                ievent = -1
                continue

        elif ievent < ds.nevents:
            try:
                evt = ds.events.next(publish=publish, init=publish)
            except:
                ievent = -1
                continue
        else:
            ievent = -1
            continue

        if len(set(eventCodes) & set(evt.Evr.eventCodes)) == 0:
            continue
       
        dtime = evt.EventId
        if dtime is None:
            continue
        
        igood += 1
        if igood+1 == nevents:
            break

        istep = ds._istep
        xbase['step'][igood] = istep
        btimes.append(dtime)
        for ec in evt.Evr.eventCodes:
            if ec in eventCodes:
                xbase['ec{:}'.format(ec)][igood] = True

        for attr, codes in code_flags.items():
            if evt.Evr.present(codes):
                xbase.coords[attr][igood] = True

        for pv, pvarray in epics_pvs.items():
            try:
                val = float(ds.epicsData.getPV(pv).data()) 
                pvarray.update({dtime: val})
            except:
                print('cannot update pv', pv, dtime)

        for det in evt._attrs:
            detector = evt._dets.get(det)
            atimes[det].append(dtime)
            aievt[det] += 1 
            ievt = aievt[det]
            aievents[det].append(ievent)
            
            for attr, attr_func in det_funcs.get(det, {}).items():
                vals = getattr(detector, attr)
                alias = attr_func.get('alias')
                if vals is not None and 'event' in attr_func:
                    # Fill event
                    try:
                        
                        axdat[det][alias][ievt] = vals
                    except:
                        print('Event Error', alias, det, attr, ievent, vals)
                        print(axdat[det][alias][ievt].shape, vals.shape)
                        vals = None

    xbase = xbase.isel(time=range(len(btimes)))
    xbase['time'] =  [e.datetime64 for e in btimes]
    for attr, dtyp in ttypes.items():
        xbase.coords[attr] = (['time'], np.array([getattr(e, attr) for e in btimes],dtype=dtyp))
        
    # fill each control PV with current step value
    scan_variables = []
    if ds.configData.ScanData and ds.configData.ScanData.nsteps > 1:
        for attr, vals in ds.configData.ScanData.control_values.items():
            alias = ds.configData.ScanData.pvAliases[attr]
            scan_variables.append(alias)
            xbase.coords[alias] = (['time'], xbase.coords[alias+'_steps'][xbase.step]) 

    xbase.attrs['scan_variables'] = scan_variables
    xbase.attrs['correlation_variables'] = []
   
    if drop_unused_codes:
        for ec in eventCodes:
            print('Dropping unused eventCode', ec)
            if not xbase['ec{:}'.format(ec)].any():
                xbase = xbase.drop('ec{:}'.format(ec))

    # add in epics_attrs (assumed fixed over run)
    for pv in epics_attrs:
        try:
            xbase.attrs.update({pv: ds.epicsData.getPV(pv).data()[0]})
        except:
            print('cannot att epics_attr', pv)
            traceback.print_exc()

    # cut down size of xdat
    det_list = [det for det in axdat]
    for det in np.sort(det_list):
        nevents = len(atimes[det])
        if nevents > 0 and det in axdat:
            try:
                print('merging', det, nevents)
                xdat = axdat.pop(det)
                if 'time' in xdat.dims:
                    xdat = xdat.isel(time=range(nevents))
                    xdat['time'] = [e.datetime64 for e in atimes[det]]
                    xdat = xdat.reindex_like(xbase)
        
                xbase = xbase.merge(xdat)
            
            except:
                print('Could not merge', det)
                return xbase, xdat, axdat, atimes, btimes

    attrs = [attr for attr,item in xbase.data_vars.items()] 
    for attr in attrs:
        for a in ['unit', 'doc']:
            if a in xbase[attr].attrs and xbase[attr].attrs[a] is None:
                xbase[attr].attrs[a] = ''
    
    for pv, pvdata in epics_pvs.items():
        xdat = xr.Dataset({pv: (['time'], np.array(pvdata.values()).squeeze())}, 
                              coords={'time': [e.datetime64 for e in pvdata.keys()]} )
        xbase = xbase.merge(xdat)

    xbase = resort(xbase)

    if save:
        try:
            to_h5netcdf(xbase)
        except:
            print('Could not save to_h5netcdf')

    return xbase


def normalize_data(x, variables=[], norm_attr='PulseEnergy', name='norm', quiet=True):
    """
    Normalize a list of variables with norm_attr [default = 'PulseEnergy']
    """
    if not variables:
        variables = [a for a in get_correlations(x) if not a.endswith('_'+name)]    

    for attr in variables:
        aname = attr+'_'+name
        try:
            x[aname] = x[attr]/x[norm_attr]
            x[aname].attrs = x[attr].attrs
            try:
                x[aname].attrs['doc'] = x[aname].attrs.get('doc','')+' -- normalized to '+norm_attr
                units = x[attr].attrs.get('unit')
                norm_units = x[norm_attr].attrs.get('unit')
                if units and norm_units:
                    x[aname].attrs['unit'] = '/'.join([units, norm_units])
            except:
                if not quiet:
                    print('cannot add attrs for', aname)
        except:
            print('Cannot normalize {:} with {:}'.format(attr, norm_attr))

    return  resort(x)

#def resort(x):
#    """
#    Resort alphabitically xarray Dataset
#    """
#    coords = sorted([c for c in x.coords.keys() if c not in x.coords.dims])
#    x = x.reset_coords()
#    x = x[sorted(x.data_vars)]
#
#    for c in coords:                                                       
#        x = x.set_coords(c)
#
#    return x

def map_indexes(xx, yy, ww):                                                                      
    """
    Simplified map method from PSCalib.GeometryAccess.img_from_pixel_arrays
    
    Parameters
    ----------
    xx : array-like
        Array of x coordinates
    yy : array-like
        Array of y coordinates
    ww : array-like
        Array of weights

    Returns
    -------
    2D image array

    """
    a = np.zeros([xx.max()+1,yy.max()+1])
    a[xx,yy] = ww
    return a

def xy_ploterr(a, attr=None, xaxis=None, title='', desc=None, fmt='o', **kwargs):
    """
    Plot summary data with error bars. 
    Matplotlib plotting options supported through pandas DataFrame.plot method
    may be passed as kwargs.

    Parameters
    ----------
    a : xarray.Dataset
        input xarray summary Dataset
    attr : str
        Dataset attribute for plotting
    xaxis : str
        Dataset attribute for x-axis 
    title : str
        Plot title
    desc : str
        Plot description [Default = 'doc' attribute of plotting data]
    fmt : str
        Matplotlib option for representing data points

    Example
    -------
    xy_ploterr(x, 'MnScatter','Sample_z',logy=True)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    if not attr:
        print('Must supply plot attribute')
        return

    if 'groupby' in kwargs:
        groupby=kwargs['groupby']
    elif 'step' in a.dims:
        groupby='step'
    else:
        groupby='run'

    run = a.attrs.get('run')
    experiment = a.attrs.get('experiment', '')
    runstr = '{:} Run {:}'.format(experiment, run)
    name = a.attrs.get('name', runstr)
    if not title:
        title = '{:}: {:}'.format(name, attr)

    if not xaxis:
        xaxis = a.attrs.get('scan_variables')
        if xaxis:
            xaxis = xaxis[0]

    ylabel = kwargs.get('ylabel', '')
    if not ylabel:
        ylabel = a[attr].name
        unit = a[attr].attrs.get('unit')
        if unit:
            ylabel = '{:} [{:}]'.format(ylabel, unit)

    xlabel = kwargs.get('xlabel', '')
    if not xlabel:
        xlabel = a[xaxis].name
        unit = a[xaxis].attrs.get('unit')
        if unit:
            xlabel = '{:} [{:}]'.format(xlabel, unit)
    
    if xaxis:
        if 'stat' in a[xaxis].dims:
            xerr = a[xaxis].sel(stat='std').values
            a[xaxis+'_axis'] = ([groupby], a[xaxis].sel(stat='mean').values)
            xaxis = xaxis+'_axis'
        else:
            xerr = None

        a = a.swap_dims({groupby:xaxis})
    
    else:
        xerr = None

    if desc is None:
        desc = a[attr].attrs.get('doc', '')

    ndims = len(a[attr].dims)
    if ndims == 2:
        c = a[attr].to_pandas().T
        if xerr is not None:
            c['xerr'] = xerr
        c = c.sort_index()
        
        plt.figure()
        plt.gca().set_position((.1,.2,.8,.7))
        p = c['mean'].plot(yerr=c['std'],xerr=c.get('xerr'), title=title, fmt=fmt, **kwargs)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if desc:
            plt.text(-.1,-.2, desc, transform=p.transAxes, wrap=True)   
 
        return p 
    elif ndims == 3:
        plt.figure()
        plt.gca().set_position((.1,.2,.8,.7))
        pdim = [d for d in a[attr].dims if d not in ['stat', groupby, xaxis]][0]
        for i in range(len(a[attr].coords[pdim])):
            c = a[attr].sel(**{pdim:i}).drop(pdim).to_pandas().T.sort_index()
            p = c['mean'].plot(yerr=c['std'], fmt=fmt, **kwargs)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        p.set_title(title)
        if desc:
            plt.text(-.1,-.2, desc, transform=p.transAxes, wrap=True)   

        return p 
    else:
        print('Too many dims to plot')

def make_image(self, pixel=.11, ix0=None, iy0=None):
    """Return image from 3-dim detector DataArray."""
    import numpy as np
    base = self.name.split('_')[0]
    xx = self[base+'_indexes_x']
    yy = self[base+'_indexes_y']
    a = np.zeros([xx.max()+1,yy.max()+1])

    x = self.coords.get(base+'_ximage')
    if not x:
        if not ix0:
            ix0 = a.shape[1]/2.
        x = (np.arange(a.shape[1])-ix0)*pixel
    y = self.coords.get(base+'_yimage')
    if not y:
        if not iy0:
            iy0 = a.shape[0]/2.
        y = (np.arange(a.shape[0])-iy0)*pixel
    a[xx,yy] = self.data
    if x is not None and y is not None:
        try:
            return xr.DataArray(a, coords=[(base+'_yimage', y), (base+'_ximage', x)],
                                   attrs=self.attrs)
        except:
            pass

    return xr.DataArray(a)
    
    #return xr.DataArray(a, coords=[(base+'_yimage', y.mean(axis=1)), (base+'_ximage', x.mean(axis=0))])


