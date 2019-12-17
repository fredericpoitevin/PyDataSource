"""
Beam statistics methods
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
import six
import logging
import traceback
from IPython.core.debugger import Tracer

meta_attrs = {'units': 'EGU', 'PREC': 'PREC', 'pv': 'name'}
_transmission_pvs = {
                'FEE': {
                        'fee_trans':     'SATT:FEE1:320:RACT', 
                        'feegas_trans':     'GATT:FEE1:310:R_ACT', 
                        #'fee_trans_set': 'SATT:FEE1:320:RDES', 
                        #'feegas_trans_set': 'GATT:FEE1:310:R_DES', 
                        },
                'XPP': {
                        'xpp_trans':        'XPP:ATT:COM:R_CUR',
                        'xpp_trans3':        'XPP:ATT:COM:R3_CUR',
                        #'xpp_trans_set':    'XPP:ATT:COM:R_DES',
                        #'xpp_trans3_set':    'XPP:ATT:COM:R3_DES',
                        },
                'XCS': {
                        'xcs_trans':        'XCS:ATT:COM:R_CUR',
                        'xcs_trans3':        'XCS:ATT:COM:R3_CUR',
                        #'xcs_trans_set':    'XCS:ATT:COM:R_DES',
                        #'xcs_trans3_set':    'XCS:ATT:COM:R3_DES',
                        },
                'MFX': {
                        'mfx_trans':        'MFX:ATT:COM:R_CUR',
                        'mfx_trans3':        'MFX:ATT:COM:R3_CUR',
                        #'mfx_trans_set':    'MFX:ATT:COM:R_DES',
                        #'mfx_trans3_set':    'MFX:ATT:COM:R3_DES',
                        },
                'CXI': {
                        'cxi_trans':        'CXI:DSB:ATT:COM:R_CUR', 
                        'dsb_trans':        'XRT:DIA:ATT:COM:R_CUR',
                        'cxi_trans3':        'CXI:DSB:ATT:COM:R3_CUR', 
                        'dsb_trans3':        'XRT:DIA:ATT:COM:R3_CUR',
                        #'cxi_trans_set':    'CXI:DSB:ATT:COM:R_DES', 
                        #'dsb_trans_set':    'XRT:DIA:ATT:COM:R_DES',
                        #'cxi_trans3_set':    'CXI:DSB:ATT:COM:R3_DES', 
                        #'dsb_trans3_set':    'XRT:DIA:ATT:COM:R3_DES',
                        },
                }

def load_exp_sum(exp, instrument=None, path=None, nctype='drop_sum', save=True):
    """
    Load drop stats summary for all runs
    
    Arguments
    ---------
    exp : str
        experiment name

    Parameters
    ----------
    instrument : str
        instrument name.  default = exp[:3] 
    
    path : str
        path of run summary files

    save : bool
        Save drop_summary file 

    """
    from . import xarray_utils
    import os
    import glob
    import xarray as xr
    import numpy as np
    import pandas as pd
    if not instrument:
        instrument = exp[0:3]
    if not path:
        path = os.path.join('/reg/d/psdm/',instrument,exp,'results','nc')

    files = sorted(glob.glob('{:}/run*_{:}.nc'.format(path,nctype)))

    axstats = {}
    dvars = []
    for f in files:                                      
        x = xr.open_dataset(f, engine='h5netcdf')
        try:
            if x.data_vars:
                axstats[x.run] = x
                dvars += x.data_vars
                print('Loading Run {:}: {:}'.format(x.run, f))
            else:
                print('Skipping Run {:}: {:}'.format(x.run, f))
        except:
            print(('cannot do ', f))

    dvars = sorted(list(set(dvars)))
    ax = []
    aattrs = {}
    for attr in dvars:
        adf = {}
        for run, x in axstats.items():
            if attr in x.data_vars:
                if attr not in aattrs:
                    aattrs[attr] = {}
                adf[run] = x[attr].to_pandas()
                for a, val in x[attr].attrs.items():
                    if val:
                        aattrs[attr][a] = val 

        ax.append(xr.Dataset(adf).to_array().to_dataset(name=attr))

    x = xr.merge(ax)
    del(ax)
    x = xarray_utils.resort(x).rename({'variable': 'run'})
    x.attrs['experiment'] = exp
    x.attrs['instrument'] = instrument
    x.attrs['expNum'] = list(axstats.values())[0].attrs['expNum']
    x.coords['dvar'] = dvars
    dattrs = {a.replace('_detected',''): str(a) for a in axstats.values()[0].attrs.keys() if a.endswith('detected')}
    advar = {}
    for attr in dattrs:
        advar[attr] = {}
        for dvar in dvars:
            advar[attr][dvar] = np.zeros((len(x.run)), dtype='bool')
    
    rinds = dict(zip(x.run.values, range(x.run.size)))
    for run, xo in axstats.items():
        irun = rinds[run]
        for attr in dattrs:
            for dvar in xo.attrs.get(attr+'_detected',[]):
                advar[attr][dvar][irun] = True
    
    for attr in aattrs:
        try:
            if attr in x:
                for a in ['doc', 'unit', 'alias']:
                    x[attr].attrs[a] = aattrs[attr].get(a, '')
        except:
            print('cannot add attrs for {:}'.format(attr))

    for attr in dattrs:
        data = np.array([advar[attr][a] for a in dvars]).T
        x[attr] = (('run', 'dvar'), data)
    
    if save:
        sum_file = '{:}/drop_summary.nc'.format(path)
        x.attrs['file_name'] = sum_file
        print('Saving drop_summary file: {:}'.format(sum_file))
        from .xarray_utils import clean_dataset
        x = clean_dataset(x)
        x.To_netcdf(sum_file, engine='h5netcdf')
        #x.To_netcdf(sum_file, engine='h5netcdf', invalid_netcdf=True)

    return x

def build_drop_stats(x, min_detected=2,  
        alert=True, to_name=None, from_name=None, html_path=None,
        report_name='drop_summary', path=None, engine='h5netcdf'):
    """
    """
    import os
    from requests import post
    from os import environ
    update_url = environ.get('BATCH_UPDATE_URL')

    if isinstance(x, str):
        x = load_exp_sum(x)
   
    exp = x.attrs.get('experiment')
    instrument = x.attrs.get('instrument')
    run = x.attrs.get('run')
    expNum = x.attrs.get('expNum')
    
    if not path:
        path = os.path.join('/reg/d/psdm/',instrument,exp,'results','nc')

    if not os.path.isdir(path):
        os.mkdir(path)
    
    filename = '{:}'.format(report_name)
    h5file = os.path.join(path,report_name+'.nc')
    
    try:
        from .PyDataSource import DataSource
        ds = DataSource(exp=exp,run=run)
        config_info = str(ds.configData.show_info(show_codes=False))
        print(config_info)
        ievt = 0
        evt = next(ds.events)
        # make sure in all detectors have been seen to get full det config
        for det, detector in ds._detectors.items(): 
            while det not in evt._attrs and ievt < min([1000,ds.nevents]):
                evt.next(publish=False, init=False)
                if evt.Evr.eventCodes:
                    ievt += 1
    except:
        config_info = None
        traceback.print_exc('Cannot get data source information for {:} Run {:}'.format(exp,run))

    report_notes = ['Report includes:']
    try:
        from .build_html import Build_html
        from .h5write import runlist_to_str
        b = Build_html(x, h5file=h5file, filename=report_name, path=html_path,
                title=exp+' Drop Summary', subtitle='Drop Summary')
        dattrs = [attr for attr, a in x.data_vars.items() if 'dvar' in a.dims]
        batch_counters = {}
        webattrs = [instrument.upper(), expNum, exp, report_name, 'report.html'] 
        weblink='http://pswww.slac.stanford.edu/experiment_results/{:}/{:}-{:}/{:}/{:}'.format(*webattrs)
        for attr in dattrs:
            inds = (x[attr].to_pandas().sum() >= min_detected)
            attrs = [a for a, val in inds.items() if val]
            gattrs = {}
            for a in attrs:
                dattr = a.split('_')[0]
                if dattr not in gattrs:
                    gattrs[dattr] = []
                gattrs[dattr].append(a)
                
            if attr == 'timing_error':
                if attrs:
                    report_notes.append(' - ALERT:  off-by-one timing errors detected \n {:}'.format([str(a) for a in sorted(gattrs)]))
                    batch_str = ', '.join(['<a href={:}#{:}_data>{:}</a>'.format(weblink, a,a) \
                                            for a in sorted(gattrs)])
                    batch_attr = '<a href={:}>{:}</a>'.format(weblink,'Off-by-one detected')
                    #batch_attr = 'Off-by-one detected'
                    batch_counters[batch_attr] = [batch_str, 'red']
                #    for a in attrs:
                #        report_notes.append('      {:}'.format(a))
                else:
                    report_notes.append(' - No off-by-one timing errors detected')
                    batch_attr = 'Off-by-one'
                    batch_counters[batch_attr] = ['None Detected', 'green']
            else:
                if attrs:
                    report_notes.append(' - {:} detected \n {:}'.format(attr.replace('_',' '), [str(a) for a in sorted(gattrs)]))
                #    for a in attrs:
                #        report_notes.append('      {:}'.format(a))
                else:
                    report_notes.append(' - No {:} detected'.format(attr.replace('_',' ')))

            for alias, attrs in gattrs.items():
                if attr == 'timing_error':
                    tbl_type = alias+'_'+attr
                    attr_cat = 'Alert Off-by-one Error'
                else:
                    tbl_type = attr
                    attr_cat = alias
                doc = ['Runs with {:} for {:} attributes'.format(attr.replace('_',' '), alias)]
                df = x[attr].to_pandas()[attrs]
                for a in attrs:
                    if a in x:
                        da = df[a]
                        name = '_'.join(a.split('_')[1:])
                        aattrs = x[a].attrs
                        runstr = runlist_to_str(da.index[da.values])
                        doc.append(' - {:} runs with {:}: {:} [{:}] \n       [{:}]'.format(df[a].sum(), name, 
                                    aattrs.get('doc',''), aattrs.get('unit',''), runstr))
                howto = ['x["{:}"].to_pandas()[{:}]'.format(attr, attrs)]
                df = df.T.rename({a: '_'.join(a.split('_')[1:]) for a in attrs}).T
                b.add_table(df, attr_cat, tbl_type, tbl_type, doc=doc, howto=howto, hidden=True)
            
        try:
            if config_info:
                report_notes.append('')
                report_notes.append(config_info)
                print(report_notes)
            else:
                print('No config info')
        except:
            traceback.print_exc('Cannot add config information')
            print(config_info)

        b.to_html(h5file=h5file, report_notes=report_notes, show_event_access=True)

        if update_url:
            try:
                print('Setting batch job counter output: {:}'.format(batch_counters))
                post(update_url, json={'counters' : batch_counters})
            except:
                traceback.print_exc('Cannot update batch submission json counter info to {:}'.format(update_url))

        if not to_name:
            if isinstance(alert, list) or isinstance(alert, str):
                to_name = alert
            else:
                to_name = from_name
            
        if alert and to_name is not None:
            try:
                alert_items = {name.lstrip('Alert').lstrip(' '): item for name, item in b.results.items() \
                        if name.startswith('Alert')}
                if alert_items:
                    from .psmessage import Message

                    message = Message('Alert {:} Off-by-one Errors'.format(exp))
                    message('')
                    for alert_type, item in alert_items.items():
                        message(alert_type)
                        message('='*len(message._message[-1]))
                        for name, tbl in item.get('table',{}).items():
                            df = tbl.get('DataFrame')
                            message('')
                            message('* '+name)
                            for a in df:
                                da = df[a]
                                message('  -{:} -- {:} Errors:'.format(a,da.sum()))
                                message('       [{:}]'.format(runlist_to_str(da.index[da.values])))
                    
                    message('')
                    message('See report:')
                    message(b.weblink)
                    print('Sending message to {:} from {:}'.format(to_name, from_name))
                    print(str(message))
                    
                    if to_name:
                        message.send_mail(to_name=to_name, from_name=from_name)

                    return message
            
            except:
                traceback.print_exc('Cannot send alerts: \n {:}'.format(str(message)))

        return b
    
    except:
        traceback.print_exc('Cannot build beam drop stats')

def get_beam_stats(exp, run, default_modules={}, 
        flatten=True, refresh=True,
        drop_code='ec162', drop_attr='delta_drop', nearest=None, drop_min=3,  
        pulse=None, gasdetcut_mJ=None,
        report_name=None, path=None, engine='h5netcdf', 
        wait=None, timeout=False,
        add_stats=False,
        pvdict={},
        **kwargs):
    """
    Get drop shot statistics to detected dropped shots and beam correlated detectors.

    Parameters
    ----------
    drop_code : str
        Event code that indicates a dropped shot with no X-rays
    
    drop_attr: str
        Name for number of shots off from dropped shot.
    
    nearest : int
        Number of nearest events to dropped shot to analayze 
        [default=5 unless rate=10 Hz, then nearest=10]

    drop_min : int
        Minumum number of dropped shots to perform analysis


    """
    from .xarray_utils import set_delta_beam
    from .xarray_utils import clean_dataset
    from .xarray_utils import merge_fill
    from . import PyDataSource
    import xarray as xr
    import numpy as np
    import pandas as pd
    import time
    import os
    from requests import post
    from os import environ
    update_url = environ.get('BATCH_UPDATE_URL')
    
    time0 = time.time() 
    logger = logging.getLogger(__name__)
    logger.info(__name__)
    
    ds = PyDataSource.DataSource(exp=exp, run=run, default_modules=default_modules, wait=wait, timeout=timeout)
#    if update_url:
#        batch_counters = {'Status': ['Loading small data...','yellow']}
#        post(update_url, json={'counters' : batch_counters})
    
    xsmd = load_small_xarray(ds, refresh=refresh, path=path)
  
    try:

        configData = ds.configData
        tstart = configData._tstart
        tend = configData._tend
        fields={
                'description':            ('DESC', 'Description'), 
                'slew_speed':             ('VELO', 'Velocity (EGU/s) '),
                'acceleration':           ('ACCL', 'acceleration time'),
                'step_size':              ('RES',  'Step Size (EGU)'),
                'encoder_step':           ('ERES', 'Encoder Step Size '),
                'resolution':             ('MRES', 'Motor Step Size (EGU)'),
                'high_limit':             ('HLM',  'User High Limit'),
                'low_limit':              ('LLM',  'User Low Limit'),
                'units':                  ('EGU',  'Units'),
    #            'device_type':            ('DTYP', 'Device type'), 
    #            'record_type':            ('RTYP', 'Record Type'), 
                }

        time_last = time.time()
        trans_pvs = [a for a in _transmission_pvs.get('FEE', {}) if a.endswith('_trans')]
        trans_pvs += [a for a in _transmission_pvs.get(ds.instrument, {}) if a.endswith('_trans')]
        trans3_pvs = [a for a in _transmission_pvs.get('FEE', {}) if a.endswith('_trans3')]
        trans3_pvs += [a for a in _transmission_pvs.get(ds.instrument, {}) if a.endswith('_trans3')]
        pvdict.update(**_transmission_pvs.get('FEE',{}))
        pvdict.update(**_transmission_pvs.get(ds.instrument,{}))
        
        pvs = {alias: pv for alias, pv in pvdict.items() if configData._in_archive(pv)} 
        
        data_arrays = {} 
        data_fields = {}
 
        for alias, pv in pvs.items():
            data_fields[alias] = {}
            dat = configData._get_pv_from_arch(pv, tstart, tend)
            if not dat:
                print('WARNING:  {:} - {:} not archived'.format(alias, pv))
                continue
            
            try:
                attrs = {a: dat['meta'].get(val) for a,val in meta_attrs.items() if val in dat['meta']}
                for attr, item in fields.items():  
                    try:
                        field=item[0]
                        pv_desc = pv.split('.')[0]+'.'+field
                        if configData._in_archive(pv_desc):
                            desc = configData._get_pv_from_arch(pv_desc)
                            if desc:
                                vals = {}
                                fattrs = attrs.copy()
                                fattrs.update(**desc['meta'])
                                fattrs['doc'] = item[1]
                                val = None
                                # remove redundant data
                                for item in desc['data']:
                                    newval =  item.get('val')
                                    if not val or newval != val:
                                        val = newval
                                        if six.PY3:
                                            vt = np.datetime64(int(item['secs']*1e9+item['nanos']), 'ns')
                                        else:
                                            vt = np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns')
                                        vals[vt] = val
                               
                                data_fields[alias][attr] = xr.DataArray(list(vals.values()), 
                                                                coords=[list(vals.keys())], dims=['time'], 
                                                                name=alias+'_'+attr, attrs=fattrs) 
                                attrs[attr] = val
         
                    except:
                        traceback.print_exc()
                        print(('cannot get meta for', alias, attr))
                        pass
                vals = [item['val'] for item in dat['data']]
                if not vals:
                    print('No Data in archive for  {:} - {:}'.format(alias, pv))
                    continue

                doc = attrs.get('description','')
                units = attrs.get('units', '')
                time_next = time.time()
              
                try:
                    if isinstance(vals[0],str):
                        vals = np.array(vals, dtype=str)
                    else:
                        if six.PY3:
                            times = [np.datetime64(int(item['secs']*1e9+item['nanos']), 'ns') for item in dat['data']]
                        else:
                            times = [np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns') for item in dat['data']]
                        dfs = pd.Series(vals, times).sort_index()
                        dfs = dfs[~dfs.index.duplicated()]
                        dfs = dfs[~(dfs.diff()==0)]
                        vals = dfs.values
                        dfs = dfs.to_xarray().rename({'index': 'time'})
                        data_arrays[alias] = dfs 
                        data_arrays[alias].name = alias
                        data_arrays[alias].attrs = attrs
                
                except:
                    traceback.print_exc()
                    print(('Error loadinig', alias))

                try:
                    print('{:8.3f} {:28} {:8} {:10.3f} {:4} {:20} {:}'.format(time_next-time_last, \
                                    gias, len(vals), np.array(vals).mean(), units, doc, pv))
                except:
                    print('{:8.3f} {:28} {:8} {:>10} {:4} {:20} {:}'.format(time_next-time_last, \
                                    alias, len(vals), vals[0], units, doc, pv))
            
            except:
                traceback.print_exc()
                print(('Error loading', alias))

        xdata = xr.merge(list(data_arrays.values()))
        if trans_pvs:
            da = xdata.reset_coords()[trans_pvs].to_array() 
            xdata['trans'] = (('time'), da.prod(dim='variable'))
            xdata['trans'].attrs['doc'] = 'Total transmission: '+'*'.join(trans_pvs)

        print(xdata)
        print(xsmd)

        xdata = merge_fill(xsmd, xdata)
        print(xdata)

    except:
        traceback.print_exc()
        print('Could not load transmission pvs')
    
    try:
        #Save sources, eventCodes and ScanData xarray datasets
        ds.configData.save_configData()
        print('')
        print('Writing Scan Data')
        print(ds.configData.ScanData.dataset)
    except:
        traceback.print_exc('Cannot write scan information for {:}'.format(ds))

    nevents = ds.nevents
    try:
        _seq_evtCodes = list(range(67,99))+list(range(167,199))+list(range(201,217))
        code_stats = [int(a.lstrip('ec')) for a in xsmd.coords if a.startswith('ec') \
                and int(a.lstrip('ec')) in _seq_evtCodes \
                and xsmd[a].sum()<nevents]
        code_stats.append(162)
        code_stats.append(-162)
    except:
        code_stats = []
    
    if drop_code not in xsmd or not xsmd[drop_code].values.any():
        logger.info('Skipping beam stats analysis for {:}'.format(ds))
        logger.info('  -- No {:} present in data'.format(drop_code))
        try:
            batch_counters = {'Warning': ['No dropped shots present in {:} events'.format(nevents), 'red']}
            print('Setting batch job counter output: {:}'.format(batch_counters))
            if update_url:
                post(update_url, json={'counters' : batch_counters})
        except:
            traceback.print_exc('Cannot update batch submission json counter info to {:}'.format(update_url))

        return xsmd

    set_delta_beam(xsmd, code=drop_code, attr=drop_attr)
    if not nearest:
        # auto set number of nearest shots to include at least two on
        # each side for up to 60 Hz data rates.
        try:
            df_drop = xsmd[drop_attr].dropna(dim='time').to_pandas()
            drop_mode = df_drop.diff().mode()[0]
            if drop_mode > 3 and drop_mode < 13:
                nearest = min([int(drop_mode*3),12])
            else:
                nearest = 5
        except:
            nearest=5
        
        xdrop = xsmd.where(abs(xsmd[drop_attr]) <= nearest, drop=True)
        if xdrop[drop_code].values.all():
            nearest=12
            xdrop = xsmd.where(abs(xsmd[drop_attr]) <= nearest, drop=True)
            if xdrop[drop_code].values.all():
                try:
                    ntimes = xdrop.time.size
                    ndrop = int(np.sum(xdrop.get(drop_code)))
                except:
                    traceback.print_exc('Cannot select dropped events')
                    ndrop = 0

                logger.info('Skipping beam stats analysis for {:}'.format(ds))
                logger.info('  -- only {:} events with {:} in {:} events'.format(ndrop, drop_code, nevents))
                batch_counters = {'Warning': 
                        ['Low event rate: {:} dropped shots present in {:} events'.format(ndrop, nevents), 'red']}
                print('Setting batch job counter output: {:}'.format(batch_counters))
                if update_url:
                    post(update_url, json={'counters' : batch_counters})
                return xdrop 
    
    if not pulse:
        pulse = 'FEEGasDetEnergy_f_21_ENRC'
        if pulse not in xdrop or xdrop[pulse].sum() == 0:
            pulse = 'FEEGasDetEnergy_f_11_ENRC'
        if pulse not in xdrop or xdrop[pulse].sum() == 0:
            pulse = None
    
    if gasdetcut_mJ is not False:
        gas_attr = 'FEEGasDetEnergy_f_11_ENRC'
        try:
            # determine gasdet cut from dropped shots
            gasdet_stats = xsmd[gas_attr].where(xsmd[drop_code]).dropna(dim='time').to_pandas().describe(percentiles=[0.05,0.95])
            gasdetcut_mJ = min([gasdet_stats['max'],gasdet_stats['std']+gasdet_stats['95%']])
        except:
            gasdetcut_mJ = 0.1
            traceback.print_exc('Cannot auto set gasdetcut_mJ - default = {:} mJ'.format(gasdetcut_mJ))
       
        try:
            xsmd.attrs['gasdetcut_mJ'] = gasdetcut_mJ
            xsmd.attrs['gasdet_attr'] = gas_attr
            xsmd.coords['Gasdet_cut'] = xsmd[gas_attr] > gasdetcut_mJ
            xsmd.coords['Gasdet_cut'].attrs['doc'] = "Gas detector cut.  Gasdet_pre_atten > {:} mJ".format(gasdetcut_mJ)
            try:
                nlowshots = ((xsmd['Gasdet_cut'] == False) & (xsmd[drop_code] == False)).values.sum()
                nshots = (xsmd[drop_code] == False).values.sum()
                xsmd.attrs['lowbeam_fraction'] = nlowshots/float(nshots) 
            except:
                traceback.print_exc('Cannot calculate lowbeam_fraction for threshold gasdetcut_mJ {:} mJ'.format(gasdetcut_mJ))
        except:
            traceback.print_exc('Cannot make Gasdet_cut with threshold gasdetcut_mJ {:} mJ'.format(gasdetcut_mJ))
        
        try:
            xsmd.coords['XrayOff'] = (xsmd[drop_code] == True) | (xsmd['Gasdet_cut'] == False)
            xsmd.coords['XrayOff'].attrs['doc'] = 'Xray Off for events with {:} or lowbeam'.format(drop_code)
            xsmd.coords['XrayOn'] = (xsmd[drop_code] == False) & (xsmd['Gasdet_cut'] == True)
            xsmd.coords['XrayOn'].attrs['doc'] = 'Xray On for events without {:} and no lowbeam'.format(drop_code)
        except:
            traceback.print_exc('Cannot set XrayOn/XrayOff with threshold gasdetcut_mJ {:} mJ'.format(gasdetcut_mJ))
   
    try:
        if 'Gasdet_cut' in xsmd.coords:
            drop_select = (abs(xsmd[drop_attr]) <= nearest) & (xsmd['Gasdet_cut'] | xsmd[drop_code])
        else:
            drop_select = abs(xsmd[drop_attr]) <= nearest
    except:
        drop_select = abs(xsmd[drop_attr]) <= nearest
        traceback.print_exc('Cannot select smd using threshold gasdetcut_mJ {:} mJ'.format(gasdetcut_mJ))
            
    try:
        xdrop = xsmd.where(drop_select, drop=True)
        ntimes = xdrop.time.size
        ndrop = int(np.sum(xdrop.get(drop_code)))
    except:
        traceback.print_exc('Cannot select dropped events')
        ndrop = 0

    if ndrop < drop_min:
        logger.info('Skipping beam stats analysis for {:}'.format(ds))
        logger.info('  -- only {:} events with {:} in {:} events'.format(ndrop, drop_code, nevents))
        try:
            batch_counters = {'Warning': ['Only {:} dropped shots present in {:} events'.format(ndrop, nevents), 'red']}
            print('Setting batch job counter output: {:}'.format(batch_counters))
            if update_url:
                post(update_url, json={'counters' : batch_counters})
        except:
            traceback.print_exc('Cannot update batch submission json counter info to {:}'.format(update_url))
        
        return xsmd

    dets = {}
    flatten_list = []
    flatten_channels = {}
    area_dets = []
    wf_dets =[]
    for det, detector in ds._detectors.items():
        methods = {}
        try:
            det_info = detector._xarray_info.get('dims',{})
            if detector._pydet is None: 
                logger.info(det, 'pydet not implemented')
            elif detector._pydet.__module__ == 'Detector.AreaDetector':
                srcstr = detector._srcstr 
                srcname = srcstr.split('(')[1].split(')')[0]
                devName = srcname.split(':')[1].split('.')[0]
                # for not just use rawsum for 'Epix10ka' and  'Jungfrau' until full
                # development of gain switching in Detector module
                if devName.startswith('Opal'):
                    method = 'rawsum'
                    name = '_'.join([det, method])
                    methods[name] = method
                    xdrop[name] = (('time'), np.zeros([ntimes]))
                    xdrop[name].attrs['doc'] = '{:} sum of raw data'.format(det)
                    xdrop[name].attrs['unit'] = 'ADU'
                    xdrop[name].attrs['alias'] = det
                    area_dets.append(name)
                    if add_stats:
                        next(detector)
                        ok_stats = detector.add.stats('raw', eventCodes=code_stats)
                        print(('Adding stats for', name, ok_stats))
                
                elif devName in ['Epix10ka']:
                    method = 'rawsum'
                    name = '_'.join([det, method])
                    methods[name] = method
                    xdrop[name] = (('time'), np.zeros([ntimes]))
                    xdrop[name].attrs['doc'] = '{:} sum of raw data'.format(det)
                    xdrop[name].attrs['unit'] = 'ADU'
                    xdrop[name].attrs['alias'] = det
                    area_dets.append(name)
                    if add_stats:
                        next(detector)
                        ok_stats = detector.add.stats('calib', eventCodes=code_stats)
                        print(('Adding stats for', name, ok_stats))
                
                elif devName in ['Jungfrau']:
                    nch = detector.configData.numberOfModules
                    if nch > 1:
                        method = 'sectors'
                        name = '_'.join([det, 'rawsum'])
                        methods[name] = method
                        xdrop[name] = (('time', det+'_ch'), np.zeros([ntimes, nch]))
                        for ich in range(nch):
                            area_dets.append('{:}_ch{:}'.format(name,ich))
                    else:
                        method = 'sector'
                        name = '_'.join([det, 'rawsum'])
                        methods[name] = method
                        xdrop[name] = (('time'), np.zeros([ntimes]))
                        area_dets.append(name)
                    xdrop[name].attrs['doc'] = '{:} sum of raw data'.format(det)
                    xdrop[name].attrs['unit'] = 'ADU'
                    xdrop[name].attrs['alias'] = det
                    if add_stats:
                        next(detector)
                        ok_stats = detector.add.stats('calib', eventCodes=code_stats)
                        print(('Adding stats for', name, ok_stats))
                
                else:
                    method = 'count'
                    name = '_'.join([det, method])
                    methods[name] = method
                    xdrop[name] = (('time'), np.zeros([ntimes]))
                    xdrop[name].attrs['doc'] = '{:} sum of pedestal corrected data'.format(det)
                    xdrop[name].attrs['unit'] = 'ADU'
                    xdrop[name].attrs['alias'] = det
                    area_dets.append(name)
                    if add_stats:
                        next(detector)
                        ok_stats = detector.add.stats('corr', eventCodes=code_stats)
                        print(('Adding stats for', name, ok_stats))

            elif detector._pydet.__module__ == 'Detector.GenericWFDetector':
                srcstr = detector._srcstr 
                srcname = srcstr.split('(')[1].split(')')[0]
                devName = srcname.split(':')[1].split('.')[0]
                if devName == 'Wave8':
                    #nch = detector.configData.NChannels 
                    nch = 8
                    method = 'wave8_height'
                    name = det
                    methods[name] = method
                    xdrop[name] = (('time', det+'_ch',), np.zeros([ntimes, nch]))
                    xdrop[name].attrs['doc'] = 'BeamMonitor {:}'.format(method.replace('_',' '))
                    xdrop[name].attrs['unit'] = 'V'
                    xdrop[name].attrs['alias'] = det
                    wf_dets.append(name)

            elif detector._pydet.__module__ == 'Detector.WFDetector':
                srcstr = detector._srcstr 
                srcname = srcstr.split('(')[1].split(')')[0]
                devName = srcname.split(':')[1].split('.')[0]
                if devName == 'Acqiris':
                    nch = detector.configData.nbrChannels 
                    for method in ['peak_height', 'peak_time']:
                        name = '_'.join([det, method])
                        methods[name] = method
                        xdrop[name] = (('time', det+'_ch',), np.zeros([ntimes, nch]))
                        xdrop[name].attrs['doc'] = 'Acqiris {:}'.format(method.replace('_',' '))
                        xdrop[name].attrs['unit'] = 'V'
                        xdrop[name].attrs['alias'] = det
                        wf_dets.append(name)
                
                elif devName == 'Imp':
                    nch = 4
                    method = 'amplitudes'
                    name = '_'.join([det, method])
                    methods[name] = method
                    xdrop[name] = (('time', det+'_ch'), np.zeros([ntimes, nch]))
                    attr = 'waveform'
                    xdrop[name].attrs['doc'] = 'IMP filtered amplitudes'
                    xdrop[name].attrs['unit'] = 'V'
                    xdrop[name].attrs['alias'] = det
                    wf_dets.append(name)
 
            elif detector._pydet.__module__ == 'Detector.UsdUsbDetector':
                srcstr = detector._srcstr 
                srcname = srcstr.split('(')[1].split(')')[0]
                devName = srcname.split(':')[1].split('.')[0]
                if devName == 'USDUSB':
                    method = 'encoder_values'
                    name = '_'.join([det, method])
                    flatten_list.append(name)
            
            elif detector._pydet.__module__ == 'Detector.DdlDetector':
                srcstr = detector._srcstr 
                srcname = srcstr.split('(')[1].split(')')[0]
                print('Setting {:}'.format(srcstr))
                if srcstr.startswith('BldInfo'):
                    if srcname.endswith('BMMON'):
                        method = 'channelValue'
                        name = '_'.join([det, method])
                        flatten_list.append(name)
                        print('Flatten {:}'.format(name))
                        #flatten_channels[name] = range(8,16)
                else:
                    devName = srcname.split(':')[1].split('.')[0]
                    if devName == 'Gsc16ai':
                        method = 'channelValue'
                        name = '_'.join([det, method])
                        flatten_list.append(name)
                        print('Flatten {:}'.format(name))

            elif detector._pydet.__module__ == 'Detector.IpimbDetector':
                pass
            
            else:
                logger.info('{:} Not implemented'.format(det))
                print('{:} Not implemented'.format(det))
        
        except:
            logger.info('Error with config of {:}'.format(det))
            print('Error with config of {:}'.format(det))
        
        if methods:
            dets[det] = methods
            xdrop.attrs['waveform_detectors'] = wf_dets
            xdrop.attrs['area_detectors'] = area_dets

    print('-'*80)
    print(xdrop)
    print('-'*80)
    print(dets)
    print('-'*80)
    ds.reload()

    times = list(zip(xdrop.sec.values,xdrop.nsec.values,xdrop.fiducials.values))
    nupdate = 100
    time_last = time0
    logger.info('Loading: {:}'.format(list(dets.keys())))
    for itime, t in enumerate(times):
        if itime % nupdate == nupdate-1:
            time_next = time.time()
            dtime = time_next-time_last
            evt_info = '{:8} of {:8} -- {:8.3f} sec, {:8.3f} events/sec'.format(itime+1, 
                    ntimes, time_next-time0, nupdate/dtime)
            logger.info(evt_info)
            time_last = time_next 
#            if update_url:
#                batch_counters = {'Status': [evt_info,'yellow']}
#                post(update_url, json={'counters' : batch_counters})

        evt = ds.events.next(t)
        for det, methods in dets.items():
            detector = evt._dets.get(det)
            if detector and detector.sourceData.eventCode in evt.Evr.eventCodes_strict:
                for name, method in methods.items():
                    try:
                        xdrop[name][itime] = globals()[method](detector) 
                    except:
                        xdrop[name][itime] = np.nan 
                        print('Cannot calculate {:} for {:}'.format(method, name))
                        traceback.print_exc('Error with {:} {:}'.format(name, method))
            else:
                for name, method in methods.items():
                    xdrop[name][itime] = np.nan 

    # Flatten waveform data with channels
    if flatten:
        #for name in xdrop.data_vars:
        for det, methods in dets.items():
            for name, method in methods.items():
                if len(xdrop[name].dims) == 2:
                    flatten_list.append(name)
        
        for name in flatten_list:
            if name in xdrop and len(xdrop[name].dims) == 2:
                nch = xdrop[name].shape[1]
                fchans = flatten_channels.get(name, list(range(nch)))
                if nch <= 16 or name in flatten_channels:
                    for ich in fchans:
                        chname = '{:}_ch{:}'.format(name,ich)
                        xdrop[chname] = xdrop[name][:,ich]
                        try:
                            if not xdrop[name].attrs.get('attr'):
                                alias = xdrop[name].attrs.get('alias','')
                                attr = name.lstrip(alias).lstrip('_')
                                xdrop[name].attrs['attr'] = attr 
                        except:
                            pass
                    del xdrop[name]

    if not path:
        path = os.path.join(ds.data_source.res_dir,'nc')
    elif path == 'home':
        path = os.path.join(os.path.expanduser('~'), 'RunSummary', 'nc')

    if not os.path.isdir(path):
        os.mkdir(path)

    if add_stats:
        try:
            ds.save_stats(path=path)
        except:
            traceback.print_exc('Cannot save stats for {:}'.format(ds))

    if not report_name:
        report_name = 'run{:04}_drop_stats'.format(ds.data_source.run)

    h5file = os.path.join(path,report_name+'.nc')
   
    print('Build {:} {:} {:}'.format(exp,run, xdrop))
    b = build_beam_stats(exp=exp, run=run, xdrop=xdrop, xsmd=xsmd, 
            nearest=nearest,
            report_name=report_name, h5file=h5file, path=path, **kwargs)
    
    try:
        xdrop = clean_dataset(xdrop)
        #xdrop.to_netcdf(h5file, engine=engine, invalid_netcdf=True)
        xdrop.to_netcdf(h5file, engine=engine)
        logger.info('Saving file to {:}'.format(h5file))
        print('Saving file to {:}'.format(h5file))
    except:
        traceback.print_exc('Cannot save to {:}'.format(h5file))
   
    return xdrop

def build_beam_stats(exp=None, run=None, 
        xdrop=None, xsmd=None, 
        instrument=None,
        report_name=None, h5file=None, path=None, 
        alert=True, to_name=None, from_name=None, html_path=None,
        make_scatter=False, cut_flag=None,
        nearest=5,
        pulse=None, **kwargs):
    """
    """
    import os
    import re
    import numpy as np
    import pandas as pd
    import xarray as xr
    from requests import post
    from os import environ
    update_url = environ.get('BATCH_UPDATE_URL')
    _seq_evtCodes = list(range(67,99))+list(range(167,199))+list(range(201,217))
    batch_counters = {}

    if not path:
        if not exp:
            exp = str(xdrop.attrs['experiment'])
        if not instrument:
            if xdrop is None:
                instrument = exp[:3]
            else:
                instrument = str(xdrop.attrs['instrument'])

        path = os.path.join('/reg/d/psdm/',instrument,exp,'results','nc')

    if not run:
        if xdrop is None:
            print('Error:  Need to specify exp and run or alternatively xdrop DataSet') 
            return None
        else:
            run = int(xdrop.attrs['run'])

    if not report_name:
        report_name = 'run{:04}_drop_stats'.format(run)
    
    if not h5file:
        h5file = os.path.join(path,report_name+'.nc')
 
    if xdrop is None:
        import xarray as xr
        xdrop = xr.open_dataset(h5file, engine='h5netcdf')

    exp = str(xdrop.attrs.get('experiment'))
    instrument = str(xdrop.attrs.get('instrument'))
    run = int(xdrop.attrs.get('run'))
    expNum = int(xdrop.attrs.get('expNum'))
    webattrs = [instrument.upper(), expNum, exp, report_name, 'report.html'] 
    weblink='http://pswww.slac.stanford.edu/experiment_results/{:}/{:}-{:}/{:}/{:}'.format(*webattrs)

    try:
        from .PyDataSource import DataSource
        ds = DataSource(exp=exp,run=run)
        configData = ds.configData
        config_info = str(ds.configData.show_info(show_codes=False))
        print(config_info)
    except:
        configData = None
        config_info = None
        traceback.print_exc('Cannot get data source information for {:} Run {:}'.format(exp,run))


    try:
        from .build_html import Build_html
   
        b = Build_html(xdrop, h5file=h5file, filename=report_name, path=html_path)
        drop_attr = str(xdrop.attrs.get('drop_attr', 'ec162'))
        if 'XrayOff' not in xdrop and drop_attr in xdrop:
            xdrop.coords['XrayOff'] = (xdrop[drop_attr] == True)
            xdrop.coords['XrayOff'].attrs['doc'] = 'Xray Off for events with {:}'.format(drop_attr)
        if 'XrayOn' not in xdrop and drop_attr in xdrop:
            xdrop.coords['XrayOn'] = (xdrop[drop_attr] == False)
            xdrop.coords['XrayOn'].attrs['doc'] = 'Xray On for events without {:}'.format(drop_attr)

        report_notes = []
        # X-ray Energy and Charge information
        try:
            energy_mean= xdrop.EBeam_ebeamPhotonEnergy.where(xdrop.XrayOn, drop=True).values.mean()
            energy_std= xdrop.EBeam_ebeamPhotonEnergy.where(xdrop.XrayOn, drop=True).values.std()
            report_notes.append('Energy = {:6.1f}+={:5.1f} eV'.format(energy_mean, energy_std))
            charge_mean= xdrop.EBeam_ebeamCharge.where(xdrop.XrayOn, drop=True).values.mean()
            charge_std= xdrop.EBeam_ebeamCharge.where(xdrop.XrayOn, drop=True).values.std()
            report_notes.append('Charge = {:6.3f}+={:5.3f} mA'.format(charge_mean, charge_std))
            gas_attr = xdrop.attrs.get('gasdet_attr')
            if gas_attr and gas_attr in xdrop: 
                pulseE_mean= xdrop[gas_attr].where(xdrop.XrayOn, drop=True).values.mean()
                pulseE_std= xdrop[gas_attr].where(xdrop.XrayOn, drop=True).values.std()
                report_notes.append('PulseE = {:6.3f}+={:5.3f} mJ'.format(pulseE_mean, pulseE_std))
            if xsmd is not None:
                try:
                    nsmd = int(xsmd.time.count()) 
                    nec_smd = int(xsmd['ec162'].sum())
                    ecfrac_smd = nec_smd/float(nsmd)
                    report_notes.append('Dropped = 1/{:}'.format(round(1./ecfrac_smd)))
                    # Check for A-line Kicker rate
                    if 'ec163' in xsmd.coords:
                        nec_smd = int(xsmd['ec163'].sum())
                        if nec_smd > 0:
                            ecfrac_smd = nec_smd/float(nsmd)
                            report_notes.append('ESA Kick = 1/{:}'.format(round(1./ecfrac_smd)))
                except:
                    print('Cannot report drop fraction')
                    traceback.print_exc('Cannot report drop fraction')

            lowbeam_frac = xdrop.attrs.get('lowbeam_fraction', 0)
            gasdetcut_mJ = xdrop.attrs.get('gasdetcut_mJ')
            if lowbeam_frac > 0.01:
                if gasdetcut_mJ: 
                    report_notes.append('LowBeam = {:4.1f} % (<{:6.2} mJ)'.format(lowbeam_frac*100., gasdetcut_mJ))
                else:
                    report_notes.append('LowBeam = {:4.1f} %'.format(lowbeam_frac*100.))
            report_notes.append('')
        except:
            pass
        
        nsteps = 1 
        step_attrs = []
        try:
            if configData is not None:
                scanData = configData.ScanData
                nsteps = scanData.nsteps
                df_steps = pd.DataFrame(scanData.control_values)
                scan_attrs = list(list(df_steps.keys())[df_steps.std() > 0])
                if nsteps > 1:
                    report_notes.append('Scan Steps = {:}'.format(nsteps))
                    report_notes.append('Scan Variables = {:}'.format(scan_attrs))
                    report_notes.append('')
            
        except:
            traceback.print_exc('Cannot get scan information for {:}'.format(ds))

        # Event code flags
        try:
            code_kicker = ['ec162']
            if 'ec163' in xdrop.coords:
                code_kicker.append('ec163')
            nevents = int(xdrop.time.count()) 
            code_flags = [a for a in xdrop.coords if a.startswith('ec') \
                    and int(a.lstrip('ec')) in _seq_evtCodes \
                    and xdrop[a].sum()<nevents]
            if xsmd is not None:
                nsmd = int(xsmd.time.count()) 
                report_str = 'Event Types in All Run Events: {:} Total'.format(nsmd)
                report_notes.append(report_str)
                report_notes.append('-'*40)
                for code in code_flags+code_kicker:
                    try:
                        doc = str(xsmd[code].attrs.get('doc','')).lstrip('event code for ')
                        nec_smd = int(xsmd[code].sum())
                        if nec_smd > 0:
                            ecfrac_smd = nec_smd/float(nsmd)
                            report_str = '{:7} {:6} - {:5.1f}%  {:20}'.format(nec_smd, code, ecfrac_smd*100., doc)
                            report_notes.append(report_str)
                    except:
                        print('Error in report of {:} rate'.format(code))
                report_notes.append('')
            
            report_str = 'Event Types in Dropped Shot Analysis: {:} of {:} ({:5.1f}%)'.format(nevents, 
                            nsmd, float(nevents)/float(nsmd)*100.)
            report_notes.append(report_str)
            report_notes.append('-'*40)
            for code in code_flags+code_kicker:
                try:
                    doc = str(xdrop[code].attrs.get('doc','')).lstrip('event code for ')
                    nec = int(xdrop[code].sum())
                    if nec > 0:
                        ecfrac = nec/float(nevents)
                        report_str = '{:7} {:6} - {:5.1f}%  {:20}'.format(nec, code, ecfrac*100., doc)
                        report_notes.append(report_str)
                except:
                    print('Error in report of {:} rate'.format(code))
            report_notes.append('')
            try:
                if code_flags:
                    flag_names = []
                    for code in code_flags:
                        flag_name = str(xdrop[code].attrs.get('doc','')).lstrip('event code for ')
                        flag_name = re.sub('-|:|\.| ','_', flag_name).replace('"','').replace('__','_').replace('__','_')
                        xdrop.coords[flag_name] = xdrop[code]
                        flag_names.append(flag_name)
                    if len(code_flags) > 1:
                        df_codes = xdrop.reset_coords()[flag_names].to_dataframe().sum()
                        if df_codes.sum() == nevents:
                            cut_flag = df_codes.argmax()
                            flag_inds = list(range(len(flag_names)))
                            xdrop.coords['tag'] = (['time'], np.zeros(nevents, dtype=int))
                            #xdrop.coords['tag_name'] = (('tag'), flag_names)
                            xdrop.attrs['tag_names'] = flag_names
                            for i in flag_inds:
                                xdrop.coords['tag'][xdrop[flag_names[i]] == 1] = i
            except:
                print('Cannot process code_flags {:}'.format(code_flags))
                traceback.print_exc('Cannot process code_flags {:}'.format(code_flags))
        except:
            print('Cannot detect code_flags')
            traceback.print_exc('Cannot detect code_flags')

        # Timing and config alerts
        try:
            configCheck = ds.configData.configCheck
            report_config = False
            det_errors = [] 
            if configCheck.alerts:
                det_errors += list(configCheck.alerts.keys())
#                report_config = True
#                batch_attr = '<a href={:}#{:}>{:}</a>'.format(weblink, 'Report_Notes', 'Timing Config Alerts')
#                batch_str = ','.join(set(configCheck.alerts.keys()))
#                batch_counters[batch_attr] = [batch_str, 'red']
            if configCheck.warnings:
                det_errors += list(configCheck.warnings.keys())
#                report_config = True
#                batch_str = ','.join(set(configCheck.warnings.keys()))
#                batch_attr = '<a href={:}#{:}>{:}</a>'.format(weblink, 'Report_Notes', 'Timing Config Warnings')
#                batch_counters[batch_attr] = [batch_str, 'red']
            #if report_config:
            if det_errors:
                batch_str = ', '.join(set(det_errors))
                batch_attr = '<a href={:}#{:}>{:}</a>'.format(weblink, 'Report_Notes', 'Timing Config Errors')
                batch_counters[batch_attr] = [batch_str, 'red']
                report_notes.append('Detector Timing and Config Errors:')
                report_notes.append('----------------------------------')
                report_notes.append(str(configCheck.show_info()))
                report_notes.append('')
                report_str='Link to Confluence Detector timing settings:'
                report_link='https://confluence.slac.stanford.edu/display/PCDS/Detector+timing+settings'
                report_notes.append('<a href={:}>{:}</a>'.format(report_link, report_str))
                report_notes.append('')
                print(report_notes)
        except:
            print('Cannot check timing and config alerts')

        report_notes.append('Report includes:')
      
       
        b._xstats = b.add_delta_beam(pulse=pulse, cut=cut_flag, nearest=nearest)
        corr_attrs = [str(aa) for aa in xdrop.attrs.get('beam_corr_detected')]
        print('Beam Correlations Detected {:}'.format(corr_attrs))
        if corr_attrs:
            pulse = str(xdrop.attrs.get('beam_corr_attr'))
            if not pulse:
                pulse = 'FEEGasDetEnergy_f_21_ENRC'
                if pulse not in xdrop or xdrop[pulse].sum() == 0:
                    pulse = 'FEEGasDetEnergy_f_11_ENRC'
                if pulse not in xdrop or xdrop[pulse].sum() == 0:
                    pulse = None
            
            if pulse:
                corr_attrs = [a for a in corr_attrs if a == pulse or not a.startswith('FEEGasDetEnergy')]
                corr_attrs.append(pulse)
                if 'PhaseCavity_charge1' in corr_attrs and 'PhaseCavity_charge2' in corr_attrs:
                    corr_attrs.remove('PhaseCavity_charge2')
                print('Adding Detector Beam Correlations {:}'.format(corr_attrs))
                try:
                    b.add_detector(attrs=corr_attrs, catagory=' Beam Correlations', confidence=0.3,
                        cut='XrayOn',
                        make_timeplot=False, make_histplot=False, make_table=False, 
                        make_scatter=make_scatter)
                except:
                    traceback.print_exc('Cannot make beam correlations {:}'.format(corr_attrs))

        #batch_counters['report'] = ['See <a href={:}>Run{:04} Report</a>'.format(weblink,run),'blue']

        try:
            for alias in ds.configData._config_srcs:
                config_data = getattr(ds.configData, alias)
                b.add_textblock(str(config_data.show_info()),alias,'config',alias+'_config')
        except:
            print('Cannot add detector configurations')
 
        try:
            offbyone_detectors = list(sorted(set([str(xdrop[a].attrs.get('alias')) for a in xdrop.timing_error_detected])))
            if offbyone_detectors:
                #batch_str = ', '.join(['<a href={:}#{:}_data>{:}</a>'.format(weblink, attr,attr) \
                #                        for attr in offbyone_detectors])
                batch_str = ', '.join([attr for attr in offbyone_detectors])
                batch_attr = '<a href={:}#{:}_data>{:}</a>'.format(weblink, "%20Alert%20Timing%20Error", 'Off-by-one detected')
                batch_counters[batch_attr] = [batch_str, 'red']
                report_notes.append(' - Off-by-one detected: '+ str(offbyone_detectors))
                for a in sorted(xdrop.timing_error_detected):
                    report_notes.append('    + '+a)
        except:
            offbyone_detectors = []

        try:
            drop_detectors = list(sorted(set([str(xdrop[a].attrs.get('alias')) for a in xdrop.drop_shot_detected])))
            if drop_detectors:
                batch_str = ', '.join([attr for attr in drop_detectors])
                #batch_str = ', '.join(['<a href={:}#{:}_data>{:}</a>'.format(weblink, attr,attr) \
                #                        for attr in drop_detectors])
                batch_counters['Dropped shot detected'] = [batch_str, 'green']
                report_notes.append(' - Dropped shot detected: '+str(drop_detectors))
                for a in sorted(xdrop.drop_shot_detected):
                    report_notes.append('    + '+a)
        except:
            drop_detectors = []


        try:
            beam_detectors = list(sorted(set([str(xdrop[a].attrs.get('alias')) for a in xdrop.beam_corr_detected])))
            if beam_detectors:
                #beam_detectors = ['<a href={:}#{:}_data>{:}</a>'.format(weblink, a, a) for a in beam_detectors]
                batch_str = ', '.join([attr for attr in beam_detectors])
                #batch_str = ', '.join(['<a href={:}#{:}_data>{:}</a>'.format(weblink, attr,attr) \
                #                        for attr in beam_detectors])
                batch_attr = '<a href={:}#{:}_data>{:}</a>'.format(weblink, "%20Beam%20Correlations", 'Beam Correlated detected')
                batch_counters[batch_attr] = [batch_str, 'green']
                report_notes.append(' - Beam correlated detected: '+str(beam_detectors))
                for a in sorted(xdrop.beam_corr_detected):
                    report_notes.append('    + '+a)
        except:
            beam_detectors = []

        try:
            warning_detectors = list(sorted(set([str(xdrop[a].attrs.get('alias')) for a in xdrop.beam_warning_detected])))
            if warning_detectors:
                batch_str = ', '.join([attr for attr in warning_detectors])
                #batch_str = ', '.join(['<a href={:}#{:}_data>{:}</a>'.format(weblink, attr,attr) \
                #                        for attr in warning_detectors])
                #batch_attr = '<a href={:}#{:}_data>{:}</a>'.format(weblink, "%20Beam%20Warnings", 'Beam Warnings detected')
                #batch_counters[batch_attr] = [batch_str, 'green']
                batch_counters['Beam warnings'] = [batch_str, 'green']
                report_notes.append(' - Beam Warnings: '+str(warning_detectors))
                for a in sorted(xdrop.beam_warning_detected):
                    report_notes.append('    + '+a)
        except:
            warning_detectors = []

        try:
            if nsteps > 1:
                batch_attr = '<a href={:}#{:}_data>{:}</a>'.format(weblink, "%20Scan", '{:} Scan Steps'.format(nsteps))
                batch_str = 'Scan Variables = {:}'.format(scan_attrs)
                batch_counters[batch_attr] = [batch_str, 'purple']
                try:
                    b.add_textblock(str(configData.ScanData.show_info()), ' Scan', 'Step Info')
                except:
                    traceback.print_exc('Cannot add Scan Data')
        except:
            pass

        if config_info:
            report_notes.append('')
            report_notes.append(config_info)
            print(report_notes)
        else:
            print('No config info')

        if b.results:
            # only make reports if not empty
            b.to_html(h5file=h5file, report_notes=report_notes, show_event_access=True)

        if update_url:
            try:
                print('Setting batch job counter output: {:}'.format(batch_counters))
                post(update_url, json={'counters' : batch_counters})
            except:
                traceback.print_exc('Cannot update batch submission json counter info to {:}'.format(update_url))

        if not to_name:
            if isinstance(alert, list) or isinstance(alert, str):
                to_name = alert
            else:
                to_name = from_name
            
        alert_items = {name.lstrip(' ').lstrip('Alert').lstrip(' '): item for name, item in b.results.items() \
                if name.lstrip(' ').startswith('Alert')}
       
        if alert_items and offbyone_detectors != ['EBeam']:
            if alert and to_name is not 'None':
                message = None
                try:
                    from .psmessage import Message
                    import pandas as pd
                    message = Message('Alert {:} Run {:}: {:}'.format(exp,run, ','.join(alert_items.keys())))
                    event_times = pd.to_datetime(b._xdat.time.values)
                    begin_time = event_times.min()
                    end_time = event_times.max()
                    run_time = (end_time-begin_time).seconds
                    minutes,fracseconds = divmod(run_time,60)

                    message('- Run Start:  {:}'.format(begin_time.ctime()))
                    message('- Run End:    {:}'.format(end_time.ctime()))
                    message('- Duration:   {:} seconds ({:02.0f}:{:02.0f})'.format(run_time, minutes, fracseconds))
                    message('- Total events: {:}'.format(len(event_times) ) )
                    message('')
 
                    for alert_type, item in alert_items.items():
                        message(alert_type)
                        message('='*len(message._message[-1]))
                        for name, tbl in item.get('figure',{}).items():
                            doc = tbl.get('doc')
                            message('* '+doc[0])
                            message('   - '+doc[1])
                    
                    message('')
                    message('See report:')
                    message(b.weblink)
                    
                    print(str(message))
                    try:
                        print('Sending message to {:} from {:}'.format(to_name, from_name))
                        print(str(message))
                        message.send_mail(to_name=to_name, from_name=from_name)
                    except:
                        traceback.print_exc('ERROR Sending message to {:} from {:}'.format(to_name, from_name))

                except:
                    traceback.print_exc('Cannot send alerts: \n {:}'.format(str(message)))

    except:
        traceback.print_exc('Cannot build drop report')

    return b


def load_small_xarray(ds, path=None, filename=None, refresh=None, 
        engine='h5netcdf', **kwargs):
    """Load small xarray Dataset with PyDataSource.DataSource. 

    Parameters
    ----------
    refresh : bool
        True = reload
        False = load from file, return None if file does not exist
        None = load from file or make new if file does not exist [default] 
         
    """
    import xarray as xr
    import glob
    import os
    if not path:
        path = os.path.join(ds.data_source.res_dir,'nc')
    elif path == 'home':
        path = os.path.join(os.path.expanduser('~'), 'RunSummary', 'nc')

    if not filename:
        filename = 'run{:04}_smd.nc'.format(ds.data_source.run)

    smd_file = os.path.join(path, filename)
    file_available = os.path.isfile(smd_file)
    if refresh is None and not file_available: 
        refresh = True

    if refresh:
        return make_small_xarray(ds, path=path, filename=filename, **kwargs)
    elif file_available:
        return xr.open_dataset(os.path.join(path,filename), engine='h5netcdf')
    else:
        return None

def make_small_xarray(self, auto_update=True,
        add_dets=True, add_counts=False, add_1d=True,
        ignore_unused_codes=True,
        ignore_attrs=['timestamp','numChannels','digital_in'],
        drop_code='ec162', drop_attr='delta_drop', 
        path=None, filename=None, save=True, engine='h5netcdf', 
        make_summary=True,
        nevents=None):
    """Make Small xarray Dataset.
    Parameters
    ----------
    ignore_unused_codes : bool
        If true drop unused eventCodes except drop_code [default=True]
    """
    from .xarray_utils import set_delta_beam
    from .xarray_utils import clean_dataset
    from .xarray_utils import to_summary
    import numpy as np
    import pandas as pd
    import time
    import os
    time0 = time.time() 
    self.reload()
    if not nevents:
        nevents = self.nevents
    
    logger = logging.getLogger(__name__)
    cnames = {code: 'ec{:}'.format(code) for  code in self.configData._eventcodes}
    try:
        code = int(drop_code[2:])
        if code not in cnames:
            cnames[code] = 'ec{:}'.format(code) 
    except:
        traceback.print_exc('Cannot add drop_code = '.format(code))
    data = {cnames[code]: np.zeros(nevents, dtype=bool) for code in cnames}
    #data = {cnames[code]: np.zeros(nevents, dtype=bool) for code in self.configData._eventcodes}
    data['step'] = np.zeros(nevents, dtype=int)
    data1d = {}
    dets1d = {}
    coords = ['fiducials', 'sec', 'nsec']
    dets = {'EventId': {'names': {a:a for a in coords}}}
    coords += list(data.keys())
    for det, attrs in dets.items():
        for name, attr in attrs['names'].items():
            data.update({name: np.zeros(nevents, dtype=int)})
    
    ievt = 0
    if add_dets:
        evt = next(self.events)
        for det, detector in self._detectors.items(): 
            #srcstr = detector._srcstr 
            #srcname = srcstr.split('(')[1].split(')')[0]
            #devName = srcname.split(':')[1].split('.')[0]
            data.update({det+'_present': np.zeros(nevents, dtype=bool)})
            
            if det == 'EBeam':
                attrs = ['ebeamCharge', 'ebeamDumpCharge', 'ebeamEnergyBC1', 'ebeamEnergyBC2', 
                         'ebeamL3Energy', 'ebeamLTU250', 'ebeamLTU450', 'ebeamLTUAngX', 'ebeamLTUAngY', 
                         'ebeamLTUPosX', 'ebeamLTUPosY', 'ebeamPhotonEnergy', 
                         'ebeamPkCurrBC1', 'ebeamPkCurrBC2', 'ebeamUndAngX', 'ebeamUndAngY', 
                         'ebeamUndPosX', 'ebeamUndPosY', 'ebeamXTCAVAmpl', 'ebeamXTCAVPhase']
                attrs = {det+'_'+attr: attr for attr in attrs}
                #attrs = {'EBeam_'+attr.lstrip('ebeam'): attr for attr in attrs}
                dets.update({'EBeam': {'names': attrs}})
                for name, attr in attrs.items():
                    data.update({name: np.zeros(nevents)})
            elif det == 'FEEGasDetEnergy':
                attrs = ['f_11_ENRC', 'f_12_ENRC', 'f_21_ENRC', 'f_22_ENRC', 'f_63_ENRC', 'f_64_ENRC']
                attrs = {det+'_'+attr: attr for attr in attrs}
                #attrs = {'GasDet_'+''.join(attr.split('_')[0:2]): attr for attr in attrs}
                dets.update({'FEEGasDetEnergy': {'names': attrs}})
                for name, attr in attrs.items():
                    data.update({name: np.zeros(nevents)})
            elif det == 'PhaseCavity':
                attrs = ['charge1', 'charge2', 'fitTime1', 'fitTime2']
                #attrs = {'PhaseCavity_'+attr: attr for attr in attrs}
                attrs = {det+'_'+attr: attr for attr in attrs}
                dets.update({'PhaseCavity': {'names': attrs}})
                for name, attr in attrs.items():
                    data.update({name: np.zeros(nevents)})
            else:
                try:
                    #detector = self._detectors.get(det)
                    while det not in evt._attrs and ievt < min([1000,nevents]):
                        evt.next(publish=False, init=False)
                        ievt += 1
                    if auto_update and hasattr(detector, '_update_xarray_info'):
                        try:
                            detector._update_xarray_info()
                        except:
                            pass

                    info = detector._xarray_info
                    attr_info = info.get('dims',{})
                    attrs = {}
                    attrs1d = {}
                    for attr, vals in attr_info.items():
                        if attr in ignore_attrs:
                            #ignore some attrs, e.g., timestamp
                            continue
                        try:
                            b = list(vals[1])
                        except:
                            b = [vals[1]]
                        try:
                            if len(vals[0]) == 0:
                                name = det+'_'+attr
                                attrs[name] = attr
                                data[name]  = np.zeros(nevents)
                            elif add_1d and len(vals[0]) == 1 and b[0] <= 16:
                                nch = b[0]
                                name = det+'_'+attr
                                attrs1d[name] = attr
                                data1d[name]  = np.zeros((nevents, nch))
                        except:
                            logger.info('Cannot add {:}, {:} -- {:}'.format(det, attr, vals))
                    if add_counts and not attrs and 'corr' in detector._attrs:
                        detector.add.count('corr')
                        for attr, vals in attr_info.items():
                            if len(vals[0]) == 0:
                                name = det+'_'+attr
                                attrs[name] = attr
                                data[name]  = np.zeros(nevents)
                    if attrs:
                        dets[det] = {'names': attrs}
                    if attrs1d:
                        dets1d[det] = {'names': attrs1d}

                except:
                    logger.info('Cannot add {:}'.format(det))

        self.reload()
    
    nupdate = 100
    time_last = time0
    logger.info('Loading Scalar: {:}'.format(list(dets.keys())))
    if add_1d:
        logger.info('Loading 1D: {:}'.format(list(dets1d.keys())))
    i = 0
    for evt in self.events:       
        # Skip events without event code with are only controls cameras 
        # nevents in ds does not include controls cameras
        if not evt.Evr.eventCodes:
            continue

        # break at nevents 
        if i >= nevents:
            logger.info('Breaking event loop at event {:} with more events available in DataSource'.format(i))
            break
        if i % nupdate == nupdate-1:
            time_next = time.time()
            dtime = time_next-time_last
            logger.info('{:8} of {:8} -- {:8.3f} sec, {:8.3f} events/sec'.format(i+1, 
                    nevents, time_next-time0, nupdate/dtime))
            time_last = time_next 

        # add step
        istep = self._istep
        data['step'][i] = istep
        if i % nupdate == nupdate-1:
            print((istep, evt))

        # add eventCodes
        for code in evt.Evr.eventCodes_strict:
            if code in cnames:
                data[cnames[code]][i] = code

        # add detectors present
        for det in self._detectors: 
            data[det+'_present'][i] = True

        if add_dets:
            # Add detector attribute scalar data
            for det, attrs in dets.items():
                if det == 'EventId' or det in evt._attrs:
                    detector = getattr(evt, det)
                    for name, attr in attrs['names'].items():
                        try:
                            data[name][i] = getattr(detector, attr)
                        except:
                            data[name][i] = np.nan
            
            if add_1d:
                # optionally add detector scalar data
                for det, attrs in dets1d.items():
                    if det in evt._attrs:
                        detector = getattr(evt, det)
                        for name, attr in attrs['names'].items():
                            try:
                                data1d[name][i] = getattr(detector, attr)
                            except:
                                data[name][i] = np.nan
                                
        i += 1


    df = pd.DataFrame(data)
    x = df.to_xarray()
    x = x.rename({'index':'time'})
    x['time'] = [np.datetime64(int(sec*1e9+nsec), 'ns') for sec,nsec in zip(x.sec,x.nsec)]
    x.attrs['data_source'] = self.data_source.__str__()
    x.attrs['instrument'] = self.data_source.instrument.upper()
    x.attrs['run'] = self.data_source.run
    x.attrs['experiment'] = self.data_source.exp
    x.attrs['expNum'] = self.expNum
    # add attributes
    self.reload()
    evt = next(self.events)
    for det, item in dets.items():
        detector = self._detectors.get(det)
        # get next event with detector in event
        # needed when controls cameras and/or groups present
        if det not in ['EventId', 'Evr']:
            try:
                next(detector)
            except:
                self.reload()
                try:
                    next(detector)
                except:
                    print('Cannot load attributes for {:}'.format(det))
            
        try:
            source_info = detector._source_info
        except:
            source_info = {}
        
        try:
            x[det+'_present'].attrs.update(source_info)
        except:
            pass

        if det == 'EventId':
            continue
        try:
            detector._update_xarray_info()
            det_info = detector._xarray_info['dims']
            for name, attr in item['names'].items():
                x[name].attrs.update(source_info)
                info = det_info.get(attr, ([],(),{},))
                x[name].attrs['attr'] = attr
                x[name].attrs['alias'] = det
                if len(info) == 3:
                    attrs = info[2]
                    x[name].attrs.update(attrs)
        except:
            logger.info('Error updating scalar data for {:}: {:}'.format(det, item))
    # add 1d
    if add_1d:
        for det, item in dets1d.items():
            detector = self._detectors.get(det)
            try:
                detector._update_xarray_info()
                det_info = detector._xarray_info['dims']
                for name, attr in item['names'].items():
                    info = det_info.get(attr, ([],(),{},))
                    if len(info[0]) == 1:
                        x[name] = (('time', det+'_'+info[0][0]), data1d.get(name))
                    if len(info) == 3:
                        attrs = info[2]
                        x[name].attrs['attr'] = attr
                        x[name].attrs['alias'] = det
                        x[name].attrs.update(attrs)
            except:
                logger.info('Error updating 1D data for {:}: {:}'.format(det, item))

    x = x.set_coords(coords)
    for code, ec in cnames.items():
        try:
            ec_doc = self.configData._eventcodes.get(code,{}).get('description','')
            if ec_doc:
                ec_doc = 'event code for '+ec_doc
            x[ec].attrs['doc'] = ec_doc 
        except:
            print('Cannot add doc for {:}'.format(ec))

    if ignore_unused_codes:
        drop_codes = []
        for code, ec in cnames.items():
            if ec != drop_code and not x[ec].any():
                drop_codes.append(code)
                x = x.drop(ec)
        if drop_codes:
            logger.info('Dropping unused eventCodes: {:}'.format(drop_codes))

    if 'ec162' in x:
        set_delta_beam(x, code=drop_code, attr=drop_attr)
        x.coords['XrayOff'] = (x.ec162 == True)
        x.coords['XrayOff'].attrs['doc'] = 'Xray Off for events with ec162'
        x.coords['XrayOn'] = (x.ec162 == False)
        x.coords['XrayOn'].attrs['doc'] = 'Xray On for events without ec162'

    self.x = x

    if save:
        if not path:
            path = os.path.join(self.data_source.res_dir,'nc')
        elif path == 'home':
            path = os.path.join(os.path.expanduser('~'), 'RunSummary', 'nc')

        if not os.path.isdir(path):
            os.mkdir(path)
        
        if not filename:
            filename = 'run{:04}_smd.nc'.format(self.data_source.run)

        try:
            self.x = clean_dataset(self.x)
            #self.x.to_netcdf(os.path.join(path,filename), engine='h5netcdf', invalid_netcdf=True)
            self.x.to_netcdf(os.path.join(path,filename), engine='h5netcdf')
        except:
            traceback.print_exc('Cannot save to {:}/{:}'.format(path,filename))

        try:
            xsum = to_summary(x)
            fsumname = os.path.join(path, '_sum.'.join(filename.split('.', 1)))
            xsum.to_netcdf(fsumname, engine='h5netcdf')
        except:
            traceback.print_exc('Cannot make smd sum file')

    return self.x

def rawsum(self, attr='raw'):
    """Return raw sum of AreaDetector
    """
    import numpy as np
    return np.sum(getattr(self, attr))

def count(self, attr='corr'):
    """Return pedestal corrected count of AreaDetector
    """
    import numpy as np
    return np.sum(getattr(self, attr))

def sector(self, attr='raw'):
    """Return calib mean of each AreaDetector sector
    """
    import numpy as np
    return getattr(self, attr).sum()

def sectors(self, attr='raw'):
    """Return calib mean of each AreaDetector sector
    """
    import numpy as np
    data = getattr(self, attr)
    return np.array([data[isector].sum() for isector in range(data.shape[0])])

def wave8_height(self, bkrange=[500,600]):
    """Peaks of 8 waveforms with background subtration from mean within bkrange.
    """
    import numpy as np
    wfs = []
    peaks = []
    for ch in range(8):
        wf = self.evtData.data_u32[ch]
        if len(wf) == 0:
            wf = self.evtData.data_u16[ch]
        back = wf[bkrange[0]:bkrange[1]].mean()
        wf = -1.*wf+back
        wfs.append(wf)
        peaks.append(max(wf))

    return np.array(peaks)

def peak_height(self):
    """Max value of each waveform for WFDetector.
    """
    if hasattr(self, 'waveform') and hasattr(self.waveform, 'max'):
        return self.waveform.max(axis=1)
    else:
        return None

def peak_time(self):
    """Time of max of each waveform for WFDetector.
    """
    import numpy as np
    return np.array([self.wftime[ch][index] for ch,index in enumerate(self.waveform.argmax(axis=1))])

def amplitudes(self):
    """
    Returns an array of max values of a set of waveforms.
    """
    return filtered(self).max(axis=1)

def filtered(self, signal_width=10):
    """
    Returns filtered waveform for WFDetector.
    """
    from scipy import signal
    import numpy as np
    waveform = self.waveform
    af = []
    hw = signal_width/2
    afilter = np.array([-np.ones(hw),np.ones(hw)]).flatten()/(hw*2)
    nch = waveform.shape[0]
    for ich in range(nch):
        wf = waveform[ich]
        f = -signal.convolve(wf, afilter)
        f[0:len(afilter)+1] = 0
        f[-len(afilter)-1:] = 0
        af.append(f[hw:wf.size+hw])
    
    return  np.array(af)


def save_exp_stats(exp, instrument=None, path=None, find_corr=True):
    """
    Save drop stats summaries for all runs 
    """
    import os
    import glob
    import xarray as xr
    from .xarray_utils import find_beam_correlations
    from .xarray_utils import clean_dataset
    if not instrument:
        instrument = exp[0:3]
    if not path:
        path = os.path.join('/reg/d/psdm/',instrument,exp,'results','nc')

    files = sorted(glob.glob('{:}/run*_{:}.nc'.format(path, 'drop_stats')))

    for f in files:                                      
        x = xr.open_dataset(f, engine='h5netcdf')
        try:
            run = x.run
            save_file='{:}/run{:04}_{:}.nc'.format(path, run, 'drop_sum')
            print(save_file)
            if not find_corr and glob.glob(save_file):
                xstats = xr.open_dataset(save_file, engine='h5netcdf')
                xout = xstats.copy(deep=True) 
                xstats.close()
                del xstats
                for avar, da in x.data_vars.items():
                    try:
                        if avar in xout:
                            for a in ['doc', 'unit', 'alias']:
                                val = x[avar].attrs.get(a, '') 
                                if isinstance(val, (six.binary_type, six.text_type)):
                                    val = str(val)
                                xout[avar].attrs[a] = val 
                    except:
                        print('Cannot add attrs for {:}'.format(avar))
                
                for attr, val in x.attrs.items():
                    try:
                        if isinstance(val, list) and len(val) > 0 and isinstance(val[0], (six.binary_type, six.text_type)):
                            val = [str(v) for v in val]
                        elif isinstance(val, (six.binary_type, six.text_type)):
                            val = str(val)
                        xout.attrs[attr] = val
                    except:
                        print('Cannot add attrs for {:}'.format(attr))

                
                xout = clean_dataset(xout)
                #xout.to_netcdf(save_file, engine='h5netcdf', invalid_netcdf=True)
                xout.to_netcdf(save_file, engine='h5netcdf')

            else:
                xout = x.copy(deep=True) 
                x.close()
                del x
                xstats = find_beam_correlations(xout, groupby='ec162', cut='XrayOn', save_file=save_file)
                #xout.to_netcdf(f, engine='h5netcdf', invalid_netcdf=True)
                xout.to_netcdf(f, engine='h5netcdf')
           
        except:
            traceback.print_exc('Cannot do {:}'.format(f))

    return load_exp_sum(exp)


