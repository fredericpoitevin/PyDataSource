from __future__ import print_function

import re
import operator
import sys
import os
import time
import traceback
from pylab import *

instrument_transmission_pvs = {
        'FEE': ['GATT:FEE1:310:R_ACT','SATT:FEE1:320:RACT'],
        'XPP': ['XPP:ATT:COM:R_CUR'],
        'XCS': ['XCS:ATT:COM:R_CUR'],
        'MFX': ['MFX:ATT:COM:R_CUR'],
        'CXI': ['CXI:DSB:ATT:COM:R_CUR','XRT:DIA:ATT:COM:R_CUR'],
        }

def write_exp_summary(self, file_name=None, path=None, **kwargs):
    """Write ExperimentSummary as pickle file.
    """
    import cPickle as pickle
    if not file_name:
        file_name = 'experiment_summary.pkl'
    
    if path == 'scratch':
        path = os.path.join(self.scratch_dir,'RunSummary')
    elif not path:
        path = os.path.join(self.res_dir,'RunSummary') 

    if not os.path.isdir(path):
        os.mkdir(path)
    
    try:
        self.xpvs.to_netcdf(path=os.path.join(path,'xpvs.nc'), engine='h5netcdf')
        self.xscan.to_netcdf(path=os.path.join(path,'xscan.nc'), engine='h5netcdf')
        self.xset.to_netcdf(path=os.path.join(path,'xset.nc'), engine='h5netcdf')
    except:
        traceback.print_exc()

    try:
        with open(path+'/'+file_name, 'wb') as pickle_file:
            pickle.dump(pickle.dumps(self, protocol=-1), pickle_file, protocol=-1)
    except:
        traceback.print_exc()
        print(('Failes writing pickle file', path+'/'+file_name))

def open_epics_data(exp=None, file_name=None, path=None, run=None, **kwargs):
    import glob
    import xarray as xr
    if not file_name:
        file_name = 'xpvs.nc'
    
    if exp is None:
        print('Must supply exp Name')
        return 

    instrument = exp[0:3]
    exp_dir = "/reg/d/psdm/{:}/{:}".format(instrument, exp)
    if path == 'scratch':
        path = os.path.join(exp_dir,'scratch', 'RunSummary')
    elif os.path.isdir(os.path.join(exp_dir,'res')): 
        path = os.path.join(exp_dir,'res','RunSummary')
    else:
        path = os.path.join(exp_dir,'results','RunSummary')

    full_file = path+'/'+file_name
    if glob.glob(full_file):
        try:
            xdata = xr.open_dataset(full_file, engine='h5netcdf')
        except:
            traceback.print_exc()
            print(('Failes reading file', full_file))
            return None

    else:
        print(('Failes finding file', full_file))
        return None
    
    return xdata

def get_pv_attrs(exp, auto_load=False):
    """
    Get dict of attributes for epics pvs

    Parameters
    ----------
    exp : str
        Experiment Name
    
    auto_load : bool
        Automatically reload exp_summary if not existing

    """
    try:
        xpvs = open_epics_data(exp, 'xpvs.nc')
    except:
        print(('Failed open_epics_data for', exp))
        xpvs = None

    if xpvs is not None:
        if auto_load:
            es = get_exp_summary(exp, reload=True)
            xpvs = es.xpvs
        else:
            return {} 

    return {a: dict(item.attrs) for a, item in xpvs.data_vars.items()}

def get_scan_pvs(exp, run=None, auto_load=False):
    """
    Get scan pvs for exp.  
    If run is provided return dict of alias, pv
    Otherwise return pandas data frame for all runs
    
    Parameters
    ----------
    exp : str
        Experiment Name
    run : int
        Run number (optional)
    auto_load : bool
        Automatically reload exp_summary if not existing
    
    """
    try:
        xscan = open_epics_data(exp, 'xscan.nc')
    except:
        xscan = None

    if xscan is not None:
        if auto_load:
            es = get_exp_summary(exp, reload=True)
            xscan = es.xscan
        else:
            return None

    attrs = [a for a in xscan.data_vars.keys()]
    dfscan = xscan.sel(stat='count').reset_coords()[attrs]
    dfscan = dfscan.dropna(dim='run', how='all').to_dataframe().astype(int)
   
    if run:
        df = dfscan.T[run]
        df = df[df>0].to_dict()
        if df is None:
            df = []
        return df
    else:
        return dfscan.T

def read_exp_summary(exp=None, file_name=None, path=None, **kwargs):
    """
    Read exp_summary pickle file.  Return none if does not exist.
    """
    import glob
    import cPickle as pickle
    if not file_name:
        file_name = 'experiment_summary.pkl'
    
    if exp is None:
        print('Must supply exp Name')
        return 

    instrument = exp[0:3]
    exp_dir = "/reg/d/psdm/{:}/{:}".format(instrument, exp)
    if path == 'scratch':
        path = os.path.join(exp_dir,'scratch', 'RunSummary')
    elif os.path.isdir(os.path.join(exp_dir,'res')): 
        path = os.path.join(exp_dir,'res','RunSummary')
    else:
        path = os.path.join(exp_dir,'results','RunSummary')

    full_file = path+'/'+file_name
    if glob.glob(full_file):
        try:
            with open(full_file, 'rb') as pickle_file:
                data = pickle.load(pickle_file)
            
            return pickle.loads(data)
        
        except:
            traceback.print_exc()
            print(('Failes reading pickle file', full_file))
            return None

    else:
        return None

def get_exp_summary(exp, reload=False, path=None, save=None, **kwargs):
    """
    Use scratch for now since results unreliable
    """
    es = None
    if not reload:
        es = read_exp_summary(exp, path=path, **kwargs)
    
    if es is None:
        if save is None:
            save = True
        es = ExperimentSummary(exp, path=path, save=save, **kwargs)

    return es

def epicsArch_dict(archfile_name,file_dir):
    """
    Load text file containing aliases and pvs with daq epicsArch convention. 
    """
    import re
    import traceback
    arch_dict = {}
    file_list = [archfile_name]
    try:
        with open(file_dir+'/'+archfile_name,'r') as f:
            for line in f:
                if line.startswith('<'):
                    file_list.append(line.strip('<').strip().strip('\n').strip())

        arch_list = []
        for file in file_list:
            with open(file_dir+'/'+file,'r') as f:
                arch_list.append(f.read().split('\n'))
                    
        arch_list = [item.strip(' ') for sublist in arch_list for item in sublist]
        pvalias = None
        for item in arch_list:
            if item.startswith(('#','<')):
                pvalias = None
            elif item.startswith('*'):
                pvalias = item.strip('*').strip(' ')
            elif len(item) > 1:
                pvname = item.strip(' ')
                pvbase = pvname.split('.')[0]
                if pvalias:
                    if pvalias in arch_dict:
                        print('Warning: duplicate alias {:}'.format(pvalias))
                else:
                    pvalias = re.sub(':|\.','_',pvname) 

                components = re.split(':|\.',pvname)
                for i,item in enumerate(components):
                    if item[0].isdigit():
                         components[i] = 'n'+components[i]
     
                arch_dict[pvname] = {'name': pvname,
                                     'alias': pvalias,
                                     'base':  pvbase,
                                     'components': components} 
    
    except:
        traceback.print_exc()
        print('Error loading {:} from {:}'.format(archfile_name,file_dir))

    return arch_dict



class ExperimentSummary(object):
    """
    Experiment information from elog and epics archive data. 

    Attributes
    ----------
    xscan : xarray.Dataset
        

    """

    _fields={
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

    def __init__(self, exp=None, instrument=None, station=0, exper_id=None,
                exp_dir=None, xtc_dir=None, h5_dir=None, scratch_dir=None, 
                save=True, init=True, **kwargs):
        """
        Experiment Summary
        """
        from RegDB import experiment_info
        if exp:
            self.exp = exp
            self.exper_id = experiment_info.name2id(exp)
        
        elif exper_id:
            self.exper_id = exper_id
            self.exp = experiment_info.getexp(exper_id)
        
        elif instrument:
            self.exp = experiment_info.active_experiment(instrument, station=station)
            self.exper_id = experiment_info.name2id(exp)

        if not instrument:
            instrument = self.exp[0:3]
            
        self.station = station
        self.instrument = instrument

        if not exp_dir:
            exp_dir = "/reg/d/psdm/{:}/{:}".format(self.instrument, self.exp)

        if not xtc_dir:
            xtc_dir =  "{:}/xtc".format(exp_dir, self.exp)

        if not h5_dir:
            h5_dir =  "{:}/hdf5".format(exp_dir, self.exp)
       
        if not scratch_dir:
            scratch_dir =  "{:}/scratch".format(exp_dir)

        self.exp_dir = exp_dir
        self.h5_dir = h5_dir
        self.xtc_dir = xtc_dir
        if os.path.isdir(os.path.join(self.exp_dir,'res')): 
            self.res_dir = os.path.join(self.exp_dir,'res')
        else:
            self.res_dir = os.path.join(self.exp_dir,'results')
        
        self.scratch_dir = os.path.join(self.exp_dir,'scratch')

        if init:
            self._init(save=save, **kwargs)

    def _init(self, save=False, **kwargs):
        """
        Initialize
        """
        self._load_run_info()
        try:
            self._init_arch(**kwargs)
        except:
            traceback.print_exc('Could not load epics Archive')
        
        try:
            self._load_exp_runs(**kwargs)
        except:
            traceback.print_exc('Could not load exp_runs')
            
        if save:
            try:
                self.save(**kwargs)
            except:
                traceback.print_exc('Could not write summary')

    def save(self, file_name=None, path=None, **kwargs):
        """Save to pickle file.
        """
        write_exp_summary(self, file_name=file_name, path=path)

    def _add_runtables(self):
        """
        Add RunTables... currently does not work in conda env.
        """
        try:
            from LogBook.runtables import RunTables

            # Setup elog pswwww RunTables
            self._RunTables = RunTables(**{'web-service-url': 'https://pswww.slac.stanford.edu/ws-kerb'})

            self._user_tables = self._RunTables.usertables(exper_name=self.exp)
            #self.add_user_run_table(RunSummary='Run Summary')
        except:
            print('currently does not work in conda env.')

    def add_user_run_table(self, **kwargs):
        """Add a user User RunTable from pswww elog server.
        """
        for alias, name in kwargs.items():
            tbl = self._RunTables.findUserTable(exper_name=self.exp, table_name=alias)
            self._user_tables.update({alias: tbl}) 
            setattr(self, alias, tbl)

    def to_html(self, **kwargs):
        """
        Build experiment html report.  
        see build_html.Build_experiment
        """
        import build_html
        return build_html.Build_experiment(self, **kwargs)

    def detectors(self, run):
        """
        Return a list of detector names configured in the DAQ system for the input run number.
        """
        from RegDB import experiment_info
        return experiment_info.detectors(self.instrument, self.exp, run)

    @property
    def calibration_runs(self):
        """
        List of calibration runs for the experiment.
        """
        from RegDB import experiment_info
        return experiment_info.calibration_runs(self.instrument, self.exp)

    def show_moved(self, attrs=None):
        """
        Show motors that have moved (or devices that have changed state).

        Parameters
        ----------
        attrs : list
            List of motor aliases.
        """
        if attrs:
            if not isinstance(attrs, list):
                attrs = [attrs]
        else:
            attrs = self.xset.data_vars.keys()

        for attr in attrs:
            xvar = self.xset.reset_coords().set_coords(['begin_time','prep_time'])[attr].dropna(dim='run')
            print('')
            print('** {:} **'.format(attr))
            for a,val in xvar.attrs.items():
                if a not in ['PREC']:
                    print('{:10} {:10}'.format(a,val))
            print(xvar.to_dataframe())

    def get_run_sets(self, delta_time=60*60*24, run_min=None, run_max=None, **kwargs):
        xscan = self.xscan
        run_last = xscan.run.values[0]
        if run_min:
            xscan = xscan.where(xscan.run > run_min, drop=True)
            run_last = run_min
        if run_max:
            xscan = xscan.where(xscan.run <= run_max, drop=True)

        runs = list(xscan.where(xscan.prep_time > delta_time, drop=True).run.values)
        runs.append(xscan.run.values[-1])
        aruns = []
        for run in runs:
            run_next = int(run)
            aruns.append((run_last, run_next,))
            run_last = run_next+1
        
        return aruns

    def get_transmission_data(self, run, time=None, **kwargs):
        """Get transmission dataframe.

        Parameters
        ----------
        run : int
            Run number
        time : array
            Array of times
        """
        df = self.get_scans(**kwargs).T
        if run in df:
            df = df[run]
            attrs = df.where(df > 0).dropna().keys() 
            a = self.xruns.sel(run=run)
            t = self.xepics.time
            xepics = self.xepics[attrs].where(t>=a.begin_time).where(t<=a.end_time).dropna(dim='time',how='all')
            
        else:
            print('No scan data for run {:}'.format(run))
            return None

        if time:    
            xepics = fill_times(xepics, time)
        
        return xepics

    def get_scan_data(self, run, time=None, **kwargs):
        """Get scan dataframe.

        Parameters
        ----------
        run : int
            Run number
        time : array
            Array of times
        """
        df = self.get_scans(**kwargs).T
        if run in df:
            df = df[run]
            attrs = df.where(df > 0).dropna().keys() 
            a = self.xruns.sel(run=run)
            t = self.xepics.time
            xepics = self.xepics[attrs].where(t>=a.begin_time).where(t<=a.end_time).dropna(dim='time',how='all')
            
        else:
            print('No scan data for run {:}'.format(run))
            return None

        if time:    
            xepics = fill_times(xepics, time)
        
        return xepics

    def plot_scan(self, run, style=None, linewidth=2,
            min_steps=4, attrs=None, device=None, min_motors=1,
            figsize=(12,8), loc='best', **kwargs):
        """
        Plot motor moves vs time and run.
        """
        fig = plt.figure(figsize=figsize, **kwargs)
        ax = fig.gca()
        x = self.get_scan_data(run, attrs=attrs, min_steps=min_steps, device=device, min_motors=min_motors)
        if x is None:
            print('Run {:} scan for attrs={:} device={:} not valid'.format(run, attrs, device))
        else:
            attrs = x.data_vars.keys()
            if len(attrs) == 1:
                ylabel = attrs[0]
            else:
                ylabel = 'Motors'
     
            aunits = {a: self.xepics[a].attrs.get('units', '') for a in attrs}
            units = list(set(aunits.values()))
            if len(set(aunits.values())) == 1 and units[0] is not '':
                ylabel += ' [{:}]'.format(units[0])

            ax.set_ylabel(ylabel)
            for attr in attrs:
                lab = attr
                if aunits[attr]:
                    lab += ' [{:}]'.format(aunits[attr])
                
                df = x[attr].dropna(dim='time').to_pandas().tz_localize('UTC').tz_convert('US/Pacific')
                if style:
                    df.plot(style=style, linewidth=linewidth, label=lab, ax=ax)
                else:
                    df.plot(drawstyle='steps', linewidth=linewidth, label=lab, ax=ax)
                
            ax.set_xlabel('time (US Pacific)')
            legend = ax.legend(loc=loc)

    def get_epics(self, attrs=None, run_min=None, run_max=None):
        """Get epics motion
        """
        xepics = self.xepics
        if attrs:
            if not isinstance(attrs, list):
                attrs = [a for a in self.xepics.data_vars.keys() if a.startswith(attrs)]
            else:
                attrs = [a for a in attrs if a in self.xepics.data_vars.keys()]
            xepics = xepics.get(attrs)

        if run_min:
            xepics = xepics.where(xepics.run >= run_min, drop=True)

        if run_max:
            xepics = xepics.where(xepics.run <= run_max, drop=True)

        return xepics

    def _pvs_with_set(self):
        aout = []
        attrs = self.xset.data_vars.keys()
        for a in attrs:
            if a.endswith('_set'):
                aout.append(a.split('_set')[0])
            else:
                aout.append(a)
        return aout

    def get_moved(self, attrs=None, run_min=None, run_max=None, group=None, min_count=2):
        """Get devices that moved (between run_min and run_max if specified)
        """
        import numpy as np
        if attrs:
            group = False

        if isinstance(attrs, list):
            for attr in set(attrs):
                if not attr.endswith('_set'):
                    attrs.append(attr+'_set')
        
        xepics = self.get_epics(attrs=attrs, run_min=run_min, run_max=run_max)
        attrs = []
        df = xepics.to_dataframe()[xepics.data_vars.keys()]
        for attr, item in df.iteritems():
            if len(set(item.astype(np.float32).dropna().drop_duplicates().values)) >= min_count:
                attrs.append(attr)
        
        if attrs:
            df = xepics[attrs].to_dataframe()
        else:
            return attrs

        attrs = list(set(attrs))

        if group is not False:
            if group == 'units':
                adets = {}
                for attr in attrs:
                    det = attr.split('_')[0] 
                    try:
                        unit = xepics[attr].attrs.get('units','')
                        if unit == '%':
                            unit='percent'
                        unit = re.sub('\W+','', unit)
                    except:
                        unit = ''
                    grp = det+'__'+unit
                    if grp not in adets:
                        adets[grp] = []
                    adets[grp].append(attr)
                return adets
            else:
                dets = list(set([a.split('_')[0] for a in attrs]))
                return {det: [a for a in attrs if a.startswith(det)] for det in dets}
        else: 
            return attrs

    def plot_move(self, attrs, run_min=None, run_max=None, style=None, linewidth=2, 
            figsize=(12,8), 
            ax_pos = (.1,.2,.8,.7), box_to_anchor=(0., 1.), 
            #ax_pos = (.1,.2,.55,.7), box_to_anchor=(1.1, 1.), 
            max_points=5000, **kwargs):
        """
        Plot motor moves vs time and run.
        """
        if not isinstance(attrs, list):
            xepics = self.get_epics(run_min=run_min, run_max=run_max)
            if attrs not in xepics.data_vars:
                ylabel = attrs
                print(attrs)
                attrs = self.get_moved(attrs, run_min=run_min, run_max=run_max, group=False)
            else:
                ylabel = attrs
                attrs = [attrs]
        else:
            ylabel = 'Motors'
    
        moved_attrs = self.get_moved(attrs, run_min=run_min, run_max=run_max, group=False)
        attrs = [a for a in attrs if a in moved_attrs]
        if not attrs:
            return None

        xepics = self.get_epics(attrs=attrs, run_min=run_min, run_max=run_max)
        
        fig = plt.figure(figsize=figsize, **kwargs)
        ax = fig.gca()
        b = xepics['begin_run'].dropna(dim='time').to_pandas()
        axr = b.plot(secondary_y=True, drawstyle='steps', label='run', ax=ax)
        axr.set_ylabel('Run')
        
        aunits = {}
        for a in attrs:
            if a in xepics:
                aunits[a] =  xepics[a].attrs.get('units', '')
            elif a+'_set' in xepics:
                aunits[a] =  xepics[a+'_set'].attrs.get('units', '')

        units = list(set(aunits.values()))
        if len(set(aunits.values())) == 1 and units[0] is not '':
            ylabel += ' [{:}]'.format(units[0])
    
        ax.set_ylabel(ylabel)
        for attr in attrs:
            lab = attr
            if aunits.get(attr):
                lab += ' [{:}]'.format(aunits[attr])
            print(attr)
            try:
                df = xepics[attr].dropna(dim='time').to_pandas()
                if df.size > max_points:
                    # need to resample
                    sample_time = (df.index.max()-df.index.min()).seconds/(max_points/50)
                    df_resampled = df.resample('{:}S'.format(sample_time))
                    ax.errorbar(df_resampled.mean().index, df_resampled.mean(), yerr=df_resampled.var(), label=lab)
                    ax.set_ylabel(ylabel)
                elif style:
                    df.plot(style=style, linewidth=linewidth, label=lab, ax=ax)
                else:
                    df.plot(drawstyle='steps', linewidth=linewidth, label=lab, ax=ax)
            except:
                traceback.print_exc()
                print(('cannot plot', attr))

            try:
                if isinstance(attr, str) and attr in xepics and xepics.get(attr+'_set'):
                    xepics[attr+'_set'].dropna(dim='time')[1:].to_pandas().plot(style='+',label=attr+'_set')
            except:
                print(('cannot plot', attr+'_set'))

        ax2 = ax.twiny()
        ax2.set_xlabel('Run')
        xticks_minor = (b.index-b.index.min())/(b.index.max()-b.index.min())*(ax.get_xlim()[1]-ax.get_xlim()[0])+ax.get_xlim()[0]
        run_step = max([1,len(xticks_minor)/100*10])
        xticks = xticks_minor.where((b % run_step) == 0).dropna()
        xticklabels = np.array(b.where((b % run_step) == 0).dropna(), dtype=int)
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(xticks)
        ax2.set_xticklabels(xticklabels, rotation='vertical')
        ax2.set_xticks(xticks_minor, minor=True)
        legendr = axr.legend(loc='upper right')
        axr.set_position(ax_pos)
        ax2.set_position(ax_pos)
        ax.set_position(ax_pos)
        legend = ax.legend(loc='upper left',  bbox_to_anchor=box_to_anchor)
        return ax

    def _add_run_details(self):
        """
        Add experiment run details to xrun xarray Dataset
        """
        xruns = self.xruns
        flist = ['drop_stats_files', 'xtc_files', 'stats_files']
        for ftyp in flist:
            try:
                df = self.dfruns.T[ftyp]
                xruns[ftyp] = (('run'), np.array([len(df[run]) for run in df.index], dtype=bool))
            except:
                pass

        dlist = []
        ddict = {}
        for run in xruns.run.values:
            dets = self.detectors(run)
            ddict[run] = dets
            dlist += dets

        aliases = {re.sub('-|:|\.| |\|','_', a): a for a in set(dlist)}
        for alias, det in aliases.items():
            try:
                xruns[alias] = (('run'), np.array([det in a for a in ddict.values()], dtype=bool))
            except:
                pass

        xruns.attrs['dets'] = aliases.keys()
        xruns.attrs['detectors'] = aliases.keys()

        return xruns

    @property
    def runs_with_xtc(self):
        """
        Return list of runs with xtc files.
        """
        runs = []
        for run, file_names in self.dfruns.T['xtc_files'].iteritems():
            if file_names != []:
                runs.append(run)

        return runs

    def last_pv_fields_changed(self):
        """
        Last time fields were changed.
        """
        pass

    def last_set(self, attrs=['value', 'units', 'run', 'description'], **kwargs):
        """Show the last time devices were set.
        """
        import pandas as pd
        avals = {}
        for attr in self.xset.data_vars:
            xvar = self.xset[attr].dropna(dim='run')
            avals[attr] = xvar.attrs
            avals[attr].update({'value': float(xvar[-1].values),
                                'run': xvar[-1].run.values, 
                                'time': xvar[-1].begin_time.values})
                
        return pd.DataFrame(avals).T[attrs]

    def _submit_summary(self, run, batchqueue='psanaq', option='beam_stats'):
        """
        Make beam_stats to check for timing error and drop shots
        """
        import subprocess
        bsubproc = 'submit_summary {:} {:} {:} {:}'.format(self.exp, run, batchqueue, option)
        print(('-> ',bsubproc))

        subproc = subprocess.Popen(bsubproc, stdout=subprocess.PIPE, shell=True)

    def save_configData(self, runs='all'):
        """
        Save configData .
        
        Parameters
        ----------
        runs : list of ints or int or 'all'
            Runs to save configData - if 'all' then submit all that have not already been created
            default = 'all'
        """
        import PyDataSource
        if runs == 'all':
            drop_stats_files = self.dfruns.T['drop_stats_files']
            runs = [run for run in self.runs_with_xtc if not drop_stats_files.get(run, [])]
        for run in runs:
            ds = PyDataSource.DataSource(exp=self.exp, run=run)
            try:
                ds.configData.save_configData()
                print('Saved {:}'.format(ds))
            except:
                traceback.print_exc()
                print(('Cannot save configData: ', ds))


    def submit_beam_stats(self, runs, batchqueue='psanaq', alert='None'):
        """
        Make beam_stats to check for timing error and drop shots

        Parameters
        ----------
        runs : list of ints or int or 'all'
            Runs to submit -- if 'all' then submit all that have not already been created
        """
        import subprocess
        if runs == 'all':
            drop_stats_files = self.dfruns.T['drop_stats_files']
            runs = [run for run in self.runs_with_xtc if not drop_stats_files.get(run, [])]

        if not isinstance(runs, list):
            runs = [runs]
        for run in runs:
            bsubproc = 'submit_beam_stats {:} {:} {:} {:}'.format(self.exp, run, batchqueue, alert)
            print(('-> ',bsubproc))
            subproc = subprocess.Popen(bsubproc, stdout=subprocess.PIPE, shell=True)

    def get_beam_summary(self):
        """
        Load drop stats summary for all runs
        """
        import beam_stats
        return beam_stats.load_exp_sum(self.exp) 

    def build_beam_stats(self, alert=True, **kwargs):
        """
        Build beam statistics summary including off-by-one.

        First run submit_beam_stats('all')
        """
        import beam_stats
        b = beam_stats.build_drop_stats(self.exp, alert=alert, **kwargs)
        return b

    def _load_run_data(self, summary=False, no_create=True):
        ax = {}
        for i in self.xruns.run.values:                         
            try:                                         
                x = PyDataSource.h5write.get_config_xarray(exp='mfx11116',run=i, summary=summary, no_create=no_create)
                if x is not None:
                    ax[i] = x
            except:     
                print(('cannot load run', i))

        self.xdata = ax
        return self.xdata

    def _load_run_info(self):
        """
        Load run info from experiment_info database
        """
        from RegDB import experiment_info
        import pandas as pd
        import xarray as xr

        rns = pd.DataFrame(experiment_info.experiment_runs(self.instrument.upper(),self.exp))
        if not rns.empty:
            xruns = xr.Dataset()
            xruns.coords['run'] = (['run'], rns['num'])
            xruns.coords['run_id'] = (['run'], rns['id'])
            xruns['begin_time'] = (['run'], pd.to_datetime(rns['begin_time']))
            xruns['end_time'] = (['run'], pd.to_datetime(rns['end_time']))
            xruns['duration'] = (['run'], (rns['end_time']/1.e9 - rns['begin_time']/1.e9))
            xruns['duration'].attrs['units'] = 'sec'
            xruns['prep_time'] = (['run'], (rns['begin_time']/1.e9)-(rns['end_time'].shift(1)/1.e9)) 
            xruns['prep_time'].attrs['units'] = 'sec'
            self.xruns = xruns
        
        else:
            self.xruns = None

        return self.xruns

    def _init_arch(self, min_run=None, max_run=None, **kwargs):
        """
        Initialize archive
        """
        from epicsarchive import EpicsArchive
        if not hasattr(self, 'xruns'):
            self._load_run_info()
        self._arch = EpicsArchive()
        df = self.xruns.to_dataframe()
        if max_run:
            df = df[df.index <= max_run]
        if min_run:
            df = df[df.index >= min_run]
        dt = df[1:]['begin_time'].min()
        tstart = [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]
        self._tstart = tstart
        dt = df['end_time'].max()
        tend = [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]
        self._tend = tend

    def _get_pv_from_arch(self, pv, tstart=None, tend=None):
        if not tstart:
            tstart = self._tstart
        if not tend:
            tend = self._tend
        vals = self._arch._get_json(pv, tstart, tend, False)
        if vals:
            return vals[0]


    def _in_archive(self, pv):
        """
        Check if pv is in archive.
        """
        return self._arch.search_pvs(pv, do_print=False) != []

    def _load_eventCodes(self):
        """Currently not archived
        """
        data_codes = {}
        seq_evtCodes = range(67,99)+range(167,199)+range(201,217)
        for num in seq_evtCodes:
            pvs = {
                    'inst_num': 'ECS:SYS0:0:EC_{:}_OWNER_ID'.format(num),
                    'desc': 'EVNT:SYS0:1:NAME{:}'.format(num),
                    #'ticks': 'EVNT:SYS0:1:ECS_{:}DLY.A'.format(num),
                  }
            data_codes[num] = {name: self._get_pv_from_arch(pv) for name, pv in pvs.items()}

        self._data_codes = data_codes

    def _load_epics(self, pvs=None, quiet=False, 
                    omit=['ABSE', 'BEAM', 'SIOC', 'USEG', 'VGBA', 'GATT'],
                    #      'MIRR:EE1:M2H.RBV', 'MIRR:FEE1:M1H.RBV'],
                    omit_rbv=['XRT:M2H','MIRR:FEE1','MIRR:XRT','STEP:M1H','STEP:M2H','MOVR:DMP1'],
                    omit_endswith=['CNT'],
                    run=None,
                    set_only=None,
                    max_size=50000,
                    **kwargs):
        """
        Make xarray and pandas objects to describe which epics variables were set before
        and during runs.
        
        For now resample with max_size for monitoring pvs in order to keep file size small enough to pickle.
        In future refactor to save/load data objects separately and write with pandas/xarray to hdf5.
        """
        import operator
        import numpy as np
        import xarray as xr
        import pandas as pd

        if not hasattr(self, 'xruns'):
            self._load_run_info()

        if not pvs:
            if set_only is None:
                set_only = True
            import PyDataSource
            autorun=False
            if not run:
                autorun=True
                if self.runs_with_xtc:
                    run = max(self.runs_with_xtc)
                #run=rns['num'].max()
            
            #pvnames = {pv:ds.epicsData.alias(pv) for pv in ds.epicsData.pvNames()}
            pvnames = {}
            if run:
                try:
                    ds = PyDataSource.DataSource(exp=self.exp, run=run)
                except:
                    if autorun:
                        # go to next to last if last run fails
                        run = self.runs_with_xtc[-2]
                        ds = PyDataSource.DataSource(exp=self.exp, run=run)
                    else:
                        pass

                for pv in ds.epicsData.pvNames():
                    pvalias = ds.epicsData.alias(pv)
                    if not pvalias:
                        pvalias = re.sub(':|\.','_',pv)
                    else:
                        pvalias = re.sub(':|\.|-| ','_',pvalias)
                   
                    # remove instrument from auto alias
                    pvalias = pvalias.lstrip(self.instrument.upper()+'_')
                    pvnames[pv] = pvalias

            self._pvnames0 = pvnames.copy()

            for a in omit or a.split(':')[0] in omit: 
                for pv in pvnames.copy():
                    if pv.startswith(a):
                        pvnames.pop(pv)

            for pv in pvnames.copy():
                for omit_pv in omit_endswith:
                    if pv.endswith(omit_pv):
                        try:
                            pvnames.pop(pv)
                        except:
                            print(('Cannot pop ',pv,'ending with',omit_pv))
                for omit_pv in omit_rbv:
                    if  omit_pv in pv and pv.endswith('.RBV'):
                        try:
                            pvnames.pop(pv)
                        except:
                            print(('Cannot pop ',pv, omit_rbv))

        else:
            if set_only is None:
                set_only = False
            pvnames = pvs

        evr_pvs = {}
        for pv, alias in pvnames.copy().items():
            if 'TRIG' in pv:
                pvsplit = pv.split('.')[0].split(':')
                pvbase = ':'.join(pvsplit[0:-1])
                if pvbase not in evr_pvs:
                    evr_pvs[pvbase] = {}
                evr_pvs[pvbase].update({pv: alias})

            elif pv.endswith('.RBV'):
                pvnames[pv.rstrip('.RBV')] = alias+'_set'
                if set_only:
                    del pvnames[pv]
            elif set_only and not pv.endswith('STATE'):
                del pvnames[pv]

        machine_pvs = []
        epicsArch_dict = self._load_machine_epicsArch()
        for attr, item in epicsArch_dict.items():
            alias = item['alias']
            pv = attr
            if pv.endswith('.RBV'):
                alias += '_set'
                attr = pv.rstrip('.RBV')
            
            pvnames[pv] = alias
            machine_pvs.append(alias)

        self._pvnames = pvnames
        #pvmots = {pv.rstrip('.RBV'): alias+'_set' for pv, alias in pvnames.items() if pv.endswith('RBV')}
        #pvnames.update(**pvmots)
        #pvs = {pv: alias for pv, alias in pvnames.items() if arch.search_pvs(pv, do_print=False) != []} 
        pvs = {pv: alias for pv, alias in pvnames.items() if self._in_archive(pv)} 

        meta_attrs = {'units': 'EGU', 'PREC': 'PREC', 'pv': 'name'}

        data_arrays = {} 
        data_machine = {} 
        data_points = {}
        data_fields = {}
        import time
        time0 = time.time()
        time_last = time0
        for pv, alias in pvs.items():
            data_fields[alias] = {}
            dat = self._get_pv_from_arch(pv)
            if dat:
                try:
                    attrs = {a: dat['meta'].get(val) for a,val in meta_attrs.items() if val in dat['meta']}
                    for attr, item in self._fields.items():  
                        try:
                            field=item[0]
                            pv_desc = pv.split('.')[0]+'.'+field
                            if self._in_archive(pv_desc):
                                desc = self._get_pv_from_arch(pv_desc)
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
                                            vt = np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns')
                                            vals[vt] = val
                                   
                                    data_fields[alias][attr] = xr.DataArray(vals.values(), 
                                                                    coords=[vals.keys()], dims=['time'], 
                                                                    name=alias+'_'+attr, attrs=fattrs) 
                                    attrs[attr] = val
             
                        except:
                            traceback.print_exc()
                            print(('cannot get meta for', alias, attr))
                            
                    vals = [item['val'] for item in dat['data']]
                    doc = attrs.get('description','')
                    units = attrs.get('units', '')
                    time_next = time.time()
                  
                    try:
                        if vals:
                            if isinstance(vals[0],str):
                                if not quiet:
                                    print((alias, 'string'))
                                vals = np.array(vals, dtype=str)
                            else:
                                times = [np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns') for item in dat['data']]
                                # Simple way to remove duplicates
                                # Duplicates will break xarray merge
        #                        times_o = [a for a in times]
        #                        vals_o = [v for v in vals]
        #                        try:
        #                            aa = dict(zip(times, vals))
        #                            #aa = sorted(dict(zip(times, vals)).items(), key=operator.itemgetter(0))
        #                        except:
        #                            print 'cannot sort', alias
        #                            aa = dict(zip(times, vals))
        #                        times = np.array(aa.keys(), dtype=times[0].dtype)
        #                        vals = np.array(aa.values())
        #                        if not alias.endswith('set') and alias not in machine_pvs and len(aa) > max_size:
        #                            print '    ... resampling', alias
        #                            inds = range(0,len(aa), len(aa)/max_size)
        #                            times = times[inds]
        #                            vals = vals[inds]
        #
                                if alias in machine_pvs:
                                    dfs = pd.Series(vals, times).sort_index()
                                    dfs = dfs[~dfs.index.duplicated()]
                                    dfs = dfs[~(dfs.diff()==0)]
                                    vals = dfs.values
                                    dfs = dfs.to_xarray().rename({'index': 'time'})
                                    data_machine[alias] = dfs 
                                    data_machine[alias].name = alias
                                    data_machine[alias].attrs = attrs
                                    #data_machine[alias] = xr.DataArray(vals, coords=[times], dims=['time'], name=alias, attrs=attrs)
                                elif len(vals) > max_size:
                                    if not quiet:
                                        print('Skipping {:} due to too many data points {:}'.format(alias, len(vals)))
                                elif len(vals) > 1:
                                    dfs = pd.Series(vals, times).sort_index()
                                    dfs = dfs[~dfs.index.duplicated()]
                                    dfs = dfs[~(dfs.diff()==0)]
                                    vals = dfs.values
                                    dfs = dfs.to_xarray().rename({'index': 'time'})
                                    data_arrays[alias] = dfs 
                                    data_arrays[alias].name = alias
                                    data_arrays[alias].attrs = attrs
                                    #data_arrays[alias] = xr.DataArray(vals, coords=[times], dims=['time'], name=alias, attrs=attrs) 
                                else:
                                    data_points[alias] = {'value': vals[0], 'time': times[0], 
                                                          'pv': attrs['pv'], 'units': attrs.get('units','')}
                    
                    except:
                        traceback.print_exc()
                        if not quiet:
                            print(('Error loadinig', alias))

                    if not quiet:
                        try:
                            print('{:8.3f} {:28} {:8} {:10.3f} {:4} {:20} {:}'.format(time_next-time_last, \
                                            alias, len(vals), np.array(vals).mean(), units, doc, pv))
                        except:
                            print('{:8.3f} {:28} {:8} {:>10} {:4} {:20} {:}'.format(time_next-time_last, \
                                            alias, len(vals), vals[0], units, doc, pv))
                
                except:
                    traceback.print_exc()
                    if not quiet:
                        print(('Error loading', alias))

        if not quiet:
            print('... Merging')
       
        self._data_machine = data_machine
        self._data_arrays = data_arrays
        self._data_points = data_points
        self._data_fields = data_fields
        self._evr_pvs = evr_pvs

        time_last = time.time()
        xepics = xr.merge(data_arrays.values())
        xepics = xepics[sorted(xepics.data_vars)]
        xmachine = xr.merge(data_machine.values())
        xmachine = xmachine[sorted(xmachine.data_vars)]

        try:
            dfend = xr.DataArray(self.xruns.run.astype('int16'), coords=[self.xruns.end_time.values],dims=['time']).rename('end_run')
            dfbegin = xr.DataArray(self.xruns.run.astype('int16'), coords=[self.xruns.begin_time.values],dims=['time']).rename('begin_run')
            xepics = xr.merge([xepics, dfend, dfbegin])
            dfend = xepics.end_run.to_dataframe().bfill().rename(columns={'end_run': 'run'})
            dfbegin = xepics.begin_run.to_dataframe().ffill().rename(columns={'begin_run': 'brun'})
            # Need to make sure there are no duplicates -- 'cxim7216' odd case where there are using this method.
            dfend = dfend[~dfend.index.duplicated()]
            dfbegin = dfbegin[~dfbegin.index.duplicated()]
            xepics = xr.merge([xepics, dfend, dfbegin])
            xepics['prep'] = xepics.run != xepics.brun 
            xepics['run'] = xepics.run.astype('int')
            xepics['end_run'] = xepics.end_run
            xepics['begin_run'] = xepics.begin_run
            xepics = xepics.drop('brun').set_coords(['run', 'begin_run', 'end_run', 'prep'])

        except:
            traceback.print_exc()
            print('Could not merge run begin_time')
            return xepics, dfend, dfbegin

        try:
            dfend = xr.DataArray(self.xruns.run.astype('int16'), coords=[self.xruns.end_time.values],dims=['time']).rename('end_run')
            dfbegin = xr.DataArray(self.xruns.run.astype('int16'), coords=[self.xruns.begin_time.values],dims=['time']).rename('begin_run')
            xmachine = xr.merge([xmachine, dfend, dfbegin])
            dfend = xmachine.end_run.to_dataframe().bfill().rename(columns={'end_run': 'run'})
            dfbegin = xmachine.begin_run.to_dataframe().ffill().rename(columns={'begin_run': 'brun'})
            # Need to make sure there are no duplicates -- 'cxim7216' odd case where there are using this method.
            dfend = dfend[~dfend.index.duplicated()]
            dfbegin = dfbegin[~dfbegin.index.duplicated()]
            xmachine = xr.merge([xmachine, dfend, dfbegin])
            xmachine['prep'] = xmachine.run != xmachine.brun 
            xmachine['run'] = xmachine.run.astype('int')
            xmachine['end_run'] = xmachine.end_run
            xmachine['begin_run'] = xmachine.begin_run
            xmachine = xmachine.drop('brun').set_coords(['run', 'begin_run', 'end_run', 'prep'])

        except:
            traceback.print_exc()
            print('Could not merge run begin_time')
            return xmachine, dfend, dfbegin

        # add in machine info
        #phot_attrs = ['photon_current','photon_beam_energy']
        # Need to add this smarter to avoid huge arrays.
        # Either interpolate or ffill and drop, e.g., photon_current.to_dataframe().ffill()
        phot_attrs = []
        self._machine_coords = phot_attrs
        if set(phot_attrs) & set(xmachine.keys()) == set(phot_attrs):
            xepics = xepics.reset_coords().merge(xmachine[phot_attrs].reset_coords()).set_coords(phot_attrs)
            xepics = xepics.set_coords(['run', 'begin_run', 'end_run', 'prep'])
        
        self.xepics = xepics
        self.xmachine = xmachine
        time_next = time_next = time.time()
        if not quiet:
            print('{:8.3f} To merge epics data'.format(time_next-time_last))
            print('{:8.3f} Total Time to load epics data'.format(time_next-time0))
       
#        if save:
#            if not path:
#                subfolder = 'results/epicsArch/'
#                path = os.path.join('/reg/d/psdm/',self.instrument,self.experiment, subfolder)
#
#            file_name = os.paht.join(path, 
#
#            #xepics.to_netcdf( 

    def _load_machine_epicsArch(self, epics_dir=None, epics_file='epicsArch_machine.txt'):
        import os
        if not epics_dir:
            # scons needs updating to copy .txt file so use fixed path
            #epics_dir = os.path.dirname(__file__)
            epics_dir = '/reg/neh/home/koglin/conda/PyDataSource/src/'
        
        self._epicsArch_dict = epicsArch_dict(epics_file,epics_dir)
        return self._epicsArch_dict

    def _load_instrument_epicsArch(self, instrument=None, epics_dir=None, epics_file=None, 
                             quiet=False, **kwargs):
        """Load epicsArch file to define aliases.
        """
        self._epics_dict = {}
        self._epics_camdict = {}
        self._devices = {}
        self._sets = {}
        if not epics_file:
            epics_file = 'epicsArch.txt'

        if not epics_dir:
            if not instrument:
                instrument = self.instrument
            
            epics_dir = '/reg/g/pcds/dist/pds/'+instrument+'/misc/'

        if epics_dir:
            if not quiet:
                print('instrument: {:}'.format(instrument))
                print(('loading epics pvs from', epics_file, ' in', epics_dir))
            
            self._epics_dict.update(epicsArch_dict(epics_file,epics_dir))
            for item in self._epics_dict.values():
                try:
                    alias, attr = item['alias'].split('_',1)
                except:
                    alias = attr = item['alias']

                if alias not in self._sets:
                    self._sets[alias] = {}
                
                self._sets[alias].update(**{attr: item['base']})
                #self._devices.update(**{item['base']: {'alias': item['alias'], 'records': {}}})

        return self._epics_dict

    def _load_exp_runs(self, stats=['mean', 'std', 'min', 'max', 'count'], quiet=False, **kwargs):
        import numpy as np
        import xarray as xr
        import pandas as pd
        from h5write import resort
        import time

        time0 = time.time()
        if not hasattr(self, 'xepics'):
            self._load_epics(**kwargs)

        if not quiet:
            print('... processing {:} epics data'.format(self.exp))
        coords = ['begin_time','end_time','duration','run_id','prep_time']
        if kwargs.get('pvs'):
            attrs = kwargs.get('pvs').values()
        else:
            attrs = [a for a,b in self.xepics.data_vars.items() \
                     if b.attrs.get('pv','').endswith('STATE') or a.endswith('_set') 
                  or 'TRIG' in b.attrs.get('pv','')] 
                 #or a.endswith('_calc') or a.endswith('_code')]
                  #or b.attrs.get('pv','').split(':')[-1] in ['RATE','TCTL','EC_RBV','TPOL','TWID','TDES','TEC','EVENTCTRL.ENM']]
        #attrs = [a for a in self.xepics.data_vars.keys() \
                 #if a.startswith('dia_s') or (a.endswith('_set') and not a.startswith('dia_a'))]
        #xpvs = self.xepics[attrs].dropna(dim='time',how='all')
        xpvs = self.xepics[attrs]
        xcut = xpvs.where(xpvs.prep == False, drop=True)
        ca = xcut.count().to_array()
        xcut = xcut.drop(ca['variable'][ca == 0].values)
        recoords = list(set(self._machine_coords) & set(xcut.coords.keys()))
        if recoords:
            xcut = xcut.reset_coords(recoords) 
        dgrp = xcut.groupby('run')
        dsets = [getattr(dgrp, func)(dim='time') for func in stats]
        mvs = xr.concat(dsets, stats).rename({'concat_dim': 'stat'})
        xscan = resort(self.xruns.merge(mvs)).set_coords(coords)
        if recoords:
            xscan = xscan.set_coords(recoords)
        for attr in xscan.data_vars:
            xscan[attr].attrs = self.xepics[attr].attrs
 
        # find the last value set in prep before a run starts
        # if only set once during run, then also consider it to be a set value
        # since run start time occurs sometimes before set is completed.
        # in future can make this a little smarter, but generally should be valid assumption

#        xcut = xpvs.where(xpvs.run_prep>0,drop=True).drop('run').rename({'run_prep': 'run'}).set_coords('run')
        xcut = xpvs.where(xpvs.prep == True, drop=True)
        ca = xcut.count().to_array()
        xcut = xcut.drop(ca['variable'][ca == 0].values)
        xset = self.xruns.copy().set_coords(coords)
        for attr in xcut.data_vars.keys():
            xpv = xcut[attr].where(xcut.run >= 0).dropna(dim='time')
            runs = set(xpv.run.values) 
            if attr in xscan:
                setonce = xscan[attr].sel(stat='count').to_pandas() == 1
                runs = runs | set(setonce.index[setonce])
            #data = {rn: xpv.where(xpv.run==rn,drop=True).values[-1] for rn in runs}
            data = {}
            for rn in runs:
                vals = xpv.where(xpv.run==rn,drop=True).values
                if len(vals) > 0:
                    val = vals[-1]
                elif setonce[rn]:
                     val = xscan[attr].sel(run=rn).sel(stat='mean')
                data[rn] = val

            xset[attr] = xr.DataArray(data.values(), coords=[data.keys()], 
                                dims=['run'], attrs=xpvs[attr].attrs)

        self.xscan = xscan
        self.xset = xset
        self.xpvs = xpvs

        attrs = [a for a in self.xscan.data_vars.keys()]
        if xscan.data_vars.keys():
            try:
                self.dfscan = self.xscan.sel(stat='count').reset_coords()[attrs]
                self.dfscan = self.dfscan.dropna(dim='run', how='all').to_dataframe().astype(int)
            except:
                traceback.print_exc()
        else:
            self.dfscan = None

        attrs = [a for a in self.xset.data_vars.keys()]
        self.dfset = self.xset.reset_coords()[attrs]
        try:
            self.dfset = self.dfset.dropna(dim='run', how='all').to_dataframe()
        except:
            traceback.print_exc()

        if not quiet:
            time_next = time.time()
            print('{:8.3f} Total Time to load and process epics data'.format(time_next-time0))
            print('... Done loading {:} epics data'.format(self.exp))
      
    def get_scan_series(self, attr, quiet=True, min_runs=5):
        """Infer sets or runs that are parametric series based on change in attr setting for each run.
        """
        df = self.xset.to_dataframe()[attr].dropna(how='all')
        dfp = df.diff()
        dfpp = dfp.diff()
        asets = {}
        nrns = 0 
        rn0 = 0
        rn_last = 0
        for run,val in dfp.iteritems():
            if abs(dfpp[run]) > 0.1 or run-rn_last > 2 or run == df.index.values[-1]:
                if run == df.index.values[-1]:
                    nrns += 1
                    rn_last = run
                if nrns >= min_runs:
                    da = self.xset[attr]
                    da = da.where((da.run>=rn0) & (da.run<=rn_last), drop=True).dropna(dim='run').to_dataset(name='val')
                    da = da.to_dataframe()[['begin_time','duration','val']]
                    asets[(rn0,rn_last)] = da
                rn0 = run-1
                rn_last = run
                nrns = 2
                if not quiet:
                    print('- {:4} {:8.3f} {:8.3f} {:8.3f}'.format(run, df[run], dfp[run], dfpp[run]))
            else:
                nrns += 1
                rn_last = run
                if not quiet:
                    print('  {:4} {:8.3f} {:8.3f} {:8.3f}'.format(run, df[run], dfp[run], dfpp[run]))
        
        return asets

    def get_scans(self, min_steps=4, attrs=None, device=None, min_motors=1, **kwargs):
        """
        Get epics scans from dfscan information.

        Parameters
        ----------
        min_steps : int
            Minimum number of steps in scan [default = 1]
        min_motors : int
            Minimum number of motors involved in scan [default = 1]
        device : str
            Base name of device for which all pvs share 
            (e.g., devide='Sample' will include 'Sample_x', 'Sample_y', 'Sample_z')
        attrs : list, optional
            List of pvs involved
        """
        if device:
            attrs = [a for a in self.dfscan.keys() if a.startswith(device)]

        if self.dfscan is None:
            return None

        dfscan = self.dfscan.copy()
        if attrs:
            dfscan = dfscan[attrs]

        # determine attributes that make miniumum step cut
        attr_cut = (dfscan > min_steps).sum() > 0
        try:
            scan_attrs = dfscan.keys()[attr_cut]
            # Make cut on runs with minimum steps
            cut = (dfscan[scan_attrs].T > min_steps).sum() >= min_motors
            # Return reduced Dataframe with number of steps
            return dfscan[scan_attrs][cut].astype(int) 
        except:
            return None

    @property
    def run_pv_set_count(self):
        """
        Number of epics PVs set before each run.
        """
        return self.dfset.T.count()

    @property
    def pv_set_count(self):
        """
        Number of times PV was set before a run.
        """
        return self.dfset.count()

    @property
    def dfruns(self):
        """
        pandas DataFrame of Experiment run information. 
        """
        import pandas as pd
        return pd.DataFrame({a.get('num'): a for a in self.runs})

    @property
    def runs(self):
        """Experiment run information from MySQL database and xtc directory.
        """
        import glob
        from RegDB import experiment_info
        if experiment_info.name2id(self.exp):
            runs_list =  experiment_info.experiment_runs(self.instrument.upper(),self.exp)
            xtc_files = glob.glob(self.xtc_dir+'/*.xtc')
            h5_files = glob.glob(self.h5_dir+'/*.h5')
            stats_files = glob.glob(self.res_dir+'/nc/*_stats.nc')
            drop_stats_files = glob.glob(self.res_dir+'/nc/*_drop_stats.nc')
            runs_dict = {}
            for item in runs_list:
                run = item['num']
                runs_dict[run] = item
                runs_dict[run].update({'xtc_files': [], 'h5_files': [], 
                        'stats_files': [], 'drop_stats_files': [], 'scratch_nc_files': []})
            for f in xtc_files:
                try:
                    run = int(f.split('-r')[1].split('-s')[0])
                    runs_dict[run]['xtc_files'].append(f)
                except:
                    pass
            for f in stats_files:
                try:
                    run = int(f.split('run')[1].split('_stats')[0])
                    runs_dict[run]['stats_files'].append(f)
                except:
                    pass
            for f in stats_files:
                try:
                    run = int(f.split('run')[1].split('_drop_stats')[0])
                    runs_dict[run]['drop_stats_files'].append(f)
                except:
                    pass
            for f in stats_files:
                try:
                    run = int(f.split('run')[1].split('.nc')[0])
                    runs_dict[run]['scratch_nc_files'].append(f)
                except:
                    pass
            for f in h5_files:
                try:
                    run = int(f.split('-r')[1].split('-s')[0])
                    if run in runs_dict:
                        runs_dict[run]['h5_files'].append(f)
                except:
                    pass
            runs_list = runs_dict.values()

        else:
            runs_list = []

        return runs_list

    @property
    def open_files(self, run=None):
        """Return a list of files created (by the DAQ system).  
           Current run if no run is specified.
        """
        from RegDB import experiment_info
        return experiment_info.get_open_files(self.exper_id,run)

    def load_run_summary(self):
        """
        Load MySQL database experiment run summary information into a dictionary.
        
        Slow process.  Generally better information is available from epics 
        information gathered in load_exp_runs method.
        """
        from RegDB import experiment_info
        vrun_attrs = {}
        print('Loading summary of {:} runs for {:} from SQL database'.format( \
                len(self.runs),self.exp))
        print('Estimate loading time ~{:} sec'.format(len(self.runs)/4))
        for run in range(1,self.runs[-1]['num']+1):
            run_attr = experiment_info.run_attributes(self.instrument,self.exp,run)
            for a in run_attr:
                if a['name'] not in vrun_attrs:
                    vrun_attrs[a['name']] = {'class': a['class'], 'desc': a['descr'],
                                             'type': a['type'], 'val':
                                             [None for i in range(1,run)]}
                vrun_attrs[a['name']]['val'].append(a['val'])
        self.run_summary = vrun_attrs


    def show_runs(self,start=0,end=99999999,csv=False):
        """Show run summary for current experiment.
        """
        if csv:
            print('{:>7}, {:>10}, {:>8}, {:>10}, {:3}, {:2}'.format('Run',
                                'Day', 'Time', 'Length', 'xtc', 'h5'))

        else:
            print('='*72)
            print('Experiment {:}'.format(self.exp))
            print('  xtc dir {:}'.format(self.xtc_dir))
            print('  hdf5 dir {:}'.format(self.h5_dir))
            print('-'*72)
            print('{:>7} {:>10} {:>8} {:>10} {:3} {:2}'.format('Run', 'Day', 'Time',
                                                  'Length', 'xtc', 'h5'))
            print('-'*72)

        for item in self.runs:
            run = item['num']
            if run >= start and run <= end:
                datestr = time.strftime('%Y-%m-%d',
                                        time.localtime(item['begin_time_unix']))
                timestr = time.strftime('%H:%M:%S',
                                        time.localtime(item['begin_time_unix']))
                if len(item['xtc_files']) > 0:
                    xtc = 'xtc'
                else:
                    xtc = ''

                if len(item['h5_files']) > 0:
                    h5 = 'h5'
                else:
                    h5 = ''

                begin_time = item['begin_time_unix']
                end_time = item['end_time_unix']
                if end_time:
                    dtime = end_time - begin_time
                    flag = ' '
                else:
                    dtime = time.time() - begin_time
                    flag = '*'

                dmin = int(dtime/60)
                dsec = int(dtime % 60)
                if dmin > 0:
                    dtstr = '{:4}m {:02}s'.format(dmin,dsec)
                else:
                    dtstr = '{:02}s'.format(dsec)

                if csv:
                    print('{:7}, {:10}, {:8}, {:>10}, {:3}, {:2}'.format(run,
                                        datestr, timestr, dtstr, xtc, h5))
                else:
                    print('{:7} {:10} {:8} {:>10} {:3} {:2}'.format(run,
                                        datestr, timestr, dtstr, xtc, h5))

                if flag in '*':
                    print('* Currently Acquiring Data for Run {:}'.format(run))

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__,str(self))
        return '< '+repr_str+' >'
    
    def __str__(self):
        return  str('{:}'.format(self.exp))


