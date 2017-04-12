
from RegDB import experiment_info
from glob import glob
import re
import operator
import sys
import os
import time
import traceback
from pylab import *

def write_exp_summary(self, file_name=None, path=None):
    """Write
    """
    import cPickle as pickle
    if not file_name:
        file_name = 'experiment_summary.pkl'
    if not path:
        path = os.path.join(self.res_dir,'RunSummary') 

    if not os.path.isdir(path):
        os.mkdir(path)
    
    with open(path+'/'+file_name, 'wb') as pickle_file:
        pickle.dump(pickle.dumps(self, protocol=-1), pickle_file, protocol=-1)

def read_exp_summary(exp=None, file_name=None, path=None, **kwargs):
    import cPickle as pickle
    if not file_name:
        file_name = 'experiment_summary.pkl'
    if not path:
        instrument = exp[0:3]
        exp_dir = "/reg/d/psdm/{:}/{:}".format(instrument, exp)
        if os.path.isdir(os.path.join(exp_dir,'res')): 
            path = os.path.join(exp_dir,'res','RunSummary')
        else:
            path = os.path.join(exp_dir,'results','RunSummary')
    
    with open(path+'/'+file_name, 'rb') as pickle_file:
        data = pickle.load(pickle_file)

    return pickle.loads(data)

def get_exp_summary(exp, load=True, **kwargs):
    if load:
        try:
            return read_exp_summary(exp, **kwargs)
        except:
            traceback.print_exc()
            load = False

    if not load:
        return ExperimentSummary(exp, **kwargs)


class ExperimentSummary(object):
    """
    Experiment information from elog and epics archive data. 

    Attributes
    ----------
    xscan : xarray.Dataset
        

    """

    _fields={
            'description':            ('DESC', 'Description'), 
            'device_type':            ('DTYP', 'Device type'), 
            'record_type':            ('RTYP', 'Record Type'), 
            'slew_speed':             ('VELO', 'Velocity (EGU/s) '),
            'acceleration':           ('ACCL', 'acceleration time'),
            'step_size':              ('RES',  'Step Size (EGU)'),
            'encoder_step':           ('ERES', 'Encoder Step Size '),
            'resolution':             ('MRES', 'Motor Step Size (EGU)'),
            'high_limit':             ('HLM',  'User High Limit  '),
            'low_limit':              ('LLM',  'User Low Limit  '),
            }

    def __init__(self, exp=None, instrument=None, station=0, exper_id=None,
                exp_dir=None, xtc_dir=None, h5_dir=None, save=True, **kwargs):
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
        
        self.exp_dir = exp_dir
        self.h5_dir = h5_dir
        self.xtc_dir = xtc_dir
        if os.path.isdir(os.path.join(self.exp_dir,'res')): 
            self.res_dir = os.path.join(self.exp_dir,'res')
        else:
            self.res_dir = os.path.join(self.exp_dir,'results')

        try:
            self._load_exp_runs(**kwargs)
        except:
            traceback.print_exc()
            print 'could not load exp_runs'
            
        if save:
            try:
                write_exp_summary(self)
            except:
                traceback.print_exc()
                print 'could not write summary'

#    def _add_runtables(self):
#        """
#        Add RunTables... currently does not work in conda env.
#        """
#        try:
#            from LogBook.runtables import RunTables
#
#            # Setup elog pswwww RunTables
#            self._RunTables = RunTables(**{'web-service-url': 'https://pswww.slac.stanford.edu/ws-kerb'})
#
#            self._user_tables = self._RunTables.usertables(exper_name=self.exp)
#            #self.add_user_run_table(RunSummary='Run Summary')
#        except:
#            print 'currently does not work in conda env.'

#    def add_user_run_table(self, **kwargs):
#        """Add a user User RunTable from pswww elog server.
#        """
#        for alias, name in kwargs.items():
#            tbl = self._RunTables.findUserTable(exper_name=self.exp, table_name=alias)
#            self._user_tables.update({alias: tbl}) 
#            setattr(self, alias, tbl)

    def to_html(self, **kwargs):
        import build_html
        return build_html.Build_experiment(self, **kwargs)

    def detectors(self, run):
        """Return a list of detector names configured in the DAQ system for the input run number.
        """
        return experiment_info.detectors(self.instrument, self.exp, run)

    @property
    def calibration_runs(self):
        return experiment_info.calibration_runs(self.instrument, self.exp)

    def show_moved(self, attrs=None):
        if attrs:
            if not isinstance(attrs, list):
                attrs = [attrs]
        else:
            attrs = self.xset.data_vars.keys()

        for attr in attrs:
            xvar = self.xset.reset_coords().set_coords(['begin_time','prep_time'])[attr].dropna(dim='run')
            print ''
            print '** {:} **'.format(attr)
            for a,val in xvar.attrs.items():
                if a not in ['PREC']:
                    print '{:10} {:10}'.format(a,val)
            print xvar.to_dataframe()

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

    def get_scan_data(self, run, **kwargs):
        """Get scan dataframe.
        """
        df = self.get_scans(**kwargs).T
        if run in df:
            df = df[run]
            attrs = df.where(df > 0).dropna().keys() 
            a = self.xruns.sel(run=run)
            t = self.xepics.time
            return self.xepics[attrs].where(t>=a.begin_time).where(t<=a.end_time).dropna(dim='time',how='all')
        else:
            return None

    def plot_scan(self, run, style=None, linewidth=2,
            min_steps=4, attrs=None, device=None, min_motors=1,
            figsize=(12,8), loc='best', **kwargs):
        """
        Plot motor moves vs time and run.
        """
        fig = plt.figure(figsize=figsize, **kwargs)
        ax = fig.gca()
        x = self.get_scan_data(run, attrs=attrs, min_steps=min_steps, device=device, min_motors=min_motors)
        if not x:
            print 'Run {:} scan for attrs={:} device={:} not valid'.format(run, attrs, device)
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
                
                df = x[attr].dropna(dim='time').to_pandas()
                if style:
                    df.plot(style=style, linewidth=linewidth, label=lab, ax=ax)
                else:
                    df.plot(drawstyle='steps', linewidth=linewidth, label=lab, ax=ax)
                
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

    def get_moved(self, attrs=None, run_min=None, run_max=None, group=None):
        """Get devices that moved (between run_min and run_max if specified)
        """
        if attrs:
            group = False

        if isinstance(attrs, list):
            for attr in set(attrs):
                if not attr.endswith('_set'):
                    attrs.append(attr+'_set')
        
        xepics = self.get_epics(attrs=attrs, run_min=run_min, run_max=run_max)
        df = xepics.to_dataframe()

        attrs = []
        for a in df.keys()[df.count() > 1]:
            #if a not in xepics.coords and (a in self.xset.data_vars or a in self.xscan.data_vars):
            if a not in xepics.coords:
                if a.endswith('_set'):
                    attrs.append(a.split('_set')[0])
                else:
                    attrs.append(a)
        
        attrs = list(set(attrs))

        if group is not False:
            if group == 'units':
                adets = {}
                for attr in attrs:
                    det = attr.split('_')[0] 
                    try:
                        unit = xepics[attr].attrs.get('units','')
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
            figsize=(12,8), ax_pos = (.1,.2,.55,.7), box_to_anchor=(1.1, 1.), **kwargs):
        """
        Plot motor moves vs time and run.
        """
        if not isinstance(attrs, list):
            xepics = self.get_epics(run_min=run_min, run_max=run_max)
            if attrs not in xepics.data_vars:
                ylabel = attrs
                print attrs
                attrs = self.get_moved(attrs, run_min=run_min, run_max=run_max, group=False)
                #attrs = [a for a in xepics.data_vars if a.startswith(attrs)]
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
            print attr
            try:
                df = xepics[attr].dropna(dim='time').to_pandas()
                if style:
                    df.plot(style=style, linewidth=linewidth, label=lab, ax=ax)
                else:
                    df.plot(drawstyle='steps', linewidth=linewidth, label=lab, ax=ax)
            except:
                print 'cannot plot', attr

            if attr+'_set' in xepics:
                try:
                    xepics[attr+'_set'].dropna(dim='time')[1:].to_pandas().plot(style='+',label=attr+'_set')
                except:
                    print 'cannot plot', attr+'_set'

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

    def last_pv_fields_changed(self):
        """
        Last time fields were changed.
        """
        pass

    def last_set(self, attrs=['value', 'units', 'run', 'description'], **kwargs):
        """Show the last time devices were set.
        """
        import pandas as pd
        #print '{:30} {:10} {:>6}   {:30}'.format('Name','Value', 'Run', 'PV')
        #print '-'*72
        avals = {}
        for attr in self.xset.data_vars:
            xvar = self.xset[attr].dropna(dim='run')
            avals[attr] = xvar.attrs
            avals[attr].update({'value': float(xvar[-1].values),'run': xvar[-1].run.values, 'time': xvar[-1].begin_time.values})
            #print '{:30} {:10.3f} {:>6}   {:30}'.format(attr, float(xvar[-1].values), xvar[-1].run.values, xvar.attrs.get('pv'))
                
        return pd.DataFrame(avals).T[attrs]

    def _load_run_data(self, summary=False, no_create=True):
        ax = {}
        for i in self.xruns.run.values:                         
            try:                                         
                x = PyDataSource.h5write.get_config_xarray(exp='mfx11116',run=i, summary=summary, no_create=no_create)
                if x:
                    ax[i] = x
            except:     
                print 'cannot load run', i

        self.xdata = ax
        return self.xdata


    def _load_epics(self, pvs=None, quiet=False, 
                    omit=['ABSE', 'SIOC', 'USEG', 'VGBA', 'GATT', 
                          'MIRR:FEE1:M2H.RBV', 'MIRR:FEE1:M1H.RBV'],
                    run=None,
                    **kwargs):
        """
        Make xarray and pandas objects to describe which epics variables were set before
        and during runs.
        """
        import numpy as np
        import xarray as xr
        import pandas as pd
        from epicsarchive import EpicsArchive
        arch = EpicsArchive()
        rns = pd.DataFrame(experiment_info.experiment_runs(self.instrument.upper(),self.exp))
        
        if not pvs:
            import PyDataSource
            if not run:
                run=rns['num'].max()
            ds = PyDataSource.DataSource(exp=self.exp, run=run)
            pvnames = {pv:ds.epicsData.alias(pv) for pv in ds.epicsData.pvNames()}
            for a in omit:
                for pv in pvnames.copy():
                    if pv.startswith(a):
                        pvnames.pop(pv)
        else:
            pvnames = pvs

        pvmots = {pv.rstrip('.RBV'): alias+'_set' for pv, alias in pvnames.items() if pv.endswith('RBV')}
        pvnames.update(**pvmots)
        pvs = {pv: alias for pv, alias in pvnames.items() if arch.search_pvs(pv, do_print=False) != []} 
    
        xruns = xr.Dataset()
        xruns.coords['run'] = (['run'], rns['num'])
        xruns.coords['run_id'] = (['run'], rns['id'])
        xruns['begin_time'] = (['run'], pd.to_datetime(rns['begin_time']))
        xruns['end_time'] = (['run'], pd.to_datetime(rns['end_time']))
        xruns['duration'] = (['run'], (rns['end_time']/1.e9 - rns['begin_time']/1.e9))
        xruns['duration'].attrs['units'] = 'sec'
        xruns['prep_time'] = (['run'], (rns['begin_time']/1.e9)-(rns['end_time'].shift(1)/1.e9)) 
        xruns['prep_time'].attrs['units'] = 'sec'

        df = xruns.to_dataframe()
        dt = df['begin_time'].min()
        tstart = [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]
        dt = df['end_time'].max()
        tend = [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]

        meta_attrs = {'units': 'EGU', 'PREC': 'PREC', 'pv': 'name'}

        data_arrays = {} 
        data_points = {}
        data_fields = {}
        import time
        time0 = time.time()
        time_last = time0
        for pv, alias in pvs.items():
            data_fields[alias] = {}
            try:
                dat = arch._get_json(pv, tstart, tend, False)[0]
                attrs = {a: dat['meta'].get(val) for a,val in meta_attrs.items() if val in dat['meta']}
                for attr, item in self._fields.items():  
                    try:
                        field=item[0]
                        pv_desc = pv.split('.')[0]+'.'+field
                        desc = arch._get_json(pv_desc, tstart, tend, False)[0]
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
                           
                            data_fields[alias][attr] = xr.DataArray(vals.values(), coords=[vals.keys()], dims=['time'], 
                                                                    name=alias+'_'+attr, attrs=fattrs) 
                            attrs[attr] = val
     
                    except:
                        pass
                vals = [item['val'] for item in dat['data']]
                doc = attrs.get('description','')
                units = attrs.get('units', '')
                time_next = time.time()
                if not quiet:
                    try:
                        print '{:8.3f} {:28} {:8} {:10.3f} {:4} {:20} {:}'.format(time_next-time_last, \
                                        alias, len(vals), np.array(vals).mean(), units, doc, pv)
                    except:
                        print '{:8.3f} {:28} {:8} {:>10} {:4} {:20} {:}'.format(time_next-time_last, \
                                        alias, len(vals), vals[0], units, doc, pv)
               
                try:
                    if isinstance(vals[0],str):
                        if not quiet:
                            print alias, 'string'
                        vals = np.array(vals, dtype=str)
                    else:
                        times = [np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns') for item in dat['data']]
                        # Simple way to remove duplicates
                        # Duplicates will break xarray merge
                        aa = dict(zip(times, vals))
                        times = np.array(aa.keys(), dtype=times[0].dtype)
                        vals = np.array(aa.values())
                        if len(vals) > 1:
                            data_arrays[alias] = xr.DataArray(vals, coords=[times], dims=['time'], name=alias, attrs=attrs) 
                        else:
                            data_points[alias] = {'value': vals[0], 'time': times[0], 'pv': attrs['pv'], 'units': attrs.get('units','')}
                except:
                    if not quiet:
                        print 'Error loadinig', alias
            except:
                if not quiet:
                    print 'Error loading', alias

        if not quiet:
            print '... Merging'
       
        self.xruns = xruns
        self._data_arrays = data_arrays
        self._data_points = data_points
        self._data_fields = data_fields

        time_last = time.time()
        xepics = xr.merge(data_arrays.values())
        xepics = xepics[sorted(xepics.data_vars)]
        try:
            dfend = xr.DataArray(self.xruns.run.astype('int16'), coords=[self.xruns.end_time.values],dims=['time']).rename('end_run')
            dfbegin = xr.DataArray(self.xruns.run.astype('int16'), coords=[self.xruns.begin_time.values],dims=['time']).rename('begin_run')
            xepics = xr.merge([xepics, dfend, dfbegin])
            dfend = xepics.end_run.to_dataframe().bfill().rename(columns={'end_run': 'run'})
            dfbegin = xepics.begin_run.to_dataframe().ffill().rename(columns={'begin_run': 'brun'})
            xepics = xr.merge([xepics, dfend, dfbegin])
            xepics['prep'] = xepics.run != xepics.brun 
            xepics['run'] = xepics.run.astype('int')
            xepics['end_run'] = xepics.end_run
            xepics['begin_run'] = xepics.begin_run
            xepics = xepics.drop('brun').set_coords(['run', 'begin_run', 'end_run', 'prep'])

        except:
            traceback.print_exc()
            print 'Could not merge run begin_time'

        self.xepics = xepics

        time_next = time_next = time.time()
        if not quiet:
            print '{:8.3f} To merge epics data'.format(time_next-time_last)
            print '{:8.3f} Total Time to load epics data'.format(time_next-time0)
       
#        if save:
#            if not path:
#                subfolder = 'results/epicsArch/'
#                path = os.path.join('/reg/d/psdm/',self.instrument,self.experiment, subfolder)
#
#            file_name = os.paht.join(path, 
#
#            #xepics.to_netcdf( 


    def _load_exp_runs(self, quiet=False, **kwargs):
        import numpy as np
        import xarray as xr
        import pandas as pd
        from h5write import resort

        if not hasattr(self, 'xepics'):
            self._load_epics(**kwargs)

        stats=['mean', 'std', 'min', 'max', 'count']
        attrs = [a for a,b in self.xepics.data_vars.items() \
                 if b.attrs.get('pv','').endswith('STATE') or a.endswith('_set')]
        #attrs = [a for a in self.xepics.data_vars.keys() \
                 #if a.startswith('dia_s') or (a.endswith('_set') and not a.startswith('dia_a'))]
        #xpvs = self.xepics[attrs].dropna(dim='time',how='all')
        xpvs = self.xepics[attrs]
        xcut = xpvs.where(xpvs.prep == False, drop=True)
        ca = xcut.count().to_array()
        xcut = xcut.drop(ca['variable'][ca == 0].values)
        dgrp = xcut.groupby('run')
        dsets = [getattr(dgrp, func)(dim='time') for func in stats]
        mvs = xr.concat(dsets, stats).rename({'concat_dim': 'stat'})
        coords = ['begin_time','end_time','duration','run_id','prep_time']
        xscan = resort(self.xruns.merge(mvs)).set_coords(coords)
        for attr in xscan.data_vars:
            xscan[attr].attrs = self.xepics[attr].attrs
 
#        xcut = xpvs.where(xpvs.run_prep>0,drop=True).drop('run').rename({'run_prep': 'run'}).set_coords('run')
        xcut = xpvs.where(xpvs.prep == True, drop=True)
        ca = xcut.count().to_array()
        xcut = xcut.drop(ca['variable'][ca == 0].values)
        xset = self.xruns.copy().set_coords(coords)
        for attr in xcut.data_vars.keys():
            xpv = xcut[attr].where(xcut.run >= 0).dropna(dim='time')
            runs = set(xpv.run.values)
            data = {rn: xpv.where(xpv.run==rn,drop=True).values[-1] for rn in runs}
            xset[attr] = xr.DataArray(data.values(), coords=[data.keys()], dims=['run'], attrs=xpvs[attr].attrs)

        self.xscan = xscan
        self.xset = xset
        self.xpvs = xpvs

        attrs = [a for a in self.xscan.data_vars.keys()]
        self.dfscan = self.xscan.sel(stat='count').reset_coords()[attrs].dropna(dim='run', how='all').to_dataframe().astype(int)
        
        attrs = [a for a in self.xset.data_vars.keys()]
        self.dfset = self.xset.reset_coords()[attrs].dropna(dim='run', how='all').to_dataframe()

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
    def runs(self):
        """Experiment run information from MySQL database and xtc directory.
        """
        if experiment_info.name2id(self.exp):
            runs_list =  experiment_info.experiment_runs(self.instrument.upper(),self.exp)
            for item in runs_list:
                runnum = item['num']
                item['xtc_files'] = glob('{:}/*-r{:04d}*.xtc'.format(
                                        self.xtc_dir,runnum))
                item['h5_files'] = glob('{:}/*-r{:04d}*.h5'.format(
                                        self.h5_dir,runnum))
        else:
            runs_list = []

        return runs_list

    @property
    def open_files(self, run=None):
        """Return a list of files created (by the DAQ system).  
           Current run if no run is specified.
        """
        return experiment_info.get_open_files(self.exper_id,run)

    def load_run_summary(self):
        """
        Load MySQL database experiment run summary information into a dictionary.
        
        Slow process.  Generally better information is available from epics 
        information gathered in load_exp_runs method.
        """
        vrun_attrs = {}
        print 'Loading summary of {:} runs for {:} from SQL database'.format( \
                len(self.runs),self.exp)
        print 'Estimate loading time ~{:} sec'.format(len(self.runs)/4)
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
            print '{:>7}, {:>10}, {:>8}, {:>10}, {:3}, {:2}'.format('Run',
                                'Day', 'Time', 'Length', 'xtc', 'h5')

        else:
            print '='*72
            print 'Experiment {:}'.format(self.exp)
            print '  xtc dir {:}'.format(self.xtc_dir)
            print '  hdf5 dir {:}'.format(self.h5_dir)
            print '-'*72
            print '{:>7} {:>10} {:>8} {:>10} {:3} {:2}'.format('Run', 'Day', 'Time',
                                                  'Length', 'xtc', 'h5')
            print '-'*72

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
                    print '{:7}, {:10}, {:8}, {:>10}, {:3}, {:2}'.format(run,
                                        datestr, timestr, dtstr, xtc, h5)
                else:
                    print '{:7} {:10} {:8} {:>10} {:3} {:2}'.format(run,
                                        datestr, timestr, dtstr, xtc, h5)

                if flag in '*':
                    print '* Currently Acquiring Data for Run {:}'.format(run)

