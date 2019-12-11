
# File and Version Information:
#  $HeaderURL$
#  $Id$
#  $LastChangedDate$
#
# Description:
#  module PyDataSource
#--------------------------------------------------------------------------
"""
Python implementation of psana DataSource object.

Example
-------

    # Import the PyDataSource module
    In [1]: import PyDataSource

    # Load example run
    In [2]: ds = PyDataSource.DataSource('exp=xpptut15:run=54')

    # Access the first event
    In [3]: evt = ds.events.next()

    # Tab to see Data objects in current event
    In [4]: evt.
    evt.EBeam            evt.FEEGasDetEnergy  evt.XppEnds_Ipm0     evt.cspad            evt.next
    evt.EventId          evt.L3T              evt.XppSb2_Ipm       evt.get              evt.yag2
    evt.Evr              evt.PhaseCavity      evt.XppSb3_Ipm       evt.keys             evt.yag_lom

    # Tab to see EBeam attributes
    In [4]: evt.EBeam.
    evt.EBeam.EventId            evt.EBeam.ebeamEnergyBC1     evt.EBeam.ebeamPhotonEnergy  evt.EBeam.epicsData
    ...
    evt.EBeam.ebeamDumpCharge    evt.EBeam.ebeamLTUPosY       evt.EBeam.ebeamXTCAVPhase    

    # Print a table of the EBeam data for the current event
    In [4]: evt.EBeam.show_info()
    --------------------------------------------------------------------------------
    EBeam xpptut15, Run 54, Step -1, Event 0, 11:37:12.4517, [140, 141, 41, 40]
    --------------------------------------------------------------------------------
    damageMask                 1.0486e+06         Damage mask.
    ebeamCharge                0.00080421 nC      Beam charge in nC.
    ...
    ebeamPhotonEnergy                   0 eV      computed photon energy, in eV
    ...

    # Print summary of the cspad detector (uses PyDetector methods for creatining calib and image data)
    In [5]: evt.cspad.show_info()
    --------------------------------------------------------------------------------
    cspad xpptut15, Run 54, Step -1, Event 0, 11:37:12.4517, [140, 141, 41, 40]
    --------------------------------------------------------------------------------
    calib                <0.010653> ADU     Calibrated data
    image               <0.0081394> ADU     Reconstruced 2D image from calibStore geometry
    raw                    <1570.2> ADU     Raw data
    shape              (32, 185, 388)         Shape of raw data array
    size                  2.297e+06         Total size of raw data

    # Print summary of cspad detector calibration data (using PyDetector access methods) 
    In [6]: evt.cspad.calibData.show_info()
    areas                  <1.0077>         Pixel area correction factor
    bkgd                      <0.0>         
    ...
    shape              (32, 185, 388)         Shape of raw data array
    size                  2.297e+06         Total size of raw data
    status             <0.00069396>         

    # Print summary of cspad detector calibration data (using PyDetector access methods) 
    In [7]: evt.cspad.configData.show_info()
    activeRunMode                       3         
    asicMask                           15         
    ...
    roiMasks                   0xffffffff         
    runDelay                        58100         
    tdi                                 4         

This software was developed for the LCLS project.
If you use all or part of it, please give an appropriate acknowledgment.

@version $Id$

@author Koglin, Jason
"""
from __future__ import print_function
from __future__ import absolute_import
#------------------------------
__version__ = "$Revision$"
##-----------------------------

import os
import sys
import operator
import imp
import re
import time
import traceback
import inspect

#from pylab import *
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# psana modules
import psana

from psmon import publish
publish.client_opts.daemon = True
from psmon.plots import Image, XYPlot, MultiPlot

# PyDataSource modules
from .DataSourceInfo import *
from .psana_doc_info import * 
from .psmessage import Message
from .build_html import Build_html
from .exp_summary import ExperimentSummary
from .exp_summary import get_exp_summary 
from .arp_tools import post_report
from .h5write import *

_eventCodes_rate = {
        40: '120 Hz',
        41: '60 Hz',
        42: '30 Hz',
        43: '10 Hz',
        44: '5 Hz',
        45: '1 Hz',
        46: '0.5 Hz',
        140: 'Beam & 120 Hz',
        141: 'Beam & 60 Hz',
        142: 'Beam & 30 Hz',
        143: 'Beam & 10 Hz',
        144: 'Beam & 5 Hz',
        145: 'Beam & 1 Hz',
        146: 'Beam & 0.5 Hz',
        150: 'Burst',
        162: 'BYKIK',
        163: 'BAKIK',
        }

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

def getattr_complete(base, args):
    """
    Recursive getattr
    
    Parameters
    ----------
    base : object
        Object

    args : str
        Argument to get from object
    """
    attrs = args.split('.')
    while len(attrs) > 0:
        base = getattr(base, attrs.pop(0))

    return base

def import_module(module_name, module_path):
    """
    Import a module from a given path.

    Parameters
    ----------
    module_name : str
        Name of module
    module_path : str
        Path of module
    """
    try:
        if not isinstance(module_path,list):
            module_path = [module_path]
        file,filename,desc = imp.find_module(module_name, module_path)
        globals()[module_name] = imp.load_module(module_name, file, filename, desc)
        return
    except Exception as err:
        print(('import_module error', err))
        print('ERROR loading {:} from {:}'.format(module_name, module_path))
        traceback.print_exc()

    sys.exit()

def get_module(module_name, module_path, reload=False):
    """
    Get module from path

    Parameters
    ----------
    module_name : str
        Name of module
    module_path : str
        Path of module
    reload : bool
        Reload module
    """
    if reload or module_name not in globals():
        import_module(module_name, module_path)
    
    try:
        new_class =  getattr(globals()[module_name],module_name)
    except:
        new_class =  getattr(globals()[module_name],module_name.capitalize())

    return new_class

def get_key_info(psana_obj):
    """Get a dictionary of the (type, src, key) for the data types of each src.
    """
    key_info = {}
    for key in psana_obj.keys():
        typ = key.type()
        src = key.src()
        if typ:
            srcstr = str(src)
            if srcstr not in key_info:
                key_info[srcstr] = [] 
            key_info[srcstr].append((typ, src, key.key()))

    return key_info

def get_keys(psana_obj):
    """Get a dictionary of the (type, src, key) for the data types of each src.
    """
    key_info = {}
    _modules = {}
    for key in psana_obj.keys():
        typ = key.type()
        src = key.src()
        if typ:
            srcstr = str(src)
            if srcstr not in key_info:
                key_info[srcstr] = [] 
            
            key_info[srcstr].append((typ, src, key.key()))
            
            type_name = typ.__name__
            module = typ.__module__.lstrip('psana.')
            if module:
                if module not in _modules:
                    _modules[module] = {}
                
                if type_name not in _modules[module]:
                    _modules[module][type_name] = []

                _modules[module][type_name].append((typ, src, key.key()))

    return key_info, _modules


def _repr_value(value):
    """Represent a value for use in show_info method.
    """
    if isinstance(value,str):
        return value
    else:
        if isinstance(value, list):
            if len(value) > 4:
                return 'list'
            else:
                return ', '.join(str(a) for a in value)
                #return str(value)
        
        elif isinstance(value, int):
            return str(value)
        
        elif hasattr(value, 'mean') and value.size > 4:
            try:
                return '<{:.4}>'.format(value.mean())
            except:
                return str(value)
        
        else:
            try:
                return ', '.join(str(a) for a in value)
                #return '{:10.5g}'.format(value)
            except:
                try:
                    return str(value)
                except:
                    return value

def _is_psana_type(value):
    """True if the input is a psana data type
    """
    return hasattr(value, '__module__') and value.__module__.startswith('psana')

def _get_typ_func_attr(typ_func, attr, nolist=False):
    """Return psana functions as properties.
    """
    value = getattr(typ_func, attr)
    module = typ_func.__module__.lstrip('psana.')
    type_name = typ_func.__class__.__name__
    try: 
        info = psana_doc_info[module][type_name].get(attr, {'unit': '', 'doc': ''}).copy()
    except:
        info = {'unit': '', 'doc': ''}

    info['typ_func'] = typ_func
    info['attr'] = attr

    if info.get('func_shape'):
        nvals = info.get('func_shape')
        if isinstance(nvals, str):
            nvals = getattr(typ_func, nvals)()[0]
       
        i0 = info.get('func0',0)
        try:
            value = [value(i+i0) for i in range(nvals)]
        except:
            pass

    elif info.get('func_len_hex'):
        nvals = getattr(typ_func, info.get('func_len_hex'))()
        try:
            value = [hex(value(i)) for i in range(nvals)]
        except:
            pass

    elif info.get('func_len'):
        nvals = info.get('func_len')
        if isinstance(nvals, str):
            nvals = getattr(typ_func, nvals)()
        
        try:
            value = [value(i) for i in range(nvals)]
        except:
            pass

    elif info.get('func_index'):
        vals = getattr(typ_func, info.get('func_index'))()
        try:
            value = [value(int(i)).name for i in vals]
        except:
            pass

    elif 'func_method' in info:
        info['value'] = info.get('func_method')(value())
        return info


    if hasattr(value, '_typ_func') and str(value._typ_func)[0].islower():
        # evaluate as name to avoid recursive psana functions 
        if 'name' in value._attrs and 'conjugate' in value._attrs:   
            info['value'] = value.name
  
    try:
        value = value()
        if hasattr(value, 'name'):
            info['value'] = value.name
            return info
    except:
        pass

    if isinstance(value, list):
        values = []
        is_type_list = False
        nvals = info.get('list_len', len(value))
        if isinstance(nvals, str):
            nvals = getattr(typ_func, nvals)()
       
        for i in range(nvals):
            val = value[i]
            if _is_psana_type(val):
                values.append(PsanaTypeData(val))
                is_type_list = True
            else:
                values.append(val)

        if is_type_list: # and not nolist:
            values = PsanaTypeList(values)

        info['value'] = values
        return info

    if _is_psana_type(value):
        info['value'] = PsanaTypeData(value)
    else:
        info['value'] = value

    return info

def psmon_publish(evt, quiet=True):
    """
    Publish psmon plots for an event
    """
    import numpy as np
    from psmon import publish
    from psmon.plots import Image, XYPlot, MultiPlot
    eventCodes = evt.Evr.eventCodes
    event_info = str(evt)
    for alias in evt._attrs:
        psplots = evt._ds._device_sets.get(alias, {}).get('psplot')
        if psplots:
            detector = evt._dets.get(alias)
            for name, psmon_args in psplots.items():
                eventCode = psmon_args['pubargs'].get('eventCode', None)
                nskip = psmon_args['pubargs'].get('nskip', None)
                reshape = psmon_args['pubargs'].get('reshape', None)
                transpose = psmon_args['pubargs'].get('transpose', None)
                fliplr = psmon_args['pubargs'].get('fliplr', None)
                flipud = psmon_args['pubargs'].get('flipud', None)
                flipud = psmon_args['pubargs'].get('flipud', None)
                rot90 = psmon_args['pubargs'].get('rot90', None)
                if nskip and (evt._ds._ievent % nskip != 0):
                    continue

                if not quiet:
                    print((eventCode, name, psmon_args))
                
                if eventCode is None or eventCode in eventCodes:
                    psplot_func = psmon_args['plot_function']
                    psmon_fnc = None
                    if psplot_func is 'Image':
                        image = getattr_complete(detector, psmon_args['attr'][0])
                        if reshape:
                            image = image.reshape(reshape)
                        if transpose:
                            image = image.transpose()
                            if not quiet:
                                print('transpose')
                        if flipud: 
                            image = np.flipud(image)
                            if not quiet:
                                print('flipud')
                        if fliplr: 
                            image = np.fliplr(image)
                            if not quiet:
                                print('fliplr')
                        if rot90: 
                            if isinstance(rot90, int):
                                krot = rot90
                            else: krot = 1
                            image = np.rot90(image, krot)
                            if not quiet:
                                print(('rot90', krot))
                        
                        if not quiet:
                            print((name, image, image.shape, psmon_args))
                        
                        if image is not None:
                            psmon_fnc = Image(
                                        event_info,
                                        psmon_args['title'],
                                        np.array(image, dtype='f'), 
                                        #np.array(image, dtype='f').transpose(), 
                                        **psmon_args['kwargs'])
                    
                    elif psplot_func is 'XYPlot':
                        ydata = np.array([getattr_complete(detector, attr) \
                                for attr in psmon_args['attr']], dtype='f').squeeze()
                        if ydata is not None:
                            if not quiet:
                                print(event_info)
                                print(psmon_args['title'])
                                print(psmon_args['xdata'].shape)
                                print(np.array(ydata, dtype='f'))
                            
                            psmon_fnc = XYPlot(
                                        event_info,
                                        psmon_args['title'],
                                        psmon_args['xdata'],
                                        ydata,
                                        **psmon_args['kwargs'])

                    if psmon_fnc:
                        #print 'publish', name, event_info, psmon_args
                        #print psmon_fnc
                        pub_info = publish.send(name,psmon_fnc)
                        psmon_args['psmon_fnc'] = psmon_fnc


class ScanData(object):
    """
    Scan configuration for Run
    
    Information from daq-scan, where things like motor positions are changed during a run 
    in "steps" or "calibcycles".

    Parameters
    ----------
    ds : object
        PyDataSource.DataSource object

    Attributes
    ----------
    nevents : list
        Number of events for each step

    """
    _array_attrs = ['pvControls_value', 'pvMonitors_loValue', 'pvMonitors_hiValue']
    _uses_attrs = ['uses_duration', 'uses_events', 'uses_l3t_events']
    _npv_attrs = ['npvControls', 'npvMonitors']

    def __init__(self, ds, quiet=True):
        import numpy as np
        self._ds = ds
        self._attrs = sorted(ds.configData.ControlData._all_values.keys())
        self._scanData = {attr: [] for attr in self._attrs}
        ds.reload()
        self.nsteps = ds._idx_nsteps
        ievent_start = []
        start_times = []
        start_datetimes = []
        ievent_end = []
        end_times = []
        end_datetimes = []
        self._ds._idx_istep = []
        istep = 0
        print('...Building ScanData configuration...')
        for step in ds.steps:
            okevt = False
            while not okevt:
                evt = next(step)
                ttup = (evt.EventId.sec, evt.EventId.nsec, evt.EventId.fiducials)
                okevt = ttup in ds._idx_times_tuple
                if not quiet:
                    print(('step:', istep, evt))

            ievent = ds._idx_times_tuple.index(ttup)
            ievent_start.append(ievent)
            if istep > 0:
                ievent_end.append(ievent-1)
            start_times.append(evt.EventId.timef64) 
            start_datetimes.append(evt.EventId.datetime64) 
 
            for attr in self._attrs:
                self._scanData[attr].append(ds.configData.ControlData._all_values[attr])
       
            istep += 1

        ievent_end.append(len(ds._idx_times_tuple)-1)       
        for istep, ievent in enumerate(ievent_end):
            end_times.append(ds.events.next(ievent).EventId.timef64)
            end_datetimes.append(ds.events.next(ievent).EventId.datetime64)
        
        self._scanData['ievent_start'] = np.array(ievent_start)
        self._scanData['ievent_end'] = np.array(ievent_end)
        self.nevents = np.array(ievent_end)-np.array(ievent_start)+1 
        self.start_times = np.array(start_times)
        self.end_times = np.array(end_times)
        self.step_times = np.array(end_times) - np.array(start_times)
        self.start_datetimes = np.array(start_datetimes)
        self.end_datetimes = np.array(end_datetimes)
        
        # Save lookup of step
        for istep, n in enumerate(self.nevents):
            for i in range(n):
                self._ds._idx_istep.append(istep)
 
        for attr in self._uses_attrs:
            setattr(self, attr, all(self._scanData.get(attr)))

        if (self.uses_duration or self.uses_events or self.uses_l3t_events) \
                and len(set(self._scanData['npvControls'])) == 1 \
                and len(set(self._scanData['npvMonitors'])) == 1 :
            self._is_simple = True
            for attr in self._npv_attrs:
                setattr(self, attr, self._scanData.get(attr)[0]) 

            if 'pvControls_name' in self._attrs:
                self.pvControls = self._scanData['pvControls_name'][0]
            else:
                self.pvControls = []
            
            if 'pvMonitors_name' in self._attrs:
                self.pvMonitors = self._scanData['pvMonitors_name'][0]
            else:
                self.pvMonitors = None

            self.pvLabels = self._scanData['pvLabels'][0]
            if not self.pvLabels:
                self.pvLabels = []
                for pv in self.pvControls:
                    alias = ds.epicsData.alias(pv)
                    if not alias:
                        alias = pv

                    self.pvLabels.append(alias) 
            
            self.control_values = {} 
            self.monitor_hivalues = {}
            self.monitor_lovalues = {}
            if self.pvControls is not None:
                for i, pv in enumerate(self.pvControls):
                    self.control_values[pv] = \
                            np.array([val[i] for val in self._scanData['pvControls_value']])
                
            if self.pvMonitors is not None:
                for i, pv in enumerate(self.pvMonitors):
                    self.monitor_hivalues[pv] = \
                            np.array([val[i] for val in self._scanData['pvMonitors_hiValue']])
                    self.monitor_lovalues[pv] = \
                            np.array([val[i] for val in self._scanData['pvMonitors_loValue']])

            self.pvAliases = {}
            for i, pv in enumerate(self.pvControls):
                alias = re.sub('-|:|\.| ','_', self.pvLabels[i])
                self.pvAliases[pv] = alias 
                setattr(self, alias, self.control_values[pv])

        ds.reload()

    @property
    def dataset(self):
        """
        xarray Dataset of ScanData
        """
        import xarray as xr
        nsteps = len(self.start_times)
        xscan = xr.Dataset({'step': range(nsteps)})
        xscan.coords['time'] = (('step'), self.start_datetimes) 
        xscan.coords['time_ns'] = (('step'), xscan.time.values.tolist())
        xscan['step_tstart'] = (('step'), self.start_datetimes) 
        xscan['step_tend'] = (('step'), self.end_datetimes) 
        xscan['step_istart'] = (('step'), self._scanData['ievent_start']) 
        xscan['step_iend'] = (('step'), self._scanData['ievent_end']) 
        xscan['step_events'] = (('step'), self.nevents) 
        for name, value in self.control_values.items():
            try:
                alias = self.pvAliases.get(name,  re.sub('-|:|\.| ','_', name))
                xscan[alias] = (('step'), value)
            except:
                print('Cannot add {:} values to ScanData.dataset'.format(name))

        return xscan

    def show_info(self, **kwargs):
        """Show scan information.
        """
        import numpy as np
        message = Message(quiet=True, **kwargs)
        
        attrs = { 
            'nsteps':      {'unit': '',     'desc': 'Number of steps'}, 
            'npvControls': {'unit': '',     'desc': 'Number of control PVs'},
            'npvMonitors': {'unit': '',     'desc': 'Number of monitor PVs'},
            }

        message('{:10}: Run {:}'.format(self._ds.data_source.exp, self._ds.data_source.run))
        message('-'*70)
        for attr, item in attrs.items():
            message('{:24} {:10} {:16}'.format(item.get('desc'), getattr(self, attr), attr))
       
        message('')
        message('{:24} {:40}'.format('Alias', 'PV'))
        message('-'*70)
        for name, alias in self.pvAliases.items():
            message('{:24} {:40}'.format(alias, name))
        message('')

        self._control_format = {}
        self._name_len = {}
        header1 = '{:4} {:6} {:>10}'.format('Step', 'Events', 'Time [s]')
        for i, name in enumerate(self.control_values):
            alias = self.pvAliases.get(name, name)
            name_len = len(alias)
            self._name_len[name] = name_len 
            name_format = ' {:>'+str(name_len)+'}'
            header1 += name_format.format(alias)
            vals = self.control_values[name]
            try:
                if self.nsteps > 1:
                    meanvals = np.mean(abs(vals[1:]-vals[:-1]))
                    if meanvals > 0:
                        sigdigit = int(np.floor(np.log10(meanvals)))
                    else:
                        sigdigit = -2
                else:
                    sigdigit = int(np.floor(np.log10(abs(vals))))
            except:
                sigdigit = 0

            if sigdigit < -5 or sigdigit > 5:
                self._control_format[name] = ' {:'+str(name_len)+'.3e}'
            elif sigdigit < 0:
                self._control_format[name] = ' {:'+str(name_len)+'.'+str(-sigdigit+1)+'f}'
            else:
                self._control_format[name] = ' {:'+str(name_len)+'}'

        message(header1)
        message('-'*(21+int(sum(self._name_len.values()))))
        for i, nevents in enumerate(self.nevents):
            a = '{:4} {:6} {:8.3f}'.format(i, nevents, self.step_times[i])
            for name, vals in self.control_values.items():
                a += self._control_format[name].format(vals[i])
            
            message(a)

        return message

    def __str__(self):
        return  'ScanData: '+str(self._ds.data_source)

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__,str(self))
        print('< '+repr_str+' >')
        self.show_info()
        return '< '+repr_str+' >'


class DataSource(object):
    """
    Python version of psana.DataSource with support for event and config
    data as well as PyDetector functions to access calibrated data.

    data_source string is build by DataSourceInfo object if not passed 
    directly as data_source string.
    
    Parameters 
    ----------
    data_source : str
        name of data source using psana convention (e.g., 'exp=xpptut15:run=54:smd')

    exp : str
        experiment name (e.g., exp='xpptut15')

    run : int
        run number (e.g., run=54)

    smd : bool , optional
        small data support -- default for experiments after Oct 2015 when this feature became standard 

    Attributes
    ----------
    data_source :  object
        DataSourceInfo object

    steps : object
        Steps object

    epicsData : object

    configData : object

    """

    _ds_funcs = ['end', 'env']
    _ds_attrs = ['empty']
    _env_attrs = ['calibDir', 'instrument', 'experiment','expNum']
#    _plugins = {}
#    _default_plugins = ['psplot.Psplot']
    _default_modules = {
            'path': '',
            'devName': {
                'Imp': 'impbox',
                'Acqiris': 'acqiris',
                'Tm6740': 'yag',
#                'Cspad': 'cspad',
#                'Opal1000': 'camera',
#                'Opal2000': 'camera',
#                'Opal4000': 'camera',
#                'Opal8000': 'camera',
##                'Epix': 'epix100',
##                'Cspad2x2': 'cspad2x2',
                },
             'srcname': {
#                'XrayTransportDiagnostic.0:Opal1000.0': 'xtcav_det',
#                'CxiDsu.0:Opal1000.0': 'timetool',     
                },
            }


    def __init__(self, data_source=None, default_modules=None, **kwargs):
        #self._device_sets = {'DataSource': {}}
        self._dataset = None
        self._exp_summary = None
        self._device_sets = {}
        self._current_data = {}
        path = os.path.dirname(__file__)
        if not path:
            path = '.'

        if default_modules is not None:
            self._default_modules = default_modules

        self._default_modules.update({'path': path})
        #self._default_modules.update({'path': os.path.join(path,'detectors')})
        loaded = self.load_run(data_source=data_source, **kwargs)
        if not loaded:
            return None

        if self.data_source.smd and not self.data_source.live:
            self._load_smd_config()

    def _load_smd_config(self, quiet=False):
        """Load configData of first calib cycle by going to first step.
           Reload so that steps can be used as an iterator.
        """
        if self.data_source.smd:
            try:
                step = next(self.steps)
                self.reload()
            except:
                traceback.print_exc()
                print('Could not load first step in smd data.' )


    def _get_exp_summary(self, reload=False, build_html=None, **kwargs):
        """
        Load experiment summary
        """
        if reload and build_html is None:
            build_html = True
        self._exp_summary = get_exp_summary(self.data_source.exp, reload=reload, build_html=build_html, **kwargs)
        return self._exp_summary

    @property
    def exp_summary(self):
        """
        Experiment summary.  Reload if not up to date with current data_source run.
        """
        if self._exp_summary is None:
            es = self._get_exp_summary()
        else:
            es = self._exp_summary

        if self.data_source.run not in es.xruns.run.values:
            es = self._get_exp_summary(reload=True)

        return es

    @property
    def scan_pvs(self):
        """
        Dict of pv alias and number of times set during run.
        For more detailed information on scan see for example xarray Dataset retruned from:
            xpvscan = ds.exp_summary.get_scan_data(run_number)
        """
        df = self.exp_summary.dfscan.T[self.data_source.run]
        df = df[df>0].to_dict()
        if df is None:
            df = []
        return df

    @property
    def moved_pvs(self):
        """
        Dict of pv alias and values that were set prior to run..
        """
        df = self.exp_summary.dfset.T.get(self.data_source.run)
        if df is not None:
            return df.dropna().to_dict()
        else:
            return None

    @property
    def dataset(self):
        """
        xarray Dataset of small data from beam_stats analysis,
        standard epics info for 'FEE' and instrument from _transmission_pvs,
        scan data automatically determined from exp_summary and 
        scan data from daq config of calib cycles.
        """
        if self._dataset is not None:
            return self._dataset
        else:
            return self.get_dataset()

    def get_dataset(self, pvdict={}, fields=None, 
            meta_attrs = {'units': 'EGU', 'PREC': 'PREC', 'pv': 'name'},
            load_default=True, load_scan=None, load_smd=None, load_psocake=None, 
            tstart=None, tend=None, quiet=True, **kwargs):
        """
        Get dataset from epics PVs

        Parameters
        ----------
        load_default : bool
            Load default PVs including transmission
        load_scan : bool
            Load PVs that are scanned -- uses exp_summary
        load_smd : bool
            Load beam_stats small data file used in off-by-one analysis
        load_psocake : bool
            Load psocake peak index 
        """
        import time
        import numpy as np
        import pandas as pd
        import xarray as xr
        from . import xarray_utils
        ds = self
        experiment = ds.data_source.exp
        run = ds.data_source.run
        configData = ds.configData
        xsmd = None
        if load_smd is not False:
            try:
                from . import beam_stats
                xsmd = beam_stats.load_small_xarray(ds) 
                if 'time_ns' not in xsmd.coords:
                    xsmd.coords['time_ns'] = (('time'), np.int64(xsmd.sec*1e9+xsmd.nsec))
            except:
                print('Cannot add smd data')

        if not tstart:
            tstart = configData._tstart
        if not tend:
            tend = configData._tend
        if not fields:
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
        time0 = time.time()
        time_last = time0
       
        if load_default:
            trans_pvs = [a for a in _transmission_pvs.get('FEE', {}) if a.endswith('_trans')]
            trans_pvs += [a for a in _transmission_pvs.get(ds.instrument, {}) if a.endswith('_trans')]
            trans3_pvs = [a for a in _transmission_pvs.get('FEE', {}) if a.endswith('_trans3')]
            trans3_pvs += [a for a in _transmission_pvs.get(ds.instrument, {}) if a.endswith('_trans3')]
            pvdict.update(**_transmission_pvs.get('FEE',{}))
            pvdict.update(**_transmission_pvs.get(ds.instrument,{}))
        else:
            trans_pvs = []
            trans3_pvs = []

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
                                        vt = np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns')
                                        vals[vt] = val
                               
                                data_fields[alias][attr] = xr.DataArray(vals.values(), 
                                                                coords=[vals.keys()], dims=['time'], 
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
                        if not quiet:
                            print((alias, 'string'))
                        vals = np.array(vals, dtype=str)
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
                    if not quiet:
                        print(('Error loadinig', alias))

                if not quiet:
                    try:
                        print('{:8.3f} {:28} {:8} {:10.3f} {:4} {:20} {:}'.format(time_next-time_last, \
                                        gias, len(vals), np.array(vals).mean(), units, doc, pv))
                    except:
                        print('{:8.3f} {:28} {:8} {:>10} {:4} {:20} {:}'.format(time_next-time_last, \
                                        alias, len(vals), vals[0], units, doc, pv))
            
            except:
                traceback.print_exc()
                if not quiet:
                    print(('Error loading', alias))

        xdata = xr.merge(data_arrays.values())
        if trans_pvs:
            da = xdata.reset_coords()[trans_pvs].to_array() 
            xdata['trans'] = (('time'), da.prod(dim='variable'))
            xdata['trans'].attrs['doc'] = 'Total transmission: '+'*'.join(trans_pvs)

        if load_psocake:
            from . import xarray_utils
            xpsocake = xarray_utils.open_cxi_psocake(experiment, run, 
                    load_smd=True, load_moved_pvs=False, **kwargs)
            if xpsocake is not None:
                for attr, item in xpsocake.data_vars.items():
                    if len(item.dims) > 1:
                        del xpsocake[attr]
                xsmd = xpsocake
        
        if xsmd is not None:
            xdata = xarray_utils.merge_fill(xsmd, xdata)
        
        xtimes = xr.Dataset({'time': ds._idx_datetime64})
        xdata = xarray_utils.merge_fill(xtimes, xdata) 
        
        if load_scan is not False:
            try:
                xscan = ds.exp_summary.get_scan_data(run)
                xdata = xarray_utils.merge_fill(xdata, xscan)
            except:
                print('Cannot add exp_summary scan data')
            try:
                xscan = ds.configData.ScanData.dataset
                if load_scan or xscan.step.size > 1:
                    xdata = xarray_utils.merge_fill(xdata, xscan)
            except:
                print('Cannot add exp_summary scan data')
      
        self._dataset = xarray_utils.resort(xdata)

        return self._dataset

    def __getattr__(self, attr):
        if attr in self._attrs:
            attr_dict = {key: pdict for key,pdict in self._pv_dict.items()
                         if pdict['components'][0] == attr}
            return PvData(attr_dict, self._ds, level=1)
        
        if attr in dir(self._ds.env().epicsStore()):
            return getattr(self._ds.env().epicsStore(),attr)

    def __dir__(self):
        all_attrs = set(self._attrs +
                        dir(self._ds.env().epicsStore()) +
                        self.__dict__.keys() + dir(EpicsData))
        return list(sorted(all_attrs))


    def load_run(self, data_source=None, reload=False, 
            try_idx=None, wait=True, timeout=20., 
            update_url=None, quiet=True, **kwargs):
        """
        Load a run with psana.

        Parameters
        ----------
        wait : bool
            Wait until files are available
        """
        self._evtData = None
        self._current_evt = None
        self._current_step = None
        self._current_run = None
        self._evt_keys = {}
        self._evt_modules = {}
        self._init_dets = []
        if not reload:
            self.data_source = DataSourceInfo(data_source=data_source, **kwargs)

        # do not reload shared memory
        if not (self.data_source.monshmserver and self._ds):
            # do not wait for monshmserver
            if wait and not self.data_source.monshmserver:
                import time
                from . import psutils
                time0 = time.time()
                while not psutils.run_available(self.data_source.exp, self.data_source.run) and 'tut' not in self.data_source.exp: 
                    time.sleep(5.)
                    print('Waiting for {:} to become available...'.format(self.data_source))
                    if update_url:
                        from requests import post
                        batch_counters = ['Waiting {:} sec for data'.format(time.time()-time0),'yellow']
                        post(update_url, json={'counters' : batch_counters})

                    if timeout and timeout > time.time()-time0:
                        print('Timeout loading {:}'.format(self.data_source))
                        return None

            if True:
                if wait and not self.data_source.monshmserver:
                    time0 = time.time()
                    if not quiet:
                        print('Try Loading {:}'.format(self.data_source))
                    while True:
                        try:
                            self._ds = psana.DataSource(str(self.data_source))
                            break
                        except:
                            time.sleep(5.)
                            print('...Waiting for {:} to become available...'.format(self.data_source))
                            if timeout and timeout > time.time()-time0:
                                print('...Timeout loading {:}'.format(self.data_source))
                                return None
                else:
                    if not quiet:
                        print('Loading {:}'.format(self.data_source))
                    self._ds = psana.DataSource(str(self.data_source))
                
                _key_info, _modules = get_keys(self._ds.env().configStore())
                if 'Partition' in _modules:
                    try_idx = False

                else:
                    #if not self.data_source.smd:
                    if False:
                        print('Exp {:}, run {:} is has no Partition data.'.format( \
                            self.data_source.exp, self.data_source.run))
                        print('PyDataSource requires Partition data.')
                        print('Returning psana.DataSource({:})'.format(str(self.data_source)))
                        return self._ds
                    elif try_idx is None:
                        try_idx = True
                        print('Exp {:}, run {:} smd data has no Partition data -- loading idx data instead.'.format( \
                            self.data_source.exp, self.data_source.run))
            else:
                if try_idx is None:
                    try_idx = True
                    print('Exp {:}, run {:} smd data file not available -- loading idx data instead.'.format( \
                                self.data_source.exp, self.data_source.run))

        if try_idx:
            if self.data_source.smd:
                try:
                    print('Use smldata executable to convert idx data to smd data.')
                    data_source_smd = self.data_source
                    data_source_idx = str(data_source_smd).replace('smd','idx')
                    self.data_source = DataSourceInfo(data_source=data_source_idx)
                    if not quiet:
                        print('Loading idx {:}'.format(self.data_source))
                    self._ds = psana.DataSource(str(self.data_source))
                except:
                    print('Failed to load either smd or idx data for exp {:}, run {:}'.format( \
                            self.data_source.exp, self.data_source.run))
                    print('Data can be restored from experiment data portal:  https://pswww.slac.stanford.edu')
                    return False

        self.epicsData = EpicsData(self._ds) 

        self._evt_time_last = (0,0)
        self._ievent = -1
        self._istep = -1
        self._irun = -1
        self._jump_last = False
        if self.data_source.idx:
            self.runs = Runs(self, **kwargs)
            self.events = self.runs.next().events
            if not reload:
                self._idx_nsteps = self.runs.current.nsteps
                self._idx_times = self.runs.current.times
                self.nevents = len(self._idx_times)
                self._idx_times_tuple = [(a.seconds(), a.nanoseconds(), a.fiducial()) \
                                        for a in self._idx_times]
                self._idx_datetime64 = [np.datetime64(int(sec*1e9+nsec), 'ns') \
                                        for sec,nsec,fid in self._idx_times_tuple]

            if 'BldInfo(EBeam)' not in self.configData._sources:
                try:
                    self.add_detector('BldInfo(EBeam)', alias='EBeam')
                except:
                    pass
            if 'BldInfo(FEEGasDetEnergy)' not in self.configData._sources:
                try:
                    self.add_detector('BldInfo(FEEGasDetEnergy)', alias='FEEGasDetEnergy')
                except:
                    pass

        elif self.data_source.smd and not self.data_source.live:
            self.steps = Steps(self, **kwargs)
            # SmdEvents automatically goes to next step if no events in current step.
            self.events = SmdEvents(self)
            if not reload:
                if wait:
                    import time
                    from . import psutils
                    time0 = time.time()
                    while not psutils.run_available(self.data_source.exp, self.data_source.run, idx=True) and 'tut' not in self.data_source.exp:
                        time.sleep(5.)
                        print('Waiting for {:} run{:} idx to become available...'.format(self.data_source.exp, self.data_source.run))
                        if timeout and timeout > time.time()-time0:
                            print('Timeout loading {:}'.format(self.data_source))
                            return None

                self._scanData = None
                if 'live' in str(self.data_source):
                    self.nevents = None
                else:
                    data_source_idx = str(self.data_source).replace('smd','idx')
                    if wait:
                        time0 = time.time()
                        while True:
                            if not quiet:
                                print('Additional wait loading idx {:}'.format(self.data_source))
                            try:
                                self._idx_ds = psana.DataSource(data_source_idx)
                                break
                            except:
                                time.sleep(5.)
                                print('Waiting for {:} to become available...'.format(data_source_idx))
                                if timeout and timeout > time.time()-time0:
                                    print('Timeout loading {:}'.format(self.data_source))
                                return None
                    else:
                        if not quiet:
                            print('Additional loading idx {:}'.format(self.data_source))
                        self._idx_ds = psana.DataSource(data_source_idx)
                    
                    self._idx_run = next(self._idx_ds.runs())
                    self._idx_nsteps = self._idx_run.nsteps()
                    self._idx_times = self._idx_run.times()
                    self.nevents = len(self._idx_times)
                    self._idx_times_tuple = [(a.seconds(), a.nanoseconds(), a.fiducial()) \
                                            for a in self._idx_times]
                    self._idx_datetime64 = [np.datetime64(int(sec*1e9+nsec), 'ns') \
                                            for sec,nsec,fid in self._idx_times_tuple]

        else:
            # For live data or data_source without idx or smd
            self.events = Events(self)
            self.nevents = None

        if not reload:
            self._init_evr()

        return str(self.data_source)

    def reload(self, reset_stats=True):
        """Reload the current run.
        """
        self.load_run(reload=True)
        if reset_stats:
            self.reset_stats()

    def reset_stats(self, attrs=None):
        """
        Reset Welford stats objects.
        
        Parameters
        ----------
        attrs:  Optional list of stats objects to reset (e.g., [CsPad.calib_stats])
        """
        if attrs and not isinstance(attrs, list):
            attrs = [attrs]

        for alias, det_config in self._device_sets.items():
            if 'stats' in det_config:
                for attr, item in det_config['stats'].items():
                    if not attrs or '.'.join([alias,attr]) in attrs:
                        item['funcs'] = {}
    
    def _add_default_stats(self, attrs=[], **kwargs):
        """
        Add default statistics for each Image and Waveform Detectors
        if not stats already created
        """
        if not attrs:
            attrs = self._detectors.keys()
        
        attrs = list(set(attrs))

        for attr in attrs:
            det = self._detectors[attr]
            try:
                if not det._det_config['stats']:
                    if det._det_class == WaveformData:
                        OKadd = det.add.stats('waveform', **kwargs)
                    elif det._det_class == ImageData:
                        try:
                            OKadd = det.add.stats('corr', **kwargs)
                        except:
                            OKadd = False
                        print((attr, OKadd))
                        if not OKadd:
                            next(det)
                            print(det)
                            OKadd = det.add.stats('corr', **kwargs)
                    else:
                        OKadd = False

                    #elif 'data' in det._attrs:
                    #    # currently just zyla
                    #    det.add.stats('data', **kwargs)
                    if OKadd:
                        print(('Added default stats for ', attr, str(det)))
                        print(det._det_config['stats'])

            except:
                traceback.print_exc()
                print(('Cannot add stats for', attr))

    @property
    def stats(self):
        """
        xarray Dataset of all Welford stats.
        """
        return self._get_stats()

    def _get_stats(self, attrs=None, aliases={}):
        """
        xarray Dataset of Welford stats.
        
        Parameters
        ----------
        attrs:  Optional list of stats objects to save (e.g., [CsPad.calib_stats])
        """
        import xarray as xr
        import os
        if attrs and not isinstance(attrs, list):
            attrs = [attrs]

        datasets = []
        for alias, det_config in self._device_sets.items():
            if 'stats' in det_config:
                for attr, item in det_config['stats'].items():
                    if not attrs or '.'.join([alias,attr]) in attrs:
                        try:
                            datasets.append(self._detectors[alias]._get_stats(attr, alias=aliases.get(alias)))
                        except:
                            traceback.print_exc()
                            print(('Cannot add stats for', alias, attr))

        if datasets == []:
            return None

        try:
            x = xr.merge(datasets)
        except:
            print('merge failed')
            return datasets
            
        try:
            xsteps = list(x.coords.get('steps').values)
            if xsteps:
                for pv, vals in self.configData.ScanData.control_values.items():
                    alias = self.configData.ScanData.pvAliases[pv]
                    x.coords[alias+'_steps'] =  (('steps'), vals[xsteps]) 
        except:
            print('could not add steps')
            
        try:
            x.attrs['data_source'] = self.data_source.data_source
            x.attrs['instrument'] = self.instrument
            x.attrs['run'] = self.data_source.run
            x.attrs['experiment'] = self.experiment
            x.attrs['expNum'] = self.expNum
        except:
            print('could not add attrs to stats Dataset')

        return x

    def save_stats(self, attrs=None, file_name=None, path=None, 
            h5folder='scratch', subfolder='nc',
            aliases={},
            engine='h5netcdf', **kwargs):
        """
        Save Welford stats as xarray compatible hdf5 files (using h5netcdf engine).
        
        Parameters
        ----------
        attrs:  Optional list of stats objects to save (e.g., [CsPad.calib_stats])
        """
        import xarray as xr
        import os
        data_sets = self._get_stats(attrs=attrs, aliases=aliases)
        try:
            if data_sets is None:
                print('No stats to save')
                return
            elif isinstance(data_sets, list):
                print('Cannot save stats')
                return
        except:
            traceback.print_exc()
            return

        if not file_name:
            run = self.data_source.run
            if not path:
                instrument = self.data_source.instrument
                exp = self.data_source.exp
                path = '/reg/d/psdm/{:}/{:}/{:}/{:}/Run{:04}/'.format(instrument, 
                        exp,h5folder,subfolder,run)
                if not os.path.isdir(path):
                    os.mkdir(path)
        
            file_base = 'run{:04}'.format(int(run))
            file_name = os.path.join(path,'{:}_{:}.nc'.format(file_base,'stats'))
   
        self.stats.to_netcdf(file_name, engine=engine)

    def _load_ConfigData(self):
        self._ConfigData = ConfigData(self)

    @property
    def xarray_kwargs(self):
        """
        Dictionary of keywords passed to to_xarray.  
        This information is saved/loaded with save_config/load_config methods.

        Parameters
        ----------
        max_size : uint
            Maximum array size of data objects to build into xarray.
        pvs: list
            List of pvs to be loaded vs time
        epics_attrs: list
            List of epics pvs to be saved as run attributes based on inital value 
            of first event.
        code_flags : dict
            Dictionary of event code flags. 
            Default = {'XrayOff': [162], 'XrayOn': [-162]}

        """
        if 'DataSource' not in self._device_sets:
            self._device_sets.update({'DataSource': {}})
        return self._device_sets.get('DataSource') 

    def submit_summary(self, **kwargs):
        """
        Batch submission of run summary with auto report.
        """
        import subprocess
        import os
        if kwargs:
            self.xarray_kwargs.update(**kwargs)
            self.save_config()

        path = os.path.dirname(__file__)
        submit_file = os.path.join(path, '../../bin/submit_summary')
        bsubproc = '{:} {:} {:}'.format(submit_file, self.data_source.exp, self.data_source.run)
        try:
            print('Submitting: {:}'.format(bsubproc))
            subproc = subprocess.Popen(bsubproc, stdout=subprocess.PIPE, shell=True)
        except:
            print('Submit from ipython not yet available')
            print('Setup conda environment the same was as for ipython then:')
            print('-> kinit')
            print('-> {:}'.fomrat(bsubproc))

    def to_xarray(self, build_html=False, **kwargs):
        """
        Build xarray object from PyDataSource.DataSource object.
        Same as to_hdf5 -- backwards compatability
        """
        from . import h5write
        xarray_kwargs = self.xarray_kwargs.copy()
        xarray_kwargs.update(**kwargs)
        x = h5write.to_hdf5(self, **xarray_kwargs)
        if build_html:
            try:
                from . import build_html
                self.html = build_html.Build_html(x, auto=True) 
            except:
                traceback.print_exc()
                print(('Could not build html run summary for', str(self)))

        return x

    def to_hdf5(self, build_html=False, **kwargs):
        """
        Write directly to hdf5 in xarray comatable netcdf4 format.
        
        Parameters
        ----------
        build_html : bool
            Build run summary if True
        max_size : uint
            Maximum array size of data objects to build into xarray.
        nchunks: int
            number of chunks when using mpi
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
        default_stats : bool
            If true include stats for waveforms and calib image data.
        min_all_save : int
            Min number of events where all 'calib' data saved.  Sets max_size = 1e10
        default_stats : bool
            If true automatically add default stats for all detectors that do not already
            have stats added
        auto_update : bool
            If true automatically update xarray info for all detectors
        auto_pvs : bool
            If true automatically add pvs that were moved during run.

        """
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)
        size = comm.Get_size()
        
        xarray_kwargs = self.xarray_kwargs.copy()
        xarray_kwargs.update(**kwargs)
        
        from . import h5write
        if size > 1:
            x = h5write.to_hdf5_mpi(self, build_html=build_html, **xarray_kwargs)
        else:
            x = h5write.to_hdf5(self, **xarray_kwargs)
            if build_html:
                try:
                    from . import build_html
                    self.html = build_html.Build_html(x, auto=True) 
                except:
                    traceback.print_exc()
                    print(('Could not build html run summary for', str(self)))

        return x

    @property
    def configData(self):
        """
        Configuration Data from ds.env().configStore().
        For effieciency only loaded at beginning of run or step unless
        working with shared memory.
        
        See Also
        --------
        ConfigData : class
            Configuation Data Access  
        
          
        """
        if self.data_source.monshmserver:
            self._load_ConfigData()
        
        return self._ConfigData

    @property
    def scanData(self):
        """
        ScanData
        
        See Also
        --------
        ScanData : class
            Scan Data Access  
        """
        return self.configData.ScanData

    @property
    def sources(self):
        """
        Detector sources
        
        See Also
        --------
        Sources : class
            Event data sources

        """
        return self.configData.Sources

    def _make_smd_file(self):
        """
        Execute script to make small data file for old files.        
 
            Usage:  ./smldata  [-f <xtc filename>] [-i <input index>] [-o <index filename>] [-s <size threshold>] [-h]
              Options:
                -h                     Show usage.
                -f <xtc filename>      Set input xtc filename
                -i <index filename>    Set input index filename
                   Note 1: -i option is used for testing. The program will parse
                     the index file to see if it is valid, and output to another
                     index file if -o option is specified
                   Note 2: -i will overwrite -f option
                -o <index filename>    Set output index filename
                -s <size threshold>    Set L1Accept xtc size threshold
        """
        import subprocess
        import glob
        import os

        exp_dir = os.path.join('/reg/d/psdm/',self.data_source.instrument,self.data_source.exp)
        xtc_dir = os.path.join(exp_dir,'xtc')
        smd_dir = os.path.join(xtc_dir,'smalldata')
        script_path = '/reg/common/package/pdsdata/8.7.3/x86_64-rhel7-opt/bin/'
        xtc_files = glob.glob(xtc_dir+'/*-r{:04}-s*.xtc'.format(self.data_source.run))
        if not os.path.isdir(smd_dir):
            os.mkdir(smd_dir)
        
        print('WARNING:  Need to have priviledge to do this...')
        print('  ask for help from CDS')
        print('Execute the following commands from psana:')
        print('------------------------------------------')
        for file_in in xtc_files:
            file_out = file_in.replace('/xtc/','/xtc/smalldata/').rstrip('xtc')+'smd.xtc'
            if not glob.glob(file_out):
                cmd = '{:}smldata -f {:} -o {:}'.format(script_path, file_in, file_out)
                print(cmd)
                #subproc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            else:
                print('# Small data file already exists: {:}'.format(file_out))


    def _init_detectors(self):
        """Initialize psana.Detector classes based on psana env information.
        """
        self._detectors = {}
        self._load_ConfigData()
        self._aliases = self.configData._aliases

        for srcstr, item in self.configData._sources.items():
            alias = item.get('alias')
            self._add_dets(**{alias: srcstr})

    def _init_evr(self):
        """
        Initialize evr. 
        Requires first event to know which evr when more than one present.
        """
        try:
            while 'EvrData' not in self._evt_modules and self._ievent < 400:
                evt = next(self.events)
        except:
            print('No EvrData in first {:} events'.format(self._ievent))
            return

        if 'EvrData' in self._evt_modules:
            evr_typ, evr_src, evr_key = self._evt_modules['EvrData'].values()[0][0]
            srcstr = str(evr_src) 
            srcname = srcstr.split('(')[1].split(')')[0]
            devName = srcname.split(':')[1].split('.')[0]
            ievr = srcname.split(':')[1].split('.')[1]
            self._evr = psana.Detector(srcstr)
            #self._evr = psana.Detector('evr{:}'.format(ievr))
            self._evr_alias = 'evr{:}'.format(ievr)
        else:
            # if no EvrData in event then try using first Evr
            ievr = 0
            while ievr is not False:
                try:
                    self._evr = psana.Detector('evr{:}'.format(ievr))
                    self._evr_alias = 'evr{:}'.format(ievr)
                    ievr = False
                except:
                    ievr += 1
                    if ievr > 10:
                        print(self.configData.keys())
                        traceback.print_exc('Erro. -- no evr present')
        
        self.reload()

    def add_detector(self, srcstr=None, alias=None, module=None, path=None, 
                     #pvs=None, desc=None, parameters={}, 
                     desc=None,
                     quiet=False,
                     **kwargs):
        """
        Add a detector 
        
        Parameters
        ----------

        srcstr : str
            Source string name or alias of detector

        alias : str
            Detector alias -- default uses srcstr

        module : str
            Name of python module that contains a user defined PyDataSource.Detector class
            with the same name as the module.  e.g., 'acqiris' loads 'acqiris.py' file 
            which must have class Acqiris(PyDatasource.Detector)

        path : str
            Name of path for python module (default is the path of PyDataSource)

        desc : str
            Description 

        """
        initialized = False
        if not alias:
            if not srcstr:
                print('Source string name or alias must be supplied')
            elif srcstr in self._aliases:
                alias = srcstr
                srcstr = self._aliases[alias]
            else:
                alias = re.sub('-|:|\.| ','_', srcstr)
        else:
            if not srcstr:
                srcstr = self._aliases.get(alias)
                if not srcstr:
                    raise Exception('{:} not a valid Detector'.format(alias))
            else:
                # update aliases with srcstr -- required for old data
                # but could overwrite aliases for new data
                if alias not in self._aliases:
                    self._aliases[alias] = srcstr
                    try:
                        psana_src = psana.Source(srcstr)
                    except:
                        psana_src = None
                    self.configData._sources[srcstr] = {'alias': alias,
                                                        'group': 0,
                                                        'src': psana_src}
                elif srcstr != self._aliases.get(alias):
                    raise Exception('{:} alias already taken'.format(alias))

        det = alias
        
        if alias not in self._device_sets:
            self._device_sets[alias] = {
                    'alias': alias, 
                    'module': {},
                    'opts': {}, 
                    'parameter': {}, 
                    'property': {},
                    'psplot': {},
                    'roi': {},
                    'count': {},
                    'peak': {},
                    'histogram': {},
                    'projection': {},
                    'stats': {},
                    'xarray': {},
                    }

            if not desc:
                desc = alias
        
        if True:
            srcname = srcstr.split('(')[1].split(')')[0]
            try:
                devName = srcname.split(':')[1].split('.')[0]
            except:
                devName = None

            det_dict = self._device_sets[alias]
            det_dict.update({'desc': desc, 'srcname': srcname, 'srcstr': srcstr, 'devName': devName})
            
            if module:
                module = module.split('.')[0]

            # First check for device configuration
            
            if det_dict.get('module'):
                module_name = det_dict['module'].get('name',None)
                module_path = det_dict['module'].get('path','')
#                module_name = det_dict['module'].get('name',None)
#                if 'path' in det_dict['module']:
#                    module_path = det_dict['module'].get('path','')
#                else:
#                    module_path = ''
            else:
                module_name = None 
                module_path = ''

            # Then use module and path keywords if applicable
            if module:
                if module_name:
                    if not quiet:
                        print('Changing {alias} detector module from {module_name} to {module}'.format( \
                               alias=alias,module=module,module_name=module_name))
                else:
                    det_dict['module'] = {}
                
                module_name = module
                det_dict['module']['name'] = module

            # Use defaults if not set by keyword or in device config
            if not module_name: 
                if srcname in self._default_modules.get('srcname', {}):
                    module_name = self._default_modules['srcname'][srcname]
                    module_path = self._default_modules.get('path','')
                elif devName and devName in self._default_modules.get('devName', {}):
                    module_name = self._default_modules['devName'][devName]
                    module_path = self._default_modules.get('path','')

            if module_name:
                is_default_class = False
            else:
                is_default_class = True

            if not is_default_class:
                if path:
                    module_path = path
        
                if module_path:
                    if not quiet:
                        print('Using the path {module_path} for {module_name}'.format( \
                               module_path=module_path, module_name=module_name))
                else:
                    module_path = self._default_modules['path']
                    
                det_dict['module']['path'] = module_path
                det_dict['module']['kwargs'] = kwargs

                new_class = get_module(module_name, module_path, reload=True)
#                import_module(module_name, module_path)
#                try:
#                    new_class =  getattr(globals()[module_name],module_name)
#                except:
#                    new_class =  getattr(globals()[module_name],module_name.capitalize())
                
                det_dict['module']['dict'] = [attr for attr in new_class.__dict__ \
                                              if not attr.startswith('_')]

                if not quiet:
                    print('Loading {alias} as {new_class} from {module_path}'.format( \
                           alias=alias,new_class=new_class,module_path=module_path))
                
                nomodule = False
                self._detectors[alias] = new_class(self, alias, **kwargs)
                initialized = True

            if is_default_class:
                if not quiet:
                    print('Loading {alias} as standard Detector class'.format(alias=alias))
                
                self._detectors[alias] = Detector(self, alias)
                initialized = True

        if initialized:
            self._init_dets.append(alias)

#    def add_plugin(self, cls, **kwargs):
#        self._plugins.update({cls.__name__: cls})

    def _add_dets(self, **kwargs):
        for alias, srcstr in kwargs.items():
            try:
                self.add_detector(srcstr, alias=alias, quiet=True)
            except Exception as err:
                print('Cannot add {:}:  {:}'.format(alias, srcstr) )
                traceback.print_exc()

    def _get_config_file(self, run=None, exp=None, instrument=None, user=None, 
            path=None):
        """
        Get the full config file name.  
        """
        if user:
            if user is True:
                path = self.data_source.user_dir
            else:
                path = self.data_source._get_user_dir(user)

        elif not path:
            exp = self.experiment
            base_path = '/reg/d/psdm/{:}/{:}/results'.format(self.instrument, exp)
            if not os.path.isdir(base_path):
                base_path = '/reg/d/psdm/{:}/{:}/res'.format(self.instrument, exp)
            path = base_path+'/summary_config'

        if not os.path.isdir(path):
            os.mkdir(path)

        if not run:
            for run in range(int(self.data_source.run), 0, -1):
                file_name = '{:}/run{:04}.config'.format(path, int(run))
                if os.path.isfile(file_name):
                    return file_name
            
        else:
            file_name = '{:}/run{:04}.config'.format(path, int(run))
            if os.path.isfile(file_name):
                return file_name
            
        err_message = 'No valid config file found for {:} Run {:}'.format(exp, run)
        print(err_message)
        return None
            

    def save_config(self, file_name=None, path=None, **kwargs):
        """
        Save DataSource configuration.
        
        Parameters
        ----------
        file_name : str
            Name of file
        path : str
            Path of file
        """
        import pandas as pd
        if not file_name:
            #file_name = self._get_config_file(run=self.data_source.run, path=path)
            if not path:
                exp = self.experiment
                base_path = '/reg/d/psdm/{:}/{:}/results'.format(self.instrument, exp)
                if not os.path.isdir(base_path):
                    base_path = '/reg/d/psdm/{:}/{:}/res'.format(self.instrument, exp)
                path = base_path+'/summary_config'
            
            if not os.path.isdir(path):
                os.mkdir(path)

            run = self.data_source.run
            file_name = '{:}/run{:04}.config'.format(path, int(run))

        if file_name:
            pd.DataFrame.from_dict(self._device_sets).to_json(file_name)
        else:
            print('No valid file_name to save config')

    def load_config(self, run=None, exp=None, user=None, file_name=None, path=None, quiet=True, **kwargs):
        """
        Load DataSource configuration.  Overwrites any existing config and reloads DataSource.
        
        Parameters
        ----------
        run : int, optional
            Run number 
        exp : str, optional
            Experiment name
        file_name : str, optional
            Name of file
        path : str, optional
            Path of file
        user : str, optional
            Name of user to load config file
        """
        import pandas as pd
        if not file_name:
            file_name = self._get_config_file(run=run, path=path, user=user)

#        attrs = ['parameter', 'property', 'psplot', 'peak', 
#                 'roi', 'projection', 'xarray']
    
        if not file_name:
            print('No config file to load')
        elif not os.path.isfile(file_name): 
            print('Config file not present: {:}'.format(file_name))
        else:
            try:
                config = pd.read_json(file_name).to_dict()
            except:
                print('Config file not present: {:}'.format(file_name))
                return

            for alias, item in config.items():
                if alias in self._device_sets:
                    if 'module' in item:
                        module_dict = item.pop('module')

                        if module_dict:
                            # Make sure start with event when loading
                            try:
                                while self._current_evt is None or alias not in self.events.current._attrs:
                                    next(self.events)
                            except:
                                traceback.print_exc()
                                print('No {:} detector object in data stream.'.format(alias))

                            #print alias, module_dict
                            kwargs = module_dict.get('kwargs', {})
                            kwargs.update({'alias': alias, 
                                           'srcstr': module_dict.get('srcstr'),
                                           'module': module_dict.get('name'),
                                           'path': module_dict.get('path'),
                                           'desc': module_dict.get('desc')})
                            print(('add_detector', kwargs))
                            try:
                                self.add_detector(**kwargs)
                                self.reload()
                            except:
                                traceback.print_exc()
                                print('Could not add detector module.')
                                print(kwargs)

                    det_config = self._device_sets[alias]
                    for attr, config_dict in item.items():
                        if False and isinstance(config_dict, dict):
                            # Does not work yet to only try updating 
                            for a, val in config_dict.items():
                                if isinstance(val, dict):
                                    if a not in det_config[attr]:
                                        det_config[attr][a] = {}
                                    for b, bval in val.items():
                                        det_config[attr][a][b] = bval
                                else:
                                    if a not in det_config[attr]:
                                        det_config[attr][a] = val
                                    else:
                                        print((alias, attr, 'No overwrite', a, val))
                        elif attr != 'stats':
                            det_config[attr] = config_dict

                    # Need to be careful with stats objects to make sure added with dims correctly
                    stats_config = item.get('stats')
                    if stats_config:
                        evt = self.events.current
                        while alias not in evt._attrs:
                            evt = next(self.events)
                        detector = getattr(evt, alias)
                        for name, stat_item in stats_config.items():
                            attr = stat_item['attr']
                            print(('adding stats for', attr, name, ))
                            detector.add.stats(attr, name=name)
                        self.reload()
        self.reload()

    def _load_small_xarray(self, path=None, filename=None, refresh=False, 
            engine='h5netcdf', **kwargs):
        """Load small xarray Dataset. 
        
        """
        from . import beam_stats
        self.x = beam_stats.load_small_xarray(self, path=path, filename=filename, 
                refresh=refresh, engine=engine, **kwargs) 
        return self.x

    def show_info(self, **kwargs):
        """
        Show DataSource information.
        """
        return self.configData.show_info(**kwargs)
    
    def __str__(self):
        return  str(self.data_source)

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__,str(self))
        if self.nevents:
            repr_str += ' {:} events'.format(self.nevents)
        print('< '+repr_str+' >')
        self.show_info()
        return '< '+repr_str+' >'

    def __getattr__(self, attr):
        if attr in self._ds_attrs:
            return getattr(self._ds, attr)()
        if attr in self._ds_funcs:
            return getattr(self._ds, attr)
        if attr in self._env_attrs:
            return getattr(self._ds.env(), attr)()
        
    def __dir__(self):
        all_attrs =  set(self._ds_attrs + 
                         self._ds_funcs + 
                         self._env_attrs +
                         self.__dict__.keys() + dir(DataSource))
        
        return list(sorted(all_attrs))


class Runs(object):
    """
    psana DataSource Run iterator from ds.runs().
    """
    def __init__(self, ds, **kwargs):
        self._ds_runs = []
        self._kwargs = kwargs
        self._ds = ds

    def __iter__(self):
        return self

    @property
    def current(self):
        return self._ds._current_run

    def next(self, **kwargs):
        """
        Returns
        -------
        iterator
            Run iterator
        
        See Also
        --------
        Run : class
            
        """
        self._ds._ds_run = next(self._ds._ds.runs())
        self._ds_runs.append(self._ds._ds_run)
        self._ds._irun +=1
        self._ds._istep = -1
        self._ds._ievent = -1
        self._ds._init_detectors()
        self._ds._current_run = Run(self._ds)

        return self._ds._current_run


class Run(object):
    """
    Python psana.Run class from psana.DataSource.runs().next().
    """
    _run_attrs = ['nsteps', 'times']
    _run_funcs = ['end', 'env']

    def __init__(self, ds, **kwargs):
        self._ds = ds

    @property
    def events(self):
        """RunEvents Iterator.

        Returns
        -------
        iterator
            Events in current run of DataSource.
        
        See Also
        --------
        RunEvents : class
        """
        return RunEvents(self._ds)

#    @property
#    def steps(self):
#        return RunSteps(self._ds)

    def __getattr__(self, attr):
        if attr in self._run_attrs:
            return getattr(self._ds._ds_run, attr)()
        if attr in self._run_funcs:
            return getattr(self._ds._ds_run, attr)
        
    def __dir__(self):
        all_attrs =  set(self._run_attrs +
                         self._run_funcs + 
                         self.__dict__.keys() + dir(Run))
        
        return list(sorted(all_attrs))


#class RunSteps(object):
#    """Step iterator from psana.DataSource.runs().steps().
#    """
#    def __init__(self, ds, **kwargs):
#        self._ds = ds
#        self._kwargs = kwargs
#        self._ds_steps = []
#        self._configSteps = []
#
#    def __iter__(self):
#        return self
#
#    def next(self):
#        try:
#            self._ds._ievent = -1
#            self._ds._istep +=1
#            self._ds._ds_step = self._ds._current_run.steps().next()
#            self._ds_steps.append(self._ds._ds_step)
#            self._ds._init_detectors()
#            return StepEvents(self._ds)
#        
#        except: 
#            raise StopIteration()


class RunEvents(object):
    """
    Event iterator from ds.runs() for indexed idx data 

    No support yet for multiple runs in a data_source
    """
    def __init__(self, ds, **kwargs):
        self._kwargs = kwargs
        self._ds = ds
        self.times = self._ds.runs.current.times 
        self._ds.nevents = len(self.times)

    def __iter__(self):
        return self

    @property
    def current(self):
        """
        Current event.
        
        See Also
        --------
        EvtDetectors
        """
        return EvtDetectors(self._ds, init=False)

    def next(self, evt_time=None, **kwargs):
        """Optionally pass either an integer for the event number in the data_source
           or a psana.EventTime time stamp to jump to an event.
        
        Parameters
        ----------
        evt_time : object or int
            psana.EventTime time stamp to jump to specific event,
            if int is supplied goto event number in DataSource (may not be exactly
            same event depending on how the data_source string is corresponding
            keywords to define the data_source is defined and also may differ
            for fast feedback and offline analysis environments.
        
        Returns
        -------
        EventDetectors : object
        """
        try:
            if evt_time is not None:
                if isinstance(evt_time, int):
                    self._ds._ievent = evt_time
                else:
                    self._ds._ievent = self.times.index(evt_time)
            else:
                self._ds._ievent += 1
            
            if self._ds._ievent >= len(self.times):
                raise StopIteration()
            else:
                evt = self._ds._ds_run.event(self.times[self._ds._ievent]) 
                self._ds._evt_keys, self._ds._evt_modules = get_keys(evt)
                self._ds._current_evt = evt
                self._ds._current_data = {}
                self._ds._current_evtData = {}
            
            #if hasattr(self._ds, '_idx_istep'):
            #if self._ds._idx_istep:
            #    self._ds._istep = self._ds._idx_istep[self._ds._ievent]

            return EvtDetectors(self._ds, **kwargs)

        except: 
            raise StopIteration()


class SmdEvents(object):
    """
    Event iterator for smd xtc data that iterates first over steps and then
    events in steps (to make sure configData is updated for each step since
    it is possible that it changes).
    """
    def __init__(self, ds, **kwargs):
        self._ds = ds

    @property
    def current(self):
        """Current event.
        """
        return EvtDetectors(self._ds, init=False)

    def __iter__(self):
        return self

    def next(self, evt_time=None, **kwargs):
        """
        Parameters
        ----------
        evt_time : object or int
            psana.EventTime time stamp (second, nanosecond, fiducial) to jump to 
            specific event,
            If int is supplied goto event number in DataSource (may not be exactly
            same event depending on how the data_source string is corresponding
            keywords to define the data_source is defined and also may differ
            for fast feedback and offline analysis environments.
        
        Returns
        -------
        EventDetectors object
            Returns next event in current step.  
            If at end of step goes to next step and returns first event.
        """
        try:
            return self._ds._current_step.next(evt_time=evt_time, **kwargs)
        except:
            try:
                next(self._ds.steps)
                return self._ds._current_step.next(**kwargs)
            except:
                raise StopIteration()


class Steps(object):
    """
    Step iterator from ds.steps().
    """
    def __init__(self, ds, **kwargs):
        self._ds = ds
        self._kwargs = kwargs
        self._ds_steps = []

    @property
    def current(self):
        """
        Current step.
        """
        return self._ds._current_step

    def __iter__(self):
        return self

    def next(self, **kwargs):
        """
        Step iteration method

        Returns
        -------
        StepEvents object
        """
        try:
            if self._ds._istep == self._ds._idx_run.nsteps()-1:
                raise StopIteration()
            else:
                self._ds._ievent = -1
                self._ds._istep +=1
                self._ds._ds_step = next(self._ds._ds.steps())
                self._ds_steps.append(self._ds._ds_step)
                self._ds._init_detectors()
                self._ds._current_step = StepEvents(self._ds)
                return self._ds._current_step

        except: 
            raise StopIteration()


class StepEvents(object):
    """
    Event iterator from ds.steps().events() 
    """
    def __init__(self, ds, **kwargs):
        self._kwargs = kwargs
        self._ds = ds

    @property
    def current(self):
        """Current event.
        """
        return EvtDetectors(self._ds, init=False)

    def __iter__(self):
        return self

    def next(self, evt_time=None, recover=False, **kwargs):
        """
        Next event in step.  If no evt_time provided, the event loop will
        procede from the last event in the step regardless of which event 
        was previously jumped to.

        Optional Parameters
        ----------
        evt_time : object or int
            psana.EventTime time stamp (second, nanosecond, fiducial) to jump to 
            specific event,
            If int is supplied goto event number in DataSource (may not be exactly
            same event depending on how the data_source string is corresponding
            keywords to define the data_source is defined and also may differ
            for fast feedback and offline analysis environments.
        recover : bool
            recover the last step before a jump

        Returns
        -------
        EventDetectors object
            Returns next event in current step.  
        """
        if evt_time is not None:
            # Jump to the specified event
            try:
                if self._ds._istep >= 0:
                    # keep trak of event index to go back to step event loop
                    self._ds._ievent_last = self._ds._ievent
                    self._ds._istep_last = self._ds._istep
                    self._ds._istep = -1
                    self._ds._jump_last = True

                if evt_time.__class__.__name__ == 'EventTime':
                    # lookup event index from time tuple
                    ttup = (evt_time.seconds(), evt_time.nanoseconds(), evt_time.fiducial())
                    self._ds._ievent = self._ds._idx_times_tuple.index(ttup)
                elif isinstance(evt_time, tuple):
                    # optionally accept a time tuple (seconds, nanoseconds, fiducial)
                    self._ds._ievent = self._ds._idx_times_tuple.index(evt_time)
                    evt_time = self._ds._idx_times[self._ds._ievent]
                else:
                    # if an integer was passed jump to the appropriate time from 
                    # the list of run times -- i.e., psana.DataSource.runs().next().times()
                    self._ds._ievent = evt_time
                    evt_time = self._ds._idx_times[evt_time]

                try:
                    if self._ds._scanData is not None:
                        self._ds._istep = self._ds._idx_istep[self._ds._ievent]
                    else:
                        #print('Warning -- must load configData.ScanData before steps can be updated when jumping to events')
                        pass
                
                except:
                    print('Error getting istep')

                #print self._ds._ievent, evt_time.seconds(), evt_time.nanoseconds()
                evt = self._ds._idx_run.event(evt_time) 
                    
                self._ds._evt_keys, self._ds._evt_modules = get_keys(evt)
                self._ds._current_evt = evt
                self._ds._current_data = {}
                self._ds._current_evtData = {}
            
            except:
                print((evt_time, 'is not a valid event time'))
        
        else:
            try:
                if self._ds._jump_last == True and recover == True:
                    # recover event and step index after previoiusly jumping to an event 
                    self._ds._ievent = self._ds._ievent_last
                    self._ds._istep = self._ds._istep_last
                
                self._ds._ievent += 1
                evt = next(self._ds._ds_step.events())
                self._ds._evt_keys, self._ds._evt_modules = get_keys(evt)
                self._ds._current_evt = evt 
                self._ds._current_data = {}
                self._ds._current_evtData = {}
                self._ds._jump_last = False
            except:
                raise StopIteration()

        return EvtDetectors(self._ds, **kwargs)


class Events(object):
    """
    Event iterator
    """

    def __init__(self, ds, **kwargs):
        self._kwargs = kwargs
        self._ds = ds
        self._ds._init_detectors()

    @property
    def current(self):
        """
        Current event
        """
        return EvtDetectors(self._ds, init=False)

    def __iter__(self):
        return self

    def next(self, **kwargs):
        """
        Returns
        -------
        EventDetectors : object
            Returns next event in DataSource.  
        """
        try:
            self._ds._ievent += 1
            evt = next(self._ds._ds.events())
            self._ds._evt_keys, self._ds._evt_modules = get_keys(evt)
            self._ds._current_evt = evt 
            self._ds._current_data = {}
            self._ds._current_evtData = {}

        except:
            raise StopIteration()

        return EvtDetectors(self._ds, **kwargs)


class PsanaTypeList(object):
    """
    Python representation for lists of psana data objects.
    """

    def __init__(self, type_list):
        import numpy as np
        self._type_list = type_list
        typ_func = type_list[0]._typ_func
        module = typ_func.__module__.lstrip('psana.')
        type_name = typ_func.__class__.__name__
        info = psana_doc_info[module][type_name].copy()
        
        self._typ_func = typ_func
        self._values = {}
        self._attr_info = {}
        for attr, item in info.items():
            item['value'] = None

        attrs = [key for key in info.keys() if not key[0].isupper()]
        for attr in attrs:
            values = [getattr(item, attr) for item in self._type_list]

            try:
                if isinstance(values[0], np.ndarray):
                    values = np.array(values)
            except:
                pass
#                print module, type_name, info
#                print values

            if hasattr(values[0], '_typ_func'):
                vals = PsanaTypeList(values)
                for name, item in vals._attr_info.copy().items():
                    alias = attr+'_'+name
                    self._values[alias] = item['value']
                    self._attr_info[alias] = item.copy()
                    self._attr_info[alias]['attr'] = alias

            else:
                self._values[attr] = values
                self._attr_info[attr] = info[attr].copy()
                self._attr_info[attr]['value'] = values
                self._attr_info[attr]['attr'] = attr

        self._attrs = self._values.keys()

    @property
    def _all_values(self):
        """All values in a flattened dictionary.
        """
        avalues = {}
        items = sorted(self._values.items(), key = operator.itemgetter(0))
        for attr, val in items:
            if hasattr(val, '_all_values'):
                for a, v in val._all_values.items():
                    avalues[attr+'_'+a] = v
            else:
                avalues[attr] = val
        return avalues

    def show_info(self, prefix='', **kwargs):
        """Show a table of the attribute, value, unit and doc information
        """
        message = Message(quiet=True, **kwargs)
        items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
        for attr, item in items:
            if attr in self._attrs:
                alias = item.get('attr')
                str_repr = _repr_value(item.get('value'))
                unit = item.get('unit')
                doc = item.get('doc')
                if prefix:
                    alias = prefix+'_'+alias
                message('{:24s} {:>12} {:7} {:}'.format(alias, str_repr, unit, doc))

        return message

    def __getattr__(self, attr):
        if attr in self._attrs:
            return self._values.get(attr)

    def __dir__(self):
        all_attrs = set(self._attrs +
                        self.__dict__.keys() + dir(PsanaTypeList))
        return list(sorted(all_attrs))


class PsanaTypeData(object):
    """
    Python representation of a psana data object (event or configStore data).
    """

    def __init__(self, typ_func, nolist=False):
        if typ_func:
            self._typ_func = typ_func
            module = typ_func.__module__.lstrip('psana.')
            type_name = typ_func.__class__.__name__
        else:
            type_name = None
        self._nolist = nolist

        if type_name in psana_doc_info[module]:
            self._info = psana_doc_info[module][type_name].copy()
            if psana_attrs.get(module,{}).get(type_name):
                self._attrs = [key for key in psana_attrs[module][type_name] if key in self._info.keys()] 
            else:
                self._attrs = [key for key in self._info.keys() if not key[0].isupper()]
            
            #print module,        
 
        else:
            self._attrs = [attr for attr in dir(typ_func) if not attr.startswith('_')]
            self._info = {}

        self._attr_info = {}
        for attr in self._attrs:
            self._attr_info[attr] = _get_typ_func_attr(typ_func, attr, nolist=nolist)
        
#        self._attr_info_new = {}
#        for attr in self._attrs_new:
#            self._attr_info_new[attr] = _get_typ_func_attr(typ_func, attr, nolist=nolist)

    @property
    def _values(self):
        """Dictionary of attributes: values. 
        """
        return {attr: self._attr_info[attr]['value'] for attr in self._attrs}

    @property
    def _all_values(self):
        """
        All values in a flattened dictionary.
        """
        avalues = {}
        items = sorted(self._values.items(), key = operator.itemgetter(0))
        for attr, val in items:
            if hasattr(val, '_all_values'):
                for a, v in val._all_values.items():
                    avalues[attr+'_'+a] = v
            else:
                avalues[attr] = val
        return avalues

    def show_info(self, prefix=None, **kwargs):
        """
        Show a table of the attribute, value, unit and doc information
        """
        message = Message(quiet=True, **kwargs)
        items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
        for attr, item in items:
            value = item.get('value')
            alias = item.get('attr')
            if prefix:
                alias = prefix+'_'+alias
            if hasattr(value, 'show_info'):
                value.show_info(prefix=alias, append=True)
            else:
                str_repr = _repr_value(item.get('value'))
                unit = item.get('unit')
                doc = item.get('doc')
                message('{:24s} {:>12} {:7} {:}'.format(alias, str_repr, unit, doc))

        return message

    def __str__(self):
        return '{:}.{:}.{:}'.format(self._typ_func.__class__.__module__,
                                    self._typ_func.__class__.__name__, 
                                    str(self._typ_func))

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__,str(self))
        return '< '+repr_str+' >'

    def __getattr__(self, attr):
        if attr in self._attrs:
            return self._attr_info[attr]['value']

    def __dir__(self):
        all_attrs = set(self._attrs +
                        self.__dict__.keys() + dir(PsanaTypeData))
        return list(sorted(all_attrs))


class PsanaSrcData(object):
    """
    Python represenation of psana data for a given detector source.
       
    Parameters
    ----------
    key_info : get_key_info(objclass) 
        for faster evt data access.
    """
    def __init__(self, objclass, srcstr, key_info=None, nolist=False):
        self._srcstr = srcstr
        if not key_info:
            key_info = get_key_info(objclass)

        self._types = {}
        self._type_attrs = {}
        self._keys = key_info.get(srcstr)
        if self._keys:
            for (typ, src, key) in self._keys:
                if key:
                    typ_func = objclass.get(*item)
                else:
                    typ_func = objclass.get(typ, src)

                if hasattr(typ_func, '__module__'):
                    module = typ_func.__module__.lstrip('psana.')
                    type_name = typ_func.__class__.__name__
                    type_alias = module+type_name+key 
                    type_data = PsanaTypeData(typ_func, nolist=nolist)
                    self._types[type_alias] = type_data 
                    self._type_attrs.update({attr: type_alias for attr in type_data._attrs})
                    #self._types[(typ,key)] = type_data 
                    #self._type_attrs.update({attr: (typ,key) for attr in type_data._attrs})
                else:
                    pass
                    # psana.EventId is in configStore keys after ana-1.3.10
                    # Not necessary so do nothing for now but avoid error here

    @property
    def _attrs(self):
        attrs = []
        for type_data in self._types.values():
            attrs.extend(type_data._attrs)

        return attrs

    @property
    def _attr_info(self):
        """Attribute information including the unit and doc information 
           and a str representation of the value for all data types.
        """
        attr_info = {}
        for type_data in self._types.values():
            attr_info.update(**type_data._attr_info)

        return attr_info

    @property
    def _values(self):
        """Dictionary of attributes: values for all data types.
        """
        values = {}
        for type_data in self._types.values():
            values.update(**type_data._values)

        return values

    @property
    def _all_values(self):
        """All values in a flattened dictionary.
        """
        values = {}
        for type_data in self._types.values():
            values.update(**type_data._all_values)

        return values

    def show_info(self, **kwargs):
        """Show a table of the attribute, value, unit and doc information
           for all data types of the given source.
        """
        message = Message(quiet=True, **kwargs)
        for type_data in self._types.values():
            type_data.show_info(append=True)

        return message

    def _get_type(self, typ):
        return self._types.get(typ)

    def __str__(self):
        return '{:}'.format(self._srcstr)

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__,str(self))
        return '< '+repr_str+' >'

    def __getattr__(self, attr):
        item = self._type_attrs.get(attr)
        if item:
            return getattr(self._types.get(item), attr)
        
        if attr in self._types:
            return self._types.get(attr)

                   #
    def __dir__(self):
        all_attrs = set(self._type_attrs.keys() +
                        self._types.keys() + 
                        self.__dict__.keys() + dir(PsanaSrcData))
        return list(sorted(all_attrs))


class ConfigData(object):
    """
    Configuration Data representation of configStore within psana.DataSource.env object

    Parameters
    ----------
    ds : DataSource object
    """
    _configStore_attrs = ['get','put','keys']
    # Alias default provides way to keep aliases consistent for controls devices like the FEE_Spec
    _alias_defaults = {
            'BldInfo(FEE-SPEC0)':       'FEE_Spec0',
            'BldInfo(NH2-SB1-IPM-01)':  'Nh2Sb1_Ipm1',
            'BldInfo(NH2-SB1-IPM-02)':  'Nh2Sb1_Ipm2',
            'BldInfo(MFX-BEAMMON-01)':  'MfxBeammon',
            }
    _seq_evtCodes = range(67,99)+range(167,199)+range(201,217)
    _lcls_evtCodes = {
            140: 'Beam & 120Hz',
            141: 'Beam & 60Hz',
            142: 'Beam & 30Hz',
            143: 'Beam & 10Hz',
            144: 'Beam & 5Hz',
            145: 'Beam & 1Hz',
            146: 'Beam & 0.5Hz',
            147: 'Full N-1',
            148: 'Full N-2',
            149: 'TCAV0',
            150: 'Burst',
            151: 'Klys Accel',
            152: 'Klys Standby',
            153: 'Klys Standby no TCAV0',
            154: 'Klys Accel at 10 Hz',
            155: 'Straight Ahead',
            156: 'TCAV3',
            157: 'Klys StdBy No TCAV3',
            158: 'Pockets Cell',
            159: 'Profile Monitors',
            160: 'TCAV3 OTR',
            161: 'BXKIK',
            162: 'BYKIK',
            163: 'A-line Kicker',
            164: 'Test Burst',
            165: 'Spare',
            40: '120 Hz',
            41: '60 Hz',
            42: '30 Hz',
            43: '10 Hz',
            44: '5 Hz',
            45: '1 Hz',
            46: '0.5 Hz',
            }

    def __init__(self, ds):
        configStore = ds.env().configStore()
        if (hasattr(ds, 'data_source') and ds.data_source.monshmserver):
            self._monshmserver = ds.data_source.monshmserver
        else:
            self._monshmserver = None 
       
        self._ds = ds
        self._configStore = configStore
        self._key_info, self._modules = get_keys(configStore)

        # Build _config dictionary for each source
        self._config = {}
        for attr, keys in self._key_info.items():
            try:
                config = PsanaSrcData(self._configStore, attr, 
                                      key_info=self._key_info, nolist=True)
                self._config[attr] = config
            except:
                print(('WARNING:  Cannot load PsanaSrcData config for ', attr, keys))

        self._sources = {}
        #Setup Partition
        if not self._modules.get('Partition'):
            #print 'ERROR:  No Partition module in configStore data.'
            self._partition = {}
            self._srcAlias = {}
            for srcstr, item in self._config.items():
                if srcstr[0:7] in ['BldInfo', 'DetInfo']:
                    alias = srcstr[8:-1]
                    alias = re.sub('-|:|\.| ','_', alias)
                    src = item._keys[0][1] 
                    self._partition[srcstr] = {
                                               #'alias': alias, 
                                               'group': 0, 
                                               'src': src}

                    self._srcAlias[alias] = (src, 0)

                self._bldMask = 0
                self._ipAddrpartition = 0 
                self._readoutGroup = {0: {'eventCodes': [], 'srcs': []}}

        elif len(self._modules['Partition']) != 1:
            print('ERROR:  More than one Partition config type in configStore data.')
            return
        else:
            #Build _partition _srcAlias _readoutGroup dictionaries based on Partition configStore data. 
            type_name = self._modules.get('Partition').keys()[0]
            if len(self._modules['Partition'][type_name]) == 1:
                typ, src, key = self._modules['Partition'][type_name][0]
                srcstr = str(src)
                config = self._config[srcstr]
                self.Partition = config
            else:
                print('ERROR:  More that one Partition module in configStore data.')
                print(('       ', self._modules['Partition'][type_name]))
                return

    # to convert ipAddr int to address 
    # import socket, struct
    # s = key.src()
    # socket.inet_ntoa(struct.pack('!L',s.ipAddr()))

            self._ipAddrPartition = src.ipAddr()
            self._bldMask = config.bldMask
            self._readoutGroup = {group: {'srcs': [], 'eventCodes': []} \
                                  for group in set(config.sources.group)}
            self._partition = {str(src): {'group': config.sources.group[i], 'src': src} \
                               for i, src in enumerate(config.sources.src)}

            self._srcAlias = {}
            if self._modules.get('Alias'):
                for type_name, keys in self._modules['Alias'].items():
                    for typ, src, key in keys:
                        srcstr = str(src)
                        config = self._config[srcstr]
                        ipAddr = src.ipAddr()
                        for i, source in enumerate(config.srcAlias.src):
                            alias = config.srcAlias.aliasName[i]
                            self._srcAlias[alias] = (source, ipAddr)

        self._aliases = {}
        for alias, item in self._srcAlias.items():
            src = item[0]
            ipAddr = item[1]
            srcstr = str(src)
            alias = re.sub('-|:|\.| ','_', alias)
            group = None
            if srcstr in self._partition:
                self._partition[srcstr]['alias'] = alias
                if srcstr.find('NoDetector') == -1:
                    self._aliases[alias] = srcstr
                
                group = self._partition[srcstr].get('group', -1)
            
            elif ipAddr != self._ipAddrPartition or self._monshmserver:
                if self._monshmserver:
                    # add data sources not in partition for live data
                    group = -2
                else:
                    # add data sources not in partition that come from recording nodes
                    group = -1

                self._partition[srcstr] = {'src': src, 'group': group, 'alias': alias}
                self._aliases[alias] = srcstr
                if group not in self._readoutGroup:
                    self._readoutGroup[group] = {'srcs': [], 'eventCodes': []}

                self._sources[srcstr] = {'group': group, 'alias': alias}
            
            if group:
                self._readoutGroup[group]['srcs'].append(srcstr)
            #else:
            #    print 'No group for', srcstr

        # Determine data sources and update aliases
        for srcstr, item in self._partition.items():
            if not item.get('alias') and 'Evr' not in srcstr:
                if srcstr in self._alias_defaults:
                    alias = self._alias_defaults.get(srcstr)
                else:
                    try:
                        alias = srcstr.split('Info(')[1].rstrip(')')
                    except:
                        alias = srcstr
                
                alias = re.sub('-|:|\.| ','_',alias)
                item['alias'] = alias
                self._aliases[alias] = srcstr

            if 'NoDetector' not in srcstr and 'NoDevice' not in srcstr:
                # sources not recorded have group None
                # only include these devices for shared memory
                if srcstr not in self._sources:
                    self._sources[srcstr] = {}
                self._sources[srcstr].update(**item)

        # Make dictionary of src: alias for sources with config objects 
        self._config_srcs = {}
        for attr, item in self._sources.items():
            config = self._config.get(attr)
            if config:
                self._config_srcs[item['alias']] = attr
    
        self._output_maps = {}
        self._evr_pulses = {}
        self._eventcodes = {}

        IOCconfig_type = None
        config_type = None
        for type_name in self._modules['EvrData'].keys():
            if type_name.startswith('IOConfig'):
                IOCconfig_type = type_name
                self._IOCconfig_type = type_name
            elif type_name.startswith('Config'):
                config_type = type_name

        if IOCconfig_type:
            # get eventcodes and combine output_map info from all EvrData config keys
            map_attrs = ['map', 'conn_id', 'module', 'value', 'source_id']
            for typ, src, key in self._modules['EvrData'][config_type]:
                srcstr = str(src)
                config = self._config[srcstr]
                for eventcode in config.eventcodes._type_list:
                    try:
                        # No archive on shared memory currently
                        if not self._monshmserver:
                            self._init_arch()
                    except:
                        traceback.print_exc('Cannot initialize archive')
                    
                    self._eventcodes.update({eventcode.code: eventcode._values})
                    try:
                        code_num = eventcode.code
                        if code_num in self._seq_evtCodes: 
                            owner_pv = 'ECS:SYS0:0:EC_{:}_OWNER_ID'.format(code_num)
                            desc_pv  = 'EVNT:SYS0:1:NAME{:}'.format(code_num)
                            owner_val = int(self._get_pv_from_arch(owner_pv)['data'][0]['val'])
                            desc_val = self._get_pv_from_arch(desc_pv)['data'][0]['val']
                            self._eventcodes[code_num]['description'] = desc_val
                            self._eventcodes[code_num]['owner'] = owner_val
                        elif code_num in self._lcls_evtCodes:
                            self._eventcodes[code_num]['description'] = self._lcls_evtCodes[code_num]
                    except:
                        pass
                    if eventcode.isReadout:
                        group = eventcode.readoutGroup
                        if group not in self._readoutGroup:
                            self._readoutGroup[group] = {'srcs': [], 'eventCodes': []}
                        self._readoutGroup[group]['eventCodes'].append(eventcode.code)

                for output_map in config.output_maps._type_list:
                    map_key = (output_map.module,output_map.conn_id)
                    if output_map.source == 'Pulse':
                        pulse_id = output_map.source_id
                        pulse = config.pulses._type_list[pulse_id]
                        evr_info = { 'evr_width': pulse.width*pulse.prescale/119.e6, 
                                     'evr_delay': pulse.delay*pulse.prescale/119.e6, 
                                     'evr_polarity': pulse.polarity}
                    else:
                        pulse_id = None
                        pulse = None
                        evr_info = {'evr_width': None, 'evr_delay': None, 'evr_polarity': None}

                    self._output_maps[map_key] = {attr: getattr(output_map,attr) for attr in map_attrs} 
                    self._output_maps[map_key].update(**evr_info) 

            # Assign evr info to the appropriate sources
            if len(self._modules['EvrData'][IOCconfig_type]) > 1:
                print('WARNING: More than one EvrData.{:} objects'.format(IOCconfig_type))

            IOCconfig_type = self._IOCconfig_type
            typ, src, key = self._modules['EvrData'][IOCconfig_type][0]
            srcstr = str(src)
            config = self._config[srcstr]
            if config._values['channels']:
                for ch in config._values['channels']._type_list:
                    map_key = (ch.output.module, ch.output.conn_id)
                    for i in range(ch.ninfo):
                        src = ch.infos[i]
                        srcstr = str(src)
                        self._sources[srcstr]['map_key'] = map_key
                        for attr in ['evr_width', 'evr_delay', 'evr_polarity']:
                            self._sources[srcstr][attr] = self._output_maps[map_key][attr]

            for group, item in self._readoutGroup.items():
                if item['eventCodes']:
                    for srcstr in item['srcs']: 
                        if srcstr in self._sources:
                            self._sources[srcstr]['eventCode'] = item['eventCodes'][0]

        # Get control data
        if self._modules.get('ControlData'):
            type_name, keys = self._modules['ControlData'].items()[0]
            typ, src, key = keys[0]
            config = self._config[str(src)]
            self._controlData = config._values
            self.ControlData = config

        if self._modules.get('SmlData'):
            type_name, keys = self._modules['SmlData'].items()[0]
            typ, src, key = keys[0]
            config = self._config[str(src)]
            self._smlData = config._values

        try:
            from . import config_check
            self.configCheck = config_check.ConfigCheck(self)
        except:
            traceback.print_exc()

    def _init_arch(self):
        """
        Epics Archive access
        """
        import pandas as pd
        from .epicsarchive import EpicsArchive
        self._arch = EpicsArchive()
        try:
            dt = pd.Timestamp(min(self._ds._idx_datetime64))
            self._tstart = [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]
            dt = pd.Timestamp(max(self._ds._idx_datetime64))
            self._tend = [dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second]
        except:
            pass

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

    def save_configData(self):
        """
        Save sources, eventCodes and ScanData xarray datasets
        """
        run = self._ds.data_source.run
        path = os.path.join(self._ds.data_source.res_dir,'nc')
        try:
            if self.ScanData.dataset.dims.get('step') > 1:
                scan_file='{:}/run{:04}_{:}.nc'.format(path, run, 'scan')
                self.ScanData.dataset.to_netcdf(scan_file, engine='h5netcdf')
        except:
            print('Cannot save scan config for {:}'.format(self._ds))
        
        try:
            eventCode_file='{:}/run{:04}_{:}.nc'.format(path, run, 'eventCodes')
            self.Sources.eventCodes.to_netcdf(eventCode_file, engine='h5netcdf')
        except:
            print('Cannot save eventCode config for {:}'.format(self._ds))
            
        try:
            eventCode_file='{:}/run{:04}_{:}.nc'.format(path, run, 'eventCodes')
            sources_file='{:}/run{:04}_{:}.nc'.format(path, run, 'sources')
            self.Sources.sources.to_netcdf(sources_file, engine='h5netcdf')
        except:
            print('Cannot save source config for {:}'.format(self._ds))

    @property
    def Sources(self):
        """
        Source information including evr config.
        
        Returns
        -------
        ConfigSources object
        
        """
        return ConfigSources(self)

    @property
    def ScanData(self):
        """
        Scan configuration from steps ControlData.  
        May take several seconds to load the first time.
        Only relevant for smd data.
        
        Returns
        -------
        ScanData object
        """
        #if self._ds.data_source.monshmserver is not None:
        if not self._ds.data_source.smd:
            return None

        if self._ds._scanData is None:
            self._ds._scanData = ScanData(self._ds)
            print(self._ds._scanData)

        return self._ds._scanData

    def show_info(self, show_codes=True, **kwargs):
        """
        Show Detector Source information.
        """
        message = Message(quiet=True, **kwargs)
        message('-'*80)
        message('Source Information:')
        message('-'*18)
        self.Sources.show_info(append=True, **kwargs)
        if show_codes:
            message('-'*80)
            message('Event Code Information:')
            message('-'*18)
            self.Sources.show_eventCodes(append=True)
 
        return message

    def get_info(self, **kwargs):
        """
        Get Detector Source Information.
        """
        return str(self.show_info(**kwargs))

    def show_all(self, **kwargs):
        """
        Show Detector Source and ScanData information.
        """
        message = Message(quiet=True, **kwargs)
        message('-'*80)
        message('Source Information:')
        message('-'*18)
        self.Sources.show_info(append=True)
        message('') 
        message('-'*80)
        message('Scan Data:')
        message('-'*18)
        self.ScanData.show_info(append=True)

        return message

    def __str__(self):
        return  'ConfigData: '+str(self._ds.data_source)

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__,str(self._ds.data_source))
        print('< '+repr_str+' >')
        print(self.Sources.show_info())
        return '< '+repr_str+' >'
    
    def __getattr__(self, attr):
        if attr in self._config_srcs:
            return self._config[self._config_srcs[attr]]

        if attr in self._configStore_attrs:
            return getattr(self._configStore, attr)
        
    def __dir__(self):
        all_attrs = set(self._configStore_attrs +
                        self._config_srcs.keys() + 
                        self.__dict__.keys() + dir(ConfigData))
        return list(sorted(all_attrs))


class EvtDetectors(object):
    """
    Psana tab accessible event detectors.
    All detectors in Partition or defined in any configStore Alias object 
    (i.e., recording nodes as well as daq) return the relevant attributes of 
    a PyDetector object for that src, but only the sources in the evt.keys()
    show up in the ipython tab accessible dir.
    
    Preserves get, keys and run method of items in psana events iterators.
    
    Parameters
    ----------
    ds : DataSource object
    """

    _init_attrs = ['get', 'keys'] #  'run' depreciated
    _event_attrs = ['EventId', 'Evr', 'L3T']

    def __init__(self, ds, publish=True, init=True, update_stats=True): 
        self._ds = ds
        if init:
            self._init(publish=publish, update_stats=update_stats)

    def _init(self, publish=True, update_stats=True):
        if publish:
            psmon_publish(self)
        if update_stats:
            _update_stats(self)
 
    @property
    def EventId(self):
        """
        EventId object
        """
        return EventId(self._ds._current_evt)

    @property
    def _attrs(self):
        """
        List of detector names in current evt data.
        """
        return [alias for alias, srcstr in self._ds._aliases.items() \
                                        if srcstr in self._ds._evt_keys]

    @property
    def _dets(self):
        """
        Dictionary of detectors.
        
        Returns
        -------
        DataSource._detectors dict
        """
        return self._ds._detectors

    @property
    def Evr(self):
        """
        Master evr from psana evt data.
        
        Returns
        -------
        EvrData object
        """
        if 'EvrData' in self._ds._evt_modules:
            return EvrData(self._ds)
        else:
            return EvrNullData(self._ds)

    @property
    def L3T(self):
        """
        L3T Level 3 trigger.
        
        Returns
        -------
        L3Tdata object
        """
        if 'L3T' in self._ds._evt_modules:
            return L3Tdata(self._ds)
        else:
            return L3Ttrue(self._ds)

    @property
    def run(self):
        """
        Run number.  For shared memory when not recording will be a large int
        
        Returns
        -------
        int
        """
        return self._ds._current_evt.run()

    def next(self, *args, **kwargs):
        """
        Returns
        -------
        EvtDetectors : object
            Next event in DataSource (behavior depends on if smd, idx options used in data_source)
        """
        evt = self._ds.events.next(*args, **kwargs)
        return evt

    def monitor(self, nevents=-1, sleep=None):
        """
        Monitor events continuously.
        """ 
        ievent = nevents
        while ievent != 0:
            try:
                next(self)
                try:
                    print(self)
                except:
                    pass
                
                if ievent < nevents and sleep:
                    time.sleep(sleep)

                ievent -= 1

            except KeyboardInterrupt:
                ievent = 0


    def __iter__(self):
        return self

    def __str__(self):
        return  '{:}, Run {:}, Step {:}, Event {:}, {:}, {:}'.format(self._ds.data_source.exp, 
                self.run, self._ds._istep, self._ds._ievent, 
                str(self.EventId), str(self.Evr))

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__, str(self))
        return '< '+repr_str+' >'

    def __getattr__(self, attr):
        if attr in self._ds._detectors:
            return self._ds._detectors[attr]
        
        if attr in self._init_attrs:
            return getattr(self._ds._current_evt, attr)

    def __dir__(self):
        all_attrs =  set(self._attrs +
                         self._init_attrs +
                         self.__dict__.keys() + dir(EvtDetectors))
        
        return list(sorted(all_attrs))


class L3Ttrue(object):
    """
    L3 Trigger default if no L3 Trigger data is in DataSource.
    Typically only used for older data where no L3 Trigger data was generated.
    """

    def __init__(self, ds):

        self._ds = ds
        self._attr_info = {'result': {'attr': 'result',
                                      'doc':  'No L3T set',
                                      'unit': '',
                                      'value': True}}

        self._attrs = self._attr_info.keys()

    @property
    def _values(self):
        """
        Dictionary of attributes: values. 
        """
        return {attr: self._attr_info[attr]['value'] for attr in self._attrs}

    def show_info(self, **kwargs):
        """
        Show a table of the attribute, value, unit and doc information
        """
        message = Message(**kwargs)
        items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
        for attr, item in items:
            value = item.get('value')
            if hasattr(value, 'show_info'):
                value.show_info(prefix=attr, append=True)
            else:
                item['str'] = _repr_value(value)
                message('{attr:24s} {str:>12} {unit:7} {doc:}'.format(**item))

        return message

    def __str__(self):
        return str(self.result)

    def __repr__(self):
        return '< {:}: {:} >'.format(self.__class__.__name__, str(self))

    def __getattr__(self, attr):
        if attr in self._attrs:
            return self._attr_info[attr]['value']

    def __dir__(self):
        all_attrs = set(self._attrs +
                        self.__dict__.keys() + dir(L3Ttrue))
        return list(sorted(all_attrs))


class L3Tdata(PsanaTypeData):
    """
    L3 Trigger data
    """
    def __init__(self, ds):

        self._typ, self._src, key = ds._evt_modules['L3T'].values()[0][0]
        typ_func = ds._current_evt.get(self._typ,self._src)
        PsanaTypeData.__init__(self, typ_func)

    def __str__(self):
        return str(self.result)

    def __repr__(self):
        return '< {:}: {:} >'.format(self.__class__.__name__, str(self))


class ConfigSources(object):
    """
    Configuration Sources

    Parameters
    ----------
    configData : ConfigData object

    """

    def __init__(self, configData):
        self._sources = configData._sources
        self._aliases = {item['alias']: src for src, item in self._sources.items()}
        self._cfg_srcs = configData._config_srcs.values()
        self._repr = str(configData._ds) 
        self._configData = configData
        self._eventcodes = configData._eventcodes

    @property
    def sources(self):
        """
        xarray of sources
        """
        import pandas as pd
        xsources = pd.DataFrame(self._sources).T.to_xarray().swap_dims({'index': 'alias'}).drop('src').rename({'index': 'src'}).drop('map_key')
        for attr in xsources.data_vars:
            xsources[attr].attrs['doc'] = SourceData._doc.get(attr,'')
            xsources[attr].attrs['unit'] = SourceData._units.get(attr,'')
            xsources.attrs['experiment'] = self._configData._ds.data_source.exp
            xsources.attrs['run'] = self._configData._ds.data_source.run
            xsources.attrs['instrument'] = self._configData._ds.data_source.instrument
            xsources.attrs['data_source'] = self._configData._ds.data_source.data_source
            xsources.attrs['expNum'] = self._configData._ds.expNum

        return xsources

    @property
    def eventCodes(self):
        """
        xarray of event codes
        """
        import pandas as pd
        xcodes = pd.DataFrame(self._eventcodes).T.to_xarray().swap_dims({'index': 'code'}).drop('index').drop('desc_shape') 
        xcodes.attrs['experiment'] = self._configData._ds.data_source.exp
        xcodes.attrs['run'] = self._configData._ds.data_source.run
        xcodes.attrs['instrument'] = self._configData._ds.data_source.instrument
        xcodes.attrs['data_source'] = self._configData._ds.data_source.data_source
        xcodes.attrs['expNum'] = self._configData._ds.expNum
        return xcodes

    def show_info(self, **kwargs):
        message = Message(quiet=True, **kwargs)
        message('*Detectors in group 0 are "BLD" data recorded at 120 Hz on event code 40')
        if self._monshmserver:
            message('*Detectors listed as Monitored are not being recorded (group -2).')
        else:
            message('*Detectors listed as Controls are controls devices with unknown event code (but likely 40).')
        message('')
        header =  '{:20} {:>3} {:>13} {:>4} {:>3} {:>11} {:>11} {:30}'.format('Alias', 'Grp', 
                 'Description', 'Code', 'Pol', 'Delay [s]', 'Width [s]', 'Source') 
        message(header)
        message('-'*(len(header)))
        data_srcs = {item['alias']: s for s,item in self._sources.items() \
                       if s in self._cfg_srcs or s.startswith('Bld')}
        
        for alias, srcstr in sorted(data_srcs.items(), key = operator.itemgetter(0)):
            item = self._sources.get(srcstr,{})

            polarity = item.get('evr_polarity', '')
            # Epics convention used -- has it always been this?
            # I was previously using the opposite convetion here
            if polarity == 0:
                polarity = 'Pos'
            elif polarity == 1:
                polarity = 'Neg'

            delay = item.get('evr_delay', '')
            if delay:
                delay = '{:11.9f}'.format(delay)

            width = item.get('evr_width', '')
            if width:
                width = '{:11.9f}'.format(width)

            group = item.get('group')
            if group == -1:
                group = 'Controls'
            elif group == -2:
                group = 'Monitor'

            if group == 0:
                eventCode = 40
            else:
                eventCode = item.get('eventCode', '')

            rate = _eventCodes_rate.get(eventCode, '')

            try:
                description = self._eventcodes.get(eventCode, {'description': ''})['description']
                if not description:
                    description = rate
            except:
                description = rate

            message('{:20} {:>3} {:>13} {:>4} {:>3} {:>11} {:>11} {:30}'.format(alias, 
                   group, description, eventCode, polarity, delay, width, srcstr))
                   #group, rate, eventCode, polarity, delay, width, srcstr))
       
        return message

    def show_eventCodes(self, **kwargs):
        """
        Show event code information
        """
        message = Message(quiet=True, **kwargs)
        header = '{:8} {:4} {:12} {:20}'.format('code',  'group', 'type', 'description') 
        message(header)
        message('-'*(len(header)+10))
        items = sorted(self._eventcodes.items(), key=operator.itemgetter(0))
        for code, item in items:
            if item.get('isReadout') == 1:
                codetype = 'Readout'
            elif item.get('isCommand') == 1:
                codetype = 'Command'
            elif item.get('isLatch') == 1:
                codetype = 'Latch'
            else:
                codetype = ''
            message('{:8} {:4} {:12} {:20}'.format(code, item.get('readoutGroup'), codetype, item.get('description')))
        
        return message

    def __str__(self):
        return  self._repr

    def __repr__(self):
        repr_str = '{:}: {:}'.format(self.__class__.__name__, self._repr)
        return '< '+repr_str+' >'
 
    def __getattr__(self, attr):
        if attr in self._aliases:
            srcstr = self._aliases[attr]
            return SourceData(self._sources[srcstr])

    def __dir__(self):
        all_attrs =  set(self._aliases.keys() +
                         self.__dict__.keys() + dir(ConfigSources))
        
        return list(sorted(all_attrs))


class SourceData(object):
    """
    Source information for a daq device.

    Parameters
    ----------
    source : str
        Source name

    Attributes
    ----------
    evr_width : float
        Detector evr trigger width ['s']
    evr_delay : float
        Detector evr trigger delay ['s']
    evr_polarity : bool
        Detector evr trigger polarity
    group: int
        Daq evr group
    map_key: int tuple
        Evr configuation map key (card,channel)
    eventCode: int
        Evr event code
    """

    _units = {'evr_width': 's',
              'evr_delay': 's'}

    _doc = {'evr_width': 'Evr trigger width',
            'evr_delay': 'Evr trigger delay',
            'evr_polarity': 'Evr trigger polarity',
            'group': 'Evr group',
            'map_key': 'Evr configuation map key (card,channel)',
            'eventCode': 'Evr event code',
            }

    def __init__(self, source):
        self._source = source

    def show_info(self, **kwargs):
        message = Message(quiet=True, **kwargs)
        for attr, val in sorted(self._source.items(), key=operator.itemgetter(0)):
#        for attr in self._source:
#            val = self._source[attr]
            if attr in self._units:
                val = '{:10.9f}'.format(val)
            item = [attr, val, self._units.get(attr, ''), self._doc.get(attr, '')]
            message('{:22} {:>12} {:3} {:40}'.format(*item))

        return message

    def __str__(self):
        return '{:} = {:}'.format(self._source.get('alias'), self._source.get('src'))

    def __repr__(self):
        return '< {:}: {:} >'.format(self.__class__.__name__, str(self))

    def __getattr__(self, attr):
        if attr in self._source:
            return self._source[attr]

    def __dir__(self):
        all_attrs =  set(self._source.keys()+
                         self.__dict__.keys() + dir(SourceData))
        
        return list(sorted(all_attrs))

class EvrData(object):

    """
    Evr eventCode information for current event in DataSource 

    Parameters
    ----------
    ds : DataSource object

    """
    _detail_attrs  = ['present', 'timestampHigh', 'timestampLow',
                      '_attr_info', '_src', '_typ', '_info', '_typ_func',
                      'show_info', 'show_table', '_present', '_values', '_all_values']

    def __init__(self, ds):
        self._ds = ds
        self._evr = None
    
    @property
    def eventCodes_strict(self):
        """
        eventcodes that were sent on precisely the fiducial corresponding to the evt.
        """
        attr = 'eventCodes'
        return getattr(self._ds._evr, attr)(self._ds._current_evt, this_fiducial_only=True)

    @property
    def eventCodes(self):
        """
        All event codes in current event
        """
        attr = 'eventCodes'
        return getattr(self._ds._evr, attr)(self._ds._current_evt)

    @property
    def _ddls(self):
        attr = '_fetch_ddls'
        return getattr(self._ds._evr, attr)(self._ds._current_evt)

    @property
    def _types(self):
        attr = '_find_types'
        return getattr(self._ds._evr, attr)(self._ds._current_evt)

    @ property
    def name(self):
        """
        Evr name
        """
        return self._ds._evr.name

    @ property
    def source(self):
        """
        Evr source
        """
        return self._ds._evr.source

    @property
    def EventId(self):
        """
        EventId object
        """
        return EventId(self._ds._current_evt)
    
    def __str__(self):
        try:
            eventCodeStr = '{:}'.format(self.eventCodes)
        except:
            eventCodeStr = ''
        
        return eventCodeStr

    def __repr__(self):
        return '< {:}: {:} >'.format(self.__class__.__name__, str(self))

    def __getattr__(self, attr):
        if attr in self._detail_attrs:
            if self._evr is None:
                self._evr = EvrDataDetails(self._ds)
            return getattr(self._evr, attr)
    
    def __dir__(self):
        all_attrs =  set(self._detail_attrs +
                         self.__dict__.keys() + dir(EvrData) + dir(EvrDataDetails))
        
        return list(sorted(all_attrs))


class EvrDataDetails(PsanaTypeData):
    """
    Evr eventCode information for current event in DataSource 

    Parameters
    ----------
    ds : DataSource object
    """

    def __init__(self, ds):

        self._typ, self._src, key = ds._evt_modules['EvrData'].values()[0][0]
        typ_func = ds._current_evt.get(self._typ,self._src)
        PsanaTypeData.__init__(self, typ_func)
        self.eventCodes = self.fifoEvents.eventCode
        self._ds = ds

    @property
    def _strict_codes(self):
        """
        Event codes strictly defined for the readout event.
        """
        fiducials = self.EventId.fiducials
        return [code for code, tsh in self.timestampHigh.items() if tsh == fiducials]

    def _present(self, eventCode, strict=True):
        """
        Return True if the eventCode is present.
        """
        try:
            if strict:
                return self.timestampHigh[eventCode] == self.EventId.fiducials
            else:
                return eventCode in self.eventCodes
        except:
            pass

        if hasattr(self._typ_func, 'present'):
            try:
                pres = self._typ_func.present(eventCode)
                if pres:
                    if strict and self.timestampHigh[eventCode] != self.EventId.fiducials:
                        return False
                    else:
                        return True
                else:
                    return False
            except:
                print('{:} {:}: error checking for eventCode {:}'.format(self.EventId, self, eventCode))
        
        try:
            return (eventCode in self.eventCodes)
        except:
            pass
            
        return False

    def present(self, *args, **kwargs):
        """
        Check if the event has specified event code.
        Multiple event codes can be tested.
        
        Parameters
        ----------
        strict: bool (default=True)
            check if code has same timestamp as EventId

        Example
        -------
        
        Assume: 
            self.eventCodes = [41, 140]
        Then:
            self.present(41, 140) = True
            self.present(42, 140) = Flase
            
        To check if an event code is not present use a negative number:
            self.present(-41) = False
            self.present(-41, 140) = False
            self.present(-42) = True
            self.present(-42, 140) = True
            self.present(-42, 41, 140) = True
        """
        if args[0] is None:
            return True

        if len(args) == 1:
            if isinstance(args[0], list):
                eventCodes = {arg for arg in args[0]}
            else:
                eventCodes = args
        else:
            eventCodes = args

        for eventCode in eventCodes:
            if (eventCode > 0 and not self._present(eventCode, **kwargs)):
                return False

            if (eventCode < 0 and self._present(abs(eventCode), **kwargs)):
                return False

        return True

    @property
    def timestampHigh(self):
        """
        timstampHigh dict from fifoEvents data
        """
        return dict(zip(self.fifoEvents.eventCode, self.fifoEvents.timestampHigh))

    @property
    def timestampLow(self):
        """
        timstampLow dict from fifoEvents data
        """
        return dict(zip(self.fifoEvents.eventCode, self.fifoEvents.timestampLow))
    
    @property
    def EventId(self):
        """
        EventId object
        """
        return EventId(self._ds._current_evt)


    def show_table(self, **kwargs):
        """Show table of event codes present'
        """
        message = Message(quiet=True, **kwargs)
        ecs = self.fifoEvents
        message('{:>6} {:>8} {:>8} {:>9}'.format('Code','TS_high', 'TS_low', 'Present'))
        for i in range(self.numFifoEvents):
            message('{:6} {:8} {:8}   {:}'.format(ecs.eventCode[i], 
                            ecs.timestampHigh[i], ecs.timestampLow[i], self.present(ecs.eventCode[i])))
        return message

    def __str__(self):
        try:
            eventCodeStr = '{:}'.format(self.eventCodes)
        except:
            eventCodeStr = ''
        
        return eventCodeStr

    def __dir__(self):
        all_attrs =  set(self.__dict__.keys() + dir(EvrDataDetails))
        
        return list(sorted(all_attrs))


class EvrNullData(object):
    """
    Evr data class when no EvrData type is in event keys.
    Occurs for controls cameras with no other daq data present.
    """

    def __init__(self, ds):
        self.eventCodes = []
        self.eventCodes_strict = []

    def __str__(self):
        return ''


class EventId(object):
    """
    Time stamp information from psana EventId. 
    """

    _attrs = ['fiducials', 'idxtime', 'run', 'ticks', 'time', 'vector']
    _properties = ['datetime64', 'EventTime', 'timef64', 'nsec', 'sec']

    def __init__(self, evt):
        import numpy as np
        self._EventId = evt.get(psana.EventId)

    @property
    def datetime64(self):
        """
        NumPy datetime64 representation of EventTime.
        
        References 
        ----------

        http://docs.scipy.org/doc/numpy/reference/arrays.datetime.html
        
        """
        return np.datetime64(int(self.sec*1e9+self.nsec), 'ns')

    @property
    def EventTime(self):
        """
        psana.EventTime for use in indexed idx xtc files.
        """
        return psana.EventTime(int((self.sec<<32)|self.nsec), self.fiducials)

    @property
    def timef64(self):
        """
        Event time represented as float64
        """
        return np.float64(self.sec)+np.float64(self.nsec)/1.e9 

    @property
    def nsec(self):
        """
        Nanosecond part of event time.
        """
        return self.time[1]

    @property
    def sec(self):
        """
        Second part of event time.
        """
        return self.time[0]

    def show_info(self, **kwargs):
        """
        Returns
        -------
        Message object
            Time stamp information for event 
        """
        message = Message(quiet=True, **kwargs)
        message(self.__repr__())
        for attr in self._attrs:
            if attr != 'idxtime': 
                message('{:18s} {:>12}'.format(attr, getattr(self, attr)))

        return message

    def __str__(self):
        try:
            EventTimeStr = time.strftime('%H:%M:%S',
                    time.localtime(self.time[0]))
            EventTimeStr += '.{:04}'.format(int(self.time[1]/1e5))
        except:
            EventTimeStr = 'NA'

        return '{:}'.format(EventTimeStr)

    def __repr__(self):
        return '< {:}: {:} >'.format(self.__class__.__name__, str(self))

    def __getattr__(self, attr):
        if attr in self._attrs:
            return getattr(self._EventId, attr)()

    def __dir__(self):
        all_attrs =  set(self._attrs+
                         self.__dict__.keys() + dir(EventId))
        
        return list(sorted(all_attrs))


class Detector(object):
    """
    Includes epicsData, configData, evrConfig info 
    Uses full ds in order to be able to access epicsData info on
    an event basis.
    """
    #import Detector as psanaDetector
    # dict of detectors in psana.DetectorTypes.detectors
   
    _tabclass = 'evtData'
    _calib_class = None
    _det_class = None
    _ds = None
    _alias = ''
    src = ''
    _pydet = None
    _pydet_name = None
    _init = False
    _xarray_init = False

    def __init__(self, ds, alias, verbose=False, **kwargs):
        """
        Initialize a psana Detector class for a given detector alias.
        Provides the attributes of the PyDetector functions for the current 
        event if applicable.  Otherwise provides the attributes from the
        raw data in the psana event keys for the given detector.
        """

        self._alias = alias
        self._ds = ds
        self.src = ds._aliases.get(alias)
        if not self._xarray_info:
            self._xarray_info.update({'coords': {}, 'dims': {}, 'attrs': {}})
        
        self._source = ds.configData._sources.get(self.src)
        #

        if self.src:
            if verbose:
                print('Adding Detector: {:20} {:40}'.format(alias, psana.Source(self.src)))
        
        else:
            print('ERROR No Detector with alias {:20}'.format(alias))
            return

        self._srcstr = str(self.src)
        #self._srcname = self._srcstr.split('(')[1].split(')')[0]
        srcname = self._det_config['srcname']
        #self._srcname = str(self.srcname)
        #self._devName = self._det_config['devName']
 
        try:
            self._pydet = psana.Detector(srcname, ds._ds.env())
            self._pydet_name = self._pydet.__class__.__name__
        except:
            pass

        _pydet_dict = {
                'AreaDetector':      {'det_class': ImageData,
                                      'calib_class': ImageCalibData},
                'GenericWFDetector': {'det_class': GenericWaveformData,
                                     'calib_class': GenericWaveformCalibData},
                'WFDetector':        {'det_class': WaveformData,
                                      'calib_class': WaveformCalibData},
                'IpimbDetector':     {'det_class': IpimbData,
                                      'calib_class': None},
                }

        item = _pydet_dict.get(self._pydet_name)
        if item:
            self._det_class = item.get('det_class')
            self._calib_class = item.get('calib_class')
            self._tabclass = 'detector'

        try:
            self._on_init(**kwargs)
            self._init = True
            # Epix10ka calib class not yet implmented
            if self._det_config['devName'] == 'Epix10ka':
                self._calib_class = None
                self._attrs.remove('image')
                self._attrs.remove('calib')
                self.shape = (352,384)
                self.size = 135168

        except:
            pass


#    # Does not quite work.  Need to invesigate how to rename robustly
#    def rename(self, alias):
#        """
#        Rename detector.
#        """
#        self._ds.rename(**{self._alias: alias})

    def _on_init(self, **kwargs):
        """
        Method for Detector instances.
        """
        pass

    @property
    def _source_info(self):
        """
        sourceData information.  Objects converted to str.
        """
        src_info = self.sourceData._source.copy()
        if 'src' in src_info:
            src_info['src'] = str(src_info['src'])
        return src_info

    @property
    def add(self):
        """
        Add data processing methods.

        See Also
        --------
        AddOn : class
            Methods to add parameters, properties, and reduction/proccesing of data 
            with roi, projection, histogram methods.
        """
        return AddOn(self._ds, self._alias)
   
    @property
    def _det_config(self):
        """
        Detector configuration

        See Also
        --------
        _device_sets : dict
            Dictionary of detector devices in DataSource object.
        """
        return self._ds._device_sets.get(self._alias)

    @property
    def _xarray_info(self):
        """
        Dictionary of information to build xarray data summaries.
        """
        return self._det_config['xarray']

    def _update_xarray_info(self):
        """
        Update default xarray information
        """
        import numpy as np
        # attrs -- not valid yet for bld but should fix this to avoid try/except here
        self._xarray_init = True
        try:
            attrs = {attr: item for attr, item in self.configData._all_values.items() \
                if np.product(np.shape(item)) <= 17}
        except:
            attrs = {}

        self._xarray_info['attrs'].update(**attrs)

        # xarray dims
        if self.src == 'BldInfo(FEE-SPEC0)':
            dims_dict = {attr: ([], ()) for attr in ['integral']}
            dims_dict['hproj'] = (['X'], self.hproj.shape)

        elif self._det_class == GenericWaveformData:
            dims_dict = {'waveform': (['ch', 't'], (16, 4096))} 

        elif self._det_class == WaveformData:
            if self._pydet.dettype == 17:
                dims_dict = {'waveform': (['ch', 't'], (4, self.configData.numberOfSamples))} 
            else:
                dims_dict = {
                        'waveform':  (['ch', 't'], 
                            (self.configData.nbrChannels, self.configData.horiz.nbrSamples)),
                        }

        elif self._det_class == IpimbData:
            attr_info = self.detector._attr_info
            dims_dict = {
                    'sum':      ([], (), attr_info.get('sum',{})),
                    'xpos':     ([], (), attr_info.get('xpos',{})),
                    'ypos':     ([], (), attr_info.get('ypos',{})),
                    'channel':  (['ch'], (4,)),
                    }

        elif self._det_class == ImageData:
            if self.calibData.ndim == 3:
                raw_dims = (['sensor', 'row', 'column'], self.calibData.shape)
                dims_dict = {
                    'calib':     raw_dims,
                    'raw':       raw_dims,
                    'corr':      raw_dims,
                    }
                try:
                    image_shape = (self.calibData.image_xaxis.size, self.calibData.image_yaxis.size)
                    image_dims = (['X', 'Y'], image_shape) 
                    dims_dict.update({'image':     image_dims})
                except:
                    pass

            else:
                if self.calibData.image_xaxis is not None and self.calibData.image_xaxis.size > 0:
                    raw_dims = (['xaxis', 'yaxis'], self.calibData.shape)
                    image_shape = (len(self.calibData.image_xaxis),len(self.calibData.image_yaxis))
                    image_dims = (['X', 'Y'], image_shape)
                    #dims_dict = {'calib':     image_dims}
                    dims_dict = {
                        'image':     image_dims,
                        'raw':       raw_dims,
                        'corr':      raw_dims,
                        'calib':     raw_dims,
                        'photons':   raw_dims,
                        }
                elif hasattr(self, 'raw') and self.raw is not None:
                    dshape = self.raw.shape
                    raw_dims = (['xaxis', 'yaxis'], dshape)
                    image_dims = (['X', 'Y'], dshape)
                    dims_dict = {
                        'image':     image_dims,
                        'raw':       raw_dims,
                        'corr':      raw_dims,
                        'calib':     raw_dims,
                        'photons':   raw_dims,
                        }
                elif hasattr(self.calibData, 'shape') and self.calibData.shape is not None:
                    # Shape is not alwasy correct as of ana-1.3.4 for roi, e.g., Timetool
                    raw_dims = (['xaxis', 'yaxis'], self.calibData.shape)
                    image_dims = (['X', 'Y'], self.calibData.shape)
                    dims_dict = {
                        'image':     image_dims,
                        'raw':       raw_dims,
                        'corr':      raw_dims,
                        'calib':     raw_dims,
                        'photons':   raw_dims,
                        }
                else:
                    try:
                        dims_dict = {'raw': (['X', 'Y'], self.evtData.data16.shape)}
                    except:
                        if hasattr(self.calibData, 'data8'):
                            dims_dict = {'raw': (['X', 'Y'], self.evtData.data8.shape)}
                        else:
                            try:
                                next(self)
                                dims_dict = {'raw': (['X', 'Y'], self.raw.shape)}
                            except:
                                print(('Error adding dims for ', str(self)))

        # temporary fix for Quartz camera not in PyDetector class
        elif self._pydet is not None and hasattr(self._pydet, 'dettype') \
                and self._pydet.dettype == 18:
            try:
                dims_dict = {'data8': (['X', 'Y'], self.data8.shape)}
            except:
                print((str(self), 'Not valid data8'))
        
        else:
            dims_dict = {}
            for attr, val in self.evtData._all_values.items():
                npval = np.array(val)
                info = self.evtData._attr_info.get(attr)
                if info:
                    xattrs = {a: b for a, b in info.items() if a in ['doc','unit']}
                else:
                    xattrs = {}
                if npval.size > 1:
                    dims_dict[attr] = (['d{:}_{:}'.format(i,a) for i,a in enumerate(npval.shape)], npval.shape, xattrs)
                else:
                    dims_dict[attr] = ([], (), xattrs)

#            dims_dict = {attr: ([], ()) for attr in self.evtData._all_values}
                    
        self._xarray_info['dims'].update(**dims_dict)

        # coords
        if self._det_class == WaveformData:
            if self._pydet.dettype == 17:
                coords_dict = {
                        't': np.arange(self.configData.numberOfSamples) 
                        }
            else:
                coords_dict = {
                        't': np.arange(self.configData.horiz.nbrSamples) \
                                *self.configData.horiz.sampInterval \
                                +self.configData.horiz.delayTime
                        }

        # temporary fix for Opal2000, Opal4000 and Opa8000
        elif self._det_class == ImageData:
            if self.calibData.ndim == 3:
                raw_dims = (['sensor', 'row', 'column'], self.calibData.shape)
                attrs = ['areas', 'coords_x', 'coords_y', 'coords_z', 
                        # 'image_xaxis', 'image_yaxis',
                         'gain', 'indexes_x', 'indexes_y', 'pedestals', 'rms']
                coords_dict = {
                        'X': self.calibData.image_xaxis,
                        'Y': self.calibData.image_yaxis
                        }

            elif self.calibData.ndim == 2 and self.calibData.shape[0] > 0 and (self.calib is None or self.calibData.image_xaxis is None):
                attrs = []
                coords_shape = self.raw.shape
                coords_dict = {'X': np.arange(coords_shape[0]), 
                               'Y': np.arange(coords_shape[1])}

            elif self.calibData.ndim == 2 and  self.calibData.shape[0] > 0:
                raw_dims = (['xaxis', 'yaxis'], self.calibData.shape)
                attrs = ['areas', 'coords_x', 'coords_y', 'coords_z', 
                        # 'image_xaxis', 'image_yaxis',
                         'gain', 'indexes_x', 'indexes_y', 'pedestals', 'rms']
                coords_dict = {
                        'X': self.calibData.image_xaxis,
                        'Y': self.calibData.image_yaxis
                        }
            else:
                coords_dict = {}
                attrs = []

            for attr in attrs:
                val = getattr(self.calibData, attr)
                if val is not None:
                    coords_dict[attr] = (raw_dims[0], val)
        
        else:
            coords_dict = {}

        self._xarray_info['coords'].update(**coords_dict)

    def _add_xarray_evtData(self, attrs=[]):
        """
        Add evtData xarray information.
        """
        import numpy as np
        dims_dict = {}
        for attr in attrs:
            if attr in self.evtData._all_values:
                val = self.evtData._all_values.get(attr)
                npval = np.array(val)
                if npval.size > 1:
                    dims_dict[attr] = (['d{:}_{:}'.format(i,a) for i,a in enumerate(npval.shape)], npval.shape)
                else:
                    info = self.evtData._attr_info.get(attr)
                    if info:
                        xattrs = {a: b for a, b in info.items() if a in ['doc','unit']}
                    else:
                        xattrs = {}
                    dims_dict[attr] = ([], (), xattrs)
               
                self.add.property('.'.join(['evtData',attr]), attr)

        self._xarray_info['dims'].update(**dims_dict)

    @property
    def _attrs(self):
        """
        Attributes of psana.Detector functions if relevant, and otherwise
        attributes of raw psana event keys for the given detector.
        """
        if self._tabclass:
            tabclass = getattr(self, self._tabclass)
            if hasattr(tabclass, '_attrs'):
                return tabclass._attrs

        return []

    def next(self, **kwargs):
        """
        Return next event that contains the Detector in the event data.
        """
        in_keys = False
        while not in_keys:
            evt = self._ds.events.next(init=False, **kwargs) 
            in_keys = self._alias in evt._attrs

        evt._init()
        return getattr(evt, self._alias)
 
    def __iter__(self):
        return self

    def monitor(self, attrs=['raw','calib'], nevents=-1, sleep=0.2):
        """
        Monitor detector attributes continuously with show_info function.
        """ 
        ievent = nevents
        while ievent != 0:
            try:
                next(self)
                try:
                    print(self.show_info(attrs=attrs))
                except:
                    pass
                
                if ievent < nevents and sleep:
                    time.sleep(sleep)

                ievent -= 1

            except KeyboardInterrupt:
                ievent = 0

    def _show_user_info(self, **kwargs):
        """
        Show user defined information from AddOn methods.
        """
        message = Message(quiet=True, **kwargs)
#        if self._det_config.get('module') and self._det_config['module'].get('dict'):
#            message('-'*80)
#            message('Class Properties:')
#            message('-'*18)
#            for attr in self._det_config['module'].get('dict', []):
#                val = getattr(self, attr)
#                try:
#                    val = val()
#                except:
#                    pass
#                
#                strval = _repr_value(val)
#                fdict = {'attr': attr, 'str': strval, 'unit': '', 'doc': ''}
#                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))

        if self._det_config['projection']:
            message('-'*80)
            message('Projections:')
            message('-'*18)
            for attr, item in self._det_config['projection'].items():
                fdict = item.copy()
                val = getattr(self, attr)
                try:
                    val = val()
                except:
                    pass
 
                strval = _repr_value(val)
                fdict.update({'attr': attr, 'str': strval})
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
 
        if self._det_config['histogram']:
            message('-'*80)
            message('Histograms:')
            message('-'*18)
            for attr, item in self._det_config['histogram'].items():
                fdict = item.copy()
                val = getattr(self, attr)
                try:
                    val = val()
                except:
                    pass
 
                strval = _repr_value(val)
                fdict.update({'attr': attr, 'str': strval})
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
 
        if self._det_config['roi']:
            message('-'*80)
            message('Region of Interests (roi):')
            message('-'*18)
            for attr, item in self._det_config['roi'].items():
                fdict = item.copy()
                val = getattr(self, attr)
                strval = _repr_value(val)
                fdict.update({'attr': attr, 'str': strval})
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
 
        if self._det_config['count']:
            message('-'*80)
            message('Detector Counts:')
            message('-'*18)
            for attr, item in self._det_config['count'].items():
                fdict = item.copy()
                val = getattr(self, attr)
                strval = _repr_value(val)
                fdict.update({'attr': attr, 'str': strval})
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
 
        if self._det_config['parameter']:
            message('-'*80)
            message('User Defined Parameters:')
            message('-'*18)
            for attr, val in self._det_config['parameter'].items():
                strval = _repr_value(val)
                fdict = {'attr': attr, 'str': strval, 'unit': '', 'doc': ''}
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
 
        if self._det_config['peak']:
            message('-'*80)
            message('Peak Result:')
            message('-'*18)
            items = sorted(self._det_config['peak'].items(), key=operator.itemgetter(0))
            for attr, item in items:
                val = getattr(self, attr)
                if hasattr(val, 'mean') and val.size > 1:
                    strval = '<{:}>'.format(val.mean())
                else:
                    try:
                        strval = '{:10.5g}'.format(val)
                    except:
                        strval = str(val)

                doc = item.get('doc','')
                unit = item.get('unit','')
                fdict = {'attr': attr, 'str': strval, 'unit': unit, 'doc': doc}
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
 
        if self._det_config['property']:
            message('-'*80)
            message('User Defined Properties:')
            message('-'*18)
            for attr, func_name in self._det_config['property'].items():
                val = self._get_property(attr)
                strval = _repr_value(val)
                fdict = {'attr': attr, 'str': strval, 'unit': '', 'doc': ''}
                try:
                    adims = self._xarray_info['dims'].get(attr)
                    if adims and len(adims) == 3:
                        fdict.update(**adims[2])
                except:
                    print(adims)
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))

        return message

    def _get_info(self, attr):
        info = self.detector._attr_info.get(attr)
        if info:
            info.update({'type': 'detector'})
            return info

        for typ in ['count','histogram','parameter','peak','projection','property','stats']:
            info = self._det_config.get(typ).get(attr)
            if info:
                info.update({'type': typ})
                return info

        return {}


    def show_all(self, **kwargs):
        """
        Show detailed detector information for current event.
        """
        message = Message(quiet=True, **kwargs)
        message('-'*80)
        message(str(self))
        message('-'*80)
        message('Event Data:')
        message('-'*18)
        self.evtData.show_info(append=True)

        self._show_user_info(append=True)

        if self._tabclass == 'detector':
            message('-'*80)
            message('Processed Data:')
            message('-'*18)
            self.detector.show_info(append=True)
            if self._calib_class:
                message('-'*80)
                message('Calibration Data:')
                message('-'*18)
                self.calibData.show_info(append=True)

        if self.configData:
            message('-'*80)
            message('Configuration Data:')
            message('-'*18)
            self.configData.show_info(append=True)
        
        if self.epicsData:
            message('-'*80)
            message('Epics Data:')
            message('-'*18)
            self.epicsData.show_info(append=True)
               
        return message

    def show_info(self, attrs=None, **kwargs):
        """
        Show basic detector information, including from user defined AddOn methods, 
        for current event.
        """
        message = Message(quiet=True, **kwargs)
        message('-'*80)
        message(str(self))
        message('-'*80)
        getattr(self, self._tabclass).show_info(attrs=attrs, append=True)
        self._show_user_info(append=True)

        return message

    @property
    def sourceData(self):
        """
        Source information for detector.
        """
        return getattr(self._ds.configData.Sources, self._alias)

    @property
    def configData(self):
        """
        Configuration data for detector.
        """
        return getattr(self._ds.configData, self._alias)

    @property
    def evtData(self):
        """
        Tab accessible raw data from psana event keys.
        """
        if self._alias not in self._ds._current_evtData:
            self._ds._current_evtData.update({self._alias: 
                PsanaSrcData(self._ds._current_evt, self._srcstr, key_info=self._ds._evt_keys)})

        return self._ds._current_evtData.get(self._alias)

    @property
    def epicsData(self):
        """
        Epics information for current detector.
        """
        return getattr(self._ds.epicsData, self._alias)

    @property
    def detector(self):
        """
        Raw, calib and image data using psana.Detector class.
        Improved speed with data cashing when accessing the same object multiple times for
        the same event (e.g., multiple roi or other access for same data).
        
        See Also
        --------
        ImageData object
        ImageCalibData object
        WaveformData object
        WaveformCalibData object
        IpimbData object
        
        """
        if self._pydet:
            if self._alias not in self._ds._current_data:
                #opts = self._det_config.get('opts', {})
                #self._ds._current_data.update({self._alias: self._det_class(self._pydet, self._ds._current_evt, opts=opts)})
                self._ds._current_data.update({self._alias: self._det_class(self._pydet, self._ds._current_evt)})
             
            return self._ds._current_data.get(self._alias)

        else:
            return None

    @property
    def calibData(self):
        """
        Calibration data using psana.Detector class
        """
        if self._pydet:
            return self._calib_class(self._pydet, self._ds._current_evt)
        else:
            return None

    def set_cmpars(self, cmpars):
        """
        Set common mode.
        """
        if 'calib' not in self._det_config['opts']:
            self._det_config['opts']['calib'] = {}
        
        # reset current data
        self._ds._current_data = {}
        # save opts for calib object
        self._det_config['opts']['calib'].update({'cmpars': cmpars})
        self._pydet.calib(self._ds._current_evt, cmpars=cmpars)

    def _get_roi(self, attr):
        """
        Get roi from roi_name as defined by AddOn.
        """
        if attr in self._det_config['roi']:
            item = self._det_config['roi'][attr]
            #img = getattr(self, item['attr'])
            #img = self._getattr(item['attr'])
            img = getattr_complete(self, item['attr'])
            if img is not None:
                roi = item['roi']
                if not roi:
                    sensor = item.get('sensor')
                    if sensor is not None:
                        return img[sensor]
                    else:
                        return img
                elif roi and len(img.shape) == 1:
                    img = img[roi[0]:roi[1]]
                    if img.size == 1:
                        try:
                            return img[0]
                        except:
                            return img
                    else:
                        return img
                
                sensor = item.get('sensor')
                if sensor is not None:
                    img = img[sensor]
                    if img.size == 1:
                        try:
                            return img[0]
                        except:
                            return img

                return img[roi[0][0]:roi[0][1],roi[1][0]:roi[1][1]]
            else:
                return None
        else:
            return None

    @property
    def psplots(self):
        """
        To kill a plot that has been created with self.add.psplot, del self.psplots[name]
        and close the plot window.  
           
        If you only close the plot window it will automatically reopen on the next event.
        If the plot is not updating as expected, simply close the window and it will refresh
        on the next event.

        Example
        -------
        evt.DscCsPad.add.projection('calib','r',publish=True)
        del evt.DscCsPad.psplots['DscCsPad_calib_r']
        
        evt.DscCsPad.add.psplot('image')                
        del evt.DscCsPad.psplots['DscCsPad_image']

        """
        return self._det_config['psplot']

    def _get_count(self, attr):
        """
        Get counts from count_name as defined by AddOn.
        
        Parameters
        ----------
        attr : str
            Attribute name
        """
        if attr in self._det_config['count']:
            item = self._det_config['count'][attr]
            #img = getattr(self, item['attr'])
            #img = self._getattr(item['attr'])
            img = getattr_complete(self, item['attr'])
            gain = item.get('gain', 1.)
            limits = item.get('limits')
            if limits:
                return img[(img >= limits[0]) & (img < limits[1])].sum()*gain
            else:
                return img.sum()*gain

        else:
            return None

    def _get_histogram(self, attr):
        """
        Returns histogram as defined by AddOn.
        """
        import numpy as np
        if attr in self._det_config['histogram']:
            item = self._det_config['histogram'][attr]
            img = getattr_complete(self, item['attr'])
            gain = item.get('gain', 1.)
            hst, hbins = np.histogram(img*gain, item.get('bins'), 
                    weights=item.get('weights'), density=item.get('density'))

            return hst
        else:
            return None

    def _get_peak(self, attr):
        """
        Returns peak information as defined in AddOn class.
        
        Parameters
        ----------
        attr : str
            Attribute name
        """
        if attr in self._det_config['peak']:
            item = self._det_config['peak'][attr]
            wf = getattr(self, item['attr'])
            ichannel = item.get('ichannel')
            if ichannel is not None:
                wf = wf[ichannel]

            background = item.get('background')
            if background:
                wf -= background

            scale = item.get('scale')
            if scale:
                wf *= scale

            method = item.get('method')
            if method == 'waveform':
                return wf
            
            roi = item.get('roi')
            if roi:
                wf = wf[roi[0]:roi[1]]

            wf = list(wf)
            maxval = max(wf)
            if method == 'max':
                return maxval
            
            if max > item.get('threshold'):
                idx = wf.index(maxval)
                if roi:
                    idx += roi[0] 
            else:
                idx = 0        
    
            if method == 'index':
                return idx

            if method in ['time', 'pos']:
                xaxis = item.get('xaxis')
                if xaxis is not None:
                    return xaxis[idx]

        return None

    def _get_stats(self, attr, alias=None, stats=['mean', 'std', 'var', 'min', 'max']):
        """
        Get xarray.DataArray of Welford stats as defined by stats AddOn.
        """
        stat_info = self._det_config['stats'].get(attr)
        import numpy as np
        import xarray as xr
        
        # Need to fix handling of coords to be robust
        #coords = stat_info.get('coords').copy()
        coords = {}
        dims = stat_info.get('dims')
        eventCodes = stat_info['funcs'].keys()
        attrs = stat_info['attrs']
        steps = []
        for ec, item in stat_info['funcs'].items():
            for step, fec in item.items():
                steps.append(step)

        steps = list(set(steps))
        nsteps = len(steps)
        neventCodes = len(eventCodes)
        nstats = len(stats)
        aevents = np.zeros(shape=(nsteps,neventCodes))
        if not alias:
            alias = self._alias
        
        dim_names = ['_'.join([alias,a]) for a in dims[0]]
        dim_shape = list(dims[1])
        dim_names.insert(0,'codes')
        dim_shape.insert(0,neventCodes)
        dim_names.insert(0,'steps')
        dim_shape.insert(0,nsteps)
        dim_names.insert(0,'stat')
        dim_shape.insert(0,nstats)

        ddims = (dim_names, tuple(dim_shape))
        coords.update({'codes': eventCodes, 'stat': stats, 'steps': steps})
        asums = np.zeros(shape=dim_shape)
        for ec, item in stat_info['funcs'].items():
            for step, fec in item.items():
                istep = steps.index(step)
                iec = eventCodes.index(ec)
                aevents[istep, iec] = fec.n
                for istat, stat in enumerate(stats):
                    vals = getattr(fec, stat)()
                    if np.any(vals):
                        asums[istat, istep, iec] = vals         

        da = xr.DataArray(asums, dims=ddims[0], coords=coords, 
                attrs=attrs, name='_'.join([alias, attr]))
        da.coords['_'.join([alias,'events'])]= (['steps', 'codes'], aevents)

        return da

    def _get_projection(self, attr):
        """
        Get projection as defined by AddOn.
        
        Parameters
        ----------
        attr : str
            Attribute name
        """
        item = self._det_config['projection'].get(attr)
        if item is None:
            return None

        #img = getattr(self, item['attr'])
        img = getattr_complete(self, item['attr'])
        #img = self._getattr(item['attr'])
        if img is None:
            return None

        axis = item['axis']
        if axis in ['r', 'az']:
            # make a radial histogram
            coord_compressed = item['coord_compressed']
            mask = item['mask']
            bins = item['bins']
            norm = item['norm']
            # make masked image with mask compressed
            img_compressed = np.ma.masked_array(img, mask).compressed()
            hst, hbins  = np.histogram(coord_compressed, bins=bins, weights=img_compressed)
            return hst/norm
            
        else:
            # perform method on oposite axis where psana convention is images have coordinates (x, y)
            iaxis = {'x': 0, 'y': 1}.get(axis,0)
            method = item.get('method', 'sum')
            return getattr(img, method)(axis=iaxis)

    def _get_property(self, attr):
        """
        Get property as defined by AddOn.

        Parameters
        ----------
        attr : str
            Attribute name
        """
        import six
        func_name = self._det_config['property'].get(attr)
        if hasattr(func_name,'__call__'):
            return func_name(self)
        elif isinstance(func_name, six.string_types)  and hasattr(self, func_name):
            return getattr_complete(self, func_name)
        else:
            return None

    def __str__(self):
        return '{:} {:}'.format(self._alias, str(self._ds.events.current))

    def __repr__(self):
        return '< {:}: {:} >'.format(self.__class__.__name__, str(self))
   
    def __getattr__(self, attr):
        if attr in self._attrs:
            return getattr(getattr(self, self._tabclass), attr)

        if attr in self._det_config['parameter']:
            return self._det_config['parameter'].get(attr)

        if attr in self._det_config['property']:
            return self._get_property(attr)
            
        if attr in self._det_config['count']:
            return self._get_count(attr)

        if attr in self._det_config['histogram']:
            return self._get_histogram(attr)

        if attr in self._det_config['roi']:
            return self._get_roi(attr)
       
        if attr in self._det_config['peak']:
            return self._get_peak(attr)
       
        if attr in self._det_config['projection']:
            return self._get_projection(attr)

        if attr in self._det_config['stats']:
            return self._get_stats(attr)

        if attr in self._ds.events.current._event_attrs:
            return getattr(self._ds.events.current, attr)

    def __dir__(self):
        all_attrs =  set(self._attrs+
                         self._det_config['parameter'].keys() +
                         self._det_config['property'].keys() +
                         self._det_config['roi'].keys() + 
                         self._det_config['count'].keys() + 
                         self._det_config['histogram'].keys() + 
                         self._det_config['peak'].keys() + 
                         self._det_config['projection'].keys() + 
                         self._det_config['stats'].keys() + 
                         self._det_config['module'].get('dict', []) +
                         self._ds.events.current._event_attrs +
                         self.__dict__.keys() + dir(Detector))
        
        return list(sorted(all_attrs))


class AddOn(object):
    """
    Collection of methods to add parameters, properties, and reduction/proccesing of data 
    with roi, projection, histogram methods.
    """

    _plugins = {}
    _init_attrs = ['_ds', '_alias']

    _attrs = ['parameter', 'property', 'peak', 
              'histogram', 'roi', 'projection', 'count', 'stats']

    def __init__(self, ds, alias):
        self._ds = ds
        self._alias = alias

    @property
    def _det_config(self):
        return self._ds._device_sets.get(self._alias)

    @property
    def _evt(self):
        return self._ds.events.current

    @property
    def _det(self):
        return getattr(self._evt, self._alias)

    def module(self, module=None, **kwargs):
        """Add detector module

        Parameters
        ----------
        module : str
            Name of python module that contains a user defined PyDataSource.Detector class
            with the same name as the module.  e.g., 'acqiris' loads 'acqiris.py' file 
            which must have class Acqiris(PyDatasource.Detector)
        path : str
            Name of path for python module (default is the path of PyDataSource)
        desc : str
            Description 

        Example
        -------
        >>> evt.acqiris.add.module('acqiris')

        See Also
        --------
        DataSource.add_detector

        """
        if module:
            self._ds.add_detector(self._alias, module=module, **kwargs)

    def parameter(self, **kwargs):
        """
        Add parameters as keword arguments.  
        Will be perserved if module is added or DataSource is reloaded
        and is saved along with other AddOn information with save_config.

        Examples
        --------
        
        >>> evt.acqiris.add.parameter(foo='bar', backlevel=1.)
        
        or alternatively

        >>> params = {foo='bar', backlevel=1.}
        >>> evt.acqiris.add.parameter(**params)

        """
        for param, value in kwargs.items():
            self._det_config['parameter'][param] = value 

    #def property(self, func_name, *args, **kwargs):
    def property(self, func_name, attr=None, **kwargs):
        """
        Add a property that operates on this detector object.
        
        The result will be added as an attribute to the detecor with the
        name of the function unless attr is provided.
       
        Examples
        --------
        Import or create a method acting on self as the Detector object and add it as a property
        
        >>> def myfunc(self):       
                return self.ebeamL3Energy
        
        >>> evt.EBeam.add.property(myfunc, 'energy')

        Or alternatively using lambda:
        
        >>> evt.EBeam.add.property(lambda self: self.ebeamL3Energy, 'energy')
        
        """
        if not attr:
            attr = func_name.__name__
        
        self._det_config['property'][attr] = func_name 
        
        doc = kwargs.get('doc', '')
        unit= kwargs.get('unit', '')
        xattrs = {}
        xattrs.update({'doc': doc, 'unit': unit})
        self._det_config['xarray']['dims'].update(
                    {attr: ([], (), xattrs)})

    def projection(self, attr=None, axis='x', 
            roi=None, sensor=None, name=None, 
            axis_name=None, method=None,  
            coords_x=None,coords_y=None,
            mask=None, bins=None, bin_size=None, rmin=None, rmax=None, 
            unit='ADU', doc='', publish=False,
            **kwargs):
        """
        Make a projection along an axis.
        
        Options implemented for x/y projections (i.e., axis='x' or axis='y')
        But radial (axis='r') and azimuth (axis='az') projections are valid 
        in limited circumstances.
        
        Parameters
        ----------
        attr : str
            Name of data object on which to make projection
        axis : str
            Projection axis
                'x': x-axis 
                'y': y-axis
                'r': radial
                'az': azimuth
                True: both 'x' and 'y' projections
        name : str
            Name of projection object [Default is to append '_'+axis and sequential number
            to input attr name.
        axis_name : str
            name of axis [default = axis+name]
        coords_x : array
            Array of X coordinates for input data object defined by attr
        coords_y : array
            Array of X coordinates for input data object defined by attr
        method : str
            by default the method is 'sum' but other numpy methods are
            also possible (e.g., 'mean', 'std', 'min', 'max', 'var').
        sensor : int, optional
            First array element for 3 dim objects
            If no roi given, then roi created for entire sensor
        roi : tuple, optional
            Region of interest (roi) parameters passed to roi method.
            roi method then acts on this reduced data
        mask : str, optional
            Name of mask to apply before making projection.
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width
            bins in the given range (10, by default). If `bins` is a sequence,
            it defines the bin edges, including the rightmost edge, allowing
            for non-uniform bin widths. (used with np.histogram)
        bin_size : float, optional
            Size of bins for polar coordinate projections.
            Default = pixelsize
        rmin : float, optional
            Minimum radius for az projection only 
            in units of calibData coords (e.g., calibData.coords_x, typically um)
        rmax : float, optional
            Maximum radius for az projection only
            in units of calibData coords (e.g., calibData.coords_x, typically um)
        doc : str, optional
            Doc string for resulting projection data
        unit : str, optional
            Units of roi data [Default assumes same as attr data]
        publish : bool, optional
            Make psplot of roi data
 
        """
        import numpy as np
        _methods = ['sum', 'mean', 'std', 'min', 'max', 'var']
        _cartesian_axes = ['x', 'y']
        _polar_axes = ['r', 'az']

        axis = axis.lower()
        xunit = ''

        if not attr:
            if axis in _polar_axes:
                attr = 'calib'

            else:
                if self._det.image is not None:
                    attr = 'image'
                elif self._det.calib is not None:
                    attr = 'calib'
                elif self._det.raw is not None:
                    attr = 'raw'
                elif self._det.evtData.data16 is not None:
                    attr = 'evtData.data16'
                else:
                    attr = 'evtData.data8'

        try:
            img = getattr_complete(self._det,attr)
        except:
            return
 
        if sensor is not None:
            img = img[sensor]
            if not roi:
                roi = ((0,img.shape[0]), (0,img.shape[1]))
            print((sensor, roi))

        if axis in _polar_axes:
            if not method:
                method = 'norm'
            elif method not in ['norm']:
                print('ERROR: {:} is not a valid method'.format(method))
                print(' - only norm method is valid for radial projections')

        elif axis in _cartesian_axes:
            if not method:
                method = 'sum'

            elif method not in _methods:
                print('ERROR: {:} is not a valid method'.format(method))
                print(' - method must be in {:}'.format(_methods))
            
            if len(img.shape) != 2:
                print('{:} is not a 2D image -- cannot make projection'.format(attr))
                return
        
        else:
            print('ERROR:  {:} is not a valid axis'.format(axis))
            print(' - Valid Cartesian Options = {:}'.format(_cartesian_axes))
            print(' - Valid Polar Options = {:}'.format(_polar_axes))
            return 

        if roi:
            if name:
                roi_name = 'roi_'+name
            else:
                roi_name = None

            roi_name = self.roi(attr, roi=roi, sensor=sensor, unit=unit, doc=doc, name=roi_name, **kwargs)
            xattrs = self._det_config['xarray']['dims'][roi_name][2]
            if not doc:
                if sensor:
                    doc = '{:}-axis projection of {:} within roi={:}'.format(axis, attr, roi)
                else:
                    doc = '{:}-axis projection of {:} sensor {:} within roi={:}'.format(axis, attr, sensor, roi)

        else:
            roi_name = attr
            if doc == '':
                doc = "{:}-axis projection of {:} data".format(axis, attr)
            
        if not name:
            name = roi_name+'_'+axis
       
        if not axis_name:
            axis_name = axis+roi_name
 
        calibData = self._det.calibData
        if axis in _polar_axes:
            animg = img.shape
            nx = animg[1]
            ny = animg[0]
            if coords_x is not None:
                if len(coords_x.shape) == 1:
                    coords_x = np.ones((ny,nx)) * coords_x

            elif calibData.coords_x is not None:
                coords_x = calibData.coords_x
            else:
                coords_x = np.ones((ny,nx)) * np.arange(nx)-(nx-1)/2.
                
            if coords_y is not None:
                if len(coords_y.shape) == 1:
                    coords_y = (np.ones((nx,ny)) * coords_y).T
            
            elif calibData.coords_y is not None:
                coords_y = calibData.coords_y
            else:
                coords_y = (np.ones((nx,ny)) * np.arange(ny)-(ny-1)/2.).T

            if not mask:
                mask = ~np.array(calibData.mask(), dtype=bool) 
                if mask.shape != img.shape:
                    mask = np.zeros_like(img, dtype=bool)

            coords_r = np.sqrt(coords_y**2+coords_x**2)
            coords_az = np.degrees(np.arctan2(coords_y, coords_x))

            if axis == 'r':
                if not xunit:
                    xunit = 'um'
                coord_hist = np.ma.masked_array(coords_r, mask)
                if not bins:
                    if not bin_size:
                        bin_size = calibData.pixel_size
                    if not bin_size:
                        bin_size = 1

                    bins = np.arange(coord_hist.compressed().min()+bin_size*10, 
                                     coord_hist.compressed().max()-bin_size*10., 
                                     bin_size)
     
            else:
                if rmin and rmax and np.array(rmin).size > 1 and np.array(rmax).size > 1:
                    # in future
                    print('Vectors for rmin and rmax not yet supported')
                    return

                else:
                    if rmin:
                        # add logical or for rmin
                        mask |= np.ma.masked_less(coords_r, rmin).mask

                    if rmax:
                        # add logical or for rmax
                        mask |= np.ma.masked_greater_equal(coords_r, rmax).mask

                    coord_hist = np.ma.masked_array(coords_az, mask)
                
                if not bins:
                    if not bin_size:
                        bin_size = 2.

                    bins = np.arange(-180., 180., bin_size)

            coord_compressed = coord_hist.compressed()
            norm, hbins = np.histogram(coord_compressed, bins=bins)
            if method != 'norm':
                norm = 1.
            
            projaxis = (hbins[1:]+hbins[0:-1])/2.
            self._det_config['xarray']['coords'].update({axis_name: projaxis})

            self._det_config['xarray']['dims'].update(
                    {name: ([axis_name], (projaxis.size))})

            self._det_config['projection'].update(
                    {name: {'attr': roi_name, 
                            'axis': axis, 
                            'method': method,
                            'axis_name': axis_name,
                            'coord_compressed': coord_compressed,
                            'mask': mask,
                            'bins': bins,
                            'rmin': rmin,
                            'rmax': rmax,
                            'norm': norm,
                            'xunit': xunit,
                            'unit': unit,
                            'doc': doc,
                            'xdata': projaxis}})

        else:
            projaxis = self._det_config['xarray']['coords'].get(axis_name) 
           
            if projaxis is None:
                # Need to add in auto axis for Data Arrays
                if not xunit:
                    xunit = 'pixel'
                
                iaxis = {'x': 1, 'y': 0}.get(axis)
                projaxis = np.arange(img.shape[iaxis])    
                self._det_config['xarray']['coords'].update({axis_name: projaxis})


            self._det_config['xarray']['dims'].update(
                    {name: ([axis_name], (projaxis.size))})

            self._det_config['projection'].update(
                    {name: {'attr': roi_name, 
                            'axis': axis, 
                            'xunit': xunit, 
                            'method': method,
                            'axis_name': axis_name, 
                            'unit': unit,
                            'doc': doc,
                            'xdata': projaxis}})

        if publish:

            self.psplot(name)

    def roi(self, attr=None, sensor=None, roi=None, name=None, 
                xaxis=None, yaxis=None, 
                xaxis_name=None, yaxis_name=None,
                doc='', unit='ADU', publish=False, projection=None, 
                graphical=None, quiet=True, **kwargs):       
        """
        Make roi for given attribute, by default this is given the name img.

        Parameters
        ----------
        attr : str
            Name of data object on which to make roi
        sensor : int
            First array element for 3 dim objects
        roi : tuple
            For 1 dim objects, roi is a tuple (xstart, xend)
            For 2 dim objects, roi is a tuple ((ystart, yend), (xstart, xend))
            For 3 dim objects (e.g., cspad raw, calib), use sensor keyword to 
            sepcify the sensor and roi as ((ystart, yend), (xstart, xend)).
        name : str
            Name of roi object [Default is to append '_roi' and sequential number
            to input attr name.
        xaxis: array_like
            X-axis coordinate of roi
        yaxis: array_like
            Y-axis coordinate of roi
        xaxis_name : str
            Name of xaxis (used for xarray Datasets - name and axis values must be unique)
        yaxis_name : str
            Name of yaxis (used for xarray Datasets - name and axis values must be unique)
        doc : str
            Doc string for resulting roi data
        unit : str
            Units of roi data [Default assumes same as attr data]
        publish : bool
            Make psplot of roi data
        projection : bool
            Make projection(s) of data.  See projection method. 
        graphical : bool
            Select roi interactively

        """
        import matplotlib.pyplot as plt
        import numpy as np
        if not attr:
            if sensor is not None:
                attr = 'calib'
            else:
                attr = 'image'

        try:
            img = self._getattr(attr)
        except:
            if not quiet:
                print('Not valid roi {:} {:}'.format(self._alias, attr))
            return
        
        if sensor is not None:
            img = img[sensor]

        if not roi and sensor is None and graphical is not False:
            try:
                plotMax = np.percentile(img, 99.5)
                plotMin = np.percentile(img, 5)
                print(('using the 5/99.5% as plot min/max: (',plotMin,',',plotMax,')'))
                fig=plt.figure(figsize=(20,10))
                gs=plt.matplotlib.gridspec.GridSpec(1,2,width_ratios=[2,1])
                plt.subplot(gs[0]).imshow(img,clim=[plotMin,plotMax],interpolation='None')

                print('Select two points to form ROI to zoom in on target location.')
                p = np.array(ginput(2))
                roi = ([int(p[:,1].min()),int(p[:,1].max())],
                       [int(p[:,0].min()),int(p[:,0].max())])
                print(('Selected ROI [y, x] =', roi))
            except:
                print('Cannot get roi')
                return None

        #if len(roi) == 3:
        #    xroi = roi[2]
        #    yroi = roi[1]
        #else:
        if not name:
            if not roi:
                name = 'sensor'+str(sensor)
            else:
                nroi = len(self._det_config['roi'])
                if nroi == 0:
                    name = 'roi'
                else:
                    name = 'roi'+str(nroi+1)

        if sensor is not None and not roi:
            if doc == '':
                doc = "{:} ROI of {:} data".format(name, attr)

            xattrs = {}
            xattrs.update({'doc': doc, 'unit': unit})
            xattrs.update({'sensor': sensor})
           
            if img.size == 1:
                self._det_config['xarray']['dims'].update(
                        {name: ([], (), xattrs)})

            elif len(img.shape) == 1:
                self._det_config['xarray']['dims'].update(
                        {name: (['x'+name], img.size, xattrs)})
            
            else:
                self._det_config['xarray']['dims'].update(
                        {name: (['y'+name, 'x'+name], img.size, xattrs)})

            self._det_config['roi'].update({name: {'attr': attr, 
                                                   'roi': None,
                                                   'sensor': sensor,
                                                   'xaxis': None,
                                                   'doc': doc,
                                                   'unit': unit}})
 
        elif len(img.shape) == 1:
            if not xaxis_name:
                xaxis_name = 'x'+name
            xroi = roi
            if xroi[0] > xroi[1]:
                raise Exception('Invalid roi {:}'.format(roi))
            xroi = (max([xroi[0], 0]), min([xroi[1], img.shape[0]]))

            if not xaxis or len(xaxis) != xroi[1]-xroi[0]:
                xaxis = np.arange(xroi[0],xroi[1])    
            
            if doc == '':
                doc = "{:} ROI of {:} data".format(name, attr)

            xattrs = {}
            xattrs.update({'doc': doc, 'unit': unit, 'roi': xroi})
            if sensor is not None:
                xattrs.update({'sensor': sensor})
            
            self._det_config['xarray']['coords'].update({xaxis_name: xaxis})
            self._det_config['xarray']['dims'].update(
                    {name: ([xaxis_name], (xaxis.size), xattrs)})

            self._det_config['roi'].update({name: {'attr': attr, 
                                                   'roi': xroi,
                                                   'sensor': sensor,
                                                   'xaxis': xaxis,
                                                   'doc': doc,
                                                   'unit': unit}})
 
        else:
            xroi = roi[1]
            yroi = roi[0]
            if xroi[0] > xroi[1] or yroi[0] > yroi[1]:
                raise Exception('Invalid roi {:}'.format(roi))
            
            xroi = (max([xroi[0], 0]), min([xroi[1], img.shape[1]]))
            yroi = (max([yroi[0], 0]), min([yroi[1], img.shape[0]]))

            if not xaxis_name:
                xaxis_name = 'x'+name
            if not xaxis or len(xaxis) != xroi[1]-xroi[0]:
                xaxis = np.arange(xroi[0],xroi[1])    
            
            if not yaxis_name:
                yaxis_name = 'y'+name
            if not yaxis or len(yaxis) != yroi[1]-yroi[0]:
                yaxis = np.arange(yroi[0],yroi[1])    
            
            if doc == '':
                doc = "{:} ROI of {:} data".format(name, attr)

            xattrs = {}
            xattrs.update({'doc': doc, 'unit': unit, 'roi': (yroi, xroi)})
            if sensor is not None:
                xattrs.update({'sensor': sensor})
            
            self._det_config['xarray']['coords'].update({xaxis_name: xaxis, 
                                                yaxis_name: yaxis})
            self._det_config['xarray']['dims'].update(
                    {name: ([yaxis_name, xaxis_name], (yaxis.size, xaxis.size), xattrs)})

            self._det_config['roi'].update({name: {'attr': attr, 
                                                   'roi': (yroi, xroi),
                                                   'sensor': sensor,
                                                   'xaxis': xaxis,
                                                   'yaxis': yaxis, 
                                                   'doc': doc,
                                                   'unit': unit}})

        if projection:
            if projection in [True, 'x']:
                self.projection(name, axis='x', publish=publish, **kwargs)
            if projection in [True, 'y']:
                self.projection(name, axis='y', publish=publish, **kwargs)
        
        elif publish:
            self.psplot(name)

        return name

    def stats(self, attr=None, doc=None, 
            name=None, eventCodes=None, attrs={}, **kwargs):
        """
        Calculate running statistics (mean, std, min, max, count) of detector data attribute 
        during event iteration using Welford algorithm.

        Parameters
        ----------
        attr : str
            Name of data object in detector object on which to act
        """
        if not name:
            name = attr+'_stats'

        try:
            img = self._getattr(attr)
        except:
            print('Stats add Not valid for {:} {:}'.format(self._alias, attr))
            return False

        if not eventCodes:
            # try first source eventCode to make sure get 140 instead of 40 
            # for example in CsPad where 40 is given in configData.eventCode
            srcstr = self._det._srcstr
            code0 = self._det._ds.configData._sources.get(srcstr,{}).get('eventCode')
            #print attr, srcstr, code0
            
            if not code0:
                # some devices like Timetool with roi do not give eventCode in configData
                code0 = self._det.configData.eventCode
            
            if not code0:
                code0 = 40
                if code0 not in self._ds.configData._eventcodes:
                    code0 = self._ds.configData._eventcodes.keys()

            eventCodes = code0

        if not isinstance(eventCodes, list):
            eventCodes = [eventCodes]
        # Need to update dims always
        if True or not self._det_config['xarray'].get('coords'):
            self._det._update_xarray_info()

        dims = self._det_config['xarray'].get('dims', {}).get(attr, {})
        if not dims:
            dims = self._det_config['xarray'].get('dims', {}).get('raw', {})
        try:
            if dims[1] != img.shape:
                print(('Fixing stat dims for ', attr, img.shape, dims))
                dims[1] = img.shape
        except:
            print(('Cannot update stat dims for ', attr, dims))

        all_coords = self._det_config['xarray'].get('coords', {})
        #coords = {dim: coord for dim, coord in all_coords.items() if dim in dims[0]}
        #coords = {dim: coord for dim, coord in all_coords.items() if not isinstance(coord[0], list) or len(set(coord[0]) & set(dims[0])) > 0}
        #coords = {dim: coord for dim, coord in all_coords.items() if coord[0] == dims[0]}
        coords = {}
        #attrs = self._det_config['xarray'].get('attrs', {})
        stat_attrs = self._det._source_info
        stat_attrs.update(**attrs)

#        try:
#            dims = self._det_config['xarray']['dims'].get(attr)
#            coords = self._det_config['xarray']['coords'].get(attr)
#        except:
#            dims = []
#            coords = {}
#
#        dims.insert(0, 'codes')
#        dims.insert(0, 'steps')
        
        self._det_config['stats'].update({name: 
                {
                'attr': attr,
                'dims': dims,
                'coords': coords,
                'funcs': {},
                'eventCodes': eventCodes,
                'attrs': stat_attrs,
                'stats': ['mean', 'std', 'min', 'max'],
                }})

        return True


    def count(self, attr=None, gain=None, unit=None, doc=None, 
            roi=None, name=None, limits=None, **kwargs):
        """
        Count (i.e., sum) of detector within optional roi.

        Parameters
        ----------
        attr : str
            Name of data object in detector object on which to act
        gain: float
            optional converion from ADU to for example X-rays.
        limits: tuple
            (low, high) values to be counted.
        doc : str
            Doc string for resulting roi data
        unit : str
            Units of histogramed data.  Default ADU (times gain factor if supplied) 
        roi : tuple, optional
            Region of interest (roi) parameters passed to roi method.
            roi method then acts on this reduced data
        
        """
        if gain:
            if not unit:
                unit = 'ADUx{:.2g}'.format(gain)

        else:
            gain = 1
            if not unit:
                unit = 'ADU'

        if roi:
            if name:
                roi_name = 'roi_'+name
            else:
                roi_name = None

            roi_name = self.roi(attr, roi=roi, unit=unit, doc=doc, name=roi_name, **kwargs)
            xattrs = self._det_config['xarray']['dims'][roi_name][2]
            if not doc:
                doc = 'Sum of {:} within roi={:}'.format(attr, roi)
       
        else:
            roi_name = attr
            xattrs = {}
            if not doc:
                doc = 'Sum of {:}'.format(attr)


        if not name:
            ncount = len(self._det_config['count'])
            if ncount == 0:
                name = roi_name+'_count'
            
            else:
                name = roi_name+'_count'+str(ncount+1)


        self._det_config['count'].update({name: {'attr': roi_name, 
                                                 'gain': gain, 
                                                 'unit': unit,
                                                 'limits': limits,
                                                 'doc': doc}})

        xattrs.update({'doc': doc, 'unit': unit})
        self._det_config['xarray']['dims'].update(
                {name: ([], (), xattrs)})

    def histogram(self, attr=None, bins=None, gain=None, 
            unit=None, doc=None, roi=None, name=None, 
            weights=None, density=None, publish=None,
            **kwargs):
        """
        Make a histogram.
        
        Parameters
        ----------
        
        attr : str
            Name of data object in detector object on which to act
        bins : int or sequence of scalars, optional
            If `bins` is an int, it defines the number of equal-width
            bins in the given range (10, by default). If `bins` is a sequence,
            it defines the bin edges, including the rightmost edge, allowing
            for non-uniform bin widths. (used with np.histogram)
        gain : float, optional
            converion from ADU to for example X-rays.
        doc : str, optional
            Doc string for resulting roi data
        unit : str, optional
            Units of histogramed data.  Default ADU (times gain factor if supplied) 
        roi : tuple, optional
            Region of interest (roi) parameters passed to roi method.
            Histogram acts on this reduced data
        name : str, optional
            Name of histogram data object
        weights : array_like, optional
            An array of weights, of the same shape as `a`.  Each value in `a`
            only contributes its associated weight towards the bin count
            (instead of 1).  If `normed` is True, the weights are normalized,
            so that the integral of the density over the range remains 1
            (used with np.histogram)
        density : bool, optional
            If False, the result will contain the number of samples
            in each bin.  If True, the result is the value of the
            probability *density* function at the bin, normalized such that
            the *integral* over the range is 1. Note that the sum of the
            histogram values will not be equal to 1 unless bins of unity
            width are chosen; it is not a probability *mass* function.
            Overrides the `normed` keyword if given.
            (used with np.histogram)
        publish : bool
            Make psplot of histogram
        
        See Also
        --------
        np.histogram

        """
#        range : (float, float), optional
#            The lower and upper range of the bins.  If not provided, range
#            is simply ``(a.min(), a.max())``.  Values outside the range are
#            An array of weights, of the same shape as `a`.  Each value in `a`
#            only contributes its associated weight towards the bin count
#            (instead of 1).  If `normed` is True, the weights are normalized,
#            so that the integral of the density over the range remains 1
#            (used with np.histogram)

        import numpy as np
        if gain:
            if not unit:
                unit = 'ADUx{:}'.format(gain)

        else:
            gain = 1
            if not unit:
                unit = 'ADU'

        if roi:
            if name:
                roi_name = 'roi_'+name
            else:
                roi_name = None

            roi_name = self.roi(attr, roi=roi, unit=unit, doc=doc, name=roi_name, **kwargs)
            xattrs = self._det_config['xarray']['dims'][roi_name][2]
            if not doc:
                doc = 'Histogram of {:} within roi={:}'.format(attr, roi)
       
        else:
            roi_name = attr
            xattrs = {}
            if not doc:
                doc = 'Histogram of {:}'.format(attr)


        if not name:
            ncount = len(self._det_config['histogram'])
            if ncount == 0:
                name = roi_name+'_hist'
            
            else:
                name = roi_name+'_hist'+str(ncount+1)
            
        if not hasattr(bins, 'shape') or bins.shape[0] < 3:
            img = self._getattr(roi_name)
            hst, bins = np.histogram(img, bins, 
                    weights=weights, density=density)

        xaxis = (bins[1:]+bins[:-1])/2.

        self._det_config['histogram'].update({name: {'attr': roi_name, 
                                                     'gain': gain, 
                                                     'unit': unit,
                                                     'bins': bins,
                                                     'xaxis': xaxis, 
                                                     'weights': weights,
                                                     'density': density,
                                                     'doc': doc}})

        xattrs.update({'doc': doc, 'unit': unit})
        self._det_config['xarray']['coords'].update({name+'_xaxis': xaxis})
        self._det_config['xarray']['dims'].update(
                {name: ([name+'_xaxis'], (xaxis.shape), xattrs)})

        if publish:
            xlabel = '[{:}]'.format(unit)
            self.psplot(name, 
                        xdata=xaxis, 
                        xlabel=xlabel)


    def peak(self, attr=None, ichannel=None, name=None,
            xaxis=None, roi=None, scale=1, baseline=None,
            methods=None,
            threshold=None, 
            docs={}, units={}):
        """
        Simple 1D peak locator.
           
        Parameters
        ----------
        
        attr : str
            Name of data object on which to find peaks
        ichannel : int
            If not present, peak locators will be created for all chanels
            [Note: acqiris waveform default names start with ch1 following AMI convention]  
        name : str
            Name of peak object [Default is to append '_peak' and sequential number
            to input attr name.
        xaxis: array_like
            X-axis coordinate of data object
            Default is to use wftime for attr='waveform'
        threshold : float
            Peak threshold.  default threshold is 1.5% of full range for Acqiris
        baseline : float
            Baseline level to subtract
        scale : float
            Peak scale factor
        method : str
            Peak finding method options include:
            'index': index at max channel.
            'max': maximum value.
        roi : tuple, optional
            Region of interest 
        docs : dict
            Doc strings for resulting peaks
        units : dict
            Units of peaks data [Default assumes same as attr data]
 
        """
        
        if not attr:
            if self._det._pydet_name == 'WFDetector':
                if self._det.__class__.__name__ == 'Detector':
                    self._det.add.module('acqiris')

                attr = 'waveform'
                if ichannel is None and not name:
                    for i in range(self._det.waveform.shape[0]):
                        print(('adding ichannel', i))
                        self.peak(ichannel=i, xaxis=xaxis, roi=roi, scale=scale, 
                                  theshold=threshold, baseline=baseline, methods=None,
                                  docs=docs,units=units)
                    
                    return

        if not methods:
            methods = {'index': 'index', 'max': 'max'}
            if attr == 'waveform':
                methods.update({'time': 'time', 'waveform': 'waveform'})
            else:
                methods.update({'pos': 'pos'})

        attr_info = self._det._get_info(attr)

        if not units:
            if attr == 'waveform':
                units = {'time': 's', 'index': 'channel', 'max': 'V', 'waveform': 'V'}
            else:
                units = {'max': '', 'index': ''}
                #units = {'max': attr_info.get('unit'), 'index': 'bin'}

        if not name:
            if attr == 'waveform':
                name = 'peak'
            else:
                name = attr+'_peak'
            
            if ichannel is not None:
                name += '_ch{:}'.format(ichannel+1)

        if not xaxis:
            if attr == 'waveform' and ichannel is not None:
                xaxis = self._det.wftime[ichannel]
            elif 'xaxis' in attr_info:
                xaxis = attr_info.get('xaxis')
            elif 'bins' in attr_info:
                xaxis = attr_info.get('bins')
            else:
                try:
                    if attr in self._det_config['xarray']['dims']:
                        adims = self._det_config['xarray']['dims'][attr][0]
                        if len(adims) == 1:
                            xaxis = self._det_config['xarray']['coords'].get(adims[0])
                except:
                    pass

        if not threshold:
            if self._det._pydet.dettype == 16:
                threshold = self._det.configData.vert.fullScale[ichannel]*0.015
            else:
                threshold = 0.

        for mname, method in methods.items():
            # make entry for each peak method
            pname = '_'.join([name, mname])
            doc = docs.get(mname, '')
            if not doc:
                doc = '{:} {:} {:}'.format(self._alias, name, method)
                if method not in ['waveform']:
                    doc+=' of peak'

            unit = units.get(mname, '')

            self._det_config['peak'].update({pname:   
                    {'attr': attr,
                     'ichannel': ichannel,
                     'roi': roi,
                     'baseline': baseline,
                     'threshold': threshold, 
                     'xaxis': xaxis,
                     'scale': scale,
                     'doc': doc,
                     'unit': unit,
                     'method': method,
                     }})

            if method in ['waveform']:
                axis_name = '_'.join([name,'t'])
                projaxis = self._det.wftime[ichannel]
                self._det_config['xarray']['coords'].update({axis_name: projaxis})
                self._det_config['xarray']['dims'].update(
                    {pname: ([axis_name], (projaxis.size))})
                
            else:
                self._det_config['xarray']['dims'].update(
                        {pname: ([], (), {'doc': doc, 'unit': unit})}
                        )



#    def mask(self, attr):
#        img = self._getattr(attr)
#        plotMax = np.percentile(img, 99.5)
#        plotMin = np.percentile(img, 5)
#        print 'using the 5/99.5% as plot min/max: (',plotMin,',',plotMax,')'
#
##        image = self.__dict__[detname].image(self.lda.run, img)
##        det = self.__dict__[detname]
##        x = self.__dict__[detname+'_x']
##        y = self.__dict__[detname+'_y']
##        iX = self.__dict__[detname+'_iX']
##        iY = self.__dict__[detname+'_iY']
#        x = self._getattr()
#        extent=[x.min(), x.max(), y.min(), y.max()]
#
#        fig=plt.figure(figsize=(10,6))
#        from matplotlib import gridspec
#        gs=gridspec.GridSpec(1,2,width_ratios=[2,1])
#        
#        mask=None
#        mask_r_nda=None
#        select=True
#        while select:
#            plt.subplot(gs[0]).imshow(image,clim=[plotMin,plotMax],interpolation='None')
#
#            shape = raw_input("rectangle(r), circle(c) or polygon(p)?:\n")
#            if shape=='r':
#                print 'select two corners: '
#                p =np.array(ginput(2))
#                mask_roi=np.zeros_like(image)
#                mask_roi[p[:,1].min():p[:,1].max(),p[:,0].min():p[:,0].max()]=1
#                mask_r_nda = np.array( [mask_roi[ix, iy] for ix, iy in zip(iX,iY)] )
#                plt.subplot(gs[1]).imshow(det.image(self.lda.run,mask_r_nda))
#                print 'mask from rectangle (shape):',mask_r_nda.shape
#                if raw_input("Done?\n") in ["y","Y"]:
#                    select = False
#            elif shape=='c':
#                plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
#                if raw_input("Select center by mouse?\n") in ["y","Y"]:
#                    c=ginput(1)
#                    cx=c[0][0];cy=c[0][1]
#                    print 'center: ',cx,' ',cy
#                else:
#                    ctot = raw_input("center (x y)?\n")
#                    c = ctot.split(' ');cx=float(c[0]);cy=float(c[1]);
#                if raw_input("Select outer radius by mouse?\n") in ["y","Y"]: 
#                    r=ginput(1)
#                    rox=r[0][0];roy=r[0][1]
#                    ro=np.sqrt((rox-cx)**2+(roy-cy)**2)
#                    if raw_input("Select inner radius by mouse?\n") in ["y","Y"]:
#                        r=ginput(1)
#                        rix=r[0][0];riy=r[0][1]
#                        ri=np.sqrt((rix-cx)**2+(riy-cy)**2)
#                    else:
#                        ri=0
#                    print 'radii: ',ro,' ',ri
#                else:
#                    rtot = raw_input("radii (r_outer r_inner)?\n")
#                    r = rtot.split(' ');ro=float(r[0]);ri=max(0.,float(r[1]));        
#                mask_router_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ro**2 for ix, iy in zip(x,y)] )
#                mask_rinner_nda = np.array( [(ix-cx)**2+(iy-cy)**2<ri**2 for ix, iy in zip(x,y)] )
#                mask_r_nda = mask_router_nda&~mask_rinner_nda
#                print 'mask from circle (shape):',mask_r_nda.shape
#                plt.subplot(gs[1]).imshow(det.image(self.lda.run,mask_r_nda))
#                if raw_input("Done?\n") in ["y","Y"]:
#                    select = False
#            elif shape=='p':
#                plt.subplot(gs[0]).imshow(np.rot90(image),clim=[plotMin,plotMax],interpolation='None',extent=(x.min(),x.max(),y.min(),y.max()))
#                nPoints = int(raw_input("Number of Points (-1 until right click)?\n"))
#                p=np.array(ginput(nPoints))
#                print p
#                mpath=path.Path(p)
#                all_p = np.array([ (ix,iy) for ix,iy in zip(x.flatten(),y.flatten()) ] )
#                mask_r_nda = np.array([mpath.contains_points(all_p)]).reshape(x.shape)
#                plt.subplot(gs[1]).imshow(det.image(self.lda.run,mask_r_nda))
#                print 'mask from polygon (shape):',mask_r_nda.shape
#                print 'not implemented yet....'
#                if raw_input("Done?\n") in ["y","Y"]:
#                    select = False
#
#            if mask_r_nda is not None:
#                print 'created a mask....'
#                if mask is None:
#                    mask = mask_r_nda.astype(bool).copy()
#                else:
#                    mask = np.logical_or(mask,mask_r_nda)
#            print 'masked now: ',np.ones_like(x)[mask_r_nda.astype(bool)].sum()
#            print 'masked tot: ',np.ones_like(x)[mask.astype(bool)].sum()
#
#            fig=plt.figure(figsize=(6,6))
#            plt.show()
#            image_mask = img.copy(); image_mask[mask]=0;
#            plt.imshow(det.image(self.lda.run,image_mask),clim=[plotMin,plotMax])
#
#        if det.is_cspad2x2():
#            mask=mask.reshape(2,185*388).transpose(1,0)
#        else:
#            mask=mask.reshape(32*185,388)
#        #2x2 save as 71780 lines, 2 entries
#        #cspad save as 5920 lines, 388 entries
#        if raw_input("Save to calibdir?\n") in ["y","Y"]:
#            if raw_input("Invert?\n") in ["n","N"]:
#                mask = (~(mask.astype(bool))).astype(int)
#            srcStr=det.source.__str__().replace('Source("DetInfo(','').replace(')")','')
#            if det.is_cspad2x2():
#                dirname='/reg/d/psdm/%s/%s/calib/CsPad2x2::CalibV1/%s/pixel_mask/'%(self.lda.expname[:3],self.lda.expname,srcStr)
#            else:
#                dirname='/reg/d/psdm/%s/%s/calib/CsPad::CalibV1/%s/pixel_mask/'%(self.lda.expname[:3],self.lda.expname,srcStr)        
#            if not os.path.exists(dirname):
#                os.makedirs(dirname)
#            fname='%s-end.data'%self.lda.run
#            np.savetxt(dirname+fname,mask)
#        elif raw_input("Save to local?\n") in ["y","Y"]:
#            if raw_input("Invert?\n") in ["n","N"]:
#                mask = (~(mask.astype(bool))).astype(int)
#            np.savetxt('%s_mask_run%s.data'%(self.lda.expname,self.lda.run),mask)
#        return mask

    def psplot(self, *attrs, **kwargs):
        """
        Add psplot.  Automatically reshape AreaDetector 3D data to 2D for plotting.

        Parameters
        ----------
              
        local: bool
            open psplot locally (default)
        eventCode: int
            check if event code(s) are in data 
            (or alternatively not in date with - sign)
            see is_eventCodePresent
        nskip : int
            number of events to skip between plot updates

        Keywords
        --------
        title : str
            Plot title
        eventCode : list
            Plot if eventCode is present
        transpose : bool
            Transpose 2D image data before plotting 
        rot90 : bool
            Rotate 2D image data before plotting [Default = True for Cspad otherwise False]
        flipud : bool
            Flip 2D image up/down before plotting [Default = True for Cspad otherwise False]
        fliplr : bool
            Flip 2D image left/right before plotting

        """
        import numpy as np
        plot_error = '' 
        alias = self._alias

        if len(attrs) == 0:
            print('A string value specifying the plot object must be suppied')
            return

        if isinstance(attrs[0], list):
            attrs = attrs[0]

        attr_name = '_and_'.join(attrs)
        attr = attrs[0]

        # by default 
        if kwargs.get('local') is False:
            local = False
        else:
            local = True
        
        if 'eventCode' in kwargs:
            ecstrs = []
            for ec in kwargs.get('eventCode'):
                if ec > 0:
                    ecstrs.append(str(ec))
                else:
                    ecstrs.append('not'+str(-ec))
            ecname = '_'+'_and_'.join(ecstrs)
            ectitle = ' '+' and '.join(ecstrs)
        else:
            ecname = ''
            ectitle = ''

        if 'name' in kwargs:
            name = kwargs['name']
        else:
            name = alias+'_'+attr_name+ecname

        if 'title' in kwargs:
            title = kwargs['title']
        else:
            title = alias+' '+attr_name+ectitle
        
        if 'ts' in kwargs:
            ts = kwargs['ts']
        else:
            ts = self._ds._ievent

        if 'plot_type' in kwargs:
            plot_type = kwargs['plot_type']
        else:
            plot_type = None

        pub_opts = ['eventCode', 'nskip', 'transpose', 'fliplr', 'flipud', 'rot90']
        pub_kwargs = {key: item for key, item in kwargs.items() \
                      if key in pub_opts}

        if not plot_error and plot_type not in ['Image','XYPlot']:
            try:
                ndim = self._getattr(attr).ndim
                if ndim in [2,3]:
                    plot_type = 'Image'
                elif ndim == 1:
                    plot_type = 'XYPlot'
                elif ndim == 5:
                    # place holder for stats plots
                    if self._getattr(attr).shape[3] <= 8:
                        plot_type = 'xXYPlot'
                    else:
                        plot_type = 'xImage'

                else:
                    plot_error = 'Data with ndim = {:} not valid'.format(ndim)
            except:
                plot_error = 'Data must be numpy array of one or two dimensions.\n'               
        
        if not plot_error:
            if plot_type is 'Image':
                plt_opts = ['xlabel', 'ylabel', 'aspect_ratio', 'aspect_lock', 'scale', 'pos', 'xrange', 'yrange']
                plt_kwargs = {key: item for key, item in kwargs.items() \
                              if key in plt_opts}
        
                if ndim == 3:
                    calibData = self._getattr('calibData')
                    pub_kwargs['reshape'] = [calibData.shape[0]*calibData.shape[1], calibData.shape[2]]

                if attr == 'image':
                    calibData = self._getattr('calibData')
#                    scale = calibData.pixel_size/1000.
#                    if not plt_kwargs.get('pos'):
#                        plt_kwargs['pos'] = (calibData.image_xaxis[0]/1000., calibData.image_yaxis[0]/1000.)
#                    if not plt_kwargs.get('scale'):
#                        plt_kwargs['scale'] = (scale, scale)
                    if not plt_kwargs.get('xlabel'):
                        plt_kwargs['xlabel'] = 'X'
                    if not plt_kwargs.get('ylabel'):
                        plt_kwargs['ylabel'] = 'Y'

                    if self._det._det_config['devName'] in ['Cspad']:
                        if 'rot90' not in pub_kwargs:
                            pub_kwargs.update({'rot90': 1})
#                        if 'flipud' not in pub_kwargs:
#                            pub_kwargs.update({'flipud': True})

                plt_args = {'det': alias,
                            'attr': attrs,  
                            'name': name,
                            'plot_function': plot_type,
                            'ts': ts,
                            'title': title,
                            'kwargs': plt_kwargs,
                            'pubargs': pub_kwargs}
            
            elif plot_type is 'XYPlot':
                plt_opts = ['xlabel','ylabel','formats']
                plt_kwargs = {key: item for key, item in kwargs.items() \
                              if key in plt_opts}
                if 'xdata' in kwargs:
                    xdata = kwargs['xdata']
                else:
                    if attr in self._det_config['projection']:
                        xdata = self._det_config['projection'][attr]['xdata']
                        axis_name = self._det_config['projection'][attr].get('axis_name')
                        xunit = self._det_config['projection'][attr].get('xunit', '')
                        unit = self._det_config['projection'][attr].get('unit')
                        xlabel = '{:} [{:}]'.format(axis_name, xunit)
                        if not plt_kwargs.get('xlabel'):
                            plt_kwargs.update({'xlabel': xlabel})
                        if not plt_kwargs.get('ylabel'):
                            plt_kwargs.update({'ylabel': '[{:}]'.format(unit)})

                    else:
                        xdata = np.arange(len(self._getattr(attrs[0])))
                
                xdata = np.array(xdata, dtype='f')
                plt_args = {'det': alias,
                            'attr': attrs,
                            'xdata': xdata,
                            'name': name,
                            'plot_function': plot_type,
                            'ts': ts,
                            'title': title,
                            'kwargs': plt_kwargs,
                            'pubargs': pub_kwargs}
            else: 
                plot_error = 'Unknown plot type {:} \n'.format(plot_type)

        if plot_error:
            print('Error adding psplot:' )
            print(plot_error)
            return None
        else:
            print('psmon plot added -- use the following to view: ')
            print('--> psplot -s {:} -p 12301 {:}'.format(os.uname()[1], name))
            print('WARNING -- see notice when adding for -p PORT specification')
            print('           if default PORT=12301 not available')
            if 'psplot' not in self._det_config:
                self._det_config['psplot'] = {}

            self._det_config['psplot'][name] = plt_args
            if not publish.initialized:
                publish.init(local=local)
            
            if local:
                psmon_publish(self._evt)  

#    def __getattr__(self, attr):
#        if attr in self._plugins:
#            plugin = self._plugins.get(attr)
#            return plugin(self._ds)._det_add(attr)
#            if plugin:
#                return getattr(plugin, 'add')

#    def __setattr__(self, attr, val):
#        if attr not in self._init_attrs:
#            self._plugins.update({attr: val})

    def _getattr(self, attr):
        """
        Get detector attribute.
        """
        return getattr_complete(self._det,attr)

    def __dir__(self):
        all_attrs =  set(self._plugins.keys() +
                         self.__dict__.keys() + dir(AddOn))
        
        return list(sorted(all_attrs))


class IpimbData(object):
    """Tab accessibile dictified psana.Detector object.
       
       Attributes come from psana.Detector 
       with low level implementation done in C++ or python.  
       Boost is used for the C++.
    """

    _attrs = ['channel', 'sum', 'xpos', 'ypos'] 

    _attr_info = {
            'channel':     {'doc': 'Array of 4 channel values',
                            'unit': 'V'},
            'sum':         {'doc': 'Sum of all 4 channels',
                            'unit': 'V'},
            'xpos':        {'doc': 'Calulated X beam position',
                            'unit': 'mm'},
            'ypos':        {'doc': 'Calulated Y beam position',
                            'unit': 'mm'},
            } 

    def __init__(self, det, evt):
        self._evt = evt
        self._det = det
        self._data = {}
    
    @property
    def instrument(self):
        """
        Instrument to which this detector belongs.
        """
        return self._det.instrument()

    def show_info(self, **kwargs):
        """
        Show information for relevant detector attributes.
        """
        message = Message(quiet=True, **kwargs)
        try:
            items = sorted(self._attr_info.items(), key=operator.itemgetter(0))
            for attr, item in items:
                fdict = {'attr': attr, 'unit': '', 'doc': ''}
                fdict.update(**item)
                value = getattr(self, attr)
                if isinstance(value, str):
                    fdict['str'] = value
                elif isinstance(value, list):
                    if len(value) < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = 'list'
                elif hasattr(value,'mean'):
                    if value.size < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = '<{:.5}>'.format(value.mean())
                else:
                    try:
                        fdict['str'] = '{:12.5g}'.format(value)
                    except:
                        fdict['str'] = str(value)

                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
        except:
            message('No Event')

        return message

    def __getattr__(self, attr):
        if attr in self._attrs:
            if attr not in self._data:
                self._data.update({attr: getattr(self._det, attr)(self._evt)})
             
            return self._data.get(attr)


    def __dir__(self):
        all_attrs =  set(self._attrs +
                         self.__dict__.keys() + dir(IpimbData))
        
        return list(sorted(all_attrs))

class GenericWaveformData(object):
    """Tab accessibile dictified psana.Detector object.
       
       Attributes come from psana.Detector 
       with low level implementation done in C++ or python.  
       Boost is used for the C++.
    """

    _attrs = ['waveform'] 

    _attr_info = {
            'waveform':    {'doc': 'Waveform array',
                            'unit': 'V'},
#            'wftime':      {'doc': 'Waveform sample time',
#                            'unit': 's'},
            } 

    def __init__(self, det, evt):
        self._evt = evt
        self._det = det

    def show_info(self, **kwargs):
        """Show information for relevant detector attributes.
        """
        message = Message(quiet=True, **kwargs)
        try:
            items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
            for attr, item in items:
                fdict = {'attr': attr, 'unit': '', 'doc': ''}
                fdict.update(**item)
                value = getattr(self, attr)
                if isinstance(value, str):
                    fdict['str'] = value
                elif isinstance(value, list):
                    if len(value) < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = 'list'
                elif hasattr(value,'mean'):
                    if value.size < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = '<{:.5}>'.format(value.mean())
                else:
                    try:
                        fdict['str'] = '{:12.5g}'.format(value)
                    except:
                        fdict['str'] = str(value)

                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
        except:
            message('No Event')

        return message

    @property
    def waveform(self):
        """
        List of waveforms
        """
        return self._det(self._evt)

#    def __getattr__(self, attr):
#        if attr in self._attrs:
#            return getattr(self._det, attr)(self._evt)
#
    def __dir__(self):
        all_attrs =  set(self._attrs +
                         self.__dict__.keys() + dir(GenericWaveformData))
        
        return list(sorted(all_attrs))


class WaveformData(object):
    """
    Tab accessibile dictified psana.Detector object.
       
    Attributes come from psana.Detector with low level implementation 
    done in C++ or python.  
    Boost is used for the C++.
    """

    _attrs = ['raw', 'waveform', 'wftime'] 

    _attr_info = {
            'waveform':    {'doc': 'Waveform array',
                            'unit': 'V'},
            'wftime':      {'doc': 'Waveform sample time',
                            'unit': 's'},
            } 

    def __init__(self, det, evt):
        self._evt = evt
        self._det = det
        self._data = {}

    @property
    def instrument(self):
        """
        Instrument to which this detector belongs.
        """
        return self._det.instrument()

    def show_info(self, **kwargs):
        """
        Show information for relevant detector attributes.
        """
        message = Message(quiet=True, **kwargs)
        try:
            items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
            for attr, item in items:
                fdict = {'attr': attr, 'unit': '', 'doc': ''}
                fdict.update(**item)
                value = getattr(self, attr)
                if isinstance(value, str):
                    fdict['str'] = value
                elif isinstance(value, list):
                    if len(value) < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = 'list'
                elif hasattr(value,'mean'):
                    if value.size < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = '<{:.5}>'.format(value.mean())
                else:
                    try:
                        fdict['str'] = '{:12.5g}'.format(value)
                    except:
                        fdict['str'] = str(value)

                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
        except:
            message('No Event')

        return message

    def __getattr__(self, attr):
        """
        Only access psana.Detector data once.
        """
        if attr in self._attrs:
            if attr not in self._data:
                self._data.update({attr: getattr(self._det, attr)(self._evt)})
             
            return self._data.get(attr)

    def __dir__(self):
        all_attrs =  set(self._attrs +
                         self.__dict__.keys() + dir(WaveformData))
        
        return list(sorted(all_attrs))

class GenericWaveformCalibData(object):
    """
    Generic Waveform Calib Data -- not implented
    """
    def __init__(self, det, evt):
        self._evt = evt
        self._det = det


class WaveformCalibData(object):
    """
    Waveform Calibration data using psana.Detector access.
    """

    _attrs = ['runnum'] 

    _attr_info = {
            'runnum':      {'doc': 'Run number',
                            'unit': ''}
            }

    def __init__(self, det, evt):
        self._evt = evt
        self._det = det

    @property
    def instrument(self):
        """
        Instrument to which this detector belongs.
        """
        return self._det.instrument()

    def print_attributes(self):
        """
        Print detector attributes.
        """
        self._det.print_attributes()

    def set_calibration(self):
        """
        On/off correction of time.'
        """
        if self._det.dettype == 16:
            self._det.set_correct_acqiris_time()
        elif self._det.dettype == 17:
            self._det.set_calib_imp()

    def show_info(self, **kwargs):
        """
        Show information for relevant detector attributes.
        """
        message = Message(quiet=True, **kwargs)
        try:
            items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
            for attr, item in items:
                fdict = {'attr': attr, 'unit': '', 'doc': ''}
                fdict.update(**item)
                value = getattr(self, attr)
                if isinstance(value, str):
                    fdict['str'] = value
                elif isinstance(value, list):
                    if len(value) < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = 'list'
                elif hasattr(value,'mean'):
                    if value.size < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = '<{:.5}>'.format(value.mean())
                else:
                    try:
                        fdict['str'] = '{:12.5g}'.format(value)
                    except:
                        fdict['str'] = str(value)

                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
        except:
            message('No Event')

        return message

    def __getattr__(self, attr):
        if attr in self._attrs:
            return getattr(self._det, attr)(self._evt)

    def __dir__(self):
        all_attrs =  set(self._attrs +
                         self.__dict__.keys() + dir(WaveformCalibData))
        
        return list(sorted(all_attrs))


class ImageData(object):
    """
    Tab accessibile dictified psana Detector object.
       
    Attributes come from psana.Detector with low level implementation 
    done in C++ or python.  Boost is used for the C++.
    """
    _attrs = ['image', 'raw', 'calib', 'corr', 'photons', 'shape', 'size'] 
    _attr_info = {
            'shape':       {'doc': 'Shape of raw data array', 
                            'unit': ''},
            'size':        {'doc': 'Total size of raw data', 
                            'unit': ''},
            'raw':         {'doc': 'Raw data', 
                            'unit': 'ADU'},
            'calib':       {'doc': 'Calibrated data',
                            'unit': 'ADU'},
            'image':       {'doc': 'Reconstruced 2D image from calibStore geometry',
                            'unit': 'ADU'},
            'photons':     {'doc': '2-d or 3-d array of integer number of merged photons',
                            'unit': 'ADU'},
            'corr':        {'doc': 'Pedestal corrected data (no common mode correction)',
                            'unit': 'ADU'},
            } 

    #def __init__(self, det, evt, opts={}):
    def __init__(self, det, evt):
        self._evt = evt
        self._det = det
        self._data = {}
        #self._opts = opts

    @property
    def instrument(self):
        """
        Instrument to which this detector belongs.
        """
        return self._det.instrument()

    @property
    def photons(self):
        """
        Returns photons with correction for split photons between neighboring pixels.
        adu_per_photon 
          - photon conversion factor parameter 
          - can be updated with set_adu_per_photon.
          - default = 30
        thr_fraction 
          - fraction of the merged intensity which gets converted to one photon
          - can be updated with set_thr_fraction.
          - default = 0.9
        see:  https://confluence.slac.stanford.edu/display/PSDM/Hit+and+Peak+Finding+Algorithms#HitandPeakFindingAlgorithms-Photoncounting
        """
        if self.adu_per_photon:
            adu_per_photon = self.adu_per_photon
        else:
            adu_per_photon = 30
        if self.thr_fraction:
            thr_fraction = self.thr_fraction
        else:
            thr_fraction = 0.9 

        return self._det.photons(self._evt, 
                    adu_per_photon=adu_per_photon,
                    thr_fraction=thr_fraction)

    @property
    def corr(self):
        """
        Pedestal corrected raw data.
        """
        ped = self._det.pedestals(self._evt)
        if self.raw is not None and ped is not None:
            return self.raw-ped
        else:
            return self.raw

    def set_thr_fraction(self, thr_fraction):
        """
        Sets threshold fraction for photons calculation
        """
        self.add.parameter(thr_fraction=thr_fraction)

    def set_adu_per_photon(self, adu_per_photon):
        """
        Sets adu per photon for photons calculation
        """
        self.add.parameter(adu_per_photon=adu_per_photon)

    def make_image(self, nda):
        """
        Make an image from the input numpy array based on the 
        geometry in the calib directory for this event.
        
        Parameters
        ----------
        nda : np.array
            input array
        """
        return self._det.image(self._evt, nda)

    def common_mode_correction(self, nda):
        """
        Return the common mode correction for the input numpy 
        array (pedestal-subtracted). 
        """
        return self._det.common_mode_correction(self._evt, nda)
        
    def common_mode_apply(self, nda):
        """
        Apply in place the common mode correction for the input 
        numpy array (pedestal-subtracted). 
        """
        self._det.common_mode_apply(self._evt, nda)

    def show_info(self, attrs=None, **kwargs):
        """
        Show information for relevant detector attributes.
        """
        if attrs:
            if not isinstance(attrs, list):
                attrs = ['raw', 'corr']
        message = Message(quiet=True, **kwargs)
        if self.size > 0 or self.raw is not None:
            items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
            for attr, item in items:
                if attrs and attr not in attrs:
                    continue
                value = getattr(self, attr)
                strval = _repr_value(value)
                fdict = {'attr': attr, 'str': strval, 'unit': '', 'doc': ''}
                fdict.update(**item)
                
#                if isinstance(value, str):
#                    fdict['str'] = value
#                elif isinstance(value, list):
#                    if len(value) < 5:
#                        fdict['str'] = str(value)
#                    else:
#                        fdict['str'] = 'list'
#                elif hasattr(value,'mean'):
#                    if value.size < 5:
#                        fdict['str'] = str(value)
#                    else:
#                        fdict['str'] = '<{:.5}>'.format(value.mean())
#                else:
#                    try:
#                        fdict['str'] = '{:12.5g}'.format(value)
#                    except:
#                        fdict['str'] = str(value)
#
                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
        else:
            message('No Event')

        return message

    def __getattr__(self, attr):
        """
        Only access psana.Detector data once.
        """
        if attr in self._attrs:
            if attr not in self._data:
#                opts = self._opts.get(attr, {})
#                self._data.update({attr: getattr(self._det, attr)(self._evt, **opts)})
                try:
                    val = getattr(self._det, attr)(self._evt)
                except:
                    if attr == 'size':
                        val =  self.raw.size
                    if attr == 'shape':
                        val = self.raw.shape

                self._data.update({attr: val})
            
            else:
                val = self._data.get(attr)

            return val

    def __dir__(self):
        all_attrs =  set(self._attrs +
                         self.__dict__.keys() + dir(ImageData))
        
        return list(sorted(all_attrs))


class ImageCalibData(object):
    """
    Calibration Data from psana Detector object.
    """

    _attrs = ['shape', 'size', 'ndim', 'pedestals', 'rms', 'gain', 'bkgd', 'status',
              'common_mode', 'runnum',
              'areas', 'indexes_x', 'indexes_y', 'pixel_size',
              'coords_x', 'coords_y', 'coords_z', 
              'image_xaxis', 'image_yaxis',
              ] 
    _attr_info = {
            'runnum':      {'doc': 'Run number',
                            'unit': ''},
            'shape':       {'doc': 'Shape of raw data array', 
                            'unit': ''},
            'size':        {'doc': 'Total size of raw data', 
                            'unit': ''},
            'ndim':        {'doc': 'Number of dimensions of raw data', 
                            'unit': ''},
            'pedestals':   {'doc': 'Pedestals from calibStore', 
                            'unit': 'ADU'},
            'rms':         {'doc': '', 
                            'unit': 'ADU'},
            'gain':        {'doc': 'Pixel Gain factor from calibStore', 
                            'unit': ''},
            'bkgd':        {'doc': '', 
                            'unit': ''},
            'status':      {'doc': '',
                            'unit': ''},
            'common_mode': {'doc': 'Common mode parameters', 
                            'unit': ''},
            'areas':       {'doc': 'Pixel area correction factor', 
                            'unit': ''},
            'indexes_x':   {'doc': 'Pixel X index', 
                            'unit': ''},
            'indexes_y':   {'doc': 'Pixel Y index', 
                            'unit': ''},
            'pixel_size':  {'doc': 'Pixel Size',
                            'unit': 'um'},
            'coords_x':    {'doc': 'Pixel X coordinate', 
                            'unit': 'um'},
            'coords_y':    {'doc': 'Pixel Y coordinate', 
                            'unit': 'um'},
            'coords_z':    {'doc': 'Pixel Z coordinate', 
                            'unit': 'um'},
            'image_xaxis': {'doc': 'Image X coordinate', 
                            'unit': 'um'},
            'image_yaxis': {'doc': 'Image Y coordinate', 
                            'unit': 'um'},
            } 

    def __init__(self, det, evt):
        self._evt = evt
        self._det = det
        self._info = {}

    @property
    def instrument(self):
        """
        Instrument to which this detector belongs.
        """
        return self._det.instrument()

    def _get_origin(self, **kwargs):
        """
        Get image origin indecies
        """
        ixo, iyo = self._det.point_indexes(self._evt, **kwargs)

        return ixo, iyo

    def _get_point_indexes(self, **kwargs):
        """
        """
        import numpy as np
        ixo, iyo = self._det.point_indexes(self._evt, **kwargs)
        if self.indexes_y is not None:
            ny = self.indexes_y.max()+1
            yaxis = (np.arange(ny)-iyo)*self.pixel_size
        else:
            if len(self.shape) == 2:
                yaxis = np.arange(self.shape[1])
            else:
                yaxis = None

        if self.indexes_x is not None:
            nx = self.indexes_x.max()+1
            xaxis = (np.arange(nx)-ixo)*self.pixel_size
        else:
            if len(self.shape) == 2:
                xaxis = np.arange(self.shape[0])
            else:
                xaxis = None
        
        self._info.update({'ixo': {'attr': 'ixo',
                                   'doc':  'Image x index of origin',
                                   'unit': '',
                                   'value': ixo}})
        if xaxis is not None:
            self._info.update({'xaxis': {'attr': 'xaxis',
                                         'doc':  'Reconstructed image xaxis',
                                         'unit': '',
                                         'value': xaxis}})

        self._info.update({'iyo': {'attr': 'iyo',
                                   'doc':  'Image y index of origin',
                                   'unit': '',
                                   'value': iyo}})
        if yaxis is not None:
            self._info.update({'yaxis': {'attr': 'yaxis',
                                         'doc':  'Reconstructed image yaxis',
                                         'unit': '',
                                         'value': yaxis}})

        return xaxis, yaxis

    @property
    def xaxis(self):
        """
        Reconstructed image x axis.
        """
        return self.image_xaxis
#        item = self._info.get('xaxis')
#        if item:
#            xaxis = item.get('value')
#        else:
#            xaxis, yaxis = self._get_point_indexes()
#            
#        return xaxis

    @property
    def yaxis(self):
        """
        Reconstruced image y axis.
        """
        return self.image_yaxis
#        item = self._info.get('yaxis')
#        if item:
#            yaxis = item.get('value')
#        else:
#            xaxis, yaxis = self._get_point_indexes()
#            
#        return yaxis

    def set_gain_mask_factor(factor=6.85):
        """
        Set Gain mask factor.  Default=6.85
        Passed to gain_mask(...) in the calib and image methods
        """
        self._det.set_gain_mask_factor(factor=factor)

    def set_do_offset(do_offset=True):
        """
        Switch mode of the Camera type of detector.
        Control parameter to turn on/off Camera intensity offset
        """
        self._det.set_do_offset(do_offset=do_offset)

    def mask(self, calib=True, status=True, 
                   edges=True, central=True, 
                   unbond=True, unbondnbrs=True):
        """
        Generate image mask.

        Parameters
        ----------
        calib: bool
            mask from file in calib directory.
        status: bool
            pixel status from file in calib director.
        edges: bool
            mask detector module edge pixels (mbit +1 in mask_geo).
        central: bool
            mask wide central columns (mbit +2 in mask_geo).
        unbond: bool
            mask unbonded pixels (mbit +4 in mask_geo).
        unbondnbrs: bool
            mask unbonded neighbour pixels (mbit +8 in mask_geo).
        
        Returns
        -------
        combined mask: array-like
        """
        return self._det.mask(self._evt, calib=calib, status=status, edges=edges, 
                              central=central, unbond=unbond, unbondnbrs=unbondnbrs)

    def mask_geo(self, mbits=15): 
        """
        Return geometry mask for given mbits keyword.
        Default is mbits=15 to mask edges, wide central columns,
        non-bo pixels and their neighbors

        Parameters
        ----------
        mbits: int 
            +1-edges; 
            +2-wide central cols; 
            +4 unbonded pixel; 
            +8-unbonded neighbour pixels;
        
        """
        return self._det.mask_geo(self._evt, mbits=mbits)

    def print_attributes(self):
        """
        Print detector attributes.
        """
        self._det.print_attributes()

    def show_info(self, **kwargs):
        """
        Show information for relevant detector attributes.
        """
        message = Message(quiet=True, **kwargs)
        if self.size > 0:
            items = sorted(self._attr_info.items(), key = operator.itemgetter(0))
            for attr, item in items:
                fdict = {'attr': attr, 'unit': '', 'doc': ''}
                fdict.update(**item)
                value = getattr(self, attr)
                if isinstance(value, str):
                    fdict['str'] = value
                elif isinstance(value, list):
                    if len(value) < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = 'list'
                elif hasattr(value,'mean'):
                    if value.size < 5:
                        fdict['str'] = str(value)
                    else:
                        fdict['str'] = '<{:.5}>'.format(value.mean())
                else:
                    try:
                        fdict['str'] = '{:12.5g}'.format(value)
                    except:
                        fdict['str'] = str(value)

                message('{attr:18s} {str:>12} {unit:7} {doc:}'.format(**fdict))
        else:
            message('No Event')

        return message

    def __getattr__(self, attr):
        if attr in self._attrs:
            return (getattr(self._det, attr)(self._evt))
        
    def __dir__(self):
        all_attrs =  set(self._attrs +
                         self.__dict__.keys() + dir(ImageCalibData))
        
        return list(sorted(all_attrs))


class EpicsConfig(object):
    """
    Tab Accessible configStore Epics information.
    Currently relatively simple, but expect this to be expanded
    at some point with more PV config info with daq update.
    """

    _pv_attrs = ['description', 'interval', 'pvId']

    def __init__(self, configStore):

        # move to PsanaSrcData objects
        self._pvs = {}
        for key in configStore.keys():
            if key.type() and key.type().__module__ == 'psana.Epics':
                a = configStore.get(key.type(),key.src())
                for pv in a.getPvConfig():
                    pvdict = {attr: getattr(pv, attr)() for attr in self._pv_attrs} 
                    self._pvs[pv.description()] = pvdict

    def show_info(self, **kwargs):
        message = Message(quiet=True, **kwargs)
        #for alias, items in self._pvs.items():
        for alias, items in sorted(self._pvs.items(), key=operator.itemgetter(0)):
            message('{:18s} {:}'.format(alias, item.pvId))

        return message

    def __getattr__(self, attr):
        if attr in self._pvs:
            return self._pvs.get(attr)

    def __dir__(self):
        all_attrs =  set(self._pvs.keys() +
                         self.__dict__.keys() + dir(EpicsConfig))
        
        return list(sorted(all_attrs))


class EpicsData(object):
    """
    Epics data from psana epicsStore.
    
    Parameters
    ----------
    ds : object
        PyDataSource.DataSource object
    
    Returns
    -------
    Dictified representation of ds.env().epicsStore()

    Example
    -------
    epicsStore = EpicsData(ds)
    """

    def __init__(self, ds):

        self._ds = ds

        pv_dict = {}
        epicsStore = self._ds.env().epicsStore()
        self.epicsConfig = EpicsConfig(self._ds.env().configStore())

        for pv in  epicsStore.names():
            name = re.sub(':|\.','_',pv)
            #check if valid -- some old data had aliases generated from comments in epicsArch files.
            if re.match("[_A-Za-z][_a-zA-Z0-9]*$", name) and not ' ' in name and not '-' in name:
                pvname = epicsStore.pvName(pv)
                if pvname:
                    pvalias = pv
                else:
                    pvalias = epicsStore.alias(pv)
                    pvname = pv

                pvalias = re.sub(':|\.|-| ','_',pvalias)
                components = re.split(':|\.|-| ',pv)
                if len(components) == 1:
                    components = re.split('_',pv,1)
                
                # check if alias has 2 components -- if not fix
                if len(components) == 1:
                    pv = '_'.join([components[0], components[0]])
                    components = re.split('_',pv,1)

                for i,item in enumerate(components):
                    try:
                        if item[0].isdigit():
                            components[i] = 'n'+components[i]
                    except:
                        pass

                pv_dict[name] =  { 'pv': pvname,
                                   'alias': pvalias,
                                   'components': components,
                                 }
        self._pv_dict = pv_dict
        self._attrs = list(set([val['components'][0] for val in self._pv_dict.values()]))

    def __getattr__(self, attr):
        if attr in self._attrs:
            attr_dict = {key: pdict for key,pdict in self._pv_dict.items()
                         if pdict['components'][0] == attr}
            return PvData(attr_dict, self._ds, level=1)
        
        if attr in dir(self._ds.env().epicsStore()):
            return getattr(self._ds.env().epicsStore(),attr)

    def __dir__(self):
        all_attrs = set(self._attrs +
                        dir(self._ds.env().epicsStore()) +
                        self.__dict__.keys() + dir(EpicsData))
        return list(sorted(all_attrs))



class PvData(object):
    """
    Epics PV Data.
    """

    def __init__(self, attr_dict, ds, level=0):
        self._attr_dict = attr_dict
        self._ds = ds
        self._level = int(level)
        self._attrs = list(set([pdict['components'][level]
                                for key,pdict in attr_dict.items()]))

    def _get_pv(self, pv):
        return EpicsStorePV(self._ds.env().epicsStore(), pv)

    def show_info(self, **kwargs):
        """
        Show information from PVdictionary for all PV's starting with 
        the specified dictified base.
        (i.e. ':' replaced by '.' to make them tab accessible in python)
        """
        message = Message(self.get_info(), quiet=True, **kwargs)
        return message

    def get_info(self):
        """
        Return string representation of all PV's starting with 
        the specified dictified base.
        (i.e. ':' replaced by '.' to make them tab accessible in python)
        """
        info = ''
        items = sorted(self._attr_dict.items(), key=operator.itemgetter(0))
        for key,pdict in items:
            alias = pdict['alias']
            if alias:
                name = alias
                pv = pdict['pv']
            else:
                name = pdict['pv']
                pv = ''

            pvfunc = self._get_pv(pdict['pv'])
            value = pvfunc.value
            if pvfunc.isCtrl:
                comment = 'isCtrl'
            else:
                comment = ''

            try:
                info += '{:30s} {:12.4g} -- {:30s} {:10}\n'.format( \
                        name, value, pv, comment)
            except:
                info += '{:30s} {:>12} -- {:30s} {:10}\n'.format( \
                        name, value, pv, comment)
        return info

    def __getattr__(self, attr):
        if attr in self._attrs:
            attr_dict = {key: pdict for key,pdict in self._attr_dict.items()
                         if pdict['components'][self._level] == attr}
            if len(attr_dict) == 1:
                key = attr_dict.keys()[0]
                if len(self._attr_dict[key]['components']) == (self._level+1):
                    pv = self._attr_dict[key]['pv']
                    return self._get_pv(pv)
            if len(attr_dict) > 0:
                return PvData(attr_dict, self._ds, level=self._level+1)

    def __repr__(self):
        return self.get_info()

    def __dir__(self):
        all_attrs = set(self._attrs +
                        self.__dict__.keys() + dir(PvData))
        return list(sorted(all_attrs))


class EpicsStorePV(object):
    """
    Epics PV access from epicsStore. 
    """

    def __init__(self, epicsStore, pv):
        self._epicsStore = epicsStore
        self._pvname = pv
        self._store = epicsStore.getPV(pv)
        self._attrs = [attr for attr in dir(self._store) \
                if not attr.startswith('_')]
        self._show_attrs = [attr for attr in self._attrs \
                if attr not in ['dbr','stamp']]

    def get_info(self):
        info = '-'*80+'\n'
        info += '{:} = {:} -- {:}\n'.format(self._pvname, \
                                self.value, self.stamp)
        info += '-'*80+'\n'
        for attr in self._show_attrs:
            val = self.get(attr)
            info += '{:20} {:12}\n'.format(attr, val)
        
        return info

    def show_info(self, **kwargs):
        message = Message(self.get_info(), quiet=True, **kwargs)
        return message

    def get(self, attr):
        if attr in self._attrs:
            if attr is 'value':
                return self._epicsStore.value(self._pvname)
            else:
                val = getattr(self._store,attr)
                try:
                    if attr is 'stamp':
                        return TimeStamp(val())
                    else:
                        return val() 
                except:
                    return val
        else:
            return None

    def __str__(self):
        import numpy as np
        if len(self.data) > 1 and isinstance(self.data, np.ndarray):
            value = '<{:}>'.format(mean(self.data))
        else:
            value = self.value

        return '{:}'.format(value)

    def __repr__(self):
        return '< {:} = {:}, {:} -- {:} >'.format(self._pvname, \
                str(self), self.stamp.time, \
                self.__class__.__name__)

    def __getattr__(self, attr):
        if attr in self._attrs:
            return self.get(attr)

    def __dir__(self):
        all_attrs = set(self._attrs +
                        self.__dict__.keys() + dir(EpicsStorePV))
        return list(sorted(all_attrs))


class TimeStamp(object):
    """
    Class to represent time stamp objects
    """

    def __init__(self, stamp):
        self.sec = stamp.sec()
        self.nsec = stamp.nsec()

    @property
    def date(self):
        """
        Time stamp date representation.
        """
        return time.strftime('%Y-%m-%d', 
                time.localtime(self.sec))

    @property
    def time(self): 
        """
        Time stamp time representation.
        """
        EventTimeStr = time.strftime('%H:%M:%S',
                time.localtime(self.sec))
        EventTimeStr += '.{:04}'.format(int(self.nsec/1e5))
        return EventTimeStr

    def __str__(self):
        return '{:}.{:} sec'.format(self.sec, self.nsec)

    def __repr__(self):
        return '< {:}: {:} >'.format(self.__class__.__name_, _self.__str__)


def _update_stats(evt):
    """
    Update Welford statistics.
    """
    from .welford import Welford
    istep = evt._ds._istep
    for alias, det in evt._dets.items():
        for name, item in det._det_config['stats'].items():
            attr = item['attr']
            vals = getattr_complete(det, attr)
            eventCodes = item['eventCodes']
            if vals is not None:
                for ec in eventCodes:
                    # Note that no Evr for Controls cameras
                    if hasattr(det.Evr, 'present') and det.Evr.present(ec):
                        if ec not in item['funcs']:
                            item['funcs'].update({ec: {}})
                        funcs = item['funcs'].get(ec) 
                        if istep not in funcs:
                            funcs.update({istep: Welford()})
                        fec = funcs.get(istep)
                        try:
                            if not fec.shape or fec.shape == vals.shape:
                                fec(vals)
                            else:
                                print(('stats update error', det._alias, name, attr, ec, istep, vals))
                        except:
                            print(('stats update error', det._alias, name, attr, ec, istep, vals))
#            else:
#                print alias, name, attr, vals


