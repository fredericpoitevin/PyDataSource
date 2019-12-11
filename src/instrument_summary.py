from __future__ import print_function
from __future__ import absolute_import

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

def get_experiment_list(instrument):
    import glob
    exp_folders = glob.glob('/reg/d/psdm/{:}/*'.format(instrument))
    experiments = [a.split('/')[5] for a in exp_folders if a.split('/')[5].startswith(instrument)]
    return experiments

def get_instrument_xarray(experiments=None, instrument=None):
    """
    get experiment xarray
    """
    from .PyDataSource import get_exp_summary
    if not experiments:
        if instrument:
            experiments = get_experiment_list(instrument)
        else:
            raise Exception('Must supply either experiments list or instrument name')       
 
    dets = []
    detectors = []
    aexperiments = {}
    for exp in experiments:
        instrument = exp[0:3]
        try:
            es = get_exp_summary(exp=exp)
            xruns = es._add_run_details()
            xruns.attrs['exp'] = exp
            for attr in ['instrument', 'station', 'exp_dir']:
                try:
                    xruns.attrs[attr] = getattr(es, attr)
                except:
                    xruns.attrs[attr] = ''

            # load run summary pv dict
            #rundict = es.load_run_summary()
            dets += xruns.attrs.get('dets', [])
            detectors += xruns.attrs.get('detectors', [])
            aexperiments[exp] = xruns 
        except:
            print('cannot add {:}'.format(exp))

    return aexperiments            


class InstrumentSummary(object):
    """
    Instrument Summary
    """
    def __init__(self, instrument):
        self.instrument = instrument
        self.experiments = get_experiment_list(instrument)
        self._xinst = get_instrument_xarray(experiments=self.experiments)
        self.xruns = xr.merge([x.swap_dims({'run':'run_id'}) for rn, x in self._xinst.items()])

class LCLSsummary(object):
    """
    LCLS Run summary
    """
    def __init__(self, instrument):
        self._instruments = []


