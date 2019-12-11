from __future__ import absolute_import
__all__ = []

#from PyDataSource import DataSource

#import matplotlib as mpl
#mpl.use('Agg')
from .PyDataSource import *
#from psxarray import * 
from . import h5write
from . import xarray_utils
from . import filter_methods
#import examples
from .h5write import *
from .exp_summary import get_exp_summary 
from .epicsarchive import EpicsArchive
from .build_html import Build_html, Build_experiment
from .arp_tools import post_report
from .psutils import get_run_from_id
from .xarray_utils import open_cxi_psocake

__version__ = '00.01.01'

import logging

logger = logging.getLogger('PyDataSource')
logger.setLevel(logging.DEBUG)

#fh = logging.FileHandler('data_summary.log')
#fh.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
#fh.setFormatter(formatter)
ch.setFormatter(formatter)

#logger.addHandler(fh)
logger.addHandler(ch)


def set_logger_level(lvl):
    logger.setLevel( getattr(logging,lvl) )
#    fh.setLevel( getattr(logging,lvl) )
    ch.setLevel( getattr(logging,lvl) )
    return

def logger_flush():
#    fh.flush()
    return


