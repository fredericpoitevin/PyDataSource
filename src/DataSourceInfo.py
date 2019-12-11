from __future__ import print_function
import sys, os, re, time
import operator
import traceback
import psana

def live_source(monshmserver='psana', **kwargs):
    """Returns psana source string for live data from shared memory on the current node.
       The standard convention is to have the shared memry servers either named 'psana'
       or the instrument name in all caps.  This will return the source string for 
       the valid one based on searching the local file system.  The process will be
       found at:

          '/dev/shm/PdsMonitorSharedMemory_'+monshmserver
    """
    from glob import glob
    import os

    shm_srvs = glob('/dev/shm/PdsMonitorSharedMemory_'+monshmserver)
    if shm_srvs == []:
        hostsplit = os.uname()[1].split('-')
        instrument = hostsplit[1]
        monshmserver = instrument.upper()
        shm_srvs = glob('/dev/shm/PdsMonitorSharedMemory_'+monshmserver)
    
    if shm_srvs != []:
        try:
            MPI_RANK = 0
            source_str = 'shmem={:}.0:stop=no'.format(monshmserver)
        except:
            print('Exception in finding shared memory server: ',shm_srvs)
            source_str = None
    else:
        source_str = None

    return source_str


class DataSourceInfo(object):
    """-------------------------------------------------------------------------------
       data_source class built from keyword arguments that can be accessed as attributes

            data_source = 'exp=CXI/cxic0115:run=10'

        The following are equivalent:

            data_source = DataSourceInfo(exp='cxic0115', run=10).data_source
            data_source = str(DataSourceInfo(exp='cxic0115', run=10))

        You can also specify the run_id instead of the experiment run number 
        (run numbers over 100000 are assumed to be a run_id as an aid in automatic run processing)

        You can also load a data source with keyword options:

            smd:  small data support -- has become standard for experiments after Oct 2015
            h5:   loads hdf5 data instead of xtc
            ffb:  use a special set of disks (fast-feedback, or "FFB") reserved for the running experiment
            live: live ffb data - waits for additional data in case the analysis "catches up" to the DAQ

        The shared memory data_source can be loaded with the monshmemsrver keyword:

            data_source = str(DataSource(monshmemsrver='psana'))

        But shared memory should alse be automatically detected if no arguments are
        supplied and you are on a shared memery server.
            
            data_source = str(DataSource())

    """
    _exp_defaults = {'instrument':    None, 
                     'exp':           None, 
                     'h5':            None,
                     'run':           0,
                     'stream':        None,
                     'smd':           None, 
                     'station':       0,
                     'idx':           None,
                     'ffb':           None,
                     'live':          None,
                     'monshmserver':  None,
                     'dir':           None,
                     'cfg':           None}

    def __init__(self, data_source=None, **kwargs):
        self.data_source = self._set_data_source(data_source=data_source, **kwargs)

        self.exp_dir = os.path.join('/reg/d/psdm/',self.instrument,self.exp)
        self.xtc_dir = os.path.join(self.exp_dir,'xtc')
        if os.path.isdir(os.path.join(self.exp_dir,'res')): 
            self.res_dir = os.path.join(self.exp_dir,'res')
        else:
            self.res_dir = os.path.join(self.exp_dir,'results')
       
    def _set_exp_defaults(self, **kwargs):
        """Sets experiment defaults based on kwargs and defaults.
        """
        for key, val in self._exp_defaults.items():
            setattr(self, key, kwargs.get(key, val))

        if self.exp is not None:
            self.exp = str(self.exp)
            self.instrument = self.exp[0:3]
        
        run_id = kwargs.get('run_id')
        if not run_id and self.run > 100000:
            run_id = self.run

        if run_id:
            import psutils
            self.run = psutils.get_run_from_id(run_id, self.exp) 
            
#        self.run = int(self.run)

#        inst_id = '{:}:{:}'.format(self.instrument.upper(), self.station)

    def _set_data_source(self, data_source=None, valid_streams=True, **kwargs): 
        import psutils
        self._set_exp_defaults(**kwargs)

        if self.monshmserver:
            self.idx = False
            if not data_source:
                data_source = live_source(**kwargs)

        if data_source:
            opts = data_source.split(':')
            for opt in opts:
                items = opt.split('=')
                key = items[0]
                if key not in self._exp_defaults:
                    self._exp_defaults.update({key: None})
                    setattr(self, key, None)

                if len(items) == 2:
                    value = items[1]
                    if key in ['run', 'station']:
                        value = int(value)
                    setattr(self, key, value)
                else:
                    setattr(self, key, True)

        else:

            if self.exp and self.run > 0:
                self.instrument = self.exp[0:3]
     
                data_source = "exp={exp}:run={run}".format(exp=self.exp,run=self.run)
               
                stream = None
                if self.stream:
                    if isinstance(self.stream, list):
                        stream = ','.join([str(a) for a in self.stream]) 
                    else:
                        stream = str(self.stream)
                elif valid_streams:
                    # specify streams if any are not valid
                    ok_stream = psutils.run_available(exp=self.exp,run=self.run, valid_streams=True)
                    all_stream = psutils.run_available(exp=self.exp,run=self.run, all_streams=True)
                    if ok_stream and ok_stream != all_stream:
                        stream = ','.join([str(a) for a in ok_stream]) 

                if self.h5:
                    data_source += ":h5"
                else:
                    if stream and not self.ffb:
                        data_source += ":stream={:}".format(stream)

                    if self.idx:
                        data_source += ":idx"
                        self.smd = False
                    else:
                        data_source += ":smd"
                        self.smd = True
                
                    if self.ffb:
                        # Set dir to standard ffb dir if not specifically passed in an unusual 
                        # circumstance such as testing.
                        if not self.dir:
                            self.dir = os.path.join('/reg/d/ffb/',self.instrument,self.exp,'xtc')
                        data_source += ":dir={:}".format(self.dir)
                        if self.live:
                            data_source += ":live"

            else:
                print('No data source specified, so assume this is shared memory.')
                data_source = live_source(**kwargs)
                self.monshmserver = data_source
                self.idx = False

        if not self.instrument:
            if self.exp is None:
                try:
                    if self.monshmserver:
                        self.exp = psutils.active_experiment()
                except:
                    print('Cannot determine active experiment for shared memory')
            if self.exp is not None:
                self.instrument = self.exp[0:3]

        if self.exp and self.monshmserver:
            calibDir = '/reg/d/psdm/{:}/{:}/calib'.format(self.instrument,self.exp)
            print('setting calibDir', self.exp, calibDir)
            psana.setOption('psana.calib-dir', calibDir)

        self._set_user_dir()
        self._set_xarray_dir()

        try:
            self.info = psutils.exp_info(self.exp)
        except:
            print('Cannot get exp_info for {:}'.format(self.exp))
            self.info = {}

        return data_source

    def _get_user_dir(self, user=None, base_path=None):
        """
        Get user dir. Default is in results (or res for older experiments) of 
        experiment folder.
        """
        import psutils

        if not base_path:
            base_path = '/reg/d/psdm/{:}/{:}/results'.format(self.instrument,self.exp)
            if not os.path.isdir(base_path):
                base_path = '/reg/d/psdm/{:}/{:}/res'.format(self.instrument,self.exp)

        if not user:
            try:
                user = psutils.get_user()
            except:
                user = 'default'

        return os.path.join(base_path, user)

    def _set_user_dir(self, user=None, base_path=None):
        """
        Set the path of the PyDataSource configuration directory
        """
        self.user_dir = self._get_user_dir(user=user, base_path=base_path)

        return self.user_dir

    def _set_xarray_dir(self, path=None):
        """
        Set the path of the xarray data directory
        """
        if not path:
            path = '/reg/d/psdm/{:}/{:}/scratch/nc'.format(self.instrument,self.exp)
        
        if not os.path.isdir(path):
            try:
                os.mkdir(path)
            except:
                return 

        self.xarray_dir = path

        return path

    def show_info(self):
        print('< {:}: {:} >'.format(self.__class__.__name__, self.data_source))
        for attr in sorted(self._exp_defaults):
            print('{:20} {:}'.format(attr, getattr(self, attr)))

    def __str__(self):
        return self.data_source

    def __repr__(self):
        self.show_info()
        return '< {:}: {:} >'.format(self.__class__.__name__, self.data_source)


