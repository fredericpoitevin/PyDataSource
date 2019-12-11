from __future__ import print_function
from __future__ import absolute_import
import time
import argparse

def initArgs():
    """Initialize argparse arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("attr", help='Input')
    parser.add_argument("option", nargs='?', default=None,
                        help='Optional attribute')
    parser.add_argument("-b", "--build", type=str,
                        help='Build html report')
    parser.add_argument("-e", "--exp", type=str,
                        help='Experiment')
    parser.add_argument("-r", "--run", type=str,
                        help='Run')
    parser.add_argument("-i", "--instrument", type=str,
                        help='Instrument')
    parser.add_argument("--nchunks", type=int,  
                        help='total number of chunks')
    parser.add_argument("--ichunk", type=int,  
                        help='chunk index')
    parser.add_argument("--det", type=str,
                        help='Cmma separated detector alias names used for plot Input method')
    parser.add_argument("--batchuser", type=str,
                        help='Batch job user name')
    parser.add_argument("--alert", type=str,
                        help='Alert emails (comma separated no spaces)')
    parser.add_argument("-c", "--config", type=str,
                        help='Config File')
    parser.add_argument("-s", "--station", type=int,
                        help='Station')
    parser.add_argument("-n", "--nevents", type=int,
                        help='Number of Events to analyze')
    parser.add_argument("-M", "--max_size", type=int,
                        help='Maximum data array size')
    parser.add_argument("--ffb", action="store_true", default=False,
                        help='Use ffb.')
    parser.add_argument("--keep_chunks", action="store_true", default=False,
                        help='Keep individual chunked files after merging.')
    #parser.add_argument("--make_summary", action="store_true", default=False,
    #                    help='Make summary for array data.')
    return parser.parse_args()


if __name__ == "__main__":
    # Make sure to use Agg matplotlib 
    import matplotlib as mpl
    mpl.use('Agg')
    time0 = time.time()
    args = initArgs()
    attr = args.attr
    exp = args.exp
    run = int(args.run.split('-')[0])
    ffb = args.ffb
    #print '{:} Run {:}'.format(exp, run)
    if attr in ['build']:
        from .build_html import Build_html
        from . import h5write
        x = h5write.open_h5netcdf(exp=exp,run=run)
        print(x)
        b = Build_html(x, auto=True)
   
    elif attr in ['beam_stats']:
        # Build beam stats to check for dropped shots and off-by-one
        from .beam_stats import get_beam_stats
        if run > 100000:
            from . import psutils
            run = psutils.get_run_from_id(run, exp) 
        format_str = "get_beam_stats(exp={:},run={:}, from_name={:}, to_name={:}, wait=True, timeout=False"
        if ffb:
            format_str += ", ffb=True"
        print(format_str.format(exp,run,args.batchuser,args.alert))
        xdrop = get_beam_stats(exp=exp,run=run, ffb=ffb,
                from_name=args.batchuser, to_name=args.alert, wait=True, timeout=False)

    else:
        from . import PyDataSource
        ds = PyDataSource.DataSource(exp=exp,run=run)
        if args.build:
            build_html=args.build
        else:
            build_html='basic'

        if attr == 'summary':
            print(ds.configData)
            es = ds.exp_summary
            if es:
                try:
                    b = es.to_html()
                except:
                    print('Cannot build epics run summary report')
               
                try:
                    b = es.build_drop_stats()
                except:
                    print('Cannot build beam ststs run summary report')

            else:
                print('Failed to load or generate exp_summary')

        if attr == 'epics':
            print(ds.configData)
            es = ds.exp_summary
            if es:
                es.to_html()
            else:
                print('Failed to load or generate exp_summary')
        if attr == 'config':
            print(ds._get_config_file())
        if attr == 'run':
            print(ds.data_source.run)
        if attr == 'runstr':
            print('Run{:04}'.format(ds.data_source.run))
        if attr == 'configData':
            print(ds.configData.show_info())
        if attr in ['steps','nsteps']:
            print(ds.configData.ScanData.nsteps)        
        if attr in ['events','nevents']:
            print(ds.nevents)        
        if attr in ['scan']:
            print(ds.configData.ScanData.show_info())
        if attr in ['plot']:
            dets = ds.configData._config_srcs.keys()
            detstr = args.det
            if detstr:
                try:
                    dets = detstr.split(',')
                except:
                    print('Cannot parse args.det')

            if dets:
                evt = next(ds.events)
                for alias, detector in ds._detectors.items():
                    if '.' in alias:
                        alias, attr = alias.split('.')
                    else:
                        attr = None
                    if alias not in dets:
                        continue
                    try:
                        while alias not in evt._attrs:
                            next(evt)
                        try:
                            detector._update_xarray_info()
                        except:
                            pass
                        if detector._pydet.__module__ == 'Detector.AreaDetector':
                            if not attr:
                                attr = 'image'
                            try:
                                detector.add.psplot(attr)
                            except:
                                print('Cannot make {:} {:} plot -- trying raw data'.format(alias, attr))
                                attr = 'raw'
                                detector.add.psplot(raw)
                        elif detector._pydet.__module__ in ['Detector.WFDetector','Detector.GenericWFDetector']:
                            try:
                                if not attr:
                                    if detector._det_config['property'].get('filtered'):
                                        attr = 'filtered'
                                    else:
                                        attr = 'waveform'
                                nwfs = detector.waveform.shape[0]
                                for sensor in range(nwfs):
                                    detector.add.roi(attr,sensor=sensor,name='{:}_Ch{:}'.format(attr,sensor+1),publish=True)
                            except:
                                print('Cannot make {:} waveform plot'.format(alias))
                        else:
                            print('No plot available for {:} {:}'.format(alias, detector._pydet.__module__))

                    except:
                        print('Cannot make {:} plot for {:}'.format(alias, detector._pydet.__module__))
                if args.nevents:
                    nevents = args.nevents
                else:
                    nevents = -1
                time.sleep(2)
                
                if nevents > 0:
                    get_new = 'y'
                    while get_new in ["y","Y"]:
                        evt.monitor(nevents)
                        get_new = raw_input("Display more events -- 'y' for {:} events or enter number of events?\n".format(nevents)) 
                        try:
                            get_int = int(get_new)
                        except:
                            get_int = False
                        if get_int:
                            nevents = get_int
                            get_new = "y"
                else:
                    evt.monitor(nevents)
                    get_new = raw_input("Run Finished - Reload ?\n".format(nevents))
                    while get_new in ["y","Y"]:
                        ds.reload()
                        evt.monitor(nevents)
                        get_new = raw_input("Run Finished - Reload ?\n".format(nevents))

        if attr in ['mpi']:
            from .h5write import to_hdf5_mpi
            print(ds.configData)
            print('to hdf5 with mpi {:}'.format(args))
            if args.config and args.config not in ['auto', 'default']:
                print('Loading config: {:}'.format(args.config))
                ds.load_config(file_name=args.config)
            else:
                print('Auto config')
                ds.load_config()

            if not args.keep_chunks:
                print('Cleanup chunked files')
                cleanup = True
            else:
                print('DO NOT Cleanup chunked files after merging')
                cleanup = False

            x = to_hdf5_mpi(ds, nevents=args.nevents, nchunks=args.nchunks, 
                    build_html=build_html, cleanup=cleanup)

        if attr in ['batch']:
            from .h5write import write_hdf5
            if args.config:
                if args.config in ['auto', 'default']:
                    config = ds._get_config_file()
                    print('Loading default config: {:}'.format(config))
                    ds.load_config()
                else:
                    print('Loading config: {:}'.format(args.config))
                    ds.load_config(file_name=args.config)
            
            #write_hdf5(ds)
            x = ds.to_hdf5(build_html=build_html) 
            print(x)
        
        if attr in ['xarray']:
            print('to_xarray')
            if args.config:
                if args.config in ['auto', 'default']:
                    config = ds._get_config_file()
                    print('Loading default config: {:}'.format(config))
                    ds.load_config()
                else:
                    print('Loading config: {:}'.format(args.config))
                    ds.load_config(file_name=args.config)
            
            print('to_xarray...')
            x = ds.to_xarray(save=True, 
                             nevents=args.nevents, 
                             nchunks=args.nchunks, 
                             ichunk=args.ichunk,
                             build=args.build)
            print(x)

        if attr in ['test']:
            from .h5write import write_hdf5
            print(ds.configData.__repr__())
            print('-'*80)
            print(ds._device_sets)
            print('-'*80)
            if args.config:
                if args.config in ['auto', 'default']:
                    config = ds._get_config_file()
                    print('Loading default config: {:}'.format(config))
                    ds.load_config()
                else:
                    print('Loading config: {:}'.format(args.config))
                    ds.load_config(file_name=args.config)
            
            print('-'*80)
            print(ds._device_sets)
            print('-'*80)
            #write_hdf5(ds)
            x = ds.to_hdf5(nevents=100) 
            print(x)
 

    #print 'Total time = {:8.3f}'.format(time.time()-time0)

