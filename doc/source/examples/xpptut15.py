from __future__ import print_function
import sys
import time

exp = 'xpptut15'

def DataSource(exp=exp, run=310, 
        save_config=True, publish=False, **kwargs):
    """Load event data source and configure for cxitut13 data.
    """
    import PyDataSource
    import numpy as np
    ds = PyDataSource.DataSource(exp=exp,run=run, **kwargs)
    evt = next(ds.events)
    
    if run in [310]:
        next(evt.OPAL1)
        evt.OPAL1.add.projection('corr','x')
        evt.OPAL1.add.projection('corr','y')
        evt.OPAL1.add.count('corr')
    elif run in [300]:
        next(evt.OPAL1)
        evt.OPAL1.add.projection('corr','x')
        if publish:
            evt.OPAL1.add.psplot('corr')
        evt.OPAL1.add.count('corr')
        evt.ACQ4.add.module('acqiris')
        evt.pnccdFront.add.psplot('image')
    elif run in [260]:
        next(evt.cspad2x2_diff)
        if publish:
            evt.cspad2x2_diff.add.psplot('image')
        evt.cspad2x2_diff.add.histogram('calib',bins=range(-15,275),publish=publish)
        evt.cspad2x2_diff.add.count('calib')
        evt.cspad2x2_diff.add.projection('image', 'x')
        evt.cspad2x2_diff.add.projection('image', 'y')
        next(evt.epix100a_diff)
        if publish:
            evt.epix100a_diff.add.psplot('image')
        evt.epix100a_diff.add.histogram('calib',bins=range(-15,275),publish=publish)
        evt.epix100a_diff.add.count('calib')
        evt.epix100a_diff.add.projection('image', 'x')
        evt.epix100a_diff.add.projection('image', 'y')

    elif run in [250]:
        next(evt.opal_1)
        evt.opal_1.add.psplot('corr')
        evt.opal_1.add.projection('corr','x')
        evt.opal_1.add.count('corr')
    elif run in [240]:
        next(evt.rayonix)
        evt.rayonix.add.psplot('corr')
        evt.rayonix.add.projection('corr','x')
        evt.rayonix.add.count('corr')
    elif run in [220]:
       next(evt.opal_1)
       evt.opal_1.add.psplot('corr')
       evt.opal_1.add.projection('corr','x')
       evt.opal_1.add.count('corr')
       evt.XppMon_Pim0.add.psplot('corr')
    elif run in [230]:
        next(evt.cspad)
        evt.cspad.add.psplot('image')
        evt.cspad.add.projection('corr','x')
        evt.cs140_0.add.psplot('corr')
        evt.cs140_0.add.projection('corr','x')
    if save_config:
        ds.save_config(path='/reg/neh/home/mkarra/RunSummary/xpptut15/')

    return ds 

def to_xarray(ds, build_html=True, default_stats=True, **kwargs):
    """
    Make xarray object
    """
    x = ds.to_hdf5(default_stats=default_stats, **kwargs) 
    if build_html:
        try:
            b = build_run(x, **kwargs)
        except:
            print('Cannot build_run')
    return x

def open_hdf5(exp=exp, run=30):
    """
    Open hdf5 file with netcdf4 convention using builtin xarray engine h5netcdf.
    """
    import PyDataSource
    return PyDataSource.open_h5netcdf(exp=exp,run=run)

def make_cuts(x, **kwargs):
    """
    Make experiment specific cuts
    """
    return x

def build_run(x, ioff=True, auto=False, **kwargs):
    """Make html web page from xarray summary dataset.
    """
    # do not perform interactive plotting
    import matplotlib.pyplot as plt
    from PyDataSource.build_html import Build_html
    if ioff:
        plt.ioff()

    x = make_cuts(x)
    run = int(x.run.values[0])

    b = Build_html(x, auto=auto, **kwargs)
    if not auto:
        try:
            b.add_detstats()
        except:
            pass

        try:
            b.add_correlations()
        except:
            pass

        attrs = [a for a in x.keys() if (a.endswith('_count') or a.endswith('_sum')) 
                    and len(x[a].dims) == 1 and 'time' in x[a].dims]
        if 'PhotonEnergy' in x:
            attrs.append('PhotonEnergy')
        if 'Gasdet_post_atten' in x:
            attrs.append('Gasdet_post_atten')
        b.add_detector(attrs=attrs, catagory='Detector Count', confidence=0.1)
 
    if run in [310]:
        b.add_detector('OPAL1')

    else:
        b.add_all(**kwargs)

    b.to_html()
    
    return b

def main():
    """Main script to create run summary.
    """
    import matplotlib as mpl
    mpl.use('Agg')
    from PyDataSource import initArgs
    time0 = time.time()
    print(time0)
    args = initArgs.initArgs()
    print(args)
    attr = args.attr
    exp = args.exp
    run = int(args.run.split('-')[0])
    print(exp, run)
    if attr in ['build']:
        from PyDataSource import h5write
        x = h5write.open_h5netcdf(exp=exp,run=run)
        print(x)
        b = build_run(x)
 
    else:
        ds = DataSource(exp=exp,run=run)
        if attr == 'epics':
            print(ds.configData)
            es = ds.exp_summary
            if es:
                es.to_html()
            else:
                print('Failed to load or generate exp_summary')

        elif attr in ['batch']:
            from PyDataSource import write_hdf5
            x = to_xarray(ds) 
            print(x)
            b = build_run(x)
        else:
            print(ds.configData)

if __name__ == "__main__":
    sys.exit(main())



