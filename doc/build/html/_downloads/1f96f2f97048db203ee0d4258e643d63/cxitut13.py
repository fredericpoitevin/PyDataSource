from __future__ import print_function
import sys
import time

exp = 'cxitut13'

def DataSource(exp=exp, run=30, 
        save_config=True, publish=False, **kwargs):
    """Load event data source and configure for cxitut13 data.
    """
    import PyDataSource
    import numpy as np
    ds = PyDataSource.DataSource(exp=exp,run=run, **kwargs)
    evt = ds.events.next()
    
    if run in [10,11]:
        evt.DscCsPad.add.parameter(gain=1/23.)
        bins = np.arange(-20,150,1)*evt.DscCsPad.gain
        evt.DscCsPad.add.histogram('calib', bins=bins, gain=evt.DscCsPad.gain, 
                name='photon_hist', unit='photons', 
                doc='Front CsPad Gain corrected histogram', 
                publish=publish)
        evt.DscCsPad.add.count('calib', limits=(12,10000), 
                gain=evt.DscCsPad.gain, 
                name='photon_count', unit='photons', 
                doc='Front CsPad Photon Count')
        evt.DscCsPad.add.stats('calib')
        for i in range(32):
            evt.DscCsPad.add.roi('calib',sensor=i,graphical=False)
            sensor = 'sensor{:}'.format(i)
            evt.DscCsPad.add.histogram(sensor, bins=bins, gain=evt.DscCsPad.gain, 
                    name=sensor+'_photon_hist', unit='photons', 
                    doc='CsPad Gain corrected histogram', 
                    publish=publish)


    elif run == 20:
        evt.Sc2Questar.next().add.stats('raw')
        evt.DsaCsPad.next().add.stats('calib')
        evt.DsaCsPad.add.parameter(gain=1/27.7)
        bins = np.arange(-20,150,1)*evt.DsaCsPad.gain
        evt.DsaCsPad.add.histogram('calib', bins=bins, gain=evt.DsaCsPad.gain, 
                name='photon_hist', unit='photons', 
                doc='CsPad Gain corrected histogram', 
                publish=publish)
        evt.DsaCsPad.add.count('calib', limits=(12,10000), 
                gain=evt.DsaCsPad.gain, 
                name='photon_count', unit='photons', 
                doc='CsPad Photon Count')
        for i in range(32):
            evt.DsaCsPad.add.roi('calib',sensor=i,graphical=False)
            sensor = 'sensor{:}'.format(i)
            evt.DsaCsPad.add.histogram(sensor, bins=bins, gain=evt.DsaCsPad.gain, 
                    name=sensor+'_photon_hist', unit='photons', 
                    doc='CsPad Gain corrected histogram', 
                    publish=publish)


    elif run == 30:
        evt.Sc2Imp.add.module('impbox')
        evt.Sc2Imp.next()
        attrs={ 'filter': evt.Sc2Imp.filter}
        evt.Sc2Imp.add.stats('filtered', doc='Filted waveforms', attrs=attrs) 
        evt.Acqiris.add.module('acqiris')
        evt.Acqiris.add.stats('waveform')

    if save_config:
        ds.save_config()

    return ds

def to_xarray(ds, build_html=True, default_stats=False, **kwargs):
    """
    Make xarray object
    """
    x = ds.to_hdf5(default_stats=default_stats, **kwargs) 
    if build_html:
        b = build_run(x, **kwargs)
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
        b.add_detstats()
        b.add_correlations()
        attrs = [a for a in x.keys() if (a.endswith('_count') or a.endswith('_sum')) 
                    and len(x[a].dims) == 1 and 'time' in x[a].dims]
        if 'PhotonEnergy' in x:
            attrs.append('PhotonEnergy')
        if 'Gasdet_post_atten' in x:
            attrs.append('Gasdet_post_atten')
        b.add_detector(attrs=attrs, catagory='Detector Count', confidence=0.1)
 
    if run in [10,11]:
        variables = ['DscCsPad_photon_hist']
        b.add_summary(variables)
        b.add_detector('DscCsPad')
    
    elif run == 20:
        variables = ['DsaCsPad_photon_hist']
        b.add_summary(variables)
        b.add_detector('DsaCsPad')
        if 'DsaCsPad_corr_stats' in x:
            b.add_stats('DsaCsPad_corr_stats')
        if 'Sc2Questar_corr_stats' in x:
            b.add_stats('Sc2Questar_corr_stats')

    else:
        if auto:
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



