from __future__ import print_function
from __future__ import division
import time
import sys, os
import traceback

import PyDataSource
from PyDataSource import models, h5write, plotting, psxarray
from pylab import *
import seaborn as sns


atten_thicks=[0,33,31,55,45,140,315,1280,2560,5120]
sicoeff=107

exp='mfx11116'
instrument=exp[0:3]
eventCodes=[40]
runs = [8, 9, 10, 11] 
#data_path = '/reg/d/psdm/{:}/{:}/'.format(instrument,exp)
data_path = os.path.join('/reg/d/psdm/',instrument,exp)
results_path = os.path.join(data_path,'results')
nc_path = os.path.join(data_path,'scratch/nc/')
#exp_path = os.path.join('/reg/neh/operator/',instrument+'opr','experiments',exp)
exp_path = results_path 

def get_atten(x, thicks=atten_thicks,sicoeff=sicoeff): 
    thick=0.
    for i in range(1,7):
        attr = 'att{:02}_in'.format(i)
        vals = x[attr].values
        thick += vals*thicks[i]
    
    return np.exp(-thick*1e-4*sicoeff)

def fit_attenuators(x, attr='BeamMonitor_intensity', name=None,  
            cut='fit_cut', sicoeff=107, PhotonEnergy0=9000., 
            step0=None, fit_intercept=False):
    if not name:
        name = attr+'_thick_calc'

    attr_corr = attr+'_corr'
    if not step0:
        step0 = int(x.step.min())
    x[attr_corr] = (['time'], x[attr]/x.Gasdet_post_atten) 
    x[attr_corr].attrs['norm'] = float(x.where(x[cut]).where(x.step == step0)[attr_corr].dropna('time').mean())
    x[attr_corr].attrs['doc'] = 'BeamMonotor_intensity normalized to Gasdet_post_atten'

    x[name] = (['time'], 
            -np.log(x[attr_corr]/x[attr_corr].norm)/(x.PhotonEnergy/PhotonEnergy0)**2/1e-4/sicoeff)
#            -np.log(x.BeamMonitor_corr/x.BeamMonitor_corr.norm)/(x.PhotonEnergy/9000.)**2/1e-4/107.)
    
    x[name].attrs['doc'] = 'Attanuator thickness calculated from Gasdet corrected BeamMonitor_intensity'
    x[name].attrs['unit'] = 'um'
    
    print('Fitting attenuators based on BeamMonitor data')
    xvars = ['att01_in', 'att02_in', 'att03_in', 'att04_in', 'att05_in', 'att06_in']
    #x = RunSummary.psxarray.linear_model(x, 'BeamMonitor_thick_calc', xvars, 
    x = models.linear_model(x, name, xvars, 
            fit_intercept=fit_intercept, cut=cut, quiet=False)
    x = h5write.resort(x)

    name_model = name+'_model'
    xft = x[name_model]
    thick = 0
    for avar, val in zip(xft.variables,xft.coef_):
        x[avar].attrs['thick'] = val
        x[avar].attrs['unit'] = 'um'
        thick += x[avar]*val
    
    x['atten_thick'] = (['time'], thick)
    x['atten_thick'].attrs['unit'] = 'um'
    x['atten_thick'].attrs['doc'] = 'Corrected total Si attenuator thickness'
    x['atten_corr'] = (['time'], 
            np.exp(-x.atten_thick.values*(x.PhotonEnergy/PhotonEnergy0)**2*1e-4*sicoeff))
    x['atten_corr'].attrs['doc'] = 'Approx Corrected total Si attenuation'

    return x

def model_plot(x):
    a = psxarray.to_summary(x)
    plotting.xy_ploterr(a, 'BeamMonitor_thick_norm', 'step')
    plt.show()

def make_all(runs=runs, **kwargs):
    for run in runs[1:]:
        make_summary(run=run, **kwargs)

def open(run=None, exp=exp, summary=False, chunk=False, **kwargs):
    """Load run.
    """
    if summary:
        try:
            x = psxarray.open_h5netcdf(exp=exp, file_base='run{:04}_sum'.format(run), chunk=False)
            return x
        except:
            pass

    try:
        x = psxarray.open_h5netcdf(exp=exp, run=run, chunk=chunk)
        try:
            x = make_cuts(x)
        except:
            print('Cannot make cuts')
            return x
        
        x = cleanup(x)
        if summary:
            x = psxarray.to_summary(x, save_summary=True)
        return x
    except:
        print('Cannot load run {:} for {:}'.format(run, exp))
        traceback.print_exc()

def to_hdf5(x, **kwargs):
    try:
        psxarray.to_h5netcdf(x, **kwargs)
    except:
        print('Cannot save to h5')
        traceback.print_exc()

def DataSource(run=None, exp=exp, publish=False, projection=False, save_config=False, **kwargs):
    """Make xarray summary dataset.
    """
    ds = PyDataSource.DataSource(exp=exp, run=run, **kwargs)
    evt = next(ds.events)
   
    if 'BeamMonitor' in ds._detectors:
        next(evt.BeamMonitor)
        evt.BeamMonitor.add.module('wave8')
        evt.BeamMonitor.add.stats('waveforms')
        #evt.BeamMonitor.add.module('wave8',path='/reg/neh/home/koglin/psana/current/PyDataSource/src/')

    if 'CsPad' in ds._detectors:
        next(evt.CsPad)
        evt.CsPad.add.histogram('corr',bins=list(range(-15,300)))
        evt.CsPad.add.count('corr')  
        evt.CsPad.add.stats('corr')
#        evt.CsPad.add.projection('corr','r') 

    pvs = [
           'atten_trans1',
           'atten_trans3',
           'atten_s01',
           'atten_s02',
           'atten_s03',
           'atten_s04',
           'atten_s05',
           'atten_s06',
           'atten_s07',
           'atten_s08',
           'atten_s09',
          ]
    
    epics_attrs = [
                   'BeamMonitor_diode_x',
                   'BeamMonitor_diode_y',
                   'BeamMonitor_diode_target',
#                   'Dg2Pim_focus',
#                   'Dg2Pim_y',
#                   'Dg2Pim_zoom',
                   ]

    ds.xarray_kwargs.update(pvs=pvs, epics_attrs=epics_attrs, max_size=100001)
    ds.save_config()

    return ds

def make_plots(ds):
    evt = next(ds.events)
   
    next(evt.BeamMonitor)
    evt.BeamMonitor.add.module('wave8')
    #evt.BeamMonitor.add.module('wave8',path='/reg/neh/home/koglin/psana/current/PyDataSource/src/')
    
    evt.InlineCamera.add.roi('raw', roi=([935, 1075], [1125, 1400]), projection=True, name='proj', publish=True)


def to_xarray(ds=None, max_size=100001, **kwargs):
    """Experiment specific setup for xarray data.
    """
    if not ds:
        ds = DataSource(**kwargs)

    code_flags={
#            'Opal': [42], 
            }
    
    #pvs = ['pi1_x','pi1_fine_x'] 
    x = psxarray.to_xarray(ds, max_size=max_size, 
            pvs=pvs, epics_attrs=epics_attrs, code_flags=code_flags, **kwargs)
    
    #x.attrs['correlation_variables'] = ['Acqiris_LaserDiode_time']
    try:
        x = make_cuts(x)
    except:
        print('Cannot make cuts')
        #traceback.print_exc()
    
    try:
        x = cleanup(x)
    except:
        print('Cannot cleanup data')
        #traceback.print_exc()

    return x

def summary_html(a, ioff=True, exp_path=exp_path):
    """RunSummary built from statistics of each step in the run. 
    """
    import pandas as pd
    if 'stat' not in a.dims:
        a = psxarray.to_summary(a)

    if ioff:
        plt.ioff()

    if not a.attrs.get('scan_variables'):
        a.attrs['scan_variables'] = ['atten_corr']

    path=os.path.join(exp_path, 'RunSummary')
    print(path)
    b = Build_html(a, path=os.path.join(exp_path, 'RunSummary'))


#    b.add_xy_ploterr('CsPad_calib_count_thick_calc_residuals', xaxis='step', logx=False, logy=False, catagory='Summary')
    b.add_xy_ploterr('BeamMonitor_intensity_thick_calc_residuals', xaxis='step', logx=False, logy=False, catagory='Summary', 
            table=pd.DataFrame(index=a.BeamMonitor_intensity_thick_calc_model.variables, data={'thick':a.BeamMonitor_intensity_thick_calc_model.coef_}), 
            text='Attenuator fit results using BeamMonitor')

    xvars = ['att01_in', 'att02_in', 'att03_in', 'att04_in', 'att05_in', 'att06_in']
    df = a.reset_coords()[xvars].sel(stat='mean').drop('stat').to_dataframe()
    add_atten(b, xaxis='step',logx=False, logy=False, table=df, text='Attenuator sequence')


#    b.add_xy_ploterr('BeamMonitor_intensity_renorm', xaxis='step', logx=False, logy=False, catagory='Summary')
#b.add_xy_ploterr('BeamMonitor_peaks_renorm', xaxis='step', logx=False, logy=True, catagory='Summary')
#    b.add_xy_ploterr('BeamMonitor_intensity', logx=True, logy=True)
#    b.add_xy_ploterr('BeamMonitor_peaks', logx=True, logy=True, catagory='Summary')

    b.add_xy_ploterr('BeamMonitor_intensity', xaxis='atten_corr', logx=True, logy=True, catagory='Corrected')
#    b.add_xy_ploterr('BeamMonitor_peaks', xaxis='atten_corr',logx=True, logy=True, catagory='Corrected')

    if 'atten_trans1' in a:
        b.add_xy_ploterr('BeamMonitor_intensity', xaxis='atten_trans1', logx=True, logy=True)
#    b.add_xy_ploterr('BeamMonitor_peaks', xaxis='atten_trans1', logx=True, logy=True)
#    b.add_xy_ploterr('BeamMonitor_intensity_norm', xaxis='atten_trans1', logx=True, logy=True)
#    b.add_xy_ploterr('BeamMonitor_peaks_norm', xaxis='atten_trans1', logx=True, logy=True)

    b.add_xy_ploterr('Gasdet_post_atten', xaxis='step', catagory='Systematics', logx=False, logy=False)
    b.add_xy_ploterr('EBeam_ebeamCharge', xaxis='step', catagory='Systematics', logx=False, logy=False)
    b.add_xy_ploterr('EBeam_ebeamDumpCharge', xaxis='step', catagory='Systematics', logx=False, logy=False)
    b.add_xy_ploterr('EBeam_ebeamPhotonEnergy', xaxis='step', catagory='Systematics', logx=False, logy=False)
    
#    add_atten(b, xaxis='step', logx=False, catagory='Step')
#@    b.add_xy_ploterr('BeamMonitor_intensity_renorm', xaxis='step', logx=False, logy=False, catagory='Step')
#    b.add_xy_ploterr('atten_corr', xaxis='step', logx=False, logy=True, catagory='Step')
#    b.add_xy_ploterr('atten_thick', xaxis='step', logx=False, logy=False, catagory='Step')
#    b.add_xy_ploterr('BeamMonitor_intensity', xaxis='step', logx=False, logy=True, catagory='Step')
#    b.add_xy_ploterr('BeamMonitor_intensity_norm', xaxis='step', logx=False, logy=True, catagory='Step')
#    b.add_xy_ploterr('CsPad_calib_count', xaxis='step', logx=False, logy=True, catagory='Step')
    b.to_html()

    return b

def add_atten(self, xaxis='atten_trans1_axis', logx=True, logy=False, howto=None, catagory=None, text=None, table=None, **kwargs):
    x = self._xdat
    if not howto:
        howto = []
    if not catagory:
        catagory = 'Attenuators'

    x = x.swap_dims({'step':xaxis})
    plt_type = '{:} vs {:}'.format('Attenuator States', xaxis)   
    if logy:
        if logx:
            plt_type+=' log-log'
        else:
            plt_type+=' log scale'
    elif logx:
        plt_type+=' lin-log'
   
    plt.figure()
    plt.gca().set_position((.1,.2,.8,.7))
    for i in range(1,7):
        #c = x['atten_s{:02}'.format(i)].to_pandas().T
        c = x['att{:02}_in'.format(i)].to_pandas().T
        c = c.sort_index()
        p = ((c['mean'])*float(i)).plot(logx=logx,**kwargs)
        #p = ((2-x['atten_s{:02}'.format(i)].sel(stat='mean'))*float(i)).plot(**kwargs)
    
    print(catagory, plt_type)
    self.add_plot(catagory, plt_type, howto=howto)
    if table is not None:
        self.results[catagory]['table'].update({'Attenuators':{'DataFrame': table, 
                                                        'name': 'df_tbl',
                                                        'howto': [], 
                                                        'doc': text}})



    return p


def make_summary(run=None, exp=exp, ds=None, nevents=None, 
        ichunk=None, nchunks=24, eventCodes=[140], max_size=100001, 
        load=False, build=True, save=True, fit=True, **kwargs):
    """Make xarray summary dataset.
    """
    if not ds:
        ds = DataSource(run=run, exp=exp, **kwargs)

    x = None
    if load:
        x = open(run=run, exp=exp, **kwargs)

    if x is None:
        print('to_xarray', nevents, nchunks, ichunk,eventCodes, max_size)
        x = to_xarray(ds, nevents=nevents, nchunks=nchunks, ichunk=ichunk, 
                eventCodes=eventCodes, max_size=max_size, **kwargs)
        
        try:
            if ichunk is None:
                fit = False
            else:
                fit = fit
            x = make_cuts(x, fit=fit)
        except:
            print('Cannot make cuts')
            traceback.print_exc()
        
        try:
            x = cleanup(x)
        except:
            print('Cannot cleanup data')
            traceback.print_exc()

    if save:
        try:
            psxarray.to_h5netcdf(x)
        except:
            psxarray.write_file(x)
            print('Cannot save to h5 -- writing picke file instead')
            traceback.print_exc()

    if build and ichunk is None:
        try:
            b = build_html(x)
        except:
            print('Cannot build html')

    return x

def cleanup(x):
    for attr in ['BeamMonitor_data_f32','BeamMonitor_data_f64','BeamMonitor_data_size','BeamMonitor_data_u8']:
        if attr in x:
            del x[attr]

    return x

def make_cuts(x, fit=None, stepmax=55):
    """Make cuts and add calculated info.
    """
#    ltuycut = (x.EBeam_ebeamLTUPosY.values < 0.) & (x.EBeam_ebeamLTUAngY.values > 0.011)
#    ltuxcut = (x.EBeam_ebeamLTUPosX.values > -0.11) & (x.EBeam_ebeamLTUAngX.values > 0.03)
#    undxcut = (x.EBeam_ebeamUndAngX.values > 0)

    pkcurrcut = (x.EBeam_ebeamPkCurrBC1.values < 240.) & (x.EBeam_ebeamPkCurrBC2.values > 4000.)
    
    phasecut = np.array(
            (x.PhaseCavity_charge1.values > 150.) & \
            (x.PhaseCavity_charge2.values > 150.) & \
            (abs(x.PhaseCavity_fitTime1.values - float(x.PhaseCavity_fitTime1.mean())) < 2.) & \
            (abs(x.PhaseCavity_fitTime2.values - float(x.PhaseCavity_fitTime2.mean())) < 2.), dtype=byte)
    
    x.coords['PhaseCavity_cut'] = (['time'], phasecut)
  
    attr = 'FEEGasDetEnergy_f_11_ENRC'
    if attr in x:
        x['Gasdet_pre_atten'] = (['time'], x[attr].values)
        x['Gasdet_pre_atten'].attrs['doc'] = "Energy measurement before attenuation ({:})".format(attr)
        for a in ['unit', 'alias']: 
            try:
                x['Gasdet_pre_atten'].attrs[a] = x[attr].attrs[a]  
            except:
                pass


    attr = 'FEEGasDetEnergy_f_21_ENRC'
    if attr in x:
        x['Gasdet_post_atten'] = (['time'], x[attr].values)
        x['Gasdet_post_atten'].attrs['doc'] = "Energy measurement afeter attenuation ({:})".format(attr)
        for a in ['unit', 'alias']: 
            try:
                x['Gasdet_post_atten'].attrs[a] = x[attr].attrs[a]  
            except:
                pass

    gasdetcut_mJ = 0.5
    gasdetcut =  np.array(x.Gasdet_pre_atten.values > gasdetcut_mJ, dtype=byte)
    x.coords['Gasdet_cut'] = (['time'], gasdetcut)
    x.coords['Gasdet_cut'].attrs['doc'] = "Gas detector cut.  Gasdet_pre_atten > {:} mJ".format(gasdetcut_mJ)
    
    #damagecut = np.array(phasecut & gasdetcut & (x.EBeam_damageMask.values == 0), dtype=byte)
    damagecut = np.array(gasdetcut, dtype=byte)
    x.coords['Damage_cut'] = (['time'], damagecut)
    x.coords['Damage_cut'].attrs['doc'] = "Combined Gas detector, Phase cavity and EBeam damage cut"

    #fitcut = np.array(damagecut & np.array(x.step < 45) & np.array(x.step > 12) )
    fitcut = np.array(damagecut & np.array(x.step <= stepmax) )
    x.coords['fit_cut'] = (['time'], fitcut)
    x.coords['fit_cut'].attrs['doc'] = "Damage_cut and only up to step {:}".format(stepmax) 

    x['BeamMonitor_sum'] = x.BeamMonitor_peaks[:,0:4].sum(axis=1)
    x['BeamMonitor_sum'].attrs.update(x.BeamMonitor_peaks.attrs)
    x['BeamMonitor_sum'].attrs['doc'] = 'Sum of first 4 BeamMonotor diode peak intensities'
#    x['BeamMonitor_normalized'] = x.BeamMonitor_sum/x.Gasdet_post_atten
#    x['BeamMonitor_normalized'].attrs.update(x.BeamMonitor_peaks.attrs)
#    x['BeamMonitor_normalized'].attrs['doc'] = 'BeamMonotor_sum normalized to Gasdet_post_atten'

    add_index(x, 'Gasdet_post_atten')
    #if 'atten_trans1' in x and 'atten_trans1' not in x.coords:
    #    x = x.set_coords('atten_trans1')

    try:
        x['PhotonEnergy'] = x['EBeam_ebeamPhotonEnergy']
    except:
        pass

    # Change attenuator to in/out
    for i in range(1,10):
        attr = 'att{:02}_in'.format(i)
        pvalias = 'atten_s{:02}'.format(i)
        if pvalias in x:
            x[attr] = (['time'], np.array(2-x[pvalias], dtype=byte))
            x[attr].attrs['doc'] = 'Attenuator {:} In/Out'.format(pvalias)
            del x[pvalias]

    if 'atten_trans1' not in x:
        try:
            atten_trans1 = get_atten(x)
            x['atten_trans1'] = (['time'], atten_trans1)
            if len(set(atten_trans1)) > 2:
                x.attrs['scan_variables'] += ['atten_trans1']
        except:
            pass
    
    if 'atten_trans1' in x:
        x['PulseEnergy'] = x['Gasdet_post_atten']*x['atten_trans1']
        x['PulseEnergy'].attrs = x['Gasdet_post_atten'].attrs
        x['PulseEnergy'].attrs['doc'] = "Energy measurement normalized by attenuators"


#   x = thick_corr(x)
    if fit:
        # fit attenuators -- creates atten_corr corrected transmission
        x = fit_attenuators(x, 'BeamMonitor_intensity')
        x = fit_attenuators(x, 'CsPad_calib_count')        
 
        x['PulseEnergy_corr'] = x['Gasdet_post_atten']*x['atten_corr']
        x['PulseEnergy_corr'].attrs = x['Gasdet_post_atten'].attrs
        x['PulseEnergy_corr'].attrs['doc'] = "Energy measurement normalized by attenuators"
        x.attrs['scan_variables'] = ['atten_corr']
        
        attrs = ['CsPad_calib_count', 'BeamMonitor_intensity','BeamMonitor_peaks','BeamMonitor_sum']
        x = psxarray.normalize_data(x, attrs)
        x = psxarray.normalize_data(x, attrs, norm_attr='PulseEnergy_corr', name='renorm')
#        for attr in attrs:
#            x[attr+'_renorm'].attrs['doc'] = '{:} with atten thicks {:}'.format(attr, atten_thicks)
        
    return x

# stat math example:
#            num = x[attr].to_pandas().T
#            den = x[norm_attr].to_pandas().T
#            rat = num/den
#            rat['std'] = np.sqrt(num['std']**2/num['mean']**2+den['std']**2/den['mean']**2)*rat['mean']


def build_html(x, add_scatter=False, do_inline=False):
    """Make html web page from xarray summary dataset.
    """
    from PyDataSource.build_html import Build_html
    b = Build_html(x)

    #b.add_detector('EBeam',cut='Damage_cut')
    #b.add_detector('FEEGasDetector',cut='Damage_cut')
    #b.add_detector('PhaseCavity',cut='Damage_cut')
    b.add_summary(variables=['BeamMonitor_waveforms'],groupby=['Damage_cut'])
    if 'atten_trans1' in x:
        b.add_summary(variables=['BeamMonitor_waveforms'],groupby=['atten_trans1'])

    attr_groups = {
#                   'Undulator X': ['EBeam_ebeamUndAngX', 'EBeam_ebeamUndPosX'],
#                   'Undulator Y': ['EBeam_ebeamUndAngY', 'EBeam_ebeamUndPosY'],
#                   'LTU X': ['EBeam_ebeamLTUAngX', 'EBeam_ebeamLTUPosX'],
#                   'LTU Y': ['EBeam_ebeamLTUAngY', 'EBeam_ebeamLTUPosY'],
#                   'LTU 250 and 450': ['EBeam_ebeamLTU250', 'EBeam_ebeamLTU450'],
                   'Energy BC': ['EBeam_ebeamEnergyBC1', 'EBeam_ebeamEnergyBC2'],
                   'Charge': ['EBeam_ebeamCharge', 'EBeam_ebeamDumpCharge'],
                   'Peak Current': ['EBeam_ebeamPkCurrBC1', 'EBeam_ebeamPkCurrBC2'],
                   'XTCAV': ['EBeam_ebeamXTCAVAmpl', 'EBeam_ebeamXTCAVPhase'],
                   'Phase Cavity Charge': ['PhaseCavity_charge1', 'PhaseCavity_charge2'],
#                   'Phase Cavity Fit': ['PhaseCavity_fitTime1', 'PhaseCavity_fitTime2'],
                   'Gasdet': ['FEEGasDetEnergy_f_63_ENRC', 'Gasdet_pre_atten'],
                   }

    if add_scatter:
        attrs = [attr for attr,item in x.data_vars.items() if item.dims == ('time',)]
        df = x[attrs].where(x.Damage_cut == 1).to_array().to_pandas().T.dropna()
        group = 'Gasdet_post_atten_index'
        for catagory, grp_attrs in attr_groups.items():
            attrs = ['EBeam_ebeamPhotonEnergy', 'Gasdet_post_atten', 'BeamMonitor_sum']
            print(group, catagory)
            for attr in grp_attrs:
                attrs.append(attr)
            howto = []
            b.add_scatter(df, catagory=catagory, attrs=attrs, howto=howto, group=group)

    b.to_html()

    return b
 
def add_index(x, attr, name=None, nbins=8, bins=None, percentiles=None):

    if not bins:
        if not percentiles:
            percentiles = (arange(nbins+1))/float(nbins)*100.

        bins = np.percentile(x[attr].to_pandas().dropna(), percentiles)

    nbins = len(bins)

    if not name:
        name = attr+'_index'
  
    per = [percentiles[i-1] for i in np.digitize(x[attr].values, bins)]
    x[name] = (['time'], per)
   
  
def main():
    """Main script to create run summary.
    """
    time0 = time.time()
    print(time0)
    args = RunSummary.initArgs()
    print(args)
    
    import sys
    print(sys.path)

    import psana
    print(psana.__path__)
    import xarray
    print(xarray.__path__)
    import pylab
    print(pylab.__file__)

    x = make_summary(build=False, fit=False, **vars(args))


    return x

if __name__ == "__main__":
    sys.exit(main())


#def add_scatter(self, df, catagory='scatter', doc=None, group=None, attrs=None, howto=[]):
#    if not attrs:
#        attrs = df.keys()
#   
#    sns.set()
#    plt.rcParams['axes.labelsize'] = 10 
#    pltattrs = [attr for attr in attrs if attr not in [group]]
#
#    howto.append("sns.set()")
#    howto.append("pltattrs = {:}".format(pltattrs))
#    
#    if group:
#        plt_type = 'correlation with {:}'.format(group)
#        g = sns.pairplot(df, hue=group,
#                x_vars=pltattrs,y_vars=pltattrs,
#                size=2.5) 
#        howto.append("g = sns.pairplot(df, hue={:}, x_vars=pltattrs,y_vars=pltattrs,size=2.5)")
#    else:
#        plt_type = 'correlation'
#        g = sns.PairGrid(df, x_vars=pltattrs,y_vars=pltattrs, size=2.5) 
#        g = g.map_upper(plt.scatter)
#        g = g.map_lower(sns.kdeplot, cmap="Blues_d")
#        g = g.map_diag(sns.kdeplot, lw=3, legend=False)
#
#        howto.append("g = sns.pairplot(df, hue={:}, x_vars=pltattrs,y_vars=pltattrs,size=2.5)")
#
#
#    self.add_plot(catagory, plt_type, howto=howto)


#def get_specmax(x):
#    for i in range(4):
#        attr = 'VonHamos_spec{:}_y'.format(i+1)
#        wf = x[attr].mean(axis=0).values
#        maxval = max(wf)
#        idx = list(wf).index(maxval)
#        print i+1, attr, idx, maxval
 
#def make_df_cut(x, stdcut=0., dropna=True, mincut=True, maxcut=True,
#        icount=3, percentiles=[0.05,0.50,0.95],
#        group='CsPadKalpha_roi_count_index',
#        runmin=None,
#        runmax=None,
#        in_alias=['EBeam', 'PhaseCavity', 'FEEGasDetEnergy'], 
#        not_attrs=['EBeam_damageMask', 'FEEGasDetEnergy_f_21_ENRC', 'FEEGasDetEnergy_f_22_ENRC', 'FEEGasDetEnergy_f_63_ENRC', 'FEEGasDetEnergy_f_64_ENRC']):
#    """
#    """
#    cut_attrs = [attr for attr,item in x.data_vars.items() \
#            if (item.attrs.get('alias') in in_alias and attr not in not_attrs)]
#    attrs = [attr for attr,item in x.data_vars.items() if item.dims == ('time',)]
#
#    df = x[attrs].where(x.Damage_cut == 1).where(x.CsPadKalpha_roi_count_index>icount).to_array().to_pandas().T.dropna()
#    if runmin:
#        df = df.where(df.runnum >= runmin).dropna()
#    if runmax:
#        df = df.where(df.runnum <= runmax).dropna()
#    
#    df_tbl = df.describe(percentiles=percentiles).T.round({'count':0})
#    df_cuts = df_tbl.T
#    df = x[attrs].where(x.Damage_cut == 1).to_array().to_pandas().T.dropna()
#    print df_tbl
#
#    perlow = '{:}%'.format(int(min(percentiles)*100))
#    perhigh = '{:}%'.format(int(max(percentiles)*100))
#    for attr in cut_attrs:
#        if mincut:
#            df = df.where(df[attr] >= df_cuts[attr][perlow]-stdcut*df_cuts[attr]['std'])
#        if maxcut:
#            df = df.where(df[attr] <= df_cuts[attr][perhigh]+stdcut*df_cuts[attr]['std'])
#   
#        print df.dropna().shape[0], 'events after cut on ', attr
#
#    if dropna:
#        df = df.dropna()
# 
#    return df


#def xcorr(x):
#    """Correct 
#    """
#    xvars =  ['atten_s01','atten_s02','atten_s03','atten_s04','atten_s05','atten_s06']
#    x = x.reset_coords()
#    for attr in  [a for a in x.coords if a != 'time']:
#        x = x.drop(attr)
#    x['BeamMonitor_corr'] = (['time'], x.BeamMonitor_intensity/x.Gasdet_post_atten) 
#    x.BeamMonitor_corr.attrs['norm'] = x.where(x.atten_thick == 0)['BeamMonitor_corr'].dropna('time').mean()
#    x['BeamMonitor_thick_calc'] = (['time'], 
#            -np.log(x.BeamMonitor_corr/x.BeamMonitor_corr.norm)/(x.PhotonEnergy/9000.)**2/1e-4/107.)
#    x = x.where(x.Damage_cut).dropna('time')
##    x.attrs['attens'] =  ['atten_s01','atten_s02','atten_s03','atten_s04','atten_s05','atten_s06']
##    x.attrs['fit_attens'] =  ['atten_s01','atten_s02','atten_s03','atten_s04','atten_s05']
#    for attr, val in xattrs.items():
#        x.attrs[attr] = val
#    print xattrs
#    return resort(x)


#def linear_model(x, attr, xvars, fit_intercept=None, name=None, 
#        cut=None, residuals=True,
#        model='LinearRegression'):
#    """Make a linear model for attr based on xvars as free parameters.
#       Currently only model='LinearRegression' implmented.
#       Uses scikit-learn.
#    """
#    if model is not 'LinearRegression':
#        raise Exception("Currently only model='LinearRegression' implmented.")
#
#    from sklearn.linear_model import LinearRegression
#    import xarray as xr
#    import pandas as pd
#    lm = LinearRegression()  
#    if not name:
#        name = '{:}_model'.format(attr)
#
#    df_xvars0 = (2-x.reset_coords()[xvars]).to_dataframe()
#    if cut:
#        df_xvars = (2-x.reset_coords()[xvars]).where(x.reset_coords()[cut] == 1, drop=True).to_dataframe()
#        xdata = x[attr].where(x[cut] == 1, drop=True).data
#    else:
#        df_xvars = df_vars0
#        xdata = x[attr].data
#    
#    if fit_intercept is not None:
#        lm.fit_intercept = fit_intercept
#   
#    lm.fit(df_xvars, xdata)
#    ft = pd.DataFrame(zip(df_xvars.columns,lm.coef_), columns=['params','estimatedCoefficients'])
#
#    x[name] = xr.DataArray(lm.predict(df_xvars0), coords=[('time', df_xvars0.index)])
#    x[name].attrs['doc'] = 'LinearRegression scikit-learn model for {:} training data'.format(attr)
#    x[name].attrs['model'] = model
#    x[name].attrs['params'] = lm.get_params()
#    x[name].attrs['variables'] = xvars
#    x[name].attrs['coef_'] = lm.coef_
#    x[name].attrs['intercept_'] = lm.intercept_
#    x[name].attrs['score'] = lm.score(df_xvars, xdata)
#    if residuals:
#        if not isinstance(residuals, str):
#            residuals = '{:}_residuals'.format(attr)
#
#        x[residuals] = (['time'], x[attr]-x[name])
#
#    return x

#def xfit(x):
#    from sklearn.linear_model import LinearRegression
#    import pandas as pd
#    lm = LinearRegression()  
#    df_attens = (2-x.reset_coords()[x.attens]).to_dataframe()
#    lm.fit_intercept = False
#    lm.fit(df_attens, x.BeamMonitor_thick_calc.data) 
#    ft = pd.DataFrame(zip(df_attens.columns,lm.coef_), columns=['attens','estimatedCoefficients'])
#    x.attrs['coeff'] = ft
#    x['BeamMonitor_thick_model'] = (['time'], lm.predict(df_attens))
#    x['BeamMonitor_thick_norm'] = (['time'], x.BeamMonitor_thick_calc-x.BeamMonitor_thick_model)
#    return x


#def to_summary(x, dim='time', groupby='step', 
#        save_summary=False,
#        normby=None,
#        omit_list=['run', 'sec', 'nsec', 'fiducials', 'ticks'],
#        stats=['mean', 'std', 'min', 'max', 'count']):
#    """Summarize a run.
#    """
#    import xarray as xr
#    xattrs = x.attrs
#    data_attrs = {attr: x[attr].attrs for attr in x}
#
#    if 'Damage_cut' in x:
#        x = x.where(x.Damage_cut).dropna(dim)
#   
#    coords = [c for c in x.coords if c != dim and c not in omit_list and dim in x.coords[c].dims] 
#    x = x.reset_coords(coords)
#
#    if isinstance(normby, dict):
#        for norm_attr, attrs in normby.items():
#            x = normalize_data(x, attrs, norm_attr=norm_attr)
#    elif isinstance(normby, list):
#        x = normalize_data(x, normby)
#
#    if groupby:
#        x = x.groupby(groupby)
#
#    dsets = [getattr(x, func)(dim=dim) for func in stats]
#    x = xr.concat(dsets, stats).rename({'concat_dim': 'stat'})
#    for attr,val in xattrs.items():
#        x.attrs[attr] = val
#    for attr,item in data_attrs.items():
#        if attr in x:
#            x[attr].attrs.update(item)
#
#    x = resort(x)
#    if save_summary:
#        to_hdf5(x)
#
#    return x
#
#def resort(x):
#    coords = sorted([c for c in x.coords.keys() if c not in x.coords.dims])
#    x = x.reset_coords()
#    x = x[sorted(x.data_vars)]
#
#    for c in coords:                                                       
#        x = x.set_coords(c)
#
#    return x
#

#    bm = x.BeamMonitor_intensity.groupby('step').mean().values

#def plot_norm(x, thicks=atten_thicks,sicoeff=sicoeff, **kwargs):
#    (x.BeamMonitor_intensity/get_corr(x,thicks=thicks,sicoeff=sicoeff)).groupby('step').mean().to_pandas().T.plot(**kwargs)
#
#def thick_corr(x, thicks=atten_thicks,sicoeff=sicoeff,dim='time'):
#    for i in range(1,10):
#        attr = 'atten_s{:02}'.format(i)
#        if attr in x:
#            x[attr].attrs['thickness'] = thicks[i]
#
#    if dim == 'step':
#        thick = x[attr].sel(stat='mean').values*0
#    else:
#        thick = x[attr].values*0
#    for i in range(1,10):
#        attr = 'atten_s{:02}'.format(i)
#        if attr in x:
#            if dim == 'step':
#                vals = x[attr].sel(stat='mean').values
#            else:
#                vals = x[attr].values
#
#            thick += (2-vals)*x[attr].attrs['thickness']
#    print dim, thick
#    x['atten_thick'] = ([dim], thick)
#    x['atten_thick'].attrs['unit'] = 'um'
#    x['atten_thick'].attrs['doc'] = 'Corrected total Si attenuator thickness'
#    x['atten_corr'] = ([dim],np.exp(-x.atten_thick.values*1e-4*sicoeff))
#    x['atten_corr'].attrs['doc'] = 'Approx Corrected total Si attenuation'
#    return x



#

#def xy_ploterr(a, attr=None, xaxis=None, title='', desc=None, **kwargs):
#    """Plot summary data with error bars, e.g.,
#        xy_ploterr(x, 'MnScatter','Sample_z',logy=True)
#    """
#    if not attr:
#        print 'Must supply plot attribute'
#        return
#
#    if 'groupby' in kwargs:
#        groupby=kwargs['groupby']
#    elif 'step' in a.dims:
#        groupby='step'
#    else:
#        groupby='run'
#
#    run = a.attrs.get('run')
#    runstr = 'Run {:}'.format(run)
#    name = a.attrs.get('name', runstr)
#    if not title:
#        title = '{:}: {:}'.format(name, attr)
#
#    if not xaxis:
#        xaxis = a.attrs.get('scan_variables')
#        if xaxis:
#            xaxis = xaxis[0]
#
#    ylabel = kwargs.get('ylabel', '')
#    if not ylabel:
#        ylabel = a[attr].name
#        unit = a[attr].attrs.get('unit')
#        if unit:
#            ylabel = '{:} [{:}]'.format(ylabel, unit)
#
#    xlabel = kwargs.get('xlabel', '')
#    if not xlabel:
#        xlabel = a[xaxis].name
#        unit = a[xaxis].attrs.get('unit')
#        if unit:
#            xlabel = '{:} [{:}]'.format(xlabel, unit)
#    
#    if xaxis:
#        if 'stat' in a[xaxis].dims:
#            xerr = a[xaxis].sel(stat='std').values
#            a[xaxis+'_axis'] = ([groupby], a[xaxis].sel(stat='mean').values)
#            xaxis = xaxis+'_axis'
#        else:
#            xerr = None
#
#        a = a.swap_dims({groupby:xaxis})
#    else:
#        xerr = None
#
#    if desc is None:
#        desc = a[attr].attrs.get('doc', '')
#
#    ndims = len(a[attr].dims)
#    if ndims == 2:
#        c = a[attr].to_pandas().T
#        if xerr is not None:
#            c['xerr'] = xerr
#        c = c.sort_index()
#        
#        plt.figure()
#        plt.gca().set_position((.1,.2,.8,.7))
#        p = c['mean'].plot(yerr=c['std'],xerr=c.get('xerr'), title=title, **kwargs)
#        plt.xlabel(xlabel)
#        plt.ylabel(ylabel)
#        if desc:
#            plt.text(-.1,-.2, desc, transform=p.transAxes, wrap=True)   
# 
#        return p 
#    elif ndims == 3:
#        plt.figure()
#        plt.gca().set_position((.1,.2,.8,.7))
#        pdim = [d for d in a[attr].dims if d not in ['stat', groupby, xaxis]][0]
#        for i in range(len(a[attr].coords[pdim])):
#            c = a[attr].sel(**{pdim:i}).drop(pdim).to_pandas().T.sort_index()
#            p = c['mean'].plot(yerr=c['std'], **kwargs)
#
#        plt.xlabel(xlabel)
#        plt.ylabel(ylabel)
#        p.set_title(title)
#        if desc:
#            plt.text(-.1,-.2, desc, transform=p.transAxes, wrap=True)   
#
#        return p 
#    else:
#        print 'Too many dims to plot'

#
#def normalize_data(x, variables=[], norm_attr='PulseEnergy', name='norm', quiet=True):
#    """Normalize a list of variables with norm_attr [default = 'PulseEnergy']
#    """
#    if not variables:
#        variables = [a for a in get_correlations(x) if not a.endswith('_'+name)]    
#
#    for attr in variables:
#        aname = attr+'_'+name
#        try:
#            x[aname] = x[attr]/x[norm_attr]
#            x[aname].attrs = x[attr].attrs
#            try:
#                x[aname].attrs['doc'] = x[aname].attrs.get('doc','')+' -- normalized to '+norm_attr
#                units = x[attr].attrs.get('unit')
#                norm_units = x[norm_attr].attrs.get('unit')
#                if units and norm_units:
#                    x[aname].attrs['unit'] = '/'.join([units, norm_units])
#            except:
#                if not quiet:
#                    print 'cannot add attrs for', aname
#        except:
#            print 'Cannot normalize {:} with {:}'.format(attr, norm_attr)
#
#    x = RunSummary.psxarray.resort(x)
#    
#    return x

#def add_xy_ploterr(self, attr, xaxis=None, howto=None, catagory=None, **kwargs):
#    x = self._xdat
#    if not howto:
#        howto = []
#    if not catagory:
#        catagory = x[attr].attrs.get('alias', 'Summary')
#    
#    howto.append("Custom plot see mfx11116.py xy_ploterr method")
#    p = xy_ploterr(x, attr, xaxis=xaxis, **kwargs)
#    if not xaxis:
#        xaxis=x.scan_variables[0]
#    plt_type = '{:} vs {:}'.format(attr, xaxis)   
#    if kwargs.get('logy'):
#        if kwargs.get('logx'):
#            plt_type+=' log-log'
#        else:
#            plt_type+=' log scale'
#    elif kwargs.get('logx'):
#        plt_type+=' lin-log'
#
#    print catagory, plt_type
#    self.add_plot(catagory, plt_type, howto=howto)
    

