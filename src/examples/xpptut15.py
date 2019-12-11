from __future__ import print_function
import time
import sys, os
import traceback
from glob import glob

import PyDataSource
from PyDataSource import psxarray

import seaborn as sns
import pandas as pd
from pylab import *

exp=__name__
instrument=exp[0:3]

data_path = '/reg/d/psdm/{:}/{:}/'.format(instrument,exp)
results_path = data_path+'results/'
nc_path = data_path+'scratch/nc/'
exp_path = '/reg/neh/operator/{:}opr/userexpts/{:}/'.format(instrument,exp)

def to_summary(x, dim='time', 
        save_summary=False,
        stats=['mean', 'std', 'min', 'max', 'count']):
    """Summarize a run.
    """
    import xarray as xr
    if 'Damage_cut' in x:
        x = x.where(x.Damage_cut).dropna('time')
    a = {func: getattr(x, func)(dim=dim) for func in stats} 
    x = xr.concat(a.values(), a.keys()).rename({'concat_dim': 'stat'})
    if save_summary:
        to_hdf5(x)

    return x

def DataSource(run=None, exp=exp, publish=False, **kwargs):
    """Make xarray summary dataset.
    """
    ds = PyDataSource.DataSource(exp=exp, run=run, **kwargs)
    evt = ds.events.next()
    if 'Sc2Imp' in ds.configData.Sources._aliases:
        evt.Sc2Imp.next()
        evt.Sc2Imp.add.module('impbox') 
    
    if 'Sc2Epix' in ds.configData.Sources._aliases:
        evt.Sc2Epix.next()
        evt.Sc2Epix.add.roi('calib', roi=((0,700),(215,250)),projection=True,name='spec', publish=publish)
        evt.Sc2Epix.add.roi('calib', roi=((0,700),(115,150)),projection=True,name='back', publish=publish)

    if 'Sc2Inline' in ds.configData.Sources._aliases:
        evt.Sc2Inline.next()
        evt.Sc2Inline.add.roi('raw', roi=((335,755),(380,830)), name='roi', publish=publish)

    return ds


def to_xarray(ds=None, max_size=100000, **kwargs):
    """Experiment specific setup for xarray data.
    """
    if not ds:
        ds = DataSource(**kwargs)

    store_data = ['Sc2Inline_roi', 'Sc2Epix_spec', 'Sc2Epix_back']

    code_flags={
#            'Opal': [42], 
            }
    
    pvs =         [
                   'Sample_x',
                   'Sample_y',
                   'dia_trans1',
                   'dia_trans3',
                   'dia_s01',
                   'dia_s02',
                   'dia_s03',
                   'dia_s04',
                   'dia_s05',
                   'dia_s06',
                   'dia_s07',
                   'dia_s08',
                   'dia_s09',
                   ]
    
    epics_attrs = [
#                   'Dg2Pim_zoom',
                    'analyzer_x',
                    'analyzer_y',
                    'crystal_theta',
                    'Sample_z',
                    'Sample_pitch',
                   # 'dia_trans1',
                   # 'dia_trans3',
                    'pi2_x',
                    'pi2_y',
                    'pi2_z',
                   ]

    x = ds.to_xarray(max_size=max_size, 
            store_data=store_data, 
            pvs=pvs, epics_attrs=epics_attrs, code_flags=code_flags, **kwargs)
    
    try:
        make_cuts(x)
    except:
        print('Cannot make cuts')
        traceback.print_exc()
    
    try:
        cleanup(x)
    except:
        print('Cannot cleanup data')
        traceback.print_exc()

    return x

def make_summary(run=None, exp=exp, ds=None, nevents=None, 
        ichunk=None, nchunks=24, eventCodes=[], max_size=100001, 
        load=False, build=True, save=True, **kwargs):
    """Make xarray summary dataset.
    """
    if not ds:
        ds = DataSource(run=run, exp=exp, **kwargs)

    x = None
    if load:
        x = open(run=run, exp=exp, **kwargs)
        if x is not None:
            save = False

    if x is None:
        x = ds.to_xarray(nevents=nevents, nchunks=nchunks, ichunk=ichunk, 
                eventCodes=eventCodes, max_size=max_size, **kwargs)
        try:
            make_cuts(x)
        except:
            print('Cannot make cuts')
            traceback.print_exc()
        
        try:
            cleanup(x)
        except:
            print('Cannot cleanup data')
            traceback.print_exc()


    if save:
        to_hdf5(x)

    if build:
        b = build_html(x)

    return x

def to_hdf5(x, **kwargs):
    try:
        psxarray.to_h5netcdf(x, **kwargs)
    except:
        print('Cannot save to h5')
        traceback.print_exc()

def cleanup(x):
    """Cleanup xarray data as needed.  i.e., remove unnecessary data for specific runs
    """
    for attr in [attr for attr in x.keys() if attr.startswith('XCS_IPM_01_')]:
        del x[attr]
    for attr in [attr for attr in x.keys() if attr.startswith('XCS_IPM_03_')]:
        del x[attr]
    
    return

def make_cuts(x, gasdetcut_mJ = 2.0):
    """Make cuts and add calculated info.  Likely needs update for CXI.
    """
    try:
        x['PhaseCavity_charge'] = x.PhaseCavity_charge1 + x.PhaseCavity_charge2
        x.PhaseCavity_charge.attrs = x.PhaseCavity_charge1 
        x.PhaseCavity_charge.attrs['doc'] = 'UND:R02:IOC:16:BAT:Charge1 + Charge2 in pico-columbs.'
    except:
        pass

    # FEEGasDetEnergy_f_12_ENRC is duplicate measurement -- can average if desired 
    try:
        attr = 'FEEGasDetEnergy_f_11_ENRC'
        x['Gasdet_pre_atten'] = (['time'], x[attr].values)
        x['Gasdet_pre_atten'].attrs['doc'] = "Energy measurement before attenuation ({:})".format(attr)
        for a in ['unit', 'alias']: 
            try:
                x['Gasdet_pre_atten'].attrs[a] = x[attr].attrs[a]  
            except:
                pass

        # FEEGasDetEnergy_f_22_ENRC is duplicate measurement -- can average if desired 
        attr = 'FEEGasDetEnergy_f_21_ENRC'
        x['Gasdet_post_atten'] = (['time'], x[attr].values)
        x['Gasdet_post_atten'].attrs['doc'] = "Energy measurement afeter attenuation ({:})".format(attr)
        for a in ['unit', 'alias']: 
            try:
                x['Gasdet_post_atten'].attrs[a] = x[attr].attrs[a]  
            except:
                pass
    except:
        pass

    try:
        x['PulseEnergy'] = (['time'], x['Gasdet_post_atten'].values*x['dia_trans1'].values)
    except:
        x['PulseEnergy'] = (['time'], x['Gasdet_post_atten'].values*x.attrs['dia_trans1'])
    
    x['PulseEnergy'].attrs = x['Gasdet_pre_atten'].attrs
    x['PulseEnergy'].attrs['doc'] = "Energy measurement normalized by attenuators"


    
    gasdetcut =  np.array(x.Gasdet_pre_atten.values > gasdetcut_mJ, dtype=byte)
    x.coords['Gasdet_cut'] = (['time'], gasdetcut)
    x.coords['Gasdet_cut'].attrs['doc'] = "Gas detector cut.  Gasdet_pre_atten > {:} mJ".format(gasdetcut_mJ)
    
#    damagecut = np.array(phasecut & gasdetcut & (x.EBeam_damageMask.values == 0), dtype=byte)
    damagecut = np.array(gasdetcut, dtype=byte)
    x.coords['Damage_cut'] = (['time'], damagecut)
    x.coords['Damage_cut'].attrs['doc'] = "Combined Gas detector, Phase cavity and EBeam damage cut"

    try:
        x['PhotonEnergy'] = x['EBeam_ebeamPhotonEnergy']
    except:
        pass

    omit_attrs = [a for a in x.data_vars if a.endswith('_index')]
    for a in omit_attrs:
        if a in x:
            del x[a]

    omit_attrs = [a for a in x.data_vars if a.startswith('FEEGasDetEnergy')]
    for a in omit_attrs:
        if a in x:
            del x[a]

    omit_attrs = ['EBeam_damageMask', 'EBeam_ebeamPhotonEnergy', 'EBeam_ebeamEnergyBC1', 'EBeam_ebeamEnergyBC2',  
                  'EBeam_ebeamLTU250', 'EBeam_ebeamLTU450', ]
    
    for a in omit_attrs:
        if a in x:
            del x[a]

    omit_attrs = [
        'PhaseCavity_charge1',
        'PhaseCavity_charge2',
        'PhaseCavity_fitTime1',
        'PhaseCavity_fitTime2',
        ]

    for a in omit_attrs:
        if a in x:
            del x[a]

    x = psxarray.normalize_data(x, attrs)

    return psxarray.resort(x)

def setup_atten(x, prefix='dia'):
    # Change attenuator to in/out
    for i in range(1,10):
        attr = 'att{:02}_in'.format(i)
        pvalias = '{:}_s{:02}'.format(prefix,i)
        if pvalias in x:
            x[attr] = (['time'], np.array(2-x[pvalias], dtype=byte))
            x[attr].attrs['doc'] = 'Attenuator {:} In/Out'.format(pvalias)
            del x[pvalias]

def setup_imp(x):
    x.coords['Sc2Imp_name'] = (['Sc2Imp_ch'], ['TiFilterBackscatter','UnfilteredBackscatter','Evr', 'MnScatter'])

    attr = 'Sc2Imp_amplitudes'
    ch=3
    x['MnScatter'] = (['time'], x[attr].values[:,ch])
    x['MnScatter'].attrs['doc'] = "Mn Scattered diode ({:})".format(attr)
    for a in ['unit', 'alias']: 
        try:
            x['MnScatter'].attrs[a] = x[attr].attrs[a]  
        except:
            pass
    ch=1
    x['UnfilteredBackscatter'] = (['time'], x[attr].values[:,ch])
    x['UnfilteredBackscatter'].attrs['doc'] = "Unfiltered Backscatter diode ({:})".format(attr)
    for a in ['unit', 'alias']: 
        try:
            x['UnfilteredBackscatter'].attrs[a] = x[attr].attrs[a]  
        except:
            pass
    ch=0
    x['TiFilteredBackscatter'] = (['time'], x[attr].values[:,ch])
    x['TiFilteredBackscatter'].attrs['doc'] = "Ti Filtered Backscatter diode ({:})".format(attr)
    for a in ['unit', 'alias']: 
        try:
            x['TiFilteredBackscatter'].attrs[a] = x[attr].attrs[a]  
        except:
            pass

    x.attrs['summary_variables'] = ['MnScatter']


def get_correlations(y, attr='PulseEnergy', confidence=0.1, method='pearson',
        omit_list=['sec', 'nsec', 'fiducials', 'ticks']):
    """Find variables that correlate with given attr.
    """
    x = y.reset_coords()
    attrs = [a for a, item in x.data_vars.items() if item.dims == ('time',)]
    cmatrix = x[attrs].to_dataframe().corr(method=method)
    cattrs = {a: item for a, item in cmatrix[attr].iteritems() if item > confidence \
            and a != attr and a not in omit_list}    
    return cattrs

def get_cov(y, attr='PulseEnergy',
        attrs=[], confidence=0.33,
        omit_list=['sec', 'fiducials', 'ticks', 'Damage_cut']):
    """Find variables that correlate with given attr.
    """
    x = y.reset_coords()
    if not attrs:
        #attrs = [a for a, item in x.data_vars.items() if item.dims == ('time',)]
        attrs = [a for a in get_correlations(x, attr, confidence=confidence).keys() if a not in omit_list \
                and not a.startswith('FEEGasDetEnergy') and not a.startswith('Gasdet')]
        attrs.append(attr)

    df = x[attrs].to_dataframe()
    cmatrix = df.cov()
    #cattrs = {a: item for a, item in cmatrix[attr].iteritems() if a != attr and a not in omit_list}    
    return cmatrix

def setup_kbencoder(x):
    try:
        x['Xmirror_pitch'] = (['time'], x.KbEncoder_encoder_values[:,0].values)
        x['Xmirror_pitch'].attrs = x.KbEncoder_encoder_values.attrs
        x['Xmirror_pitch'].attrs['doc'] = 'Corrected KB2 HFM pitch (CXI:KB2:MMS:04 encoder)'
        x['Xmirror_pitch'].attrs['unit'] = 'mrad'
        
        x['Ymirror_pitch'] = (['time'], x.KbEncoder_encoder_values[:,1].values)
        x['Ymirror_pitch'].attrs = x.KbEncoder_encoder_values.attrs
        x['Ymirror_pitch'].attrs['doc'] = 'Corrected KB2 VFM pitch (CXI:KB2:MMS:08 encoder)'
        x['Ymirror_pitch'].attrs['unit'] = 'mrad'

    except:
        pass
    
    for attr in [attr for attr in x.keys() if attr.startswith('KbEncoder_')]:
        del x[attr]

def setup_sample(x):
    """Setup sample config -- see also data from cxi12016.
    """
    #x['Sc2Epix_spec_sum'] = x.Sc2Epix_spec_y.sum(axis=1)-x.Sc2Epix_back_y.sum(axis=1)
    #x['Sc2Epix_spec_norm'] = x.Sc2Epix_spec_y-x.Sc2Epix_back_y.values
    x['Spec_sum'] = x.Sc2Epix_spec_y.sum(axis=1) 
    #- x.Sc2Epix_back_y.values.sum(axis=1)
    
    try:
        x.attrs['Sample_y0'] = 60.0422
        x.attrs['Sample_x0'] = -34.892
      
        x['Xfocus'] = (['time'], x.Sample_x.values - x.Sample_x0)
        x['Yfocus'] = (['time'], x.Sample_y.values - x.Sample_y0)
        x['Xfocus'].attrs['doc'] = "Focus X axis roughly normalized from Sample_x"
        x['Yfocus'].attrs['doc'] = "Focus Y axis roughly normalized from Sample_y"
        x['Xfocus'].attrs['Sample_x0'] = x.attrs['Sample_x0']
        x['Yfocus'].attrs['Sample_y0'] = x.attrs['Sample_y0']
        psxarray.add_steps(x, 'Xfocus')
        psxarray.add_steps(x, 'Yfocus')
        
        # Mirror Scans
        if x.attrs['run'] >= 9 and x.attrs['run'] <= 11:
            if x.Xfocus_step.max() > x.Yfocus_step.max():
                x.attrs['scan_variables'] = ['Xfocus', 'Xmirror_pitch']
                x.coords['step'] = x.Xfocus_step
            else:
                x.attrs['scan_variables'] = ['Yfocus', 'Ymirror_pitch']
                x.coords['step'] = x.Yfocus_step
       
            x.attrs['summary_variables'] = ['UnfilteredBackscatter']
        elif x.attrs['run'] >= 21 and x.attrs['run'] <= 999:
            #add_steps(x, 'Xfocus')
            #add_steps(x, 'Yfocus')
            x.attrs['scan_variables'] = []
            for a in ['Xfocus_step', 'Yfocus_step', 'Xfocus', 'Yfocus']:
                if a in x:
                    del x[a]

        for alias, item in peaks.items():
            #ka = item.get('center')-item.get('width')
            #kb = item.get('center')+item.get('width')+1
            #x[alias] = x.Sc2Epix_spec_y[:,ka:kb].sum(axis=1)
            dy = item.get('width')
            yo = item.get('center')
            x[alias] = x.Sc2Epix_spec_y.where((x.Sc2Epix_yspec-yo)<=dy).sum(axis=1)
            for attr, val in item.items():
                x[alias].attrs[attr] = val

            #- x.Sc2Epix_back_y[:,400:600].values.sum(axis=1)
            #- x.Sc2Epix_back_y[:,660:].values.sum(axis=1)

    except:
        pass

def build_html(x, dorunnum=False, ioff=True):
    """Make html web page from xarray summary dataset.
    """
    # do not perform interactive plotting
    from PyDataSource.build_htmml import Build_html
    if ioff:
        plt.ioff()

    make_cuts(x)
    attrs = [attr for attr,item in x.data_vars.items() if item.dims == ('time',)]
    df = x[attrs].where(x.Damage_cut == 1).to_array().to_pandas().T.dropna()
    
    b = Build_html(x, path='/reg/neh/operator/cxiopr/userexpts/cxi12116/RunSummary/')
     
    variables = ['Sc2Inline_roi', 'Sc2Epix_spec', 'Sc2Epix_spec_y'] 
    b.add_summary(variables)
    b.add_detector('Sc2Imp', cut='Damage_cut') 
    #run_summary(b)
    
    #b.add_summary(['Sc2Imp_filtered'],groupby='run')
    #b.add_all(cut='Damage_cut')

    attr_groups = {
#                   'Undulator X': ['EBeam_ebeamUndAngX', 'EBeam_ebeamUndPosX'],
#                   'Undulator Y': ['EBeam_ebeamUndAngY', 'EBeam_ebeamUndPosY'],
#                   'LTU X': ['EBeam_ebeamLTUAngX', 'EBeam_ebeamLTUPosX'],
#                   'LTU Y': ['EBeam_ebeamLTUAngY', 'EBeam_ebeamLTUPosY'],
#                   'LTU 250 and 450': ['EBeam_ebeamLTU250', 'EBeam_ebeamLTU450'],
#                   'Energy BC': ['EBeam_ebeamEnergyBC1', 'EBeam_ebeamEnergyBC2'],
#                   'Charge': ['EBeam_ebeamCharge', 'EBeam_ebeamDumpCharge'],
#                   'Peak Current': ['EBeam_ebeamPkCurrBC1', 'EBeam_ebeamPkCurrBC2'],
#                   'XTCAV': ['EBeam_ebeamXTCAVAmpl', 'EBeam_ebeamXTCAVPhase'],
#                   'Phase Cavity Charge': ['PhaseCavity_charge1', 'PhaseCavity_charge2'],
#                   'Phase Cavity Fit': ['PhaseCavity_fitTime1', 'PhaseCavity_fitTime2'],
#                   'Gasdet': ['FEEGasDetEnergy_f_63_ENRC', 'Gasdet_pre_atten'],
#                   'Diodes Backscatter': ['TiFilteredBackscatter','UnfilteredBackscatter'],
#                   'Diodes Filtered': ['TiFilteredBackscatter'],
                   'Summary': ['Spec_sum', 'Kalpha1'],
                   'Satellite': ['Satellite1', 'Satellite2'],
                   'Kalpha': ['Kalpha1', 'Kalpha2'],
                   'Diodes': ['TiFilteredBackscatter', 'XCS_IPM_02_sum']
                   }

    group = 'PulseEnergy_index'
    for catagory, grp_attrs in attr_groups.items():
        attrs = ['PhotonEnergy', 'PulseEnergy', 'MnScatter']
        #if 'XCS_IPM_02_sum' in x:
        #    attrs.append('XCS_IPM_02_sum')
        
        OK = False
        for attr in grp_attrs:
            if attr in x:
                attrs.append(attr)
                OK = True

        howto = []
        if OK:
            print(group, catagory)
            try:
                b.add_scatter(df, catagory=catagory, attrs=attrs, howto=howto, group=group)
            except:
                pass

    b.to_html()

    return b

def summary_html(a, ioff=True, exp_path=exp_path):
    """RunSummary built from statistics of each step in the run. 
    """
    from PyDataSource.build_htmml import Build_html
    if 'stat' not in a.dims:
        a = psxarray.to_summary(a)

    if ioff:
        plt.ioff()

    path=os.path.join(exp_path, 'RunSummary')
    print(path)
    b = Build_html(a, path=os.path.join(exp_path, 'RunSummary'))

    b.add_xy_ploterr('MnScatter_norm', xaxis='atten_thick', logx=False, logy=False, catagory='FitResults')
    b.add_xy_ploterr('MnScatter_norm', xaxis='PulseEnergy', logx=True, logy=False, catagory='FitResults')
    
    b.add_xy_ploterr('atten_corr', xaxis='run', logx=False, logy=True, catagory='RunSummary')
    b.add_xy_ploterr('atten_thick', xaxis='run', logx=False, logy=False, catagory='RunSummary')
    b.to_html()

    return b


#attrs = ['Spec_sum_norm','MnScatter_norm','TiFilteredBackscatter_norm','atten_thick']
#x.reset_coords()[attrs].where(x.Damage_cut,drop=True).groupby('atten_thick').mean().to_array().to_pandas().T

def run_summary(self, catagory='Summary', howto=[]):
   
    #if catagory not in self.results:
    #    self.results[catagory] = {'figure': {}, 'table': {}, 'text': {}}

    x = self._xdat
    attrs = x.summary_variables
    attrs.append(x.scan_variables[0])
    df = x[attrs].where(x.Damage_cut).to_array().to_pandas().T.dropna().groupby(x.scan_variables[0]).mean()
    howto.append("attrs = {:}".format(attrs))
    howto.append("df = x[attrs].where(x.Damage_cut).to_array().to_pandas().T.dropna().groupby(x.scan_variables[0]).mean()")
    #self.results[catagory]['text'].update({'setup':{'howto': howto}})
    self.add_setup(catagory, howto)

    for attr in df.columns:
        df[attr].plot()
        howto = ["df['{:}'].plot()".format(attr)]
        plt_type = '{:} vs {:} Wire Scan'.format(attr, x.scan_variables[0])
        self.add_plot(catagory, plt_type, howto=howto)

    #df.plot(subplots=True, sharex=True, layout=(3, 1), figsize=(8, 12))
    #howto.append("df.plot(subplots=True, sharex=True, layout=(3, 1), figsize=(8, 12))")


def main():
    """Main script to create run summary.
    """
    import RunSummary
    time0 = time.time()
    print(time0)
    args = RunSummary.initArgs()
    print(args)
    x = make_summary(**vars(args))
    return x

if __name__ == "__main__":
    sys.exit(main())


