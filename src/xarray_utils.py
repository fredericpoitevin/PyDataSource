from __future__ import print_function
from __future__ import absolute_import
def to_summary(x, dim='time', groupby='step', 
        save_summary=False,
        normby=None,
        omit_list=['run', 'sec', 'nsec', 'fiducials', 'ticks'],
        #stats=['mean', 'std', 'min', 'max',],
        stats=['mean', 'std', 'var', 'min', 'max', 'count'],
        cut='Damage_cut',
        **kwargs):
    """
    Summarize a run.
    
    Parameters
    ---------
    x : xarray.Dataset
        input xarray Dataset
    dim : str
        dimension to summarize over [default = 'time']
    groupby : str
        coordinate to groupby [default = 'step']
    save_summary : bool
        Save resulting xarray Dataset object
    normby : list or dict
        Normalize data by attributes
    omit_list : list
        List of Dataset attributes to omit
    stats : list
        List of statistical operations to be performed for summary.
        Default = ['mean', 'std', 'var', 'min', 'max', 'count']
    """
    import xarray as xr
    xattrs = x.attrs
    data_attrs = {attr: x[attr].attrs for attr in x.data_vars}
    if cut in x:
        x = x.where(x[cut]).dropna(dim)
   
    coords = [c for c in x.coords if c != dim and c not in omit_list and dim in x.coords[c].dims] 
    x = x.reset_coords(coords)
    if isinstance(normby, dict):
        for norm_attr, attrs in normby.items():
            x = normalize_data(x, attrs, norm_attr=norm_attr)
    elif isinstance(normby, list):
        x = normalize_data(x, normby)

    # With new xarray 0.9.1 need to make sure loaded otherwise h5py error
    x.load()
    xgroups = {}
    if groupby:
        try:
            sattrs = [a for a in x.data_vars if 'stat' in x[a].dims or groupby+'s' in x[a].dims]
            if sattrs:
                xstat = x[sattrs].rename({groupby+'s': groupby})
                x = x.drop(sattrs)
            x = x.groupby(groupby)
    #        group_vars = [attr for attr in x.keys() if groupby+'s' in x[attr].dims]
    #        if group_vars:
    #            xgroups = {attr: x[attr].rename({groupby+'s': groupby}) for attr in group_vars}
    #            for attr in group_vars:
    #                del x[attr]
            #x = x.groupby(groupby)
        except:
            print('Cannot groupby {:} -- ignoring groupby'.format(groupby))


    dsets = [getattr(x, func)(dim=dim) for func in stats]
    x = xr.concat(dsets, stats).rename({'concat_dim': 'stat'})
#    try:
#        x = xr.concat([dsets, stats], dim='stat')
#    except:
#        return dsets
    
    for attr,val in xattrs.items():
        x.attrs[attr] = val
    for attr,item in data_attrs.items():
        if attr in x:
            x[attr].attrs.update(item)

    for c in coords:                                                       
        x = x.set_coords(c)

    if groupby and sattrs:
        x = x.merge(xstat)
#    try:
#        for attr in group_vars:
#            x[attr] = xgroups[attr]
#    except:
#        return x, xgroups

    x = resort(x)
    if save_summary:
        to_h5netcdf(x)
    return x

def resort(x):
    """
    Resort alphabitically xarray Dataset
    """
    coords = sorted([c for c in x.coords.keys() if c not in x.coords.dims])
    x = x.reset_coords()
    x = x[sorted(x.data_vars)]

    for c in coords:                                                       
        x = x.set_coords(c)

    return x


def normalize_data(x, variables=[], norm_attr='PulseEnergy', name='norm', quiet=True):
    """
    Normalize a list of variables with norm_attr [default = 'PulseEnergy']
    """
    if not variables:
        variables = [a for a in get_correlations(x) if not a.endswith('_'+name)]    
    elif isinstance(variables, str):
        variables = [variables]

    for attr in variables:
        aname = attr+'_'+name
        try:
            x[aname] = x[attr]/x[norm_attr]
            x[aname].attrs = x[attr].attrs
            try:
                x[aname].attrs['doc'] = x[aname].attrs.get('doc','')+' -- normalized to '+norm_attr
                units = x[attr].attrs.get('unit')
                norm_units = x[norm_attr].attrs.get('unit')
                if units and norm_units:
                    x[aname].attrs['unit'] = '/'.join([units, norm_units])
            except:
                if not quiet:
                    print('cannot add attrs for', aname)
        except:
            print('Cannot normalize {:} with {:}'.format(attr, norm_attr))

    return  resort(x)

def add_index(x, attr, name=None, nbins=8, bins=None, percentiles=None):
    """
    Add index for attribute
    """
    import numpy as np
    if not bins:
        if not percentiles:
            percentiles = (np.arange(nbins+1))/float(nbins)*100.

        bins = np.percentile(x[attr].to_pandas().dropna(), percentiles)

    if not name:
        name = attr+'_index'

    #per = [percentiles[i-1] for i in np.digitize(x[attr].values, bins)]
    x[name] = (['time'], np.digitize(x[attr].values, bins))

def add_steps(x, attr, name=None):
    """
    Add step coordinate for attr
    """
    vals = getattr(x, attr).values
    steps = np.sort(list(set(vals)))
    asteps = np.digitize(vals, steps)
    if not name:
        name = attr+'_step'
 
    x.coords[name] = (['time'], asteps)

def add_butterworth_filter(x, attr, lowcut=None, highcut=None, order=10, 
            norm_attr=None, norm_threshold=None, 
            filt_name=None, pass_name=None, imputer_name=None,
            dim='time', drop_attr='XrayOff', 
            threshold=None, quiet=True, **kwargs):
    """
    Add Butterworth high, low or band pass filter for attr.
    
    Events with drop_attr and below threshold are considered missing 
    during filtering.  They are replaced using sklearn Imputer method

    Parameters
    ----------
    attr : str
        Name of Dataset attribute
    lowcut : float
        Low pass frequency cutoff
    highcut : float
        High pass frequency cutoff
    order : int
        Butterworth filter order [Default = 10, i.e., 10th order Butterworth filter]
    filt_name : str, optional
        Name of resutling filtered data.  [Default = attr+'_filt']
    pass_name : str, optional
        Name of resutling band pass data.
        - Lowpass Default = attr+'_lowpass'
        - Highpass Default = attr+'_highpass'
        - Bandpass Default = attr+'_bandpass'
    imputer_name : str, optional
        Name of corrected data using Imputer method for drop_attr and below 
        threshold data
    dim : str
        Dimension to filter over [Default = 'time']
    drop_attr : str
        Replace points where drop_attr is True with average of nearest points 
        during filtering.  
        Resulting filtered data where drop_attr is True is unaffected.
    threshold : float
        Minimum threshold to filter
    """
    import xarray as xr
    from .filter_methods import butter_bandpass_filter 
    if dim == 'time':
        fs = x.time.size/float(x.time.sec.max()-x.time.sec.min())
    else:
        fs = x[dim].size/float(x[dim].max()-x[dim].min())

    if lowcut and highcut:
        btype = 'band'
    elif highcut:
        btype = 'high'
    elif lowcut:
        btype = 'low'
    else:
        print('Error -- must supply lowcut, highcut or both')

    if 'time_ns' not in x:
        x.coords['time_ns'] = x['sec']*1e9+x['nsec']

    if drop_attr and drop_attr in x:
        cut = ~x[drop_attr]
        if not quiet:
            print(('drop sum', cut.sum()))
        if threshold:
            cut =  xr.ufuncs.logical_and(cut, x[attr]>threshold)
            if not quiet:
                print(('threshold cut', cut.sum()))
    
        if norm_attr and norm_attr in x:
            if norm_threshold:
                cut =  xr.ufuncs.logical_and(cut, x[norm_attr]>norm_threshold)
                if not quiet:
                    print(('norm threshold cut', cut.sum()))
     
            dfnorm = x[[norm_attr]].where(cut).reset_coords()[norm_attr].to_dataframe()
            dfnorm0 = x[[norm_attr]].reset_coords()[norm_attr].to_dataframe()
        
        else:
            dfnorm = None
            if norm_attr:
                print('{:} not available for normalization'.format(norm_attr))
    
        df = x[[attr]].where(cut).reset_coords()[attr].to_dataframe()
        df0 = x[[attr]].reset_coords()[attr].to_dataframe()
    
    else:
        cut = None
        df = x[attr].reset_coords()[attr].to_dataframe()
        
    try:
        from sklearn.preprocessing import Imputer
        mean_imputer = Imputer(missing_values='NaN', strategy='mean', axis=0)
        if dfnorm is None:
            mean_imputer = mean_imputer.fit(df)
            data = mean_imputer.transform(df.values)[:,0]
        else:
            mean_imputer = mean_imputer.fit(dfnorm)
            data = mean_imputer.transform(df.values)[:,0]
            data /= mean_imputer.transform(dfnorm.values)[:,0]
   
    except:
        mean_imputer = None
        cut = None
        data = df.values[:,0]
        print('Could not preprocess {:} using Imputation of events with {:}'.format(attr, drop_attr))

    if cut is not None:
        data0 = df0.values/dfnorm0.values

    filt = butter_bandpass_filter(data, fs, lowcut=lowcut, highcut=highcut, order=order)

    if imputer_name and mean_imputer is not None:
        x[imputer_name] = ((dim), data)
        x[imputer_name].attrs = x[attr].attrs
        x[imputer_name].attrs['doc'] = '{:} imputer corrected data'.format(attr)
        x[imputer_name].attrs['drop_attr'] = drop_attr
        if threshold:
            x[imputer_name].attrs['threshold'] = threshold 

    if filt is None:
        print('Error making filter for {:}'.format(attr))
   
    else:

        if dfnorm is None:
            if not filt_name:
                filt_name = '{:}_filt'.format(attr)
            if not pass_name:
                pass_name = '{:}_{:}pass'.format(attr,btype)
        else:
            if not filt_name:
                filt_name = '{:}_norm_filt'.format(attr)
            if not pass_name:
                pass_name = '{:}_norm_{:}pass'.format(attr,btype)

        x[filt_name] = ((dim), filt)
        x[filt_name].attrs = x[attr].attrs
        x[filt_name].attrs['doc'] = '{:}-pass butterworth filter of {:}'.format(btype, attr)
        x[filt_name].attrs['order'] = order
        x[filt_name].attrs['btype'] = btype

        x[pass_name] = ((dim), data-filt) 
#        if cut is not None:
#            try:
#                x[pass_name][~cut] = data0[~cut]
#            except:
#                print('Cannot reset cut data')
        
        x[pass_name].attrs = x[attr].attrs
        x[pass_name].attrs['doc'] = '{:} {:}-pass butterworth filtered data'.format(attr, btype)
        x[pass_name].attrs['order'] = order
        x[pass_name].attrs['btype'] = btype

        if highcut:
            x[pass_name].attrs['highcut'] = highcut
            x[filt_name].attrs['highcut'] = highcut
        if lowcut:
            x[pass_name].attrs['lowcut'] = lowcut
            x[filt_name].attrs['lowcut'] = lowcut
        
    return x

def get_correlations(y, attr, confidence=0.33, method='pearson',
        omit_list=['sec', 'nsec', 'fiducials', 'ticks', 'Damage_cut', 'EBeam_damageMask']):
    """Find variables that correlate with given attr.
    """
    x = y.reset_coords()
    attrs = [a for a, item in x.data_vars.items() if item.dims == ('time',)]
    cmatrix = x[attrs].to_dataframe().corr(method=method)
    cattrs = {a: item for a, item in cmatrix[attr].iteritems() if item > confidence \
            and a != attr and a not in omit_list}    
    return cattrs

def get_cov(y, attr, attrs=[], confidence=0.33,
        omit_list=['sec', 'nsec', 'fiducials', 'ticks', 'Damage_cut', 'EBeam_damageMask']):
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

def heatmap(df, attrs=[], method='pearson', confidence=0.33, position=(0.3,0.35,0.45,0.6), show=False):
    import matplotlib.pyplot as plt
    import seaborn as sns
    fig = plt.figure() 
    ax1 = fig.add_subplot(111)
    #plt.gca().set_position(position)
    if attrs:
        df = df[attrs]
    
    corr = df.corr(method=method)
    if confidence:
        cattrs = df.keys()[((abs(corr)>=confidence).sum()>1).values]
        corr = df[cattrs].corr(method=method)
    
    sns.heatmap(corr,annot=True)
    plt.yticks(rotation=0)
    plt.xticks(rotation=90) 
    plt.gca().set_position(position)
    if show:
        plt.show()   

    return corr

def xy_ploterr(a, attr=None, xaxis=None, title='', desc=None, 
        fmt='o', position=(.1,.2,.8,.7), **kwargs):
    """Plot summary data with error bars, e.g.,
        xy_ploterr(x, 'MnScatter','Sample_z',logy=True)
    """
    import matplotlib.pyplot as plt
    if not attr:
        print('Must supply plot attribute')
        return

    if 'groupby' in kwargs:
        groupby=kwargs['groupby']
    elif 'step' in a.dims:
        groupby='step'
    elif 'steps' in a.dims:
        groupby='steps'
    else:
        groupby='run'

    run = a.attrs.get('run')
    experiment = a.attrs.get('experiment', '')
    runstr = '{:} Run {:}'.format(experiment, run)
    name = a.attrs.get('name', runstr)
    if not title:
        title = '{:}: {:}'.format(name, attr)

    if not xaxis:
        xaxis = a.attrs.get('scan_variables')
        if xaxis:
            xaxis = xaxis[0]

    if xaxis:
        if 'stat' in a[xaxis].dims:
            xerr = a[xaxis].sel(stat='std').values
            a[xaxis+'_axis'] = ([groupby], a[xaxis].sel(stat='mean').values)
            xaxis = xaxis+'_axis'
        else:
            xerr = None

        a = a.swap_dims({groupby:xaxis})
    
    else:
        xaxis = groupby
        xerr = None

    ylabel = kwargs.get('ylabel', '')
    if not ylabel:
        ylabel = a[attr].name
        unit = a[attr].attrs.get('unit')
        if unit:
            ylabel = '{:} [{:}]'.format(ylabel, unit)

    xlabel = kwargs.get('xlabel', '')
    if not xlabel:
        xlabel = a[xaxis].name
        unit = a[xaxis].attrs.get('unit')
        if unit:
            xlabel = '{:} [{:}]'.format(xlabel, unit)
    
    if desc is None:
        desc = a[attr].attrs.get('doc', '')

    ndims = len(a[attr].dims)
    if ndims == 2:
        c = a[attr].to_pandas().T
        if xerr is not None:
            c['xerr'] = xerr
        c = c.sort_index()
        
        plt.figure()
        plt.gca().set_position(position)
        p = c['mean'].plot(yerr=c['std'],xerr=c.get('xerr'), title=title, fmt=fmt, **kwargs)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if desc:
            plt.text(-.1,-.2, desc, transform=p.transAxes, wrap=True)   
 
        return p 
    elif ndims == 3:
        plt.figure()
        plt.gca().set_position((.1,.2,.8,.7))
        pdim = [d for d in a[attr].dims if d not in ['stat', groupby, xaxis]][0]
        for i in range(len(a[attr].coords[pdim])):
            c = a[attr].sel(**{pdim:i})
            if pdim in c:
                c = c.drop(pdim)
            c = c.to_pandas().T.sort_index()
            p = c['mean'].plot(yerr=c['std'], fmt=fmt, **kwargs)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        p.set_title(title)
        if desc:
            plt.text(-.1,-.2, desc, transform=p.transAxes, wrap=True)   

        return p 
    else:
        print('Too many dims to plot')

def ttest_groupby(xo, attr, groupby='ec162', ishot=0, nearest=None, verbose=False):
    """
    nearest neighbors compare with
    """
    from scipy import stats
    try:
        # Select DataArray for attr with groupby as coord
        da = xo.reset_coords()[[attr,groupby]].set_coords(groupby)[attr].load().dropna(dim='time')
    except:
        return None
    ntime = len(da.time)
    g = da.groupby(groupby).groups
    if len(g) > 1:
        if nearest is not None:
            k0 = []
            for inear in range(-nearest,0)+range(1,nearest+1):
                for itest in g[1]:
                    a = itest+inear
                    if a>0 and a<ntime:
                        k0.append(a)
        else:
            k0 = [a+ishot for a in g[0] if (a+ishot)>0 and (a+ishot)<ntime]
        k1 = [a+ishot for a in g[1] if (a+ishot)>0 and (a+ishot)<ntime]
        df1 = da[k1].to_pandas().dropna()
        df0 = da[k0].to_pandas().dropna()
        ttest = stats.ttest_ind(df1,df0)
        return ttest
    else:
        if verbose:
            print('{:} has only one group -- cannot compare'.format(attr))
        return None

def test_correlation(x, attr0, attr1='Gasdet_post_atten', cut=None, shift=None, dim='time'):
    from scipy import stats
    import numpy as np
    xds = x[[attr0,attr1]]
    if cut and cut in xds:
        xds = xds.where(x[cut]).where(x[cut].shift(**{dim:-shift})).dropna(dim) 

    xds = xds.reset_coords()
    
    df0 = xds[attr0].to_pandas()
    if shift:
        df0 = df0.shift(-shift)
    df1 = xds[attr1].to_pandas()
    kind = np.isfinite(df1) & np.isfinite(df0)
    result = stats.pearsonr(df0[kind],df1[kind])
    return result

def set_delta_beam(x, code='ec162', attr='delta_drop'):
    """Find the number of beam codes to nearest code
       using time stamps
    """
    import traceback
    import pandas as pd
    import numpy as np
    if code not in x:
        return None
    else:
        try:
            df0 = x.reset_coords()[code].to_pandas()
            if df0.any():
                df_beam = x.reset_coords().fiducials.to_pandas()/3
                df_drop = df_beam[df0]
                x.coords[attr] = df_beam
                i = 0
                for a, b in df_beam.iteritems():
                    x.coords[attr][i] = b-df_drop.values[np.argmin(abs(df_drop.index-a))]
                    i += 1
                # Not robust way to get rid of outliers 
                #dbeam_max = df_beam.diff().max()
                #x.coords[attr][abs(x.coords[attr]) > dbeam_max] = dbeam_max
                x.coords[attr].attrs['doc'] = "number of beam codes to nearest {:}".format(code) 
                return x.coords[attr]
            else:
                print('No event code {:} to drop'.format(code))
                return None
        except:
            traceback.print_exc('Cannot set drop_code = '.format(code))


def find_beam_correlations(xo, pvalue=1e-20, pvalue0_ratio=0.1, corr_pvalue=0.0001,
            adet_pvalue0_ratio=0.33, adet_pvalue=1e-5,
            groupby='ec162', nearest=5, corr_coord='delta_drop',
            pulse=None, confidence=0.1, sigma0=5, sigma_scale=2,
            percentiles=[0.5], 
            save_file=None,
            cut=None, verbose=False, **kwargs):
    """
    """
    import traceback
    xcoor_coord = set_delta_beam(xo, code=groupby, attr=corr_coord)
    import xarray as xr
    import pandas as pd
    xstats = xr.Dataset()
    attrs = [a for a in xo.variables if xo[a].dims == ('time',)]
    xds = xo[attrs].load()
    xo.attrs['drop_shot_detected'] = []
    xo.attrs['timing_error_detected'] = []
    xo.attrs['beam_warning_detected'] = []
    xo.attrs['beam_corr_detected'] = []
    area_detectors = [a for a in xds.attrs.get('area_detectors',[])]
    wf_detectors = [a for a in xds.attrs.get('wf_detectors',[])]
    adets = list(sorted(set(area_detectors + wf_detectors)))
    for attr in adets:
        try:
            if len(xds[attr].dropna(dim='time')) < 5:
                adets.remove(attr)
        except:
            adets.remove(attr)

    if 'EBeam_damageMask' in xds:
        xds = xds.drop('EBeam_damageMask')
    if pulse not in xds:
        if 'FEEGasDetEnergy_f_21_ENRC' in xds:
            pulse = 'FEEGasDetEnergy_f_21_ENRC'
            if xds[pulse].sum() == 0:
                pulse = 'FEEGasDetEnergy_f_11_ENRC'

        else:
            pulse = 'GasDet_f21'
    xo.attrs['beam_corr_attr'] = pulse
    xo.attrs['beam_corr_confidence'] = confidence 
    xo.attrs['drop_attr'] = groupby

    print('Analyzing beam correlations for {:} Run {:}'.format(xo.experiment, xo.run))
    for attr in [a for a in xds.data_vars if xds[a].dims == ('time',)]:
        adf = xds[attr].dropna(dim='time').to_pandas()
        if len(adf) < 5:
            print('Skipping {:} -- too few events'.format(attr))
            continue
        adf_unique = adf.unique()
        if len(adf_unique) < 2:
            print('Skipping {:} -- only one unique value in data'.format(attr))
            continue

        #if verbose:
        print('*****', attr, '*******')
        attrs = [attr, groupby, pulse]
        if cut:
            if cut in xds:
                attrs.append(cut)
            else:
                print('Ommitting cut: {:} not valid cut'.format(cut))
                cut = None

        # use lower threshold for area detector
        if attr in adets:
            pvalue_thresh = adet_pvalue
            pval0_ratio = adet_pvalue0_ratio
        else:
            pvalue_thresh = pvalue
            pval0_ratio = pvalue0_ratio
        
        x = xds[attrs]
        alias = x[attr].attrs.get('alias')
        attest = {}
        actest = {}

        #ashots = range(-nearest,nearest+1)
        #for ishot in ashots:
        # Need to rework ttest_groupby and test_correlation
        # to handle non 120Hz data in a more natural way
        # This fixes timing errors but not necessarily robustly
        ashots = sorted([a for a in set(xo[corr_coord].values) if abs(a) <= nearest])
        for iashot,idelta in enumerate(ashots):
            ishot = iashot-ashots.index(0)
            if x[pulse].attrs.get('alias') == x[attr].attrs.get('alias'):
                # if FEEGasDetEnergy detector then test agains all other with shifted drops
                tnearest = None 
            else:
                tnearest = nearest
            ttest = ttest_groupby(x, attr, groupby=groupby, ishot=ishot, 
                                    nearest=tnearest)
            attest[idelta] = ttest
            if ttest is None:
                if verbose:
                    print(attr, groupby, 'Not valid test')
                continue
         
            # Do not test correlation of same attr
            ctest = test_correlation(x, attr, pulse, cut=cut, shift=ishot)
            actest[idelta] = ctest

        xstd = x[attr].groupby(groupby).std()
        xmean = x[attr].groupby(groupby).mean()
        try:
            # Statistics for delta_drop
            df = xo.reset_coords()[[attr,corr_coord]].to_dataframe()
            df_nearest = df.where(abs(df[corr_coord]) <= nearest).dropna()
            df_table = df_nearest[[attr,corr_coord]].groupby(corr_coord).describe(percentiles=percentiles)[attr]

            # T-test for the means of *two independent* samples of scores.
            # see scipy.stats.ttest_ind
            # The calculated t-statistic, The two-tailed p-value
            df_ttest = pd.DataFrame(attest,index=['t_stat','t_pvalue']).T
            # Pearson correlation coefficient and the p-value for testing non-correlation.
            # see scipy.stats.pearsonr
            # Pearson's correlation coefficient, 2-tailed p-value
            df_ctest = pd.DataFrame(actest,index=['beam_corr','c_pvalue']).T        
            df_stats = df_table.join(df_ctest).join(df_ttest)
            
            tag_shot_corr = False
            t_pvalue0 = df_stats['t_pvalue'][0]
            if abs(float((xmean[1]-xmean[0])/xstd[0])) > sigma0 \
                    or xstd[1] > xstd[0]*sigma0 \
                    or xstd[0] > xstd[1]*sigma_scale:
                # First check if dropped shot and regular shot mean values differ by > sigma0 
                # If drop_shot std (i.e., xstd[1] is < sigma/sigma_scale times other shots then timing is correct
                #    e.g., opal_1 from exp=xppls7917:run=98 which has laser on/off
                #          and the signal is strong and stable with no X-rays but
                #          varies significantly with X-rays because of laser x-ray effects
                # If drop_shot std is < sigma0 times other shots then also timing is correct
                ishot = 0
                tag_shot_corr = True
            else:
                # Otherwise choose shot with greatest ttest significance
                ishot = df_stats['t_stat'].abs().idxmax()
            
            try:
                # ishot can be nan when t_stat are nan
                t_pvalue = df_stats['t_pvalue'][ishot]
                sig_significance = (df_stats['mean']/df_stats['std'])[ishot]
                # Check pvalue valid and if not timed with drop_code then
                # check ratio of found ishot pvalue is less than pvalue on drop code
                # Ignore if mean/std for time detected is too big to be reasonable
                if sig_significance < 1.e10 and t_pvalue0 != 1 and t_pvalue <= pvalue_thresh \
                            and (t_pvalue/t_pvalue0 < pval0_ratio or ishot == 0):
                    tag_shot_corr = True

            except:
                t_pvalue = None

           
            try:
                shot_corr_detected = False
                shot_corr = df_stats['beam_corr'].abs().idxmax()
                c_pvalue = df_stats['c_pvalue'][shot_corr]
                beam_corr = df_stats['beam_corr'][shot_corr]
                if shot_corr != 0:
                    c_pvalue0 = df_stats['c_pvalue'][0]
                    beam_corr0 = df_stats['beam_corr'][0]
                    # default to correlated with X-rays if not clear correlation on timing error
                    if c_pvalue < corr_pvalue and abs(beam_corr) > 0.2:
                        shot_corr_detected = True
                        # Make sure shot correlation is not also OK on drop shot
                        if (beam_corr0 > beam_corr/2. or (c_pvalue>0 and c_pvalue0/c_pvalue >= pval0_ratio)):
                            shot_corr = 0
                            beam_corr = df_stats['beam_corr'][shot_corr]
                            c_pvalue = df_stats['c_pvalue'][shot_corr]
                else:
                    dfs = df_stats['c_pvalue'].copy()
                    dfs.pop(0)
                    c_pvalue0 = dfs.min()
                    dfs = df_stats['beam_corr'].copy()
                    dfs.pop(0)
                    beam_corr0 = dfs.max()
                    #beam_corr = df_stats['beam_corr'][shot_corr]
                    # If off-by-one make sure not beam correlated on ec162
                    if ishot != 0 and c_pvalue < corr_pvalue and \
                                (abs(beam_corr) > 0.2 or abs(beam_corr) > abs(beam_corr0)*2.):
                        tag_shot_corr = False
                    elif ishot != 0 and (c_pvalue == 1 or c_pvalue0 == 1):
                        tag_shot_corr = False

                    # Make checks to be sure beam_corr is valid 
                    if c_pvalue0 == 0 or c_pvalue/c_pvalue0 > 1./pval0_ratio:
                        shot_corr_detected = False
                    elif c_pvalue == 0 and c_pvalue0 < pvalue**2 and beam_corr0>beam_corr/2.:
                        shot_corr_detected = False
                    elif c_pvalue < corr_pvalue and (abs(beam_corr) > 0.2 or abs(beam_corr) > abs(beam_corr0)*2.):
                        shot_corr_detected = True

                # Make sure beam correlation is consistent with drop shot detection
                # If beam_corr is on ec162 then do not tag
                if shot_corr_detected:
                    if tag_shot_corr and ishot != shot_corr:
                        if shot_corr == 0:
                            tag_shot_corr = False
                            ishot = shot_corr
                            t_pvalue = df_stats['t_pvalue'][ishot]
                    
                    xo[attr].attrs['beam_corr'] = beam_corr
                    xo[attr].attrs['shot_corr'] = shot_corr
                    if attr not in xo.attrs['beam_corr_detected']:
                        xo.attrs['beam_corr_detected'].append(attr)
                    note = 'Beam-corr  of {:5.3f}'.format(beam_corr)
                    corr_note = '{:} on {:2} shot for {:} detector {:} (corr = {:}, c_pvalue = {:})'.format(note, shot_corr, 
                                alias, attr, beam_corr, c_pvalue)
                    print(corr_note)

            except:
                shot_corr = 0
                c_pvalue = None
                beam_corr = None
                shot_corr_detected = False

            # Tag shot correlation after checking beam correlation
            if tag_shot_corr: 
                xo[attr].attrs['delta_beam'] = ishot 
                xo[attr].attrs['delta_beam_pvalue'] = t_pvalue
                if ishot == 0:
                    note = 'Drop-shot detected '
                    xo[attr].attrs['drop_shot_detected'] = True
                    if attr not in xo.attrs['drop_shot_detected']:
                        xo.attrs['drop_shot_detected'].append(attr)
                elif alias in ['EBeam','PhaseCavity'] or attr.endswith('xpos') or attr.endswith('ypos'):
                    note = 'Beam-warning detected'
                    xo[attr].attrs['beam_warning_detected'] = True
                    if attr not in xo.attrs['beam_warning_detected']:
                        xo.attrs['beam_warning_detected'].append(attr)
                else:
                    note = 'Time-error detected'
                    xo[attr].attrs['timing_error_detected'] = True
                    if attr not in xo.attrs['timing_error_detected']:
                        xo.attrs['timing_error_detected'].append(attr)
                
                shot_note = '{:} on {:2} shot for {:} detector {:} (t_pvalue={:})'.format(note, ishot, 
                        alias, attr, t_pvalue)
                print(shot_note)


            #xstats[attr] = ((corr_coord,'drop_stats'), df_stats)
            xstats[attr] = df_stats

        except:
            traceback.print_exc()
            print('xmean', xmean)
            print('xstd', xstd)
            print('Cannot calc stats for', attr, ishot)
        
    for avar, da in xo.data_vars.items():
        try:
            if avar in xstats:
                for a in ['doc', 'unit', 'alias']:
                    val = xo[avar].attrs.get(a, '') 
                    if isinstance(val, unicode):
                        val = str(val)
                    xstats[avar].attrs[a] = val 
        except:
            print('Cannot add attrs for {:}'.format(avar))
    
    for attr, val in xo.attrs.items():
        try:
            if isinstance(val, list) and len(val) > 0 and isinstance(val[0], unicode):
                val = [str(v) for v in val]
            elif isinstance(val, unicode):
                val = str(val)
            xstats.attrs[attr] = val
        except:
            print('Cannot add attrs for {:}'.format(attr))

    if 'dim_1' in xstats:
        xstats = xstats.rename({'dim_1':'drop_stats'})
    else:
        print('No Drop stats')
        return None

    if save_file:
        try:
            xstats.to_netcdf(save_file, engine='h5netcdf')
            print('Save drop summary file: {:}'.format(save_file))
        except:
            print('Cannot save file: {:}'.format(save_file))

    return xstats

def clean_dict(attrs):
    """
    Replace unicode with str in dict 
    """
    import numpy as np
    for attr, item in attrs.items():
        if isinstance(attr, unicode):
            attr = str(attr)
        if isinstance(item, list):
            clean_item = []
            for a in item:
                if isinstance(a, unicode):
                    clean_item.append(str(a))
                else:
                    clean_item.append(a)
            attrs[attr] = clean_item
        elif isinstance(item, unicode):
            attrs[attr] = str(item)
        elif isinstance(item, np.ndarray):
            attrs[attr] = [a for a in item]
        else:
            attrs[attr] = item

    return attrs

def clean_dataset(xds):
    xds.attrs = clean_dict(xds.attrs)
    for attr in xds.coords.keys():
        xcoord = xds.coords[attr]
        attrs = clean_dict(xcoord.attrs)
        xcoord.attrs = attrs
        if isinstance(attr, unicode):
            del xds.coords[attr]
            attr = str(attr)
        xds.coords[attr] = xcoord
    for attr, xdata in xds.data_vars.items():
        attrs = clean_dict(xdata.attrs)
        xdata.attrs = attrs
        if isinstance(attr, unicode):
            del xds[attr]
            attr = str(attr)
        xds[attr] = xdata
    return xds

    return xds

def get_psocake_runs(exp):
    """
    Find psocake runs for experiment.
    """
    import os
    import glob
    instrument = exp[0:3]
    base_path = os.path.join('/reg/d/psdm',instrument,exp)
    
    files = glob.glob(base_path+'/*/*/psocake/r*/{:}_*.cxi'.format(exp))
    runs = [int(a.split('/')[-1].split('_')[1].split('.')[0]) for a in files]
    return sorted(runs)

def open_cxi_psocake(exp=None, run=None, folder=None, file_name=None, 
        load_summary=None, load_smd=None, load_moved_pvs=True,
        save=True, refresh=False, **kwargs):
    """
    load psocake cxidb formatted hdf5 file and return xarray
    Parameters
    ----------
    load_smd : bool
        load small data from 
    """
    import os
    import glob
    import xarray as xr
    import numpy as np
    import time
    if not exp:
        print('Experiment must be supplied as exp="xxx" or as first argument')
        return

    instrument = exp[0:3]
    base_path = os.path.join('/reg/d/psdm',instrument,exp)
    
    if not run:
        runs = get_psocake_runs(exp)
        if runs:
            print('No run specified:  Runs with psocake files include:')
            print(runs)
            print('Enter run number as run= or second argument')
        else:
            print('No psocake files exist for {:}'.format(exp))
        
        return

    file_nc=None
    if not folder:
        files = glob.glob(base_path+'/*/*/psocake/r{:04}/{:}_{:04}.nc'.format(run,exp,run))
        if files:
            file_nc = files[-1]
    else:
        file_nc = os.path.join(base_path, folder, 'r{:04}'.format(run), '{:}_{:04}.nc'.format(exp,run))
        if not os.path.isfile(file_nc):
            file_nc = None

    if not refresh and file_nc:
        try:
            xdata = xr.open_dataset(file_nc,engine='h5netcdf')
            print('Loading netcdf4 file {:}'.format(file_nc))
            return xdata
        except:
            print('Could not open netcdf4 file {:}'.format(file_nc))

    if not file_name:
        if not folder:
            files = glob.glob(base_path+'/*/*/psocake/r{:04}/{:}_{:04}.cxi'.format(run,exp,run))
            if files:
                print(files)
                file_name = files[-1]
                print('...choosing {:}'.format(file_name))
                print('use folder keyword to select a different psocake file')
            else:
                print('No psocake files available')
                return None
        else:
            file_name = os.path.join(base_path, folder, 'r{:04}'.format(run), '{:}_{:04}.cxi'.format(exp,run))

    if not file_nc:
        file_nc = file_name.rstrip('cxi')+'nc'

    xdata = open_cxi_dataset(file_name, **kwargs)
    xdata.attrs['experiment'] = exp
    xdata.attrs['run'] = run
    if load_summary is not False:
        time_last = time.time() 
        try:
            file_name = os.path.join(base_path, 'scratch', 'nc', 'run{:04}.nc'.format(run))
            print('... opening {:}'.format(file_name))
            xsmd = xr.open_dataset(file_name, engine='h5netcdf')
            if 'time_ns' not in xsmd.coords:
                xsmd.coords['time_ns'] = (('time'), np.int64(xsmd.sec*1e9+xsmd.nsec))
            # datetime64 not consistent at the us time level so need to recreate it
            # datetime64 is curriosly slow to calculate -- ~16 ms each time point
            print('... recalculating datetime64 for scratch summary data')
            #xsmd['time'] = [np.datetime64(int(sec*1e9+nsec), 'ns') for sec,nsec in zip(xsmd.sec,xsmd.nsec)]
            print('Load Time smd hdf5 = {:8.3f} sec'.format(time.time()-time_last))
            time_last = time.time() 
            xdata = xsmd.swap_dims({'time': 'time_ns'}).merge(xdata).swap_dims({'time_ns': 'time'})
            print('merge time psocake & smd = {:8.3f} sec'.format(time.time()-time_last))
            if load_smd is None:
                load_smd = False
        except:
            print('Cannot load and merge smd')

    if load_smd is not False:
        time_last = time.time() 
        try:
            file_name = os.path.join(base_path, 'results', 'nc', 'run{:04}_smd.nc'.format(run))
            print('... opening {:}'.format(file_name))
            xsmd = xr.open_dataset(file_name, engine='h5netcdf')
            if 'time_ns' not in xsmd.coords:
                xsmd.coords['time_ns'] = (('time'), np.int64(xsmd.sec*1e9+xsmd.nsec))
            # datetime64 not consistent at the us time level so need to recreate it
            # datetime64 is curriosly slow to calculate -- ~16 ms each time point
            print('... recalculating datetime64 for smd')
            #xsmd['time'] = [np.datetime64(int(sec*1e9+nsec), 'ns') for sec,nsec in zip(xsmd.sec,xsmd.nsec)]
            print('Load Time smd hdf5 = {:8.3f} sec'.format(time.time()-time_last))
            time_last = time.time() 
            xdata = xsmd.swap_dims({'time': 'time_ns'}).merge(xdata).swap_dims({'time_ns': 'time'})
            print('merge time psocake & smd = {:8.3f} sec'.format(time.time()-time_last))
        except:
            print('Cannot load and merge smd')


    if load_moved_pvs:
        try:
            xdata = add_moved_pvs(xdata)
            print('Adding moved pvs')
        except:
            print('Cannot add moved pvs')

    #if save:
    #    xdata.to_

    return xdata


def open_cxi_dataset(file_name, load_peakpos=False, add_time=False, **kwargs):
    """
    Read cxi format hdf5 file as xarray.  
    Currently skips much and is first intended to retrieve psocake results_1.nPeaks
    
    Paarameters:
    ------------
    load_peakpos : bool
        lood bragg peak positions in event images 
    add_time : bool
        add datetime64
    """
    import h5py
    import xarray as xr
    import pandas as pd
    import numpy as np
    import time
    time0 = time.time() 
    f5 = h5py.File(file_name, 'r')
    f5keys = f5.keys()

    xdata = xr.Dataset()

    attrs = {}
    try:
        attr = 'cxi_version'
        f5keys.remove(attr)
        xdata.attrs[attr] = f5.get(attr).value
    except:
        print('Error getting cxi_version')

    attr = 'psocake'
    if attr in f5.keys():
        f5keys.remove(attr)
        try:
            a = f5.get(attr+'/input')
            psocake_attrs = {b.split(' ')[0]: b.split(' ')[1] for b in a.value[0].split('\n') if b}        
        except:
            print('Error getting {:} attrs'.format(attr))
            psocake_attrs = {}

    status_attrs = {}
    attr = 'status'
    if attr in f5.keys():
        f5keys.remove(attr)
        try:
            for key in f5.get(attr).keys():
                status_attrs[attr+'_'+key] = f5.get(attr).get(key)
        except:
            print('Error getting {:} attrs'.format(attr))

    eventNumber = f5.get('LCLS').get('eventNumber').value
    nevents = len(eventNumber)
    print('Loading {:} events'.format(nevents))

    for name in f5keys:
        item = f5.get(name)
        for attr in item.keys():
            try:
                data = item.get(attr)
                if attr.startswith('data') and not load_data:
                    print('Skipping {:} {:} {:} event data'.format(name, attr))
                elif hasattr(data, 'keys'):
                    time_next = time.time()
                    for a in data.keys():
                        adata = data.get(a)
                        ashape = adata.shape
                        try:
                            if ashape[0] == nevents:
                                if len(ashape) == 2 and load_peakpos:
                                    val = adata.value
                                    xdata[attr+'_'+a] = (('time_ns','peak'), val)
                                elif len(ashape) == 1:
                                    val = adata.value
                                    if val[0].shape:
                                        print('Skipping {:} {:} {:} {:}'.format(name, attr, a, adata))
                                        continue
                                    xdata[attr+'_'+a] = (('time_ns'), val)
                                else:
                                    print('Skipping {:} {:} {:} {:}'.format(name, attr, a, adata))
                                    continue

                            else:
                                print('Skipping {:} {:} {:} {:}'.format(name, attr, a, adata))
                        except:
                            print('Error getting {:} {:} {:} event data'.format(name, attr, a))
                
                else:
                    try:
                        val = data.value
                        if isinstance(val, str) or len(val) == 1:
                            xdata.attrs[attr] = val
                        else:
                            xdata[attr] = (['time_ns',], val)
                    except:
                        print('Error getting {:} {:} {:} event data'.format(name, attr))

            except:
                print('Error getting {:} {:} event data'.format(name, attr))

    try:
        xdata.coords['sec'] = xdata.machineTime
        xdata.coords['nsec'] = xdata.machineTimeNanoSeconds
        xdata.coords['time_ns'] = np.int64(xdata.sec*1e9+xdata.nsec)
        if add_time:
            xdata['time'] = [np.datetime64(int(sec*1e9+nsec), 'ns') for sec,nsec in zip(xdata.sec,xdata.nsec)]
    except:
        print('Error making time from machintTime and machineTimeNanSeconds')

    print('Load Time psocake hdf5 = {:8.3f} sec'.format(time.time()-time0))
    return xdata

def add_transmission(xdata, exp=None, run=None):
    """
    Add transmission from epicsArch attenuation pvs

    Parameters
    ----------
    exp : str
        Experiment name
        Default use xdata.attrs['experiment'] or xdata.attrs['exp']
    run : int
        Run number
        Default use xdata.attrs['run']

    """
    from .exp_summary import get_exp_summary
    if not exp:
        if 'experiment' in xdata.attrs:
            exp = str(xdata.attrs['experiment'])
        elif 'exp' in xdata.attrs:
            exp = str(xdata.attrs['exp'])
        else:
            print('exp must be specified or in Dataset attrs')
            return
    print(exp)
    
    if not run:
        if 'run' in xdata.attrs:
            run = int(xdata.attrs['run'])
        else:
            print('run must be specified or in Dataset attrs')
            return

    es = get_exp_summary(exp) 
 
def add_moved_pvs(xdata, exp=None, run=None):
    """
    Add epics pvs that were moved during run to Dataset.
    - Moved pvs are automatically determined from epicsArch data
    loaded in PyDataSource.ExperimentSummary class.
    - merge_fill method used to match moved pvs using timestamp info

    Parameters
    ----------
    exp : str
        Experiment name
        Default use xdata.attrs['experiment'] or xdata.attrs['exp']
    run : int
        Run number
        Default use xdata.attrs['run']

    """
    from .exp_summary import get_exp_summary
    if not exp:
        if 'experiment' in xdata.attrs:
            exp = str(xdata.attrs['experiment'])
        elif 'exp' in xdata.attrs:
            exp = str(xdata.attrs['exp'])
        else:
            print('exp must be specified or in Dataset attrs')
            return
    print(exp)
    
    if not run:
        if 'run' in xdata.attrs:
            run = int(xdata.attrs['run'])
        else:
            print('run must be specified or in Dataset attrs')
            return

    es = get_exp_summary(exp) 
    xadd = es.get_scan_data(run)
    xdata = merge_fill(xdata, xadd, bfill=True)
    xdata.attrs['scan_pvs'] = xadd.data_vars.keys()
    return xdata

def merge_fill(xdata, xadd, bfill=True, 
        keep_attrs=True, keep_new_times=False): 
    """
    Merge Datasets.  Fill second Dataset forward then backward unless bfill=False
    Currently only 1D arrays with dim='time' are added.
    """
    import operator
    import xarray as xr
    if 'time_ns' not in xdata.coords:
        xdata.coords['time_ns'] = (('time'), xdata.time.values.tolist())
    da = xdata.swap_dims({'time': 'time_ns'}).reset_coords()[['time_ns']]
    if not keep_new_times:
        da.coords['_orig_data'] = (('time_ns'), xdata.time_ns > 0)

    if 'time_ns' not in xadd.coords:
        xadd.coords['time_ns'] = (('time'), xadd.time.values.tolist())

    if 'time_ns' not in xadd.dims:
        xadd = xadd.swap_dims({'time': 'time_ns'})

    xmerge = da.combine_first(xadd)
    xfill = xmerge.ffill(dim='time_ns')
    if bfill:
        xfill = xfill.bfill(dim='time_ns') 
    if not keep_new_times:
        xfill = xfill.where(xfill['_orig_data'] == 1, drop=True) 
        del xfill['_orig_data']
    xfill = xfill.swap_dims({'time_ns': 'time'})
    # merge and concat are not efficient when already same dims
    # may still be more clever way, but this is fast
    for attr, item in xfill.data_vars.items():
        if len(item.dims) == 1 and item.dims[0] == 'time':
            xdata[attr] = (('time'), item.values)
            for a, val in sorted(xadd[attr].attrs.items(),key=operator.itemgetter(0)):
                xdata[attr].attrs[a] = val

    if keep_attrs:
        for attr, val in xadd.attrs.items():
            if attr not in xdata.attrs:
                xdata.attrs[attr] = val

    return xdata

def dataset_fill(xdata, time, bfill=True):
    """
    Dataset with time and fill forward then backward
    """
    import xarray as xr
    import numpy as np
    if 'time_ns' not in xdata.coords:
        xdata.coords['time_ns'] = (('time'), xdata.time.values.tolist())
    
    if isinstance(time, np.datetime64):
        dim = 'time'
    else:
        dim = 'time_ns'
    xadd = xr.DataArray(time, dims=[dim])
    xadd = xadd.to_dataset(name='time_stamp')

    if 'time_ns' not in xadd.coords:
        xadd.coords['time_ns'] = (('time'), xadd.time.values.tolist())

    if 'time_ns' not in xadd.dims:
        xadd = xadd.swap_dims({'time': 'time_ns'})

    xmerge = xdata.swap_dims({'time': 'time_ns'}).merge(xadd).swap_dims({'time_ns': 'time'})
    xfill = xmerge.ffill(dim='time').bfill(dim='time') 

    return xfill

# Placeholder to make independent of PyDataSource.get_dataset
#def get_epics_dataset(exp=None, run=None, pvdict={}, fields=None,
#            meta_attrs = {'units': 'EGU', 'PREC': 'PREC', 'pv': 'name'},
#            tstart=None, tend=None, quiet=True, **kwargs):
#        import time
#        import numpy as np
#        import pandas as pd
#        import xarray as xr
#        import xarray_utils
#        if not fields:
#            fields={
#                    'description':            ('DESC', 'Description'), 
#                    'slew_speed':             ('VELO', 'Velocity (EGU/s) '),
#                    'acceleration':           ('ACCL', 'acceleration time'),
#                    'step_size':              ('RES',  'Step Size (EGU)'),
#                    'encoder_step':           ('ERES', 'Encoder Step Size '),
#                    'resolution':             ('MRES', 'Motor Step Size (EGU)'),
#                    'high_limit':             ('HLM',  'User High Limit'),
#                    'low_limit':              ('LLM',  'User Low Limit'),
#                    'units':                  ('EGU',  'Units'),
#        #            'device_type':            ('DTYP', 'Device type'), 
#        #            'record_type':            ('RTYP', 'Record Type'), 
#                    }
#        time0 = time.time()
#        time_last = time0
#       
#        if load_default:
#            trans_pvs = [a for a in _transmission_pvs.get('FEE', {}) if a.endswith('_trans')]
#            trans_pvs += [a for a in _transmission_pvs.get(ds.instrument, {}) if a.endswith('_trans')]
#            trans3_pvs = [a for a in _transmission_pvs.get('FEE', {}) if a.endswith('_trans3')]
#            trans3_pvs += [a for a in _transmission_pvs.get(ds.instrument, {}) if a.endswith('_trans3')]
#            pvdict.update(**_transmission_pvs.get('FEE',{}))
#            pvdict.update(**_transmission_pvs.get(ds.instrument,{}))
#        else:
#            trans_pvs = []
#            trans3_pvs = []
#
#        pvs = {alias: pv for alias, pv in pvdict.items() if configData._in_archive(pv)} 
#        
#        data_arrays = {} 
#        data_fields = {}
# 
#        for alias, pv in pvs.items():
#            data_fields[alias] = {}
#            dat = configData._get_pv_from_arch(pv, tstart, tend)
#            if not dat:
#                print('WARNING:  {:} - {:} not archived'.format(alias, pv))
#                continue
#            
#            try:
#                attrs = {a: dat['meta'].get(val) for a,val in meta_attrs.items() if val in dat['meta']}
#                for attr, item in fields.items():  
#                    try:
#                        field=item[0]
#                        pv_desc = pv.split('.')[0]+'.'+field
#                        if configData._in_archive(pv_desc):
#                            desc = configData._get_pv_from_arch(pv_desc)
#                            if desc:
#                                vals = {}
#                                fattrs = attrs.copy()
#                                fattrs.update(**desc['meta'])
#                                fattrs['doc'] = item[1]
#                                val = None
#                                # remove redundant data
#                                for item in desc['data']:
#                                    newval =  item.get('val')
#                                    if not val or newval != val:
#                                        val = newval
#                                        vt = np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns')
#                                        vals[vt] = val
#                               
#                                data_fields[alias][attr] = xr.DataArray(vals.values(), 
#                                                                coords=[vals.keys()], dims=['time'], 
#                                                                name=alias+'_'+attr, attrs=fattrs) 
#                                attrs[attr] = val
#         
#                    except:
#                        traceback.print_exc()
#                        print('cannot get meta for', alias, attr)
#                        pass
#                vals = [item['val'] for item in dat['data']]
#                if not vals:
#                    print('No Data in archive for  {:} - {:}'.format(alias, pv))
#                    continue
#
#                doc = attrs.get('description','')
#                units = attrs.get('units', '')
#                time_next = time.time()
#              
#                try:
#                    if isinstance(vals[0],str):
#                        if not quiet:
#                            print(alias, 'string')
#                        vals = np.array(vals, dtype=str)
#                    else:
#                        times = [np.datetime64(long(item['secs']*1e9+item['nanos']), 'ns') for item in dat['data']]
#                        dfs = pd.Series(vals, times).sort_index()
#                        dfs = dfs[~dfs.index.duplicated()]
#                        dfs = dfs[~(dfs.diff()==0)]
#                        vals = dfs.values
#                        dfs = dfs.to_xarray().rename({'index': 'time'})
#                        data_arrays[alias] = dfs 
#                        data_arrays[alias].name = alias
#                        data_arrays[alias].attrs = attrs
#                
#                except:
#                    traceback.print_exc()
#                    if not quiet:
#                        print('Error loadinig', alias)
#
#                if not quiet:
#                    try:
#                        print('{:8.3f} {:28} {:8} {:10.3f} {:4} {:20} {:}'.format(time_next-time_last, \
#                                        gias, len(vals), np.array(vals).mean(), units, doc, pv))
#                    except:
#                        print('{:8.3f} {:28} {:8} {:>10} {:4} {:20} {:}'.format(time_next-time_last, \
#                                        alias, len(vals), vals[0], units, doc, pv))
#            
#            except:
#                traceback.print_exc()
#                if not quiet:
#                    print('Error loading', alias)
#
#        xdata = xr.merge(data_arrays.values())
#        if trans_pvs:
#            da = xdata.reset_coords()[trans_pvs].to_array() 
#            xdata['trans'] = (('time'), da.prod(dim='variable'))
#            xdata['trans'].attrs['doc'] = 'Total transmission: '+'*'.join(trans_pvs)
#
#
#   
