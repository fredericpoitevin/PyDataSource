'''
Includes methods for creating tables for analyzing active detectors in an experiment
'''

def get_active_dict(files,print_on=False):
    ''' Creates a dictionary which maps variable names to corresponding lists of run numbers
    
    Parameters
    ----------
    files : list(str)
        List of file paths to .nc files where the variables are read
    print_on : bool, optional
        Decides whether the progress is printed on screen
    
    Returns
    -------
    datadict : dict(str, list(int))
        Dictionary mapping variable names to numbers of runs where the variable has been active
    runs : list(int)
        List of all run numbers
    '''
    import xarray as xr
    import sys
    
    datadict = {}
    runs = []
    det_alias = {}
    for i in range(len(files)):
        if print_on:
            sys.stdout.write('\r{}/{}'.format(i,len(files)))
            sys.stdout.flush()
        x = xr.open_dataset(files[i], engine='h5netcdf')
        data = list(x.data_vars.keys())
        for det in x.data_vars:
            alias = x[det].attrs.get('alias')
            if alias:
                det_alias[det] = str(alias)

        runs.append(x.attrs['run'])
        for key in data:
            datadict[key]=datadict.get(key, [])+[x.attrs['run']]
   

    if print_on:
        sys.stdout.write('\n')
        sys.stdout.flush()
    return datadict, runs, det_alias

def get_active_table(datadict,runs, det_alias={}):
    '''Takes the dictionary and list of runs and returns a pandas DataFrame,
    where rows are different runs and columns different variables,
    values 1 or 0 depending on whether the variable was active in the run
    
    Parameters
    ----------
    datadict : dict(str, list(int))
        Dictionary mapping variable names to numbers of runs where the variable has been active
    runs : list(int)
        List of all run numbers
        
    Returns
    -------
    df : DataFrame
        Table of 1s and 0s, showing whether a variable has been active in a run
    '''
    
    import pandas as pd
    import numpy as np
    
    # dictionary: maps variable names to pandas series, for creating the dataframe
    tabledict = {}
    runs.sort()

    for var in datadict:
#        if not var.endswith('present'):
#            continue
        alias = det_alias.get(var)
        if not alias or alias is 'None' or alias in list(tabledict.keys()):
            continue
        s = pd.Series(np.asarray([runs[i] in datadict[var] for i in range(len(runs))], dtype=int),
            runs)
        tabledict[alias]=s

    df = pd.DataFrame(tabledict)
    return df


def get_damage_table(files, datadict, runs, det_alias={}, print_on=False):
    ''' Creates a pandas DataFrame showing the percentages of valid, 
    damaged events for each variable in a run,
    where x-axis is variables and y-axis is runs
    
    Parameters
    ----------
    files : list(str)
        List of file paths 
    datadict : dict(str, list(int))
        Dictionary mapping variable names to numbers of runs where the variable has been active
        (this is only used to get the variable names)
    runs : list(int)
        List of all run numbers
    print_on : bool, optional
        Decides whether the progress is printed
    
    Returns
    -------
    table : DataFrame
        Table of percentages
    '''
    import xarray as xr
    import sys
    import pandas as pd
    import numpy as np
    
    tabledict = {} # Dict of run->Series, for creating the table

    for i in range(len(files)):
        if print_on:
            sys.stdout.write('\r{}/{}'.format(i,len(files)))
            sys.stdout.flush()
            
        ds = xr.open_dataset(files[i], engine='h5netcdf')
         
        vardict = {}
        tot = ds['time'].shape[0]
        # Use next line for events = events in drop_stats
        #vardict['events'] = tot 

        # Use next line for events = total events
        x = xr.open_dataset(files[i].split('_')[0]+'_smd.nc', engine='h5netcdf')
        vardict['events'] = x.dims['time']
        vardict['drop_events'] = tot
        for var in datadict.keys():
            var = str(var)
            alias = det_alias.get(var)
            if not alias or alias is 'None' or alias in list(vardict.keys()):
                continue
            try:
                arr = ds[var]
                while len(arr.shape)>1:
                    arr = arr[:,0] 
                # Remove second part of next line if accepting zero values
                validnumbers = np.asarray(1-np.isnan(arr)).sum()-(arr.shape[0]-np.count_nonzero(arr))
                
                vardict[alias] = 100*(1.-float(validnumbers)/float(tot))
            except KeyError:
                vardict[alias] = np.nan

        tabledict[ds.attrs['run']]=pd.Series(vardict)
    
    if print_on:
        sys.stdout.write('\n')
        sys.stdout.flush()
    table = pd.DataFrame(tabledict).T
    # deleting zero columns
    table = table.loc[:,table.astype(bool).sum(axis=0).astype(bool).tolist()]

    # arranging columns
    cols = list(table.columns.values)
    cols.remove('events')
    cols.insert(0,'events')
    table = table[cols]
    return table


