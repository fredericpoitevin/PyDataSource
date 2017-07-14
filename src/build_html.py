import time
import os
import logging
import operator
import traceback

import matplotlib as mpl
mpl.use('Agg')

from pylab import *
import pandas as pd
import xarray as xr

import output_html
import markup

pd.set_option('precision',2)
pd.set_option('max_colwidth', 66)
pd.set_option('display.width', 110)
plt.rcParams['axes.labelsize'] = 16 

default_path = 'stats'

def hover(hover_color="#ffff99"):
    return dict(selector="tr:hover",
                props=[("background-color", "%s" % hover_color)])

styles = [
    hover(),
    dict(selector="th", props=[("font-size", "100%"),
                               ("text-align", "center")]),
    dict(selector="caption", props=[("caption-side", "bottom")])
]

def get_run_size(run, start_path = '.'):
    tt = 0
    nn = 0
    runstr = 'r{:04.0f}'.format(run)
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            if 'xtc' in f.lower() and runstr in f.lower():
                fp = os.path.join(dirpath,f)
                tt += os.path.getsize(fp)
                nn += 1
    return nn, tt

def load_DataSource(**kwargs):
    import PyDataSource
    try:
        return PyDataSource.DataSource(**kwargs)
    except:
        print 'DataSource not available'

def to_html(x, **kwargs):
    """
    Make default html run summary
    """
    return Build_html(x, auto=True, **kwargs)


# Make tidy pandas dataframe
#a = x.dims.keys()
#a.remove('time')
#df = x.drop(a).to_dataframe()

class Build_experiment(object):
    """Class to build an html RunSummary report.
    """

    def __init__(self, exp, auto=True, **kwargs):
        import psutils 
        if psutils.interactive_mode():
            plt.ioff()
        
        self.logger = logging.getLogger(__name__+'.build_html')
        self.logger.info(__name__)
        if exp.__class__.__module__ == 'PyDataSource.exp_summary': 
            self._es = exp
            self.exp = self._es.exp
        else:
            self.exp = exp
            self._es = self._get_exp_summary()
       
        self.instrument = self.exp[0:3]
        self.results = {}
        self.title = title='{:}'.format(self.exp)
      
        self.exp_dir = os.path.join('/reg/d/psdm/',self.instrument,self.exp)
        self.xtc_dir = os.path.join(self.exp_dir,'xtc')
        if os.path.isdir(os.path.join(self.exp_dir,'res')): 
            self.res_dir = os.path.join(self.exp_dir,'res')
        else:
            self.res_dir = os.path.join(self.exp_dir,'results')
        
        self.scratch_dir = os.path.join(self.exp_dir,'scratch')
        self.stats_dir = os.path.join(self.exp_dir,'stats')

        self._init_output(**kwargs)
        if auto:
            # assume first 2% are setup
            run_min = int(self._es.xruns.run.max()*.01)
            self.add_defaults(run_min=run_min, **kwargs)
            self.to_html()
 
    def _get_exp_summary(self, **kwargs):
        from exp_summary import get_exp_summary
        return get_exp_summary(self.exp, **kwargs)

    def add_defaults(self, run_min=None, run_max=None, **kwargs):
        """
        Add defaults.
        """
        self.logger.info('Adding defaults')
        try:
            self.add_scans()
        except:
            traceback.print_exc()

        try:
            #self.add_xy_ploterr('photon_beam_energy', 'run', title='Photon Beam Energy')
            #self.add_xy_ploterr('photon_current', 'run', title='Photon Current')
            for attr in self._es._machine_coords:
                self.add_xy_ploterr(attr, 'run', title=attr.replace('_',' '), catagory='Summary')
            
            df = self._es.get_scans()
            if df is not None:
                for attr in df.keys():
                    self.add_xy_ploterr(attr, 'run', title=attr.replace('_',' '), catagory='Scan Summary')

        except:
            traceback.print_exc()

        run_sets = self._es.get_run_sets(run_min=run_min, run_max=run_max, **kwargs)
        if len(run_sets) > 3:
                self.logger.info('Adding moved ')
                self.add_epics(run_min=run_min, run_max=run_max, **kwargs)
        else:
            for run_set in run_sets:
                self.logger.info('Adding moved {:}'.format(run_set))
                self.add_epics(run_min=run_set[0], run_max=run_set[1], **kwargs)

        try:
            self.add_last_set()
        except:
            traceback.print_exc()
            self.logger.info('Failed adding last set')

    def add_table(self, df, catagory=None, tbl_type=None, name=None, howto=[], doc=[], hidden=False): 
        if not catagory:
            catagory = 'Tables'
        if not tbl_type:
            tbl_type = ', '.join(a.columns)
        if not name:
            name = '__'.join(a.columns)
        
        self._add_catagory(catagory)
        self.results[catagory]['table'].update({tbl_type: 
                                                   {'DataFrame': df, 
                                                    'name': name,
                                                    'howto': howto, 
                                                    'hidden': hidden, 
                                                    'doc': doc}})


    def add_last_set(self, **kwargs):
        es = self._es
        df = es.last_set(**kwargs) 
   
        alias = 'Epics'
        self._add_catagory(alias)
        doc = []
        doc.append('Last time Epics PV was set')
        howto = []
        if kwargs:
            howto.append("es.last_set({:})".format(kwargs))
        else:
            howto.append("es.last_set({:})".format(kwargs))
        
        self.results[alias]['table'].update({'Last Set': 
                                                   {'DataFrame': df, 
                                                    'name': 'dflastset',
                                                    'howto': howto, 
                                                    'hidden': False, 
                                                    'doc': doc}})


    def add_epics(self, dets=None, run_min=None, run_max=None, group='units', **kwargs):
        """
        Add epics.
        """
        es = self._es
        df = es.xset.to_dataframe()
        if dets:
            if not isinstance(dets, list):
                dets = [dets]
        else:
            #dets = list(set([a.split('_')[0] for a in es.xset.data_vars.keys()]))
            dets = es.get_moved(run_min=run_min, run_max=run_max, group=group)
    
        alias = 'Epics'
        if run_min:
            alias += ' from Run {:}'.format(run_min)
        if run_max:
            alias += ' to Run {:}'.format(run_max)

        self.logger.info('Adding {:}'.format(alias))
        for det in dets:
            howto=['es.plot_move({:}, run_min={:}, run_max={:})'.format(det, run_min, run_max)]
            print 'Try plot_move', det
            if isinstance(dets, dict):
                attrs = dets[det]
            else:
                attrs = det

            result = self._es.plot_move(attrs, run_min=run_min, run_max=run_max)
            
            if result is not None:
                self.logger.info('Adding {:} plot for {:}'.format(alias, det))
                self.add_plot(alias, det+' Data', howto=howto)

            dfattrs = pd.DataFrame({attr: es.xset[attr].attrs for attr in attrs if attr in es.xset})
            self.results[alias]['table'].update({det+' Configuration': 
                                                       {'DataFrame': dfattrs, 
                                                        'name': 'attrs',
                                                        'howto': None, 
                                                        'hidden': True, 
                                                        'doc': None}})

    def add_scans(self, dets=None, min_steps=4, **kwargs):
        """
        Add summary of scans.
        """
        es = self._es
        df = es.get_scans(min_steps=min_steps)
        if df is None:
            return

        if dets:
            if not isinstance(dets, list):
                dets = [dets]
        else:
            dets = list(set([a.split('_')[0] for a in df.keys()]))

        alias = 'Scans'
        for det in dets:
            attrs = [a for a in df.keys() if a.startswith(det)]
            dfscan = df.where(df[attrs].T.sum() > min_steps).dropna()
            dfscan = dfscan.T.where(dfscan.sum() > min_steps).dropna().T.astype(int)
           
            self._add_catagory(alias)
            doc = []
            doc.append('Scans involving {:} devices'.format(det))
            howto = []
            howto.append("attrs={:})".format(attrs))
            howto.append("xattrs = pd.DataFrame({attr: self._es.xscan[attr].attrs for attr in attrs})")
            howto.append("min_steps={:}".format(min_steps))
            howto.append("dfscan = df.where(df[attrs].T.sum() > min_steps).dropna()")
            howto.append("dfscan = dfscan.T.where(dfscan.sum() > min_steps).dropna().T")
            self.results[alias]['table'].update({det: 
                                                       {'DataFrame': dfscan, 
                                                        'name': 'dfscan',
                                                        'howto': howto, 
                                                        'hidden': False, 
                                                        'doc': None}})
            doc='Number of times each {:} device was moved (or state changed) during run.'.format(det)
            dfattrs = pd.DataFrame({attr: es.xscan[attr].attrs for attr in attrs if attr in es.xscan})
            self.results[alias]['table'].update({det+' Configuration': 
                                                       {'DataFrame': dfattrs, 
                                                        'name': 'attrs',
                                                        'howto': None, 
                                                        'hidden': True, 
                                                        'doc': doc}})

            for run in dfscan.index:
                self._es.plot_scan(run, attrs=attrs)
                howto=['es.plot_scan({:}, attrs={:})'.format(run, attrs)]
                self.add_plot(alias, 'Run {:03} {:}'.format(run, det), howto=howto )

    def _add_catagory(self, catagory):
        if catagory not in self.results:
            self.results[catagory] = {'figure': {}, 'table': {}, 'text': {}, 'textblock': {}}
    
    def add_plot(self, catagory, plt_type, howto=[], show=False):
        """Add a plot to RunSummary.
        """
        try:
            self._add_catagory(catagory)

            plt_file = '{:}_{:}.png'.format(catagory, plt_type).replace(' ','_') 
            self.results[catagory]['figure'].update({plt_type: {'path': self.output_dir, 
                                                             'png': plt_file,
                                                             'howto': howto}})
            plt.savefig(os.path.join(self.output_dir, plt_file))
            if show:
                plt.show()
            else:
                plt.close()

        except:
            traceback.print_exc()
            self.logger.info('Could not add Plot {:} {:}'.format(catagory, plt_type))

    def _init_output(self, path=default_path, filename='experiment', **kwargs):
        """Set output path and build appropriate directories for html
        """
        if not path or path == 'stats':
            path = os.path.join(self.stats_dir, 'summary')
        elif path == 'res':
            path = os.path.join(self.res_dir, 'RunSummary')
        elif path == 'home':
            path = os.path.join(os.path.expanduser('~'), 'RunSummary', self.exp)
        elif path == 'scratch':
            path = os.path.join(self.scratch_dir, 'RunSummary')
        
        print 'Setting output path to', path

        self.path = path
        self.filename = filename

        self.output_dir = os.path.join(self.path, self.filename) 
       
        if not os.path.isdir(self.path):
            os.mkdir(self.path)

        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
 
    def to_html(self, path=None, quiet=False, **kwargs):
        """Write out html file
        """
        self._init_html(path=path)
        self._add_html()
        self._close_html()
        if not quiet:
            print 'Writing html to: ', os.path.join(self.html.output_dir, self.html.output_file)

    def _init_html(self, path=None, **kwargs):
        """Initialize html page.
        """
        if path:
            self._init_output(path=path, **kwargs)

        self.html = output_html.report(self.exp, 'Experiment Summary', 
                 title=self.title,
                 css=('css/bootstrap.min.css','jumbotron-narrow.css','css/mine.css'),
                 script=('js/ie-emulation-modes-warning.js','js/jquery.min.js','js/toggler.js','js/sticky.js'),
                 output_dir=self.output_dir)

        self.html.start_block('Experiment Summary', id="metadata")
        self.html.start_subblock('Experiment Information',id='datatime')
        
        try:
            event_times = pd.to_datetime(self._es.xruns.begin_time.values)

            sformat = 'First Run: {:}<br/>Last Run: {:}'
            self.html.page.p('Number of Runs: {:}'.format(len(event_times) ) )
            self.html.page.p(sformat.format(event_times.min(), event_times.max()) )
        except:
            pass
            
        self.html.page.p('Report time: {:}'.format(time.ctime()))
        self.html.end_subblock()

        self._make_ExperimentInfo_html()

        self.html.end_block()         
 
    def _make_ExperimentInfo_html(self, **kwargs):
        """Make html notes for accessing ExperimentInfo.
        """
        self.html.start_subblock('Access the Experiment Information')
        self.html.page.p('Access overall experiment information:')
        text = ['import PyDataSource']
        text.append('es = PyDataSource.ExperimentInfo("{:}")'.format(self.exp))
        self.html.add_textblock(text)
        self.html.end_subblock()

        self.html.page.p('For questions and feedback contact koglin@slac.stanford.edu')
 
    def _add_html(self, table_caption='Hover to highlight.', show_attrs=['attrs'], **kwargs):
        """Add html pages for results.
        """
        cat_items = sorted(self.results.items(), key=operator.itemgetter(0))
        for catagory, item in cat_items:
            self.html.start_block('{:} Data'.format(catagory), id='{:}_data'.format(catagory))
            ahowto = []
            
            if item['figure']:
                ahowto.append('# For interactive plotting -- plt.ioff() to turn off interactive plotting.')
                ahowto.append('plt.ion()')
                ahowto.append('# Alternatively make plt.show() after each plot and close window to make next')

            datatyp = 'text'
            data_items = sorted(item[datatyp].items(), key=operator.itemgetter(0))
            for name, data in data_items:
                howto = data.get('howto')
                if howto is not None:
                    if not isinstance(howto, list):
                        howto = [howto]
                    howto_step = "# Howto {:} {:}:\n".format(name, catagory)
                    howto_step += '\n'.join(howto)
                    ahowto.append(howto_step)

            datatyp = 'table'
            data_items = sorted(item[datatyp].items(), key=operator.itemgetter(0))
            for name, data in data_items:
            #for name, data in item[datatyp].items():
                howto = data.get('howto')
                if howto is not None:
                    if not isinstance(howto, list):
                        howto = [howto]
                    howto_step = "# Howto make the {:} {:} {:}:\n".format(catagory, name, datatyp)
                    howto_step += '\n'.join(howto)
                else:
                    howto_step = ''
                
                if data.get('format'):
                    formatters = {a: b.format for a,b in data.get('format').items()}
                    formatterstr = {a: b+'.format' for a,b in data.get('format').items()}
                else:
                    formatters = None

                df = data.get('DataFrame')
                doc = data.get('doc')
                dfname = data.get('name')
                if df is not None:
                    if name in show_attrs or data.get('hidden') is False:
                        hidden = False
                    else:
                        hidden = True
                   
                    pd.set_option('display.max_rows', len(df))
                    if formatters:
                        dfstr = df.to_string(justify='left',formatters=formatters)
                        #if dfname:
                        #    howto_step += '\n# to reprsent with formatting'
                        #    howto_step += "\nprint {:}.to_string(justify='left',formatters={:})".format(dfname, formatterstr)
                    else:
                        dfstr = str(df)
                        #if dfname:
                        #    howto_step += "\n print {:}".format(dfname)
                    
                    if howto_step:
                        ahowto.append(howto_step)

                    self.html.add_textblock(dfstr, doc=doc, 
                            subblock='{:} {:} {:}'.format(catagory,name,datatyp), 
                            hidden=hidden)
                    
                    pd.reset_option('display.max_rows')

            datatyp = 'figure'
            data_items = sorted(item[datatyp].items(), key=operator.itemgetter(0))
            for name, data in data_items:
            #for name, data in item[datatyp].items():
                png = data.get('png')
                if png:
                    self.html.start_subblock('{:} {:} {:}'.format(catagory, name, datatyp), hidden=True)
                    self.html.page.a( markup.oneliner.img(src=png,style='width:100%;'), 
                            href=png )
                    self.html.end_subblock(hidden=True)
            
                howto = data.get('howto')
                if howto:
                    if not isinstance(howto, list):
                        howto = [howto]
                    howto_step = "# Howto make the {:} {:} {:}:\n".format(catagory, name, datatyp)
                    howto_step += '\n'.join(howto)
                    ahowto.append(howto_step)

            if ahowto:
                self.html.add_textblock(ahowto, 
                        subblock='HowTo make {:} tables and figures.'.format(catagory), 
                        id='howto_{:}'.format(catagory.replace(' ','_')), 
                        hidden=True)

            self.html.end_block()         

    def _close_html(self, **kwargs):
        """Close html files.
        """
        # this closes the left column
        self.html.page.div.close()
        
        self.html.mk_nav()
        
        self.html._finish_page()
        self.html.myprint(tofile=True)

    def add_xy_ploterr(self, attr, xaxis=None, howto=None, catagory=None, text=None, table=None, **kwargs):
        from plotting import xy_ploterr
        x = self._es.xscan
        if not howto:
            howto = []
        if not catagory:
            catagory = x[attr].attrs.get('alias', 'Summary')
        
        howto.append("Custom plot see PyDataSource.plotting.xy_ploterr method")
        p = xy_ploterr(x, attr, xaxis=xaxis, **kwargs)
        if not xaxis:
            xaxis=x.scan_variables[0]
        plt_type = '{:} vs {:}'.format(attr, xaxis)   
        if kwargs.get('logy'):
            if kwargs.get('logx'):
                plt_type+=' log-log'
            else:
                plt_type+=' log scale'
        elif kwargs.get('logx'):
            plt_type+=' lin-log'

        print catagory, plt_type
        self.add_plot(catagory, plt_type, howto=howto)
        if table is not None:
            self.results[catagory]['table'].update({attr:{'DataFrame': table, 
                                                        'name': 'df_tbl',
                                                        'howto': [], 
                                                        'doc': text}})


class Build_html(object):
    """Class to build an html RunSummary report.
    """

    def __init__(self, xdat=None, ds=None, auto=False, **kwargs):
       
        import psutils 
        if psutils.interactive_mode():
            plt.ioff()
        
        self.logger = logging.getLogger(__name__+'.build_html')
        self.logger.info(__name__)
       
        self.results = {}

        if xdat:
            if xdat.__class__.__name__ == 'DataSource':
                # 1st arg is actually PyDataSource.DataSource
                ds = xdat
                xdat = None

        elif not ds:
            if 'run' in kwargs and 'exp' in kwargs:
                self.logger.info('Loading PyDataSource data')
                ds = DataSource(**kwargs)
            else:
                print 'Require valid xdat=xarray object and/or ds=PyDataSource object'
                print 'or alternatively exp and run kwargs'
                return

        if xdat:
            self._xdat = xdat
        else:
            self.logger.info('Building xarray data')
            self._xdat = ds.to_xarray(**kwargs)
           
        self.exp = str(self._xdat.experiment)
        self.instrument = str(self._xdat.instrument)
        try:
            self.run = int(self._xdat.run.values[0])
        except:
            self.run = self._xdat.run

        if ds:
            self._ds = ds
        else:
            self.logger.info('Loading PyDataSource data')
            self._ds = load_DataSource(exp=self.exp, run=self.run)

        self.title = title='{:} Run {:}'.format(self.exp ,self.run)
      
        self.exp_dir = os.path.join('/reg/d/psdm/',self.instrument,self.exp)
        self.xtc_dir = os.path.join(self.exp_dir,'xtc')
        if os.path.isdir(os.path.join(self.exp_dir,'res')): 
            self.res_dir = os.path.join(self.exp_dir,'res')
        else:
            self.res_dir = os.path.join(self.exp_dir,'results')
        
        self.scratch_dir = os.path.join(self.exp_dir,'scratch')
        self.stats_dir = os.path.join(self.exp_dir,'stats')
        
        self.nc_file = '{:}/nc/run{:04}.nc'.format(self.scratch_dir, self.run)
       
        self.aliases = list(set([item.attrs.get('alias') for attr,item in self._xdat.data_vars.items() if item.attrs.get('alias')]))

        if 'steps' in self._xdat:
            self.nsteps = len(self._xdat.steps)
        else:
            self.nsteps = None

        if 'codes' in self._xdat:
            self.ncodes = len(self._xdat.codes)
        else:
            self.ncodes = None

        self._init_output(**kwargs)
        if auto:
            self.add_all(**kwargs)
            self.to_html()

    def _init_output(self, path=default_path, filename=None, **kwargs):
        """Set output path and build appropriate directories for html
        """
        if not path or path == 'stats':
            path = os.path.join(self.stats_dir, 'summary')
        elif path == 'res':
            path = os.path.join(self.res_dir, 'RunSummary')
        elif path == 'home':
            path = os.path.join(os.path.expanduser('~'), 'RunSummary', self.exp)
        elif path == 'scratch':
            path = os.path.join(self.scratch_dir, 'RunSummary')
        
        print 'Setting output path to', path

        self.path = path
        if not filename:
            filename = 'run{:04}'.format(self.run)

        self.filename = filename

        self.output_dir = os.path.join(self.path, self.filename) 
       
        if not os.path.isdir(self.path):
            os.mkdir(self.path)

        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
 
    def add_correlations(self, cut=None, confidence=0.4, **kwargs):
        from xarray_utils import get_correlations
        from xarray_utils import heatmap
        
        x = self._xdat
        catagory = 'PhotonEnergy'
        if catagory in x:
            attrs = get_correlations(x, catagory, confidence=confidence).keys()
            self.add_detector('EBeam', catagory=catagory, cut=cut,
                    make_scatter=attrs, make_timeplot=False, make_histplot=False)
            attrs.append('PhotonEnergy')
            if 'Gasdet_post_atten' in x:
                attrs.append('Gasdet_post_atten')
            df = x.reset_coords()[attrs].to_dataframe()
            heatmap(df)
            plt_type = 'correlation'
            howto = []
            self.add_plot(catagory, plt_type, howto=howto, tight=False)


        alias = 'Gasdet_post_atten'
        catagory = 'PulseEnergy'
        if alias in x:
            attrs = get_correlations(x, 'Gasdet_post_atten', confidence=confidence).keys()
            self.add_detector('FEEGasDetEnergy', catagory=catagory, cut=cut, 
                    make_scatter=attrs, make_timeplot=False, make_histplot=False)
            attrs.append('Gasdet_post_atten')
            if 'PhotonEnergy' in x:
                attrs.append('PhotonEnergy')
            df = x.reset_coords()[attrs].to_dataframe()
            heatmap(df)
            plt_type = 'correlation'
            howto = []
            self.add_plot(catagory, plt_type, howto=howto, tight=False)


    def add_all(self, quiet=True, 
                        min_step_data=5, 
                        min_in_cut=4,
                        **kwargs):
        """Add All detectors.
        """
        #self.make_default_cuts()
        
        if 'cuts' in self._xdat.attrs:
            cuts = self._xdat.attrs.get('cuts')
            if not isinstance(cuts, list):
                try:
                    cuts = list(cuts)
                except:
                    cuts = [cuts]
        else:
            cuts = []
        
        # Add RunStats
        try:
            self.add_detstats()
        except:
            print 'cannot add detstats'
        try:
            self.add_config()
        except:
            print 'cannot add config'
       
        if cuts:
            print cuts
            for cut in cuts:
                if cut in self._xdat and self._xdat[cut].sum() >= min_in_cut: 
                    print 'adding', cut
                    try:
                        self.add_correlations(cut=cut)
                    except:
                        traceback.print_exc()
                        print 'Cannot add correlelations for cut', cut
        else:
            try:
                self.add_correlations()
            except:
                traceback.print_exc()
                print 'Cannot add correlelations'

        nevents = len(self._xdat.time)

        groupby = [group for group in self._xdat.attrs.get('scan_variables') \
                            if group in self._xdat and len(set(self._xdat[group].data)) > min_step_data]
        print ''
        print '{:24} {:8}'.format('Group','points')
        print '-'*40
        for group in groupby:
            print '{:24} {:>8}'.format(group, len(set(self._xdat[group].data)))
        print ''

        for alias in self.aliases:
            if not quiet:
                print 'adding detector', alias
           
            if cuts:
                for cut in cuts:
                    if cut in self._xdat and self._xdat[cut].sum() >= min_in_cut: 
                        self.add_detector(alias, cut=cut, groupby=groupby, **kwargs)
            else:
                self.add_detector(alias, groupby=groupby, **kwargs)

            attrs = [attr for attr,item in self._xdat.data_vars.items() \
                     if item.attrs.get('alias') == alias and 'time' in item.dims and len(item.dims) > 1]
            if attrs:
                try:
                    self.add_summary(attrs, groupby=False)
                except:
                    traceback.print_exc()
                    print 'Could not add summary', attrs
                if groupby:
                    try:
                        self.add_summary(attrs, groupby=groupby)
                    except:
                        traceback.print_exc()
                        print 'Could not add summary', attrs

        # Add sum/count info
        x = self._xdat
        attrs = [a for a in x.keys() if (a.endswith('_count') or a.endswith('_sum')) 
                    and len(x[a].dims) == 1 and 'time' in x[a].dims]
        if 'PhotonEnergy' in x:
            attrs.append('PhotonEnergy')
        if 'Gasdet_post_atten' in x:
            attrs.append('Gasdet_post_atten')
        self.add_detector(attrs=attrs, catagory='Detector Count', confidence=0.1)
        
        attrs = [a for a in x.keys() if (a.endswith('_PosX') or a.endswith('_AngX') or a.endswith('xpos'))
                    and len(x[a].dims) == 1 and 'time' in x[a].dims]
        if attrs:
            if 'PhotonEnergy' in x:
                attrs.append('PhotonEnergy')
            if 'Gasdet_post_atten' in x:
                attrs.append('Gasdet_post_atten')
            self.add_detector(attrs=attrs, catagory='Xaxis', confidence=0.1)
            
        attrs = [a for a in x.keys() if (a.endswith('_PosY') or a.endswith('_AngY') or a.endswith('ypos'))
                    and len(x[a].dims) == 1 and 'time' in x[a].dims]
        if attrs:
            if 'PhotonEnergy' in x:
                attrs.append('PhotonEnergy')
            if 'Gasdet_post_atten' in x:
                attrs.append('Gasdet_post_atten')
            self.add_detector(attrs=attrs, catagory='Yaxis', confidence=0.1)

        xdat = self._xdat
        if 'codes' in xdat.coords and len(xdat['codes']) > 0:
            attrs = [a for a in self._xdat.data_vars.keys() if 'stat' in self._xdat[a].dims and len(self._xdat[a].dims) in [5,6]]
            for attr in attrs:
                print 'adding stats for ', attr
                if nevents > 4:
                    self.add_stats(attr, **kwargs)

        attrs = [a for a in self._xdat.data_vars.keys() if 'time' in self._xdat[a].dims and len(self._xdat[a].dims) == 3]
        if nevents <= 100:
            self.add_event(attrs)
        plt.cla()
        plt.close('all')

    def make_default_cuts(self, gasdetcut_mJ=0.5):
        """
        Make default cuts.
        """
        import h5write
        try:
            self._xdat = h5write.make_default_cuts(self._xdat)
        except:
            traceback.print_exc()
            print 'Cannot make default cuts'

    def make_summary(self):
        """Make summary
        """
        import PyDataSource
        self._xsummary = PyDataSource.to_summary(self._xdat)

    def add_config(self, attrs=[]):
        """
        Add detector configurations
        """
        howto = []
        doc = ''
        if not attrs:
            attrs = [a for a in self._xdat.variables.keys() if a.endswith('present')]
        for attr in attrs:
            alias = self._xdat[attr].attrs.get('alias')
            self._add_catagory(alias)
            data = []
            for key, val in self._xdat[attr].attrs.items():
                data.append('{:20} {:<50}'.format(key, val))
            self.results[alias]['textblock'].update({'config':{'text': data, 
                                                        'name': 'textblock_text',
                                                        'howto': howto, 
                                                        'doc': doc}})

    def add_detstats(self, alias='RunStats'):
        """Add detector statistics. 
        """
        attrs = [a for a in self._xdat.variables.keys() if a.endswith('events')]
        if attrs and self.nsteps == 1:
            self._add_catagory(alias)
            df_tbl = self._xdat.reset_coords()[attrs].sel(steps=0).to_array().to_pandas()
            self.eventStats = df_tbl
            doc = []
            doc.append("A summary of the mean, stdev, min and max values were created "
                      +"for each eventCode for the detectors in this table:")
            howto = []
            howto.append("attrs={:})".format(attrs))
            howto.append("df_tbl = x[attrs].sel(steps=0).to_array().to_pandas()")
            self.results[alias]['table'].update({alias:{'DataFrame': df_tbl, 
                                                        'name': 'df_tbl',
                                                        'howto': howto, 
                                                        'doc': doc}})

    def add_detector(self, alias='FEEGasDetEnergy', percentiles=[0.05,0.50,0.95], 
                           attrs=None,
                           catagory=None, 
                           figsize=None, 
                           layout=None,
                           groupby=None,
                           scat_name=None,
                           make_scatter=None,
                           make_table=True,
                           make_timeplot=None,
                           make_histplot=None,
                           make_correlation=True,
                           min_step_data=5,
                           confidence=0.5,
                           labelsize=20,
                           robust_attrs=None,
                           plot_errors=False,
                           bins=50, 
                           plt_style='.',
                           cut=None, 
                           show=False):
        """Add a single detector based on alias
        """
        import seaborn as sns
        x = self._xdat
        nevents = x.time.size
        if nevents < 10:
            if make_scatter is None:
                make_scatter = False
            if make_timeplot is None:
                make_timeplot = False
            if make_histplot is None:
                make_histplot = False
        elif nevents > 100000:
            if make_scatter is None:
                make_scatter = False
        else:
            if make_scatter is None:
                make_scatter = True
            if make_timeplot is None:
                make_timeplot = True
            if make_histplot is None:
                make_histplot = True 

        #if 'Damage_cut' in x:
        #    x = x.where(x.Damage_cut == 1)

        plt.rcParams['axes.labelsize'] = labelsize 
        desc_attrs = ['unit','doc']
        default_scatter_attrs = {
                'EBeam': ['ebeamCharge', 'ebeamDumpCharge','ebeamUndPosX', 'ebeamUndPosY']
                #'EBeam': ['ebeamCharge', 'ebeamDumpCharge','ebeamL3Energy', 'ebeamPhotonEnergy', 'ebeamUndPosX', 'ebeamUndPosY']
                #'EBeam': ['ebeamCharge', 'ebeamL3Energy', 'ebeamPhotonEnergy', 'ebeamXTCAVAmpl', 'ebeamXTCAVPhase']
                }

#        print default_scatter_attrs

        if groupby is None:
            groupby = []
            for gattr in x.attrs.get('scan_variables'):
                if gattr in x.reset_coords().data_vars and len(set(x[gattr].data)) > min_step_data:
                    groupby.append(gattr)

            if groupby is None and self.nsteps > 4:
                groupby = ['step']

        elif not isinstance(groupby, list):
            groupby = [groupby]

        #if groupby and not np.shape(groupby):
        if groupby is not None:
            groupby = list(set(groupby))

        attr_names = None
        if not attrs:
            attrs = [attr for attr,item in x.data_vars.items() \
                     if item.attrs.get('alias') == alias and item.dims == ('time',)]
            attrs0 = [attr for attr in attrs]
            for attr in x.attrs.get('correlation_variables', []):
                attrs.append(attr)
            if isinstance(make_scatter, list):
                for attr in make_scatter:
                    attrs.append(attr)
        else:
            if isinstance(attrs, dict):
                attr_names = attrs
                attrs = attr_names.values()
            
            attrs0 = [attr for attr in attrs]
       
        if not attrs0:
            print 'No scalar data for ', alias
            return

        if groupby:
            for group in groupby:
                attrs.append(group)
                if group in x.coords:
                    x = x.reset_coords(group)
        if not attr_names:
            attr_names = {attr: attr.replace(alias+'_','') for attr in attrs if attr in x}

        if not catagory:
            catagory = alias
#        if cut is None:
#            if 'Damage_cut' not in x:
#                self.make_damage_cut()
#            cut = 'Damage_cut'
        
        if cut:
            catagory = catagory+' '+cut

        if attrs:
            self._add_catagory(catagory)
            nattrs = len(attrs)
            if nattrs == 4:
                nrows = 2
                ncolumns = 2
            elif nattrs == 2:
                nrows = 2
                ncolumns = 1
            else:
                ncolumns = int(min([nattrs,3]))
                nrows = int(np.ceil(nattrs/float(ncolumns)))

            if not layout:
                layout = (nrows,ncolumns)
            if not figsize:
                figsize = (max([ncolumns*4,8]),max([nrows*3.0,6]))
            
            xselect = x[attrs]
            tselect = 'time'
            if cut:
                tselect = 'points'
                itimes = xselect.groupby(cut).groups
                if len(itimes) > 1:
                    xselect = xselect.isel_points(time=itimes[1])

            try:
                df = xselect.to_array().to_pandas().T.dropna()
                
                howto = ["plt.rcParams['axes.labelsize'] = {:}".format(labelsize)]
                howto.append("attrs = {:}".format(attrs))
                if cut:
                    howto.append("itimes = x[attrs].groupby('{:}').groups[1]".format(cut))
                    howto.append("xselect = x[attrs].isel_points(times=itimes)")
                else:
                    howto.append("xselect = x[attrs]")

                howto.append("df = xselect.to_array().to_pandas().T")

                df_attrs = pd.DataFrame({attr_names[attr]: {a: x[attr].attrs.get(a) for a in desc_attrs} \
                                         for attr in attrs}).T[desc_attrs]
                
                # Need to add in howto print formatted df_attrs, but tricky
                self.results[catagory]['table'].update({'attrs': {'DataFrame': df_attrs, 'name': 'df_attrs',
                            'format': {'unit': '{:<10s}', 'doc': '{:<69s}'}}})
                
                df.rename(inplace=True, columns=attr_names)
                howto.append("attr_names={:}".format(attr_names))
                howto.append("df.rename(inplace=True, columns=attr_names)")
                self.results[catagory]['text'].update({'setup':{'howto': howto}})
                
                #df_attrs.to_string(formatters={'unit': '{:<10s}'.format, 'doc': '{:<60s}'.format},justify='left')

                df_tbl = df.describe(percentiles=percentiles).T.round({'count':0})
                howto = ["df_tbl = df.describe(percentiles={:}).T.round({:})".format(percentiles,{'count':0})]
                howto.append("print df_tbl")

                # make table using pandas
                if make_table:
                    self.results[catagory]['table'].update({'stats':{'DataFrame': df_tbl, 'name': 'df_tbl', 'howto': howto}})
               
                # make time plots
                if make_timeplot:
                    plt_type = 'time'
                    df.plot(subplots=True, sharex=True, style=plt_style, layout=layout, figsize=figsize)
                    #plt.tight_layout()
                    howto = ["df.plot(subplots=True, sharex=True, style='{:}', layout={:}, figsize={:})".format(plt_style, layout, figsize)]
                    self.add_plot(catagory, plt_type, howto=howto)

            except:
                traceback.print_exc()
                print 'Could not build df to describe ', attrs
                print xselect

            # Make groupby plots
            if groupby:
                if not isinstance(groupby, list):
                    groupby = [groupby]
                for grp in groupby:
                    
#                    try:
                    try:
                        # not sure why this is a list sometimes
                        if isinstance(grp, list):
                            grp = grp[0]

                        group = str(grp)
                        #print group, tselect
                        if plot_errors:
                            df_group = x[attrs].to_array().to_pandas().T.groupby(group)
                            xaxis = df_group[group].mean().values
                            xlab = group
                            unit = x[group].attrs.get('unit')
                            if unit:
                                xlab = '{:} [{:}]'.format(xlab, unit)
                            
                            for attr in attrs0:
                                plt_type = '{:} vs {:}'.format(attr, group)
                                plt.figure()
                                yaxis = df_group[attr].mean().values
                                yerr = df_group[attr].std().values
                                plt.errorbar(xaxis,yaxis,yerr=yerr)
                                plt.xlabel(xlab)
                                ylab = attr
                                unit = x[attr].attrs.get('unit')
                                if unit:
                                    ylab = '{:} [{:}]'.format(ylab, unit)
                                
                                plt.ylabel(ylab)
                                #plt.tight_layout()
                                howto = ["xaxis = x['{:}'].values".format(group)]
                                howto.append("df_group = xselect.to_array().to_pandas().T.groupby('{:}')".format(group))
                                howto.append("xaxis = df_group['{:}'].mean()".format(group)) 
                                howto.append("yaxis = df_group['{:}'].mean()".format(attr)) 
                                howto.append("yerr = df_group['{:}'].std()".format(attr)) 
                                howto.append("figure()")
                                howto.append("plt.errorbar(xaxis,yaxis,yerr=yerr)")
                                howto.append("plt.xlabel('{:}')".format(xlab))
                                howto.append("plt.ylabel('{:}')".format(ylab))
                                self.add_plot(catagory, plt_type, howto=howto)

                        else: 
                            # Switch to something like 
                            # g = sns.PairGrid(df_group, y_vars=attrs0, x_vars=group, size=5, aspect=.5)
                            # g.map(sns.pointplot, color=sns.xkcd_rgb["plum"])
                            plt_type = 'vs {:}'.format(group)
                            x_group = xselect.groupby(group).mean(dim=tselect)
                            df_group = x_group.to_array().to_pandas().T
                            ndfg = df_group.columns.size
                            ncolumns = int(min([ndfg,3]))
                            nrows = int(np.ceil(ndfg/float(ncolumns)))
                            df_group.plot(subplots=True, sharex=True, layout=(nrows,ncolumns), figsize=figsize)
                            #plt.tight_layout()
                            howto = ["x_group = xselect.groupby('{:}').mean(dim={:})".format(group,tselect)]
                            howto.append("df_group = x_group.to_array().to_pandas().T")
                            howto.append("df_group.plot(subplots=True, sharex=True, layout={:}, figsize={:})".format(layout, figsize))
                            self.add_plot(catagory, plt_type, howto=howto)
                        
                        #grp_values = set(xselect[grp].values)
                        #if len(grp_values) < 9:
                        #    for attr in attrs0:
                        #        plt.cla()

                        #        xselect[attr] plot(robust=True)
                        #        howto = []
                        #        howto.append("xdat.squeeze('steps').sel(codes={:}).plot(robust=True)".format(code))
                        #        howto.append("plt.show()")
                        #        name = attr.replace(alias+'_','')
                        #        plt_type = 'summary_{:}_ec{:}'.format(name, code) 
                        #        self.add_plot(alias, plt_type, howto=howto)

                    except:
                    #except:
                        traceback.print_exc()
                        print 'Failed adding:', group
                        plt.cla()
                        plt.close()

            # make hist plots:
            if make_histplot:
                #perlow = '{:}%'.format(int(min(percentiles)*100))
                #perhigh = '{:}%'.format(int(max(percentiles)*100))
                #dfcut = df[(df > df_tbl[perlow]-2*df_tbl['std']).all(axis=1) & (df < df_tbl[perhigh]+2*df_tbl['std']).all(axis=1)]
                #dfcut = df[(df > df_tbl['5%']-2*df_tbl['std']).all(axis=1) & (df < df_tbl['95%']+2*df_tbl['std']).all(axis=1)]
                dfcut = df
                try:
                    #howto = ["dfcut = df[(df > df_tbl['5%']-2*df_tbl['std']).all(axis=1) & (df < df_tbl['95%']+2*df_tbl['std']).all(axis=1)]"] 
                    howto = []
                    plt_type = 'hist'
                    #dfcut.hist(alpha=0.2, layout=layout, figsize=figsize)
                    df.hist(bins=bins, alpha=0.2, layout=layout, figsize=figsize)
                    #plt.tight_layout()
                    howto.append("df.hist(alpha=0.2, layout={:}, figsize={:})".format(layout, figsize))
                    self.add_plot(catagory, plt_type, howto=howto)
                except:
                    plt.cla()
                    plt.close()
                    print 'make_histplot failed'
                    print howto
                    print dfcut.keys()

            if make_correlation:
                try:
                    from xarray_utils import heatmap
                    heatmap(df, confidence=confidence)
                    plt_type = 'correlation'
                    howto = ['from PyDataSource.xarray_utils import heatmap']
                    howto.append('heatmap(df, confidence={:})'.format(confidence))
                    self.add_plot(catagory, plt_type, howto=howto, tight=False)
                except:
                    print 'Could not make correlation'


            # make scatter matrix using pandas -- cut outliers
            if make_scatter:
                if not robust_attrs:
                    robust_attrs = df.keys()
                    
                robust_attrs = [a for a in [attr.replace(alias+'_','') for attr in robust_attrs] if a in df.keys()]
                dfr = df[robust_attrs]
                df_tblr = df_tbl.T[robust_attrs].T
                dfcut = df[(dfr > df_tblr['5%']-2*df_tblr['std']).all(axis=1) & (dfr < df_tblr['95%']+2*df_tblr['std']).all(axis=1)]
                if dfcut.shape[0] == 0:
                    dfcut = df

                if isinstance(make_scatter, list):
                    #scat_attrs = [attr.replace(alias+'_','') for attr in make_scatter if attr in dfcut.keys()]
                    scat_attrs = [a for a in [attr.replace(alias+'_','') for attr in make_scatter] if a in dfcut.keys()]

                else:
                    scat_attrs = default_scatter_attrs.get(alias, attr_names.values())
                
                for attr in x.attrs.get('correlation_variables', []):
                    if attr in dfcut.keys() and attr not in scat_attrs:
                        scat_attrs.append(attr)

                if len(scat_attrs) > 8:
                    print 'Too many parameters to make scatter plots:', scat_attrs
                elif len(scat_attrs)-len(groupby) > 1:
                    try:
                        dfscat = dfcut[scat_attrs]
                    except:
                        plt.cla()
                        plt.close()
                        print 'make_scatter failed'
                        return dfcut

                    try:
                        howto = ["scat_attrs = {:}".format(scat_attrs)]
                        howto.append("robust_attrs = [a for a in [attr.replace(alias+'_','') for attr in df.keys()]]")
                        howto.append("dfr = df[robust_attrs]")
                        howto.append("df_tblr = df_tbl.T[robust_attrs].T")
                        howto.append("dfcut = df[(dfr > df_tblr['5%']-2*df_tblr['std']).all(axis=1) & (dfr < df_tblr['95%']+2*df_tblr['std']).all(axis=1)]")
                        howto.append("dfscat = dfcut[scat_attrs]")
                        if groupby and self.nsteps < 20:
                            if scat_name:
                                plt_type = scat_name
                            else:
                                plt_type = 'correlation with {:}'.format(group)
                            
                            pltattrs = scat_attrs
                            for group in groupby:
                                if group in pltattrs:
                                    pltattrs.remove(group)
                                
                            sns.set()
                            plt.rcParams['axes.labelsize'] = 10 
                            g = sns.pairplot(dfcut, hue=group,
                                    x_vars=pltattrs,y_vars=pltattrs,
                                    #palette="Set1", 
                                    size=2.5) 
                                    #palette="Set1", diag_kind="kde", size=2.5) 
                                    #x_vars=pltattrs,y_vars=pltattrs)
                            #g = sns.PairGrid(dfscat, hue=group, 
                            #        palette="Set1", 
                            #        x_vars=pltattrs,y_vars=pltattrs)
                            #g = g.map_offdiag(plt.scatter)
                            #g = g.add_legend()
                            #plt.tight_layout()
                            self.add_plot(catagory, plt_type, howto=howto)
    #df = x[attrs].to_dataframe()
    #sns.pairplot(df, hue="LaserOn", vars=attrs0)
                        else:

                            plt_type = 'scatter_matrix'
                            #Axes = pd.tools.plotting.scatter_matrix(dfscat, alpha=0.2, figsize=(20, 20), diagonal='kde')
                            #plt.tight_layout()
                            #howto.append("Axes = pd.tools.plotting.scatter_matrix(dfscat, alpha=0.2, figsize=(20, 20), diagonal='kde')")
                            plt.rcParams['axes.labelsize'] = 10 
                            g = sns.PairGrid(dfscat, diag_sharey=False)
                            g.map_lower(sns.kdeplot, cmap="Blues_d")
                            g.map_upper(plt.scatter)
                            g.map_diag(sns.kdeplot, lw=3)
                            howto.append('g = sns.PairGrid(dfscat, diag_sharey=False)')
                            howto.append('g.map_lower(sns.kdeplot, cmap="Blues_d")')
                            howto.append('g.map_upper(plt.scatter)')
                            howto.append('g.map_diag(sns.kdeplot, lw=3)')
                            
                            self.add_plot(catagory, plt_type, howto=howto)
        
        #                    #y ticklabels
        #                    [plt.setp(item.yaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
        #                    howto.append("[plt.setp(item.yaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]")
        #                    #x ticklabels
        #                    [plt.setp(item.xaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
        #                    howto.append("[plt.setp(item.xaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]")
        #                    #y labels
        #                    [plt.setp(item.yaxis.get_label(), 'size', 20) for item in Axes.ravel()]
        #                    howto.append("[plt.setp(item.yaxis.get_label(), 'size', 20) for item in Axes.ravel()]")
        #                    #x labels
        #                    [plt.setp(item.xaxis.get_label(), 'size', 20) for item in Axes.ravel()]
        #                    howto.append("[plt.setp(item.xaxis.get_label(), 'size', 20) for item in Axes.ravel()]")

                    except:
                        plt.cla()
                        plt.close()
                        print dfscat
                        traceback.print_exc()
                        return dfscat

        plt.close('all')
    
    def add_xy_ploterr(self, attr, xaxis=None, howto=None, catagory=None, text=None, table=None, **kwargs):
        from plotting import xy_ploterr
        x = self._xdat
        if not howto:
            howto = []
        if not catagory:
            catagory = x[attr].attrs.get('alias', 'Summary')
        
        howto.append("Custom plot see PyDataSource.plotting.xy_ploterr method")
        p = xy_ploterr(x, attr, xaxis=xaxis, **kwargs)
        if not xaxis:
            xaxis=x.scan_variables[0]
        plt_type = '{:} vs {:}'.format(attr, xaxis)   
        if kwargs.get('logy'):
            if kwargs.get('logx'):
                plt_type+=' log-log'
            else:
                plt_type+=' log scale'
        elif kwargs.get('logx'):
            plt_type+=' lin-log'

        print catagory, plt_type
        self.add_plot(catagory, plt_type, howto=howto)
        if table is not None:
            self.results[catagory]['table'].update({attr:{'DataFrame': table, 
                                                        'name': 'df_tbl',
                                                        'howto': [], 
                                                        'doc': text}})


    def add_scatter(self, df, catagory='scatter', doc=None, group=None, attrs=None, howto=[]):
        """Add scatter plot.
        """
        import seaborn as sns
        if not attrs:
            attrs = df.keys()
       
        sns.set()
        plt.rcParams['axes.labelsize'] = 10 
        pltattrs = [attr for attr in attrs if attr not in [group]]

        howto.append("sns.set()")
        howto.append("pltattrs = {:}".format(pltattrs))
        
        if group:
            plt_type = 'correlation with {:}'.format(group)
            g = sns.pairplot(df, hue=group,
                    x_vars=pltattrs,y_vars=pltattrs,
                    size=2.5) 
            howto.append("g = sns.pairplot(df, hue='{:}', x_vars=pltattrs,y_vars=pltattrs,size=2.5)".format(group))
        else:
            plt_type = 'correlation'
            g = sns.PairGrid(df, x_vars=pltattrs,y_vars=pltattrs, size=2.5) 
            g = g.map_upper(plt.scatter)
            g = g.map_lower(sns.kdeplot, cmap="Blues_d")
            g = g.map_diag(sns.kdeplot, lw=3, legend=False)

            howto.append("g = sns.pairplot(df, x_vars=pltattrs,y_vars=pltattrs,size=2.5)")

        #plt.tight_layout()
        self.add_plot(catagory, plt_type, howto=howto)
        plt.close('all')

    def add_scatter_groups(self, attr_groups=None, group=None, howto=None, attrs=None):
        """Add Scatter Groups
        """
        x = self._xdat
        if not attrs:
            attrs=['EBeam_ebeamPhotonEnergy', 'FEEGasDetEnergy_f_21_ENRC']
            for attr in x.attrs.get('correlation_variables', []):
                attrs.append(attr)

        if not group:
            group = 'XrayOn'
        
        if group not in x:
            print group, 'is not a valid group -- ignoring group keyword'
            group = None

        if not howto:
            howto = []
        
        if not attr_groups:
            attr_groups = {'Undulator X': ['EBeam_ebeamUndAngX', 'EBeam_ebeamUndPosX'],
                           'Undulator Y': ['EBeam_ebeamUndAngY', 'EBeam_ebeamUndPosY'],
                           'LTU X': ['EBeam_ebeamLTUAngX', 'EBeam_ebeamLTUPosX'],
                           'LTU Y': ['EBeam_ebeamLTUAngY', 'EBeam_ebeamLTUPosY'],
                           'LTU 250 and 450': ['EBeam_ebeamLTU250', 'EBeam_ebeamLTU450'],
                           'Energy BC': ['EBeam_ebeamEnergyBC1', 'EBeam_ebeamEnergyBC2'],
                           'Charge': ['EBeam_ebeamCharge', 'EBeam_ebeamDumpCharge'],
                           'Peak Current': ['EBeam_ebeamPkCurrBC1', 'EBeam_ebeamPkCurrBC2'],
                           'XTCAV': ['EBeam_ebeamXTCAVAmpl', 'EBeam_ebeamXTCAVPhase'],
                           'Phase Cavity Charge': ['PhaseCavity_charge1', 'PhaseCavity_charge2'],
                           'Phase Cavity Fit': ['PhaseCavity_fitTime1', 'PhaseCavity_fitTime2'],
                           'Gasdet': ['FEEGasDetEnergy_f_63_ENRC', 'FEEGasDetEnergy_f_11_ENRC'],
                           }

        all_attrs = [attr for attr,item in x.data_vars.items() if item.dims == ('time',)]
        
        if 'Damage_cut' not in x:
            self.make_damage_cut()

        df = self._xdat[all_attrs].where(self._xdat.Damage_cut == 1).to_array().to_pandas().T
        
        if group not in df:
            df[group] = x[group].values

        df = df.dropna()

        for catagory, grp_attrs in attr_groups.items():
            gattrs = [attr for attr in attrs]
            for attr in grp_attrs:
                gattrs.append(attr)
            ghowto = howto
            print catagory, gattrs, group, ghowto
            self.add_scatter(df, catagory=catagory, attrs=gattrs, howto=ghowto, group=group)

    def make_damage_cut(self, charge_min=1., 
            fitTime1_min=0.3, fitTime1_max=1.1, 
            fitTime2_min=-0.4, fitTime2_max=0.2,
            gasdet_min=0.5):
        """Make damage cut based on EBeam and Gasdet data.
        """
        x = self._xdat
        phasecut = (x.PhaseCavity_charge1.values > charge_min) & \
                (x.PhaseCavity_charge2.values > charge_min) & \
                (x.PhaseCavity_fitTime1.values <= fitTime1_max) & \
                (x.PhaseCavity_fitTime1.values >= fitTime1_min) & \
                (x.PhaseCavity_fitTime2.values <= fitTime2_max) & \
                (x.PhaseCavity_fitTime2.values >= fitTime2_min) 
        gasdetcut =  x.FEEGasDetEnergy_f_11_ENRC.values > 0.5
        damagecut = phasecut & gasdetcut & (x.EBeam_damageMask.values == 0)  
        self._xdat.coords['Damage_cut'] = (['time'], damagecut)

#    def add_scatter(self, df, catagory='scatter', group=None, attrs=None, howto=[]):
#        if not attrs:
#            attrs = df.keys()
#       
#        pltattrs = [attr for attr in attrs if attr not in [group]]
#        if group:
#            plt_type = 'correlation with {:}'.format(group)
#        else:
#            plt_type = 'correlation'
#
#        sns.set()
#        plt.rcParams['axes.labelsize'] = 10 
#        g = sns.pairplot(df, hue=group,
#                x_vars=pltattrs,y_vars=pltattrs,
#                size=2.5) 
#
#        howto.append("sns.set()")
#        howto.append("pltattrs = {:}".format(pltattrs))
#        howto.append("g = sns.pairplot(df, hue={:}, x_vars=pltattrs,y_vars=pltattrs,size=2.5)")
#
#        self.add_plot(catagory, plt_type, howto=howto)
#
#                #palette="Set1", diag_kind="kde", size=2.5) 
#                #x_vars=pltattrs,y_vars=pltattrs)
#        #g = sns.PairGrid(dfscat, hue=group, 
#        #        palette="Set1", 
#        #        x_vars=pltattrs,y_vars=pltattrs)
#        #g = g.map_offdiag(plt.scatter)
#        #g = g.add_legend()

    def add_setup(self, catagory, howto):
        """Add text to beginning of HowTo.
        """
        self._add_catagory(catagory)
        self.results[catagory]['text'].update({'setup':{'howto': howto}})

    def _add_catagory(self, catagory):
        if catagory not in self.results:
            self.results[catagory] = {'figure': {}, 'table': {}, 'text': {}, 'textblock': {}}
    
    def add_plot(self, catagory, plt_type, howto=[], tight=True, show=False):
        """Add a plot to RunSummary.
        """
        if tight:
            try:
                plt.tight_layout()
            except:
                # seems like tight_layout does not work with add_detector 
                print 'cannot make tight', catagory, plt_type
        self._add_catagory(catagory)
        plt_file = '{:}_{:}.png'.format(catagory, plt_type).replace(' ','_') 
        self.results[catagory]['figure'].update({plt_type: {'path': self.output_dir, 
                                                         'png': plt_file,
                                                         'howto': howto}})
        plt.savefig(os.path.join(self.output_dir, plt_file))
        if show:
            plt.show()
        else:
            plt.close()

    def add_summary(self, 
                        variables=None, 
                        max_columns=10,
                        bins=50, 
                        plt_style='.',
                        show=False,
                        groupby=None,
                        codes=None,
                        wrap=False,
                        labelsize=None, 
                        figsize=None,
                        plot_axis_sums=None,
                        pltsize=8, 
                        nlines=16,
                        layout=None):
        """Add summary for variables.
        """
        import seaborn as sns
        if not variables:
            if 'scan_variables' in self._xdat.attrs:
                variables = self._xdat.scan_variables
            else:
                variables = []

        plt.close()
        plt.cla()
        variables = [v for v in variables if v in self._xdat.variables.keys()]

        if groupby is None:
            groupby = [group for group in self._xdat.attrs.get('scan_variables') \
                            if group in self._xdat and len(set(self._xdat[group].data)) > 1]

        if groupby and not variables:
            variables = [attr for attr, item in self._xdat.variables.items() \
                            if len(item.shape) == 2 and item.shape[1] > 8]

        for attr in variables:
            if 'alias' in self._xdat[attr].attrs:
                alias = self._xdat[attr].alias
            else:
                alias = attr.split('_')[0]

            catagory = alias
            self._add_catagory(catagory)
            # Make groupby plots
            if groupby:
                if not isinstance(groupby, list) and not hasattr(groupby, 'shape'):
                    groupby = [groupby]
                for group in groupby:
                    print 'Make groupby plots for', group
                    xdat = self._xdat[attr]
                    x_group = xdat.groupby(group).mean(dim='time')
                    howto = ["plt.rcParams['axes.labelsize'] = {:}".format(24)]
                    howto.append("x_group = x['{:}'].T.groupby('{:}').mean(dim='time')".format(attr, group))
                    self.results[alias]['text'].update({'setup':{'howto': howto}})
 
                    if len(x_group.shape) == 2:
                        if x_group.shape[1] > nlines:
                            x_group.plot()
                            plt_type = '{:} vs {:}'.format(attr.lstrip(alias), group) 
                            howto = []
                            howto.append('x_group.plot()')
                            self.add_plot(alias, plt_type, howto=howto)
                        
                        else:
                            pass
                            # make scatter with groupby like in add_detector 
                            #df = xdat.to_pandas().dropna()
                            #pltattrs = scat_attrs
                            #for group in groupby:
                            #    if group in pltattrs:
                            #        pltattrs.remove(group)
                            #    
                            #sns.set()
                            #plt.rcParams['axes.labelsize'] = 10 
                            #g = sns.pairplot(dfcut, hue=group,
                            #        x_vars=pltattrs,y_vars=pltattrs,
                            #        #palette="Set1", 
                            #        size=2.5) 
 
                            # Not very interesting generally to have mean plots like this:
                            #plt.cla()
                            #x_group.to_pandas().T.plot()
                            #howto = []
                            #howto.append("x_group.to_pandas().T.plot()")
                            #howto.append("plt.show()")
                            #name = attr.replace(alias+'_','')
                            #plt_dim = x_group.dims[1]
                            #plt_type = 'summary_{:}_{:}'.format(name, str(plt_dim)) 
                            #self.add_plot(alias, plt_type, howto=howto)
                            # Instead Add pairplot next ..
                            #df = xdat.to_pandas()
                            #df[group] = xdat[group].to_pandas()
                            #g = sns.PairGrid(df, hue=group, diag_sharey=False)
                            #g.map_lower(sns.kdeplot, cmap="Blues_d")
                            #g.map_upper(plt.scatter)
                            #g.map_diag(sns.kdeplot, lw=3)

                    elif len(x_group.shape) == 3:
                        if x_group.shape[1] <= nlines:
                            plt_dim = x_group.dims[1]
                            for ich, ch in enumerate(x_group[plt_dim]):
                                plt.cla()
                                x_group[:,ich].to_pandas().T.plot()
                                howto = []
                                howto.append("x_group[:,{:}].to_pandas().T.plot()".format(ich))
                                howto.append("plt.show()")
                                name = attr.replace(alias,'')
                                plt_type = 'summary_{:}_{:}_{:}'.format(name, str(plt_dim), ich) 
                                self.add_plot(alias, plt_type, howto=howto)
            
            else:
                nax = len(self._xdat[attr].shape)
                print nax, attr
                if nax == 2:
    #                self._xdat[attr].mean(axis=0).plot()
    #                plt_type = '{:} summary'.format(attr.lstrip(alias+'_')) 
    #                howto = ["x['{:}'].mean(axis=0).plot()".format(attr)]
    #                self.add_plot(alias, plt_type, howto=howto)
                   
                    if self._xdat[attr].shape[1] <= nlines:
                        try:
                            df = self._xdat[attr].to_pandas().dropna()
                            df.plot(style=plt_style)
                            plt_type = '{:} with {:}'.format(attr.lstrip(alias), 'time') 
                            howto = ["df = x['{:}'].to_pandas().dropna()".format(attr)]
                            howto.append("df.plot(style='{:}')".format(plt_style))
                            self.add_plot(alias, plt_type, howto=howto)
                            
                            df.plot.hist(bins=bins)
                            plt_type = '{:} {:}'.format(attr.lstrip(alias), 'histogram') 
                            howto = ["df = x['{:}'].to_pandas().dropna()".format(attr)]
                            howto.append("df.plot.hist(bins={:})".format(bins))
                            self.add_plot(alias, plt_type, howto=howto)

                            g = sns.PairGrid(df, diag_sharey=False)
                            g.map_lower(sns.kdeplot, cmap="Blues_d")
                            g.map_upper(plt.scatter)
                            g.map_diag(sns.kdeplot, lw=3)
                            howto = ["df = x['{:}'].to_pandas().dropna()".format(attr)]
                            howto.append('g = sns.PairGrid(df, diag_sharey=False)')
                            howto.append('g.map_lower(sns.kdeplot, cmap="Blues_d")')
                            howto.append('g.map_upper(plt.scatter)')
                            howto.append('g.map_diag(sns.kdeplot, lw=3)')
                            plt_type = '{:} {:}'.format(attr.lstrip(alias), 'scatter matrix') 
                            self.add_plot(alias, plt_type, howto=howto)
                        except:
                            traceback.print_exc()
                            print 'scatter matrix plot failed for', alias

                    else:
                        self._xdat[attr].plot()
                        plt_type = '{:} with {:}'.format(attr.lstrip(alias), 'time') 
                        howto = ["x['{:}'].plot()".format(attr)]
                        self.add_plot(alias, plt_type, howto=howto)


                if nax in [2,3]:
                    if self._xdat[attr].shape[1] > nlines:
                        if plot_axis_sums:
                            for iaxis in range(nax):
                                pltaxis = self._xdat[attr].dims[iaxis]
                                self._xdat[attr].sum(axis=iaxis).plot()
                                plt_type = '{:} sum over {:}'.format(attr.lstrip(alias), pltaxis) 
                                howto = ["x['{:}'].mean(axis=0).plot()".format(attr)]
                                self.add_plot(alias, plt_type, howto=howto)
                        else:
                            iaxis = 0
                            pltaxis = self._xdat[attr].dims[iaxis]
                            self._xdat[attr].sum(axis=iaxis).plot()
                            plt_type = '{:} sum over {:}'.format(attr.lstrip(alias), pltaxis) 
                            howto = ["x['{:}'].mean(axis=0).plot()".format(attr)]
                            self.add_plot(alias, plt_type, howto=howto)
                    
                    elif nax == 3 and codes:
                        attr_codes = [code for code in codes if code in self._xdat.codes.values]
                        for code in attr_codes:
                            plt.cla()
                            xdat = self._xdat[attr]
                            df = xdat.where(xdat['ec{:}'.format(code)] == 1).to_pandas().dropna().mean(axis=0).T
                            df.plot()
                            pltaxis = df.index.name
                            plt_type = '{:} summary for ec{:}'.format(attr.lstrip(alias), code) 
                            howto = ["xdat = x['{:}']".format(attr)]
                            howto.append("df = xdat.where(xdat[{:}] == 1).to_pandas().dropna().mean(axis=0).T".format('ec{:}'.format(code)))
                            self.add_plot(alias, plt_type, howto=howto)

    #                    plt_dim = xdat.dims[1]
    #                    for ich, ch in enumerate(xdat[plt_dim]):
    #                        plt.cla()
    #                        xdat[:,ich].to_pandas().T.plot()
    #                        howto = []
    #                        howto.append("xdat.squeeze('steps')[:,{:}].to_pandas().T.plot()".format(ich))
    #                        howto.append("plt.show()")
    #                        name = attr.replace(alias+'_','')
    #                        plt_type = 'summary_{:}_{:}_{:}'.format(name, str(plt_dim), ich) 
    #                        self.add_plot(alias, plt_type, howto=howto)
    #

#### Moved to add_stats method
#            if len(self._xdat[attr].shape) == 4 and 'codes' in self._xdat[attr]:
#                xdat = self._xdat[attr]
#                attr_codes = None
#                if codes:
#                    attr_codes = codes
#                elif hasattr(self, 'eventStats'):
#                    try:
#                        tbl = self.eventStats.T.get(alias+'_events')
#                        if tbl is not None:
#                            attr_codes = [c for c,n in tbl.to_dict().items() if n > 0]
#                    except:
#                        print 'WARNING: Need to fix auto getting attr_codes'
#                
#                if not attr_codes:
#                    if 'codes' in self._xdat:
#                        attr_codes = self._xdat.codes.values
#                    else:
#                        attr_codes = [140]
#                else:
#                    attr_codes = [code for code in attr_codes if code in self._xdat.codes.values]
#
#                print 'det_summary', attr, attr_codes
#
#                if self.nsteps < max_columns**2:
#                    ncolumns = int(sqrt(self.nsteps))
#                else:
#                    ncolumns = int(min([self.nsteps,max_columns]))
#                    
#                nrows = int(np.ceil(self.nsteps/float(ncolumns)))
#                
#                if not labelsize:
#                    if ncolumns > 9:
#                        labelsize = 12
#                    elif ncolumns > 5:
#                        labelsize = 18
#                    else:
#                        labelsize = 24
#
#                if not layout:
#                    layout = (nrows,ncolumns)
#                if not figsize:
#                    figsize = (ncolumns*4,nrows*4)
#                
#                plt.rcParams['axes.labelsize'] = 24 
#                howto = ["plt.rcParams['axes.labelsize'] = {:}".format(24)]
#                howto.append("attr = '{:}'".format(attr))
#                howto.append("xdat = x[attr]")
#                self.results[alias]['text'].update({'setup':{'howto': howto}})
#
#                # Need to fix errors calculation with std and plot with error bars
#                #fig, ax = plt.subplots()
#                #errors = sqrt(self._xdat[attr].sum(axis=(2,3)).to_pandas()**2)
#                #df.plot(yerr=errors, ax=ax)
#
#                if self.ncodes > 4:
#                    col_wrap = int(min([self.ncodes,3]))
#                    aspect = 0.9
#                elif self.ncodes == 4:
#                    col_wrap = 2
#                    aspect = 0.9
#                else:
#                    col_wrap = 1
#                    aspect = 0.9
#                
#                if self.nsteps == 1:
#                    if wrap:
#                        howto = []
#                        df = xdat.squeeze('steps')
#                        howto.append("df = xdat.squeeze('steps')")
#                        if self.ncodes == 1:
#                            df.squeeze('codes').plot(robust=True)
#                            howto.append("df.squeeze('codes').plot()")
#                            howto.append("plt.show()")
#                        else:
#                            df.plot(col='codes', col_wrap=col_wrap, size=20, aspect=aspect, robust=True)
#                            howto.append("df.plot(col='codes', col_wrap={:})".format(col_wrap))
#                            howto.append("plt.show()")
#
#                        name = attr.replace(alias+'_','')
#                        plt_type = 'summary_{:}'.format(name) 
#                        self.add_plot(alias, plt_type, howto=howto)
#                        
#                    else:
#                        if xdat.shape[2] > 8: 
#                            for code in attr_codes:
#                                plt.cla()
#                                plt.figure(figsize=(16,14))
#                                xdat.squeeze('steps').sel(codes=code).plot(robust=True,size=pltsize,aspect='equal')
#                                howto = []
#                                howto.append('plt.figure(figsize=(16,14))')
#                                howto.append("xdat.squeeze('steps').sel(codes={:}).plot(robust=True,size={:}aspect='equal')".format(code, pltsize))
#                                howto.append("plt.show()")
#                                name = attr.replace(alias+'_','')
#                                plt_type = 'summary_{:}_ec{:}'.format(name, code) 
#                                self.add_plot(alias, plt_type, howto=howto)
#                        else:
#                            for code in attr_codes:
#                                plt.cla()
#                                xdat.squeeze('steps').sel(codes=code).to_pandas().T.plot()
#                                howto = []
#                                howto.append("xdat.squeeze('steps').sel(codes={:}).to_pandas().T.plot()".format(code))
#                                howto.append("plt.show()")
#                                name = attr.replace(alias+'_','')
#                                plt_type = 'summary_{:}_ec{:}'.format(name, code) 
#                                self.add_plot(alias, plt_type, howto=howto)
#                           
#                            plt_dim = xdat.dims[2]
#                            for ich, ch in enumerate(xdat[plt_dim]):
#                                plt.cla()
#                                xdat.squeeze('steps')[:,ich].to_pandas().T.plot()
#                                howto = []
#                                howto.append("xdat.squeeze('steps')[:,{:}].to_pandas().T.plot()".format(ich))
#                                howto.append("plt.show()")
#                                name = attr.replace(alias+'_','')
#                                plt_type = 'summary_{:}_{:}_{:}'.format(name, str(plt_dim), ich) 
#                                self.add_plot(alias, plt_type, howto=howto)
#
#
#                else:
#                    # plot mean of image vs steps for each event code
#                    df = xdat.sum(axis=(2,3)).to_pandas()
#                    df.plot(subplots=True)
#                   
#                    howto = []
#                    howto.append("df = xdat.sum(axis=(2,3)).to_pandas()")
#                    howto.append("df.plot(subplots=True)")
#                    howto.append("plt.show()")
#                    name = attr.replace(alias+'_','')
#                    plt_type = 'scan_{:}_summary'.format(name) 
#                    self.add_plot(alias, plt_type, howto=howto)
#     
#                    for iaxis in [2,3]:
#                        df = xdat.mean(axis=iaxis)
#                        
#                        howto = ["plt.rcParams['axes.labelsize'] = 24"]
#                        howto.append("df = xdat.mean(axis={:}).to_pandas()".format(iaxis))
#                        if self.ncodes == 1:
#                            df.squeeze('codes').plot(robust=True)
#                            howto.append("df.squeeze('codes').plot(robust=True)")
#                            howto.append("plt.show()")
#                        else:
#                            df.plot(col='codes', col_wrap=col_wrap, size=20, aspect=aspect, robust=True)
#                            sformat = "df.plot(col='codes', col_wrap={:}, size=20, aspect={:}, robust=True)"
#                            howto.append(sformat.format(col_wrap, aspect))
#                            howto.append("plt.show()")
#
#                        name = df.dims[2].replace(alias+'_','')
#                        plt_type = 'scan_vs_{:}'.format(name) 
#                        self.add_plot(alias, plt_type, howto=howto)
#
#                    plt.rcParams['axes.labelsize'] = labelsize
#                    for code in self._xdat[attr].codes.values:
#                        df = xdat.sel(codes=code)
#                        df.plot(col='steps', col_wrap=ncolumns, robust=True)
#                        
#                        howto = ["plt.rcParams['axes.labelsize'] = {:}".format(labelsize)]
#                        howto.append("df = xdat.sel(codes={:})".format(code))
#                        howto.append("df.plot(col='steps', col_wrap={:})".format(ncolumns))
#                        name = attr.replace(alias+'_','')
#                        plt_type = 'scan_{:}_2D_ec{:}'.format(name, str(code)) 
#                        self.add_plot(alias, plt_type, howto=howto)
        plt.cla()

    def add_stats(self, attr, codes=None, stats=['mean','std','max'],
                        catagory=None, 
                        make_images=False,
                        max_columns=10,
                        show=False,
                        wrap=False,
                        labelsize=None, 
                        figsize=None,
                        pltsize=8, 
                        nlines=16,
                        layout=None):
        """
        Add image stats
        """
        plt.close('all')
        xdat = self._xdat[attr]
        if len(xdat.shape) not in [5,6]:
            print "add_stats only valid for objects with dimensions ('stat','step','codes', xdim, ydim)"
            return

        if 'alias' in xdat.attrs:
            alias = xdat.alias
        else:
            alias = attr.split('_')[0]

        if not catagory:
            catagory = alias+' Stats'
        self._add_catagory(catagory)

        attr_codes = None
        if codes:
            attr_codes = codes
        elif hasattr(self, 'eventStats'):
            try:
                tbl = self.eventStats.T.get(alias+'_events')
                if tbl is not None:
                    attr_codes = [c for c,n in tbl.to_dict().items() if n > 0]
            except:
                print 'WARNING: Need to fix auto getting attr_codes'
        
        if not attr_codes:
            if 'codes' in xdat:
                attr_codes = xdat.codes.values
            else:
                attr_codes = [140]
        else:
            attr_codes = [code for code in attr_codes if code in xdat.codes.values]
        
        plt.rcParams['axes.labelsize'] = 24 
        howto = ["plt.rcParams['axes.labelsize'] = {:}".format(24)]
        howto.append("attr = '{:}'".format(attr))
        howto.append("xdat = x[attr]")
        self.results[catagory]['text'].update({'setup':{'howto': howto}})

        for stat in stats:
            if self.nsteps == 1:
                xstat = xdat.sel(stat=stat).squeeze('steps')
                howto0 = ["xstat = xdat.sel(stat='{:}').squeeze('steps')".format(stat)]
                name = attr.replace(alias,'')+' '+stat
            else:
                if stat == 'mean':
                    xstat = xdat.sel(stat=stat).mean(dim='steps')
                    howto0 = ["xstat = xdat.sel(stat='{:}').sum(dim='steps')".format(stat)]
                    name = attr.replace(alias,'')+' sum_over_steps'
                else:
                    continue

            if len(xstat.shape) == 4:
                # make images
                from h5write import make_image
                try:
                    xcodes = xstat.codes
                    xstat = xr.concat([make_image(xstat.sel(codes=code)) for code in xstat.codes], dim='codes')
                    if 'codes' not in xstat.coords:
                        xstat.coords['codes'] = xcodes
                except:
                    print 'Cannot make image', attr, stat
                    continue

           #print 'det_summary', attr, attr_codes, stat

           # Need to fix errors calculation with std and plot with error bars
            #fig, ax = plt.subplots()
            #errors = sqrt(self._xstat[attr].sum(axis=(2,3)).to_pandas()**2)
            #df.plot(yerr=errors, ax=ax)

            if self.ncodes > 4:
                col_wrap = int(min([self.ncodes,3]))
                aspect = 0.9
            elif self.ncodes == 4:
                col_wrap = 2
                aspect = 0.9
            else:
                col_wrap = 1
                aspect = 0.9
            
            if wrap:
                howto = howto0
                if self.ncodes == 1:
                    xstat.squeeze('codes').plot(robust=True)
                    howto.append("xstat.squeeze('codes').plot()")
                    howto.append("plt.show()")
                else:
                    xstat.plot(col='codes', col_wrap=col_wrap, size=pltsize, aspect=aspect, robust=True)
                    howto.append("xstat.plot(col='codes', col_wrap={:})".format(col_wrap))
                    howto.append("plt.show()")

                #plt.tight_layout()
                plt_type = '{:}'.format(name) 
                self.add_plot(catagory, plt_type, howto=howto)
                
            else:
                if xstat.shape[1] > nlines:
                    print 'doing image', stat
                    for code in attr_codes:
                        plt.cla()
                        #plt.figure(figsize=(16,14))
                        plt.close('all')
                        try:
                            #ax = xstat.sel(codes=code).reset_coords()[attr].plot(robust=True, size=pltsize)
                            ax = xstat.sel(codes=code).plot(robust=True, size=pltsize)
                            ax.axes.set_aspect('equal')
                            #plt.tight_layout()
                            howto = howto0
                            howto.append('plt.figure(figsize=(16,14))')
                            howto.append("xstat.sel(codes={:}).plot(robust=True, size={:}, aspect=True)".format(code, pltsize))
                            howto.append("plt.show()")
                            plt_type = '{:}_ec{:}'.format(name, code) 
                            self.add_plot(catagory, plt_type, howto=howto)
                        except:
                            traceback.print_exc()
                            print xstat
                        plt.close()
                elif stat == 'mean':
                    for code in attr_codes:
                        plt.cla()
                        xstat.sel(codes=code).to_pandas().T.plot()
                        howto = howto0
                        howto.append("xstat.sel(codes={:}).to_pandas().T.plot()".format(code))
                        howto.append("plt.show()")
                        plt_type = 'summary_{:}_ec{:}'.format(name, code) 
                        self.add_plot(catagory, plt_type, howto=howto)
                  
                    if len(attr_codes) > 1:
                        plt_dim = xstat.dims[1]
                        for ich, ch in enumerate(xstat[plt_dim]):
                            print ich, ch
                            plt.cla()
                            xstat[:,ich,:].to_pandas().T.plot()
                            howto = []
                            howto.append("xstat[:,{:}].to_pandas().T.plot()".format(ich))
                            howto.append("plt.show()")
                            plt_type = '{:}_{:}_{:}_{:}'.format(name, stat, str(plt_dim), ich) 
                            self.add_plot(catagory, plt_type, howto=howto)

        if self.nsteps > 1 and xdat.shape[3] <= nlines:
            print 'Smart waveforms'

        elif self.nsteps > 1:
            if self.nsteps < max_columns**2:
                ncolumns = int(sqrt(self.nsteps))
            else:
                ncolumns = int(min([self.nsteps,max_columns]))
                
            nrows = int(np.ceil(self.nsteps/float(ncolumns)))
            
            if not labelsize:
                if ncolumns > 9:
                    labelsize = 12
                elif ncolumns > 5:
                    labelsize = 18
                else:
                    labelsize = 24

            if not layout:
                layout = (nrows,ncolumns)
            if not figsize:
                figsize = (ncolumns*4,nrows*4)
            
            plt.cla()
            # plot mean of image vs steps for each event code
            print xdat
            sumaxes = tuple(range(2,len(xdat.shape)-1))
            df = xdat.sel(stat='mean').sum(axis=sumaxes).to_pandas()
            df.plot()
            #plt.tight_layout()
            #df.plot(subplots=True)

            howto = []
            howto.append("df = xdat.sel(stat='mean').sum(axis=(2,3)).to_pandas()")
            #howto.append("df.plot(subplots=True)")
            howto.append("df.plot()")
            howto.append("plt.show()")
            name = attr.replace(alias,'')
            plt_type = 'scan_{:}_summary'.format(name) 
            self.add_plot(catagory, plt_type, howto=howto)
            plt.close()

            # need to handle CsPad calib data as images and radial/azimuth projections
            if len(xdat.shape) <= 5:
                for iaxis in [2,3]:
                    plt.cla()
                    da = xdat.sel(stat='mean').mean(axis=iaxis)
                    #plt.tight_layout()
                    howto = ["plt.rcParams['axes.labelsize'] = 24"]
                    howto.append("da = xdat.sel(stat='mean').mean(axis='{:}')".format(iaxis))
                    if self.ncodes == 1:
                        da.squeeze('codes').plot(robust=True)
                        howto.append("da.squeeze('codes').plot(robust=True)")
                        howto.append("plt.show()")
                    else:
                        da.plot(col='codes', col_wrap=col_wrap, size=pltsize, aspect=aspect, robust=True)
                        sformat = "da.plot(col='codes', col_wrap={:}, size={:}, aspect={:}, robust=True)"
                        howto.append(sformat.format(col_wrap, pltsize, aspect))
                        howto.append("plt.show()")

                    plt_type = '{:} scan_vs_{:}'.format(name, da.dims[2]) 
                    self.add_plot(catagory, plt_type, howto=howto)
                    plt.close()

                    da.sum(dim='steps').squeeze('codes').to_pandas().plot()
                    #plt.tight_layout()
                    howto.append("da.sum(dim='steps').squeeze('codes').plot()")
                    howto.append("plt.show()")

                    plt_type = '{:} sum over steps vs {:}'.format(name, da.dims[2]) 
                    self.add_plot(catagory, plt_type, howto=howto)
                    plt.close()

            if make_images and len(xdat.shape) <= 5:
                plt.rcParams['axes.labelsize'] = labelsize
                for code in xdat.codes.values:
                    plt.cla()
                    da = xdat.sel(stat='mean').sel(codes=code)
                    da.plot(col='steps', col_wrap=ncolumns, robust=True)
                    #plt.tight_layout()
                    
                    howto = ["plt.rcParams['axes.labelsize'] = {:}".format(labelsize)]
                    howto.append("da = xdat.sel(stat='mean').sel(codes={:})".format(code))
                    howto.append("da.plot(col='steps', col_wrap={:})".format(ncolumns))
                    plt_type = '{:} 2D each step for ec{:}'.format(name,  str(code)) 
                    self.add_plot(catagory, plt_type, howto=howto)
                    plt.close()

        elif len(xdat.shape) <= 5:
            xstat = xdat.sel(stat='mean').squeeze('steps')
            if xstat.shape[1] <= nlines: 
                for code in attr_codes:
                    plt.cla()
                    plt.figure(figsize=(16,14))
                    df = xstat.sel(codes=code).reset_coords()[attr].to_pandas().T
                    df.plot()
                    #plt.tight_layout()

                    howto = howto0
                    howto.append('plt.figure(figsize=(16,14))')
                    howto.append("df = xstat.sel(codes={:}).to_pandas()".format(code))
                    howto.append("df.plot()")
                    howto.append("plt.show()")
                    plt_type = '{:}_ec{:}'.format(name, code) 
                    self.add_plot(catagory, plt_type, howto=howto)
                    plt.close()
#                if xstat.shape[1] > 1 and xstat.shape[0] > 1:
#                    for ich in 
#                        plt.cla()
#                        plt.figure(figsize=(16,14))
#                        xstat.sel(codes=code).plot()
#
#                        howto = howto0
#                        howto.append('plt.figure(figsize=(16,14))')
#                        howto.append("xstat.sel(codes={:}).plot(robust=True, size={:}, aspect=True)".format(code, size))
#                        howto.append("plt.show()")
#                        plt_type = '{:}_ec{:}'.format(name, code) 
#                        self.add_plot(catagory, plt_type, howto=howto)
     
            else:
                for iaxis in [1,2]:
                    plt.close()
                    da = xstat.reset_coords()[attr].mean(axis=iaxis)
                    plt_type = '{:} scan_vs_{:}'.format(name, da.dims[1]) 
                    try:
                        howto = ["plt.rcParams['axes.labelsize'] = 24"]
                        howto.append("xstat = xdat.sel(stat='mean').squeeze('steps')")
                        howto.append("da = xstat.mean(axis={:})".format(iaxis))
                        if self.ncodes == 1:
                            da.squeeze('codes').plot(size=pltsize)
                            #plt.tight_layout()
                            howto.append("da.squeeze('codes').plot()")
                            howto.append("plt.show()")
                        else:
                            df = da.to_pandas().T.plot()
                            howto.append('df = da.to_pandas().T.plot()')
                            #df.plot(subplots=True, layout=(self.ncodes,1), figsize=(self.ncodes*6,10))
                            #sformat = "da.plot(col='codes', col_wrap={:}, size={:})"
                            #plt.tight_layout()
                            #howto.append(sformat.format(col_wrap, pltsize, aspect))
                            howto.append("plt.show()")

                        self.add_plot(catagory, plt_type, howto=howto)
                        plt.close()

                    except:
                        traceback.print_exc()
                        print 'cannot add ', catagory, plt_type
    #
#                da.squeeze('codes').to_pandas().plot()
#                howto.append("da.squeeze('codes').plot()")
#                howto.append("plt.show()")
#
#                plt_type = '{:} sum over steps vs {:}'.format(name, da.dims[1]) 
#                self.add_plot(catagory, plt_type, howto=howto)


        plt.cla()
        plt.close()

    def add_event(self, variables=None,
                        max_columns=10,
                        nmax_events=100,
                        show=False,
                        pltsize=10,
                        labelsize=None, 
                        figsize=None,
                        layout=None):
        """Add every event of variables.
        """

        if not variables:
            if 'event_variables' in self._xdat.attrs:
                variables = self._xdat.event_variables
            else:
                variables = []

        variables = [v for v in variables if v in self._xdat.variables.keys()]

        for attr in variables:
            alias = self._xdat[attr].alias
            catagory = alias
            self._add_catagory(catagory)

            if len(self._xdat[attr].shape) == 3 and 'time' in self._xdat[attr]:
                xdat = self._xdat[attr]
                itimes = xdat.groupby('ec{:}'.format(xdat.eventCode)).groups[1]
                df = xdat.isel_points(time=itimes)
                #cut = [i for i,a in enumerate(self._xdat[attr][:,0,0].values) if isfinite(a)]
               
                nevents = len(df.time)
                if nevents == 2:
                    ncolumns = 2 
                elif nevents < max_columns**2:
                    ncolumns = int(sqrt(nevents))
                else:
                    ncolumns = int(min([nevents,max_columns]))
                    
                nrows = int(np.ceil(nevents/float(ncolumns)))
                
                if not labelsize:
                    if ncolumns > 9:
                        labelsize = 12
                    elif ncolumns > 5:
                        labelsize = 18
                    else:
                        labelsize = 24

                if not layout:
                    layout = (nrows,ncolumns)
                if not figsize:
                    figsize = (ncolumns*4,nrows*4)
 
                plt.rcParams['axes.labelsize'] = 24
                
                howto = ["plt.rcParams['axes.labelsize'] = 24"]
                howto.append("attr = ['{:}']".format(attr))
                howto.append("xdat = x[attr]")
                howto.append("itimes = xdat.groupby('ec{:}').groups[1]".format(xdat.eventCode))
                howto.append("df = xdat.isel_points(time=itimes)")
                self.results[alias]['text'].update({'begin':{'howto': howto}})
                if nevents <= nmax_events:
                    #plt.rcParams['axes.labelsize'] = labelsize
                    plt.cla()
                    plt.close('all')
                    ax = df.reset_coords()[attr].plot.pcolormesh(col='points', col_wrap=ncolumns, robust=True, size=pltsize)
                    
                    #howto = ["plt.rcParams['axes.labelsize'] = {:}".format(labelsize)]
                    howto = ["df.plot(col='points', col_wrap={:}, robust=True)".format(ncolumns)]
                    howto.append("plt.show()")
                    name = attr.replace(alias+'_','')
                    plt_type = 'events_{:}_2D'.format(name) 
                    self.add_plot(alias, plt_type, howto=howto)
                    plt.close()

                for iaxis in [1,2]:
                    df.mean(axis=iaxis).plot(robust=True)
                    
                    howto = ["df.mean(axis={:}).plot(robust=True)".format(iaxis)]
                    howto.append("plt.show()")
                    name = df.dims[iaxis].replace(alias+'_','')
                    plt_type = 'events_vs_{:}'.format(name) 
                    self.add_plot(alias, plt_type, howto=howto)
                           
# Methods to build html file

    def to_html(self, path=None, quiet=False, **kwargs):
        """Write out html file
        """
        self._init_html(path=path)
        self._add_html()
        self._close_html()
        if not quiet:
            print 'Writing html to: ', os.path.join(self.html.output_dir, self.html.output_file)


    def _init_html(self, path=None, **kwargs):
        """Initialize html page.
        """
        if path:
            self._init_output(path=path, **kwargs)

        self.html = output_html.report(self.exp, 'Run {:}'.format(self.run), 
                 title=self.title,
                 css=('css/bootstrap.min.css','jumbotron-narrow.css','css/mine.css'),
                 script=('js/ie-emulation-modes-warning.js','js/jquery.min.js','js/toggler.js','js/sticky.js'),
                 output_dir=self.output_dir)

        self.html.start_block('Data Summary', id="metadata")
        self.html.start_subblock('Data Information',id='datatime')
        
        try:
            event_times = pd.to_datetime(self._xdat.time.values)
            run_time = (event_times.max()-event_times.min()).seconds
            minutes,fracseconds = divmod(run_time,60)

            sformat = 'Start Time: {:}<br/>End Time: {:}<br/>Duration: {:} seconds ({:02.0f}:{:02.0f})'
            self.html.page.p('Total events: {:}'.format(len(event_times) ) )
            self.html.page.p(sformat.format(event_times.min(), event_times.max(), run_time, minutes,fracseconds) )
            self.html.page.p('Total files: {:}<br/>Total bytes: {:} ({:0.1f} GB)<br/>'.format(
                    self._nfile,self._nbytes,self._nbytes/(1000.**3)))
        except:
            pass
            
        self.html.page.p('Report time: {:}'.format(time.ctime()))
        self.html.end_subblock()

        self._make_PyDataSource_html()

        self.html.end_block()         
 
    def _make_PyDataSource_html(self, **kwargs):
        """Make html notes for accessing data with PyDataSource.
        """
        import PyDataSource

        self.html.start_subblock('Access the Data')
        self.html.page.p('Access event data with PyDataSource python module on a psana node:')
        #ana_env = '.  /reg/g/psdm/etc/ana_env.sh'
        conda_setup = 'source /reg/g/psdm/bin/conda_setup'
        idatasource = 'idatasource --exp {:} --run {:}'.format(self.exp, self.run) 
        #self.html.add_textblock('\n'.join([ana_env, idatasource]))
        self.html.add_textblock('\n'.join([conda_setup, idatasource]))
        self.html.page.p('Where idatasource is a shortcut for starting ipython and loading the event data:')
        self.html.add_textblock('\n'.join([conda_setup, 
                                           'ipython --pylab',
                                           '...',
                                           'import PyDataSource', 
                                           'ds = PyDataSource.DataSource(exp="{:}",run={:})'.format(self.exp, self.run)]))
        self.html.end_subblock()

        try:
            pyds = PyDataSource.PyDataSource
        except:
            pyds = PyDataSource

        self.html.add_textblock(pyds.__doc__, 
                subblock='HowTo access event data with PyDataSource python package', hidden=True)
        
        self.html.page.p('Analyze run summary data on a psana node using pylab, pandas and xarray:')
#        ixarray = 'ixarray --exp {:} --run {:}'.format(self.exp, self.run)
#        self.html.add_textblock('\n'.join([conda_setup, ixarray]))
#        self.html.page.p('Where ixarray is a shortcut for starting ipython and loading the data summary:')
        self.html.add_textblock('\n'.join([conda_setup, 
                                           'ipython --pylab',
                                           '...',
                                           'import PyDataSource', 
                                           'x = PyDataSource.open_h5netcdf(exp="{:}",run={:})'.format(self.exp, self.run)]))
        self.html.add_textblock(str(self._xdat), 
                subblock='View of RunSummary data with xarray python package ', hidden=True)
        weblink='http://pswww.slac.stanford.edu/swdoc/ana/PyDataSource'
        self.html.page.p('See <a href="{:}">{:}</a> for details on how to use PyDataSource.'.format(weblink,weblink))
        weblink='http://xarray.pydata.org'
        self.html.page.p('See <a href="{:}">{:}</a> for details on how to use data using xarray.'.format(weblink,weblink))
        self.html.page.p('HDF5 summary file (using netCDF format) located at:')
        self.html.add_textblock('{:}'.format(self.nc_file))
        self.html.page.p('For questions and feedback contact koglin@slac.stanford.edu')
 
    def _add_html(self, table_caption='Hover to highlight.', show_attrs=['attrs'], **kwargs):
        """Add html pages for results.
        """
        cat_items = sorted(self.results.items(), key=operator.itemgetter(0))
        for catagory, item in cat_items:
            self.html.start_block('{:} Data'.format(catagory), id='{:}_data'.format(catagory))
            ahowto = []
            
            if item['figure']:
                ahowto.append('# For interactive plotting -- plt.ioff() to turn off interactive plotting.')
                ahowto.append('plt.ion()')
                ahowto.append('# Alternatively make plt.show() after each plot and close window to make next')

            datatyp = 'text'
            data_items = sorted(item[datatyp].items(), key=operator.itemgetter(0))
            for name, data in data_items:
                howto = data.get('howto')
                if howto is not None:
                    if not isinstance(howto, list):
                        howto = [howto]
                    howto_step = "# Howto {:} {:}:\n".format(name, catagory)
                    howto_step += '\n'.join(howto)
                    ahowto.append(howto_step)

            datatyp = 'textblock'
            data_items = sorted(item[datatyp].items(), key=operator.itemgetter(0))
            for name, data in data_items:
                text = data.get('text')
                if not isinstance(text, list):
                    text = [text]
                text_data = '\n'.join(text)
                if text_data:
                    self.html.add_textblock(text_data, 
                        subblock='{:} {:}'.format(catagory,name), 
                        hidden=True)

            datatyp = 'table'
            data_items = sorted(item[datatyp].items(), key=operator.itemgetter(0))
            for name, data in data_items:
            #for name, data in item[datatyp].items():
                howto = data.get('howto')
                if howto is not None:
                    if not isinstance(howto, list):
                        howto = [howto]
                    howto_step = "# Howto make the {:} {:} {:}:\n".format(catagory, name, datatyp)
                    howto_step += '\n'.join(howto)
                else:
                    howto_step = ''
                
                if data.get('format'):
                    formatters = {a: b.format for a,b in data.get('format').items()}
                    formatterstr = {a: b+'.format' for a,b in data.get('format').items()}
                else:
                    formatters = None

                df = data.get('DataFrame')
                doc = data.get('doc')
                dfname = data.get('name')
                if df is not None:
                    if name in show_attrs:
                        hidden = False
                    else:
                        hidden = True
                   
                    pd.set_option('display.max_rows', len(df))
                    if formatters:
                        dfstr = df.to_string(justify='left',formatters=formatters)
                        #if dfname:
                        #    howto_step += '\n# to reprsent with formatting'
                        #    howto_step += "\nprint {:}.to_string(justify='left',formatters={:})".format(dfname, formatterstr)
                    else:
                        dfstr = str(df)
                        #if dfname:
                        #    howto_step += "\n print {:}".format(dfname)
                    
                    if howto_step:
                        ahowto.append(howto_step)

                    self.html.add_textblock(dfstr, doc=doc, 
                            subblock='{:} {:} {:}'.format(catagory,name,datatyp), 
                            hidden=hidden)
                    
                    pd.reset_option('display.max_rows')

            datatyp = 'figure'
            data_items = sorted(item[datatyp].items(), key=operator.itemgetter(0))
            for name, data in data_items:
            #for name, data in item[datatyp].items():
                png = data.get('png')
                if png:
                    self.html.start_subblock('{:} {:} {:}'.format(catagory, name, datatyp), hidden=True)
                    self.html.page.a( markup.oneliner.img(src=png,style='width:100%;'), 
                            href=png )
                    self.html.end_subblock(hidden=True)
            
                howto = data.get('howto')
                if howto:
                    if not isinstance(howto, list):
                        howto = [howto]
                    howto_step = "# Howto make the {:} {:} {:}:\n".format(catagory, name, datatyp)
                    howto_step += '\n'.join(howto)
                    ahowto.append(howto_step)

            if ahowto:
                self.html.add_textblock(ahowto, 
                        subblock='HowTo make {:} tables and figures.'.format(catagory), 
                        id='howto_{:}'.format(catagory.replace(' ','_')), 
                        hidden=True)

            self.html.end_block()         

    def _close_html(self, **kwargs):
        """Close html files.
        """
        # this closes the left column
        self.html.page.div.close()
        
        self.html.mk_nav()
        
        self.html._finish_page()
        self.html.myprint(tofile=True)

#    def make_default_cuts(self, gasdetcut_mJ=0.5):
#        """
#        Make default cuts.
#        """
#        x = self._xdat
#        # FEEGasDetEnergy_f_12_ENRC is duplicate measurement -- can average if desired 
#        try:
#            attr = 'FEEGasDetEnergy_f_11_ENRC'
#            x['Gasdet_pre_atten'] = (['time'], x[attr].values)
#            x['Gasdet_pre_atten'].attrs['doc'] = "Energy measurement before attenuation ({:})".format(attr)
#            for a in ['unit', 'alias']: 
#                try:
#                    x['Gasdet_pre_atten'].attrs[a] = x[attr].attrs[a]  
#                except:
#                    pass
#
#            # FEEGasDetEnergy_f_22_ENRC is duplicate measurement -- can average if desired 
#            attr = 'FEEGasDetEnergy_f_21_ENRC'
#            x['Gasdet_post_atten'] = (['time'], x[attr].values)
#            x['Gasdet_post_atten'].attrs['doc'] = "Energy measurement afeter attenuation ({:})".format(attr)
#            for a in ['unit', 'alias']: 
#                try:
#                    x['Gasdet_post_atten'].attrs[a] = x[attr].attrs[a]  
#                except:
#                    pass
#        
#            x = x.drop(['FEEGasDetEnergy_f_11_ENRC', 'FEEGasDetEnergy_f_12_ENRC',
#                        'FEEGasDetEnergy_f_21_ENRC', 'FEEGasDetEnergy_f_22_ENRC'])
#
#        except:
#            pass
#
##      Need to add in experiment specific PulseEnergy with attenuation
##        try:
##            x['PulseEnergy'] = (['time'], x['Gasdet_post_atten'].values*x['dia_trans1'].values)
##        except:
##            x['PulseEnergy'] = (['time'], x['Gasdet_post_atten'].values*x.attrs['dia_trans1'])
##        
##        x['PulseEnergy'].attrs = x['Gasdet_pre_atten'].attrs
##        x['PulseEnergy'].attrs['doc'] = "Energy measurement normalized by attenuators"
##
#        try:
#            gasdetcut =  np.array(x.Gasdet_pre_atten.values > gasdetcut_mJ, dtype=byte)
#            x.coords['Gasdet_cut'] = (['time'], gasdetcut)
#            x.coords['Gasdet_cut'].attrs['doc'] = "Gas detector cut.  Gasdet_pre_atten > {:} mJ".format(gasdetcut_mJ)
#        except:
#            gasdetcut = np.ones(x.time.data.shape)
#
#    #    damagecut = np.array(phasecut & gasdetcut & (x.EBeam_damageMask.values == 0), dtype=byte)
#        damagecut = np.array(gasdetcut, dtype=byte)
#        x.coords['Damage_cut'] = (['time'], damagecut)
#        x.coords['Damage_cut'].attrs['doc'] = "Combined Gas detector, Phase cavity and EBeam damage cut"
#
#        try:
#            x = x.rename( {'EBeam_ebeamPhotonEnergy': 'PhotonEnergy'} )
#        except:
#            pass
#
#        if not x.attrs.get('correlation_variables'):
#            cvars = []
#            for cvar in ['PhotonEnergy','Gasdet_post_atten']:
#                if cvar in x:
#                    cvars.append(cvar)
#            x.attrs['correlation_variables'] = cvars
#
#        try:
#            sattrs = list(set(x.attrs.get('scan_variables',[])))
#            for pv in [a.replace('_steps','') for a in x.keys() if a.endswith('_steps')]:
#                if pv not in x:
#                    x.coords[pv] = (('time',), x[pv+'_steps'].data[list(x.step.data.astype(int))])
#                    x.coords[pv].attrs = x[pv+'_steps'].attrs
#                    sattrs.append(pv)
#        except:
#            print 'Cannot add _steps as events'
#
#        try:
#            for pv in x.attrs.get('pvControls', []):
#                if pv+'_control' in x:
#                    sattrs.append(pv+'_control')
#                elif pv in x:
#                    sattrs.append(pv)
#            sattrs = list(set([attr for attr in sattrs if attr in x and len(set(x[attr].data)) > 3]))
#            x.attrs['scan_variables'] = sattrs
#
#        except:
#            print 'Cannot add pvControls to scan_variables'
#
#        try:
#            if 'XrayOff' in x and x.XrayOff.data.sum() > 0 and not x.attrs.get('cuts'):
#                x.attrs['cuts'] = ['XrayOn','XrayOff']
#                print 'Setting cuts', x.attrs['cuts']
#
#        except:
#            print 'Cannot make default cuts'
#
#
#        self._xdat = x
#
#        return x




#                Axes = pd.tools.plotting.scatter_matrix(dfscat, figsize=figsize)
#                #y ticklabels
#                [plt.setp(item.yaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
#                #x ticklabels
#                [plt.setp(item.xaxis.get_majorticklabels(), 'size', 15) for item in Axes.ravel()]
#                #y labels
#                [plt.setp(item.yaxis.get_label(), 'size', 20) for item in Axes.ravel()]
#                #x labels
#                [plt.setp(item.xaxis.get_label(), 'size', 20) for item in Axes.ravel()]


