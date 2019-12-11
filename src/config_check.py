from __future__ import print_function

class ConfigCheck(object):
    """
    Check the detector config and trigger timing of PyDataSource.DataSource object.

    """
    def __init__(self, configData, 
            config_file='config_alert.json', 
            trigger_file='trigger_alert.json', 
            path='', 
            **kwargs):
        """
        Parameters:
        -----------
        config_file : str
            json file for daq alert config 
        trigger_file : str
            json file for trigger alertconfig
        path : str
            path of config_file and trigger file [Default = '/reg/g/psdm/utils/arp/config',instrument.lower()]
        """
        import os
        import pandas as pd
        self._configData = configData
        instrument = configData._ds.instrument
        if not instrument:
            calibDir = configData._ds._ds.env().calibDir()
            instrument = calibDir.split('/')[4]
            if not instrument:
                instrument = 'default'

        default_path = '/reg/g/psdm/utils/arp/config/default'
        if not path:
            path = os.path.join('/reg/g/psdm/utils/arp/config',instrument.lower())
        
        filename = os.path.join(path,trigger_file)
        try:
            self.trigger_dict = pd.read_json(filename).to_dict()
        except:
#            import inspect
#            default_path = os.path.dirname(inspect.getfile(self.__class__))
            print('ERROR: cannot read json file {:}'.format(filename))
            trigger_file = 'trigger_alert.json'
            filename = os.path.join(default_path,trigger_file)
            print('reading default {:}'.format(filename))
            self.trigger_dict = pd.read_json(filename).to_dict()
        
        filename = os.path.join(path,config_file)
        try:
            self.config_dict = pd.read_json(filename).to_dict()
        except:
#            import inspect
#            default_path = os.path.dirname(inspect.getfile(self.__class__))
            print('ERROR: cannot read json file {:}'.format(filename))
            config_file = 'config_alert.json'
            filename = os.path.join(default_path,config_file)
            print('reading default {:}'.format(filename))
            self.config_dict = pd.read_json(filename).to_dict()
        
        self._clear_alerts()
        self._clear_warnings()
        self._set_trigger_alerts()
        self._set_config_alerts()

    @property
    def det_info(self):
        dets = {}
        for alias, srcstr in self._configData._config_srcs.items():
            src_info = self._configData._sources.get(srcstr,{})
            group = src_info.get('group')
            if not group or group == -1:
                # controls camera do not have evr info
                continue
            try:
                srcname = srcstr.split('(')[1].split(')')[0]
                devName = srcname.split(':')[1].split('.')[0]
                dets[alias] = { 
                    'srcstr': srcstr,
                    'srcname': srcname,
                    'devName': devName,
                    }
            except:
                pass
        return dets

    def get_alert_info(self, typ='error'):
        """
        Get dict of alert info for each detector

        Parameters
        ----------
        typ : str
            type of alert info [default = 'error']
            options include ['error','value','valid_setting','doc','unit']
        """
        return {a[0]: {b: c.get(typ) for b,c in a[1].items()} for a in self._alert_dict.items()}

    def get_warning_info(self, typ='error'):
        """
        Get dict of warning info for each detector

        Parameters
        ----------
        typ : str
            type of warning info [default = 'error']
            options include ['error','value','valid_setting','doc','unit']
        """
        return {a[0]: {b: c.get(typ) for b,c in a[1].items()} for a in self._warning_dict.items()}
    
    @property
    def alerts(self):
        """
        List of alerts for each detector
        """
        return {a[0]: ['{:}: {:}'.format(b,c.get('error')) for b,c in a[1].items()] 
                for a in self._alert_dict.items()}

    @property
    def warnings(self):
        """
        List of warnings for each detector
        """
        return {a[0]: ['{:}: {:}'.format(b,c.get('error')) for b,c in a[1].items()] 
                for a in self._warning_dict.items()}

    def show_info(self, **kwargs):
        from psmessage import Message
        message = Message(quiet=True, **kwargs)
        if self.alerts or self.warnings:
            header = '{:15} {:14} {:8} {:22}'.format('detector',  'parameter', 'level', 'error') 
            message(header)
            message('-'*(len(header)))
        if self.alerts:
            for alias, alerts in self.get_alert_info().items():
                for attr, alert in alerts.items():
                    message('{:15} {:14} {:8} {:22}'.format(alias,attr,'alert',alert)) 
        if self.warnings:
            for alias, alerts in self.get_warning_info().items():
                for attr, alert in alerts.items():
                    message('{:15} {:14} {:8} {:22}'.format(alias,attr,'warning',alert)) 
 
        return message

    def __str__(self): 
        return '{:}: {:}'.format(self.__class__.__name__, self._configData._ds)

    def __repr__(self):
        print('< '+str(self)+' >')
        print(self.show_info())
        return '< '+str(self)+' >'

    def _clear_warnings(self):
        self._warning_dict = {}

    def _clear_alerts(self):
        self._alert_dict = {}

    def _set_alert(self, alias, attr, value, **kwargs):
        if alias not in self._alert_dict:
            self._alert_dict[alias] = {}
        self._alert_dict[alias][attr] = {'value': value}
        if kwargs:
            self._alert_dict[alias][attr].update(**kwargs)

    def _set_warning(self, alias, attr, value, **kwargs):
        if alias not in self._warning_dict:
            self._warning_dict[alias] = {}
        self._warning_dict[alias][attr] = {'value': value}
        if kwargs:
            self._warning_dict[alias][attr].update(**kwargs)

    def _set_trigger_alerts(self, types=['alert','warning']):
        """
        Set trigger alerts and warnings
        """
        import numpy as np
        units = {'delay': 'sec', 'width': 'sec'}
        for alias, info in self.det_info.items():
            devName = info.get('devName')
            sourceData = getattr(self._configData.Sources, alias)
            if not sourceData:
                continue
            for typ in types:
                alert_info = self.trigger_dict.get(devName,{}).get(typ,{})
                # json read will put nan if warning or alert empty for devName when set for others
                if not alert_info or alert_info is np.nan:
                    continue
                for attr, avalue in alert_info.items():
                    set_value = getattr(sourceData, attr)
                    if attr == 'evr_polarity':
                        set_value = {0: 'Pos', 1: 'Neg'}.get(set_value, set_value)
                    doc=None
                    unit = units.get(attr)
                    if isinstance(avalue, list):
                        if set_value < avalue[0]:
                            err = '{:} < {:}'.format(set_value, avalue[0])
                            doc = '{:} evr {:} = {:} is too low -- valid range = {:}'.format(alias, 
                                        attr, set_value, avalue)
                        elif set_value > avalue[1]:
                            err = '{:} > {:}'.format(set_value, avalue[1])
                            doc = '{:} evr {:} = {:} is too high -- valid range = {:}'.format(alias, 
                                        attr, set_value, avalue)
                    else:
                        if set_value != avalue:
                            err = '{:} != {:}'.format(set_value, avalue)
                            doc = '{:} evr {:} = {:} is incorrect -- valid setting = {:}'.format(alias, 
                                        attr, set_value, avalue)
                    
                    if doc:
                        if typ == 'alert':
                            self._set_alert(alias, attr, set_value, valid_setting=avalue, 
                                    doc=doc, error=err, unit=unit) 
                        elif typ == 'warning':
                            self._set_warning(alias, attr, set_value, valid_setting=avalue, 
                                    doc=doc, error=err, unit=unit) 
                

    def _set_config_alerts(self, types=['alert','warning']):
        """
        Set config alerts and warnings
        """
        import numpy as np
        for alias, info in self.det_info.items():
            devName = info.get('devName')
            det_config = getattr(self._configData, alias)
            if not det_config:
                continue
            attr_info = getattr(det_config, '_attr_info')
            for typ in types:
                alert_info = self.config_dict.get(devName,{}).get(typ,{})
                # json read will put nan if warning or alert empty for devName when set for others
                if not alert_info or alert_info is np.nan:
                    continue
                for attr, avalue in alert_info.items():
                    ainfo = attr_info.get(attr)
                    set_value = ainfo.get('value')
                    unit = ainfo.get('unit')
                    doc=None
                    if isinstance(avalue, list):
                        if set_value < avalue[0]:
                            err = '{:} < {:}'.format(set_value, avalue[0])
                            doc = '{:} {:} = {:} {:} is too low -- valid range = {:}'.format(alias, 
                                        attr, set_value, unit, avalue)
                        elif set_value > avalue[1]:
                            err = '{:} > {:}'.format(set_value, avalue[1])
                            doc = '{:} {:} = {:} {:} is too high -- valid range = {:}'.format(alias, 
                                        attr, set_value, unit, avalue)
                    else:
                        if set_value != avalue:
                            err = '{:} != {:}'.format(set_value, avalue)
                            doc = '{:} {:} = {:} {:} is incorrect -- valid setting = {:}'.format(alias, 
                                        attr, set_value, unit, avalue)
                    
                    if doc:
                        if typ == 'alert':
                            self._set_alert(alias, attr, set_value, 
                                    valid_setting=avalue, doc=doc, error=err, unit=unit) 
                        elif typ == 'warning':
                            self._set_warning(alias, attr, set_value, 
                                    valid_setting=avalue, doc=doc, error=err, unit=unit) 
                

