from __future__ import print_function
# Basic 

def post_report(exp, run, report_list=[], report_name=None, weblink=None):
    """
    Helper method for making ARP report posts
    
    Parameters
    ----------
    report_list : list
        e.g., [[name1, info1, color1], [name2, info2]]
    """
    from requests import post
    import os
    from os import environ
    from RegDB import experiment_info
    expNum = experiment_info.name2id(exp)
    instrument = exp[0:3]
    
    update_url = environ.get('BATCH_UPDATE_URL')
    print('Updating ARP URL dict...')
    print(update_url)
    try:
        if not weblink and report_name:
            webattrs = [instrument.upper(), expNum, exp, report_name, 'report.html'] 
            weblink='http://pswww.slac.stanford.edu/experiment_results/{:}/{:}-{:}/{:}/{:}'.format(*webattrs)

    except:
        print('Error making weblink')

    if not isinstance(report_list[0], list):
        report_list = [report_list]

    batch_counters = {}
    try:
        for item in report_list:
            name = item[0]
            info = item[1]

            if weblink:
                batch_attr = '<a href={:}>{:}</a>'.format(weblink,name)
            else:
                batch_attr = name

            if len(item) > 2:
                color = item[2]
                batch_counters[batch_attr] = [info, color]
            else:
                batch_counters[batch_attr] = [info]
        
        if update_url and batch_counters:
            post(update_url, json={'counters' : batch_counters})
            print('Done making ARP counters')
        else:
            print('No ARP counters to update')
    
    except:
        print('Error making ARP counters')


