from __future__ import print_function

def get_config_file(run=None, exp=None, instrument=None, path=None):
    """
    Get the configuration file for a run and experiment
    """
    import os
    if not path:
        if not instrument:
            instrument = exp[0:3]
        base_path = '/reg/d/psdm/{:}/{:}/results'.format(instrument,exp)
        if not os.path.isdir(base_path):
            base_path = '/reg/d/psdm/{:}/{:}/res'.format(instrument,exp)
        path = base_path+'/summary_config'

    if not run:
        run = 9999
    else:
        run = int(run)

    while run > 0:
        file_name = '{:}/run{:04}.config'.format(path, int(run))
        if os.path.isfile(file_name):
            return file_name
        run -= 1

    return ''

if __name__ == "__main__":
    from sys import argv
    exp = argv[1]
    run = argv[2]
    print(get_config_file(run, exp))

