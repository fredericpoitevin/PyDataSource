from __future__ import print_function
def get_run_from_id(run_id, exp):
    from RegDB import experiment_info
    instrument = exp[0:3]
    exp_runs = experiment_info.experiment_runs(instrument.upper(),exp)
    run_dict = {int(item['id']): int(item['num']) for item in exp_runs}

    return run_dict.get(run_id, run_id)
 
if __name__ == "__main__":
    from sys import argv
    exp = argv[1]
    run = int(argv[2])
    print(get_run_from_id(run, exp))


