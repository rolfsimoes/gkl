import gkl_files_zeta
import subprocess
import time
import os

def run(gkl, experiments, append = False):
    print 'Making experiments...\n'
    max_procs = 4
    actual_procs = 0
    procs = []
    logfiles = []
    for experiment_index in xrange(len(experiments)):
        ##File name
        file_name = ''
        for parameter in experiments[experiment_index]:
            file_name += str(parameter).replace('.', '_').replace('-', '_')
        file_name += '.txt'
        ##Open file...
        try:
            if append:
                logfile = open(file_name, 'a')
                logfile.seek(0, os.SEEK_END)
            else:
                logfile = open(file_name, 'w')
        except:
            print 'Cannot open the file', file_name, 'in write mode!\n'
            sys.exit(1)
        logfiles.append(logfile)
        ##With file opened, start a new process
        procs.append(subprocess.Popen(gkl + experiments[experiment_index], stdout = logfiles[-1], shell = False))
        actual_procs += 1
        while (actual_procs >= max_procs):
            time.sleep(0.25)
            actual_procs = 0
            for proc in xrange(len(procs)): 
                if (procs[proc].poll() == None): actual_procs += 1 
                else: logfiles[proc].close()
    while (actual_procs > 0):
        time.sleep(0.25)
        actual_procs = 0
        for proc in xrange(len(procs)): 
            if (procs[proc].poll() == None): actual_procs += 1 
            else: logfiles[proc].close()
    print 'Finished!\n'



###########################
###########################

experiments = []

experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.005', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.01', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.015', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.02', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.025', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.03', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.035', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.04', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.045', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.05', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.055', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.06', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.065', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.07', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.075', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.08', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.085', '-R', '1000', \
    '-transient', '1000', '-report', 'zeta'])

######
######

run(gkl = ['..\\Release\\gkl.exe'], experiments, append = False)
gkl_files_zeta.plot_zeta_experiments(experiments)