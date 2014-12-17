import gkl_files_psi_count
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
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.01', '-R', '1000', \
    '-x1', '0', '-z0', '0', '-transient', '1000', '-report', 'psi_count'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.06', '-R', '1000', \
    '-x1', '1', '-z0', '0', '-transient', '1000', '-report', 'psi_count'])
experiments.append(['s4', '-L', '400', '-T', '9000', '-n', '0.06', '-R', '1000', \
    '-transient', '1000', '-report', 'psi_count'])

######
######

__graph_layout_rows = 3
__graph_layout_cols = 3
__graph_cols_to_be_plotted = [1, 2, 3, 4, 5, 6, 7, 8, 9]
__graph_experiments_to_be_plotted = [0, 1, 2]
__graph_experiment_colors = ['g', 'r', 'b']
__graph_experiment_alphas = [1.0, 0.6, 0.4]

run(gkl = ['..\\Release\\gkl.exe'], experiments, append = False)
gkl_files_psi_count.plot_experiments(experiments, __graph_layout_rows, \
    __graph_layout_cols, __graph_cols_to_be_plotted, \
    __graph_experiments_to_be_plotted, __graph_experiment_colors, \
    __graph_experiment_alphas)
gkl_files_psi_count.show_psi_prob_experiments(experiments)

############################
############################