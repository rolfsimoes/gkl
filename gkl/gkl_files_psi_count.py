import matplotlib.pyplot as plt
import sys
    
def show_psi_prob_experiments(experiments):
    __report_signature = '2'
    __report_gkl_s4 = 's4'
    __from_rules = [0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, \
        0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, \
        2, 3, 3, 3, 3, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
    __to_rules = [0, 0, 0, 0, 3, 1, 3, 1, 0, 2, 0, 0, 0, 2, 0, 0, 3, 3, 3, 3, \
        3, 1, 3, 1, 2, 1, 2, 2, 2, 1, 2, 2, 3, 3, 3, 3, 3, 1, 3, 1, 2, 1, 2, 2, \
        2, 1, 2, 2, 0, 0, 0, 0, 3, 1, 3, 1, 2, 1, 2, 2, 2, 1, 2, 2]

    print '"From-To" probabilities:\n'
    for experiment_index in xrange(len(experiments)):
        mean = [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], \
            [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], \
            [0.0, 0.0, 0.0, 0.0]]
        mean_2 = [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], \
            [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], \
            [0.0, 0.0, 0.0, 0.0]]
        ##File name
        file_name = ''
        for parameter in experiments[experiment_index]:
            file_name += str(parameter).replace('.', '_').replace('-', '_')
        file_name += '.txt'
        ##Open file
        try:
            logfile = open(file_name, "r")
        except:
            print 'File', file_name, 'was not found!\n'
            sys.exit(1)
        ##With file...
        with logfile:
            row_index = 0
            for row in logfile:
                row_index += 1
                ##Report's version signature
                if row_index == 1:
                    fields = row.split()
                    if (fields[0] != '#GKL_report_version:') or (fields[1] != __report_signature):
                        print 'This report version is not supported!\n'
                        sys.exit(1)
                if row_index == 2:
                    fields = row.split()
                    if (fields[0] != '#Command_line:') or (fields[2] != __report_gkl_s4):
                        print 'This report type is not supported!\n'
                        sys.exit(1)
                ##Report's parse
                if row[0] != '#' and row[0] != '\n':
                    cols = row.split()
                    for rule_index in xrange(64):
                        mean[__from_rules[rule_index]][__to_rules[rule_index]] += float(cols[rule_index + 1])
                        mean_2[__from_rules[rule_index]][__to_rules[rule_index]] += pow(float(cols[rule_index + 1]), 2)
            ##Print experiment data
            print 'from\tprob_to_0\tprob_to_1\tprob_to_2\tprob_to_3\tsd_to_0\tsd_to_1\tsd_to_2\tsd_to_3',
            for row_index in xrange(4):
                print '\n' + str(row_index),
                for col_index in xrange(4):
                    print '\t' + str(mean[row_index][col_index] / 1000),
                for col_index in xrange(4):
                    print '\t' + str(pow(pow(mean[row_index][col_index] / 1000, 2) - mean_2[row_index][col_index] / 1000, 0.5)),
            print '\n'
    print 'Finished!\n'
    
    
def plot_psi_count_experiments(experiments, __graph_layout_rows, __graph_layout_cols, \
    __graph_cols_to_be_plotted, __graph_experiments_to_be_plotted, \
    __graph_experiment_colors, __graph_experiment_alphas):
    
    __report_signature = '2'
    print 'Ploting graphs...\n'
    graphs = []
    fig = plt.figure()
    for experiment_index in xrange(len(experiments)):
        data = []
        ##File name
        file_name = ''
        for parameter in experiments[experiment_index]:
            file_name += str(parameter).replace('.', '_').replace('-', '_')
        file_name += '.txt'
        ##Open file
        try:
            logfile = open(file_name, "r")
        except:
            print 'File', file_name, 'was not found!\n'
            sys.exit(1)
        ##With file...
        with logfile:
            row_index = 0
            for row in logfile:
                row_index += 1
                ##Report's version signature
                if row_index == 1:
                    fields = row.split()
                    if (fields[0] != '#GKL_report_version:') or (fields[1] != __report_signature):
                        print 'This report version is not supported!\n'
                        sys.exit(1)
                ##Report's parse
                if row[0] != '#' and row[0] != '\n':
                    cols = row.split()
                    if len(data) == 0:
                        while (len(data) < len(cols)):
                            data.append([])
                    for col_index in xrange(len(cols)):
                        data[col_index].append(float(cols[col_index]))
            ##Plot experiment data
            for col_index in xrange(len(__graph_cols_to_be_plotted)):
                graphs.append(fig.add_subplot(__graph_layout_rows, __graph_layout_cols, col_index + 1))
                graphs[-1].set_title('Rule: ' + str(__graph_cols_to_be_plotted[col_index]))
                graphs[-1].hist(data[__graph_cols_to_be_plotted[col_index]], 25, \
                    facecolor = __graph_experiment_colors[experiment_index], \
                    alpha = __graph_experiment_alphas[experiment_index], \
                    linewidth = 0)
    print 'Finished!\n'
    plt.show()