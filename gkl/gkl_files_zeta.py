from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
import sys

def plot_zeta_experiments(experiments):
    
    __report_signature = '2'
    print 'Ploting graphs...\n'
    x = []
    y = []
    z = []
    fig = plt.figure()
    graph = fig.add_subplot(1, 1, 1, projection = '3d')
    for experiment_index in xrange(len(experiments)):
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
                    x.append(float(cols[0]))
                    try:
                        noise_index = experiments.index('-n')
                    except:
                        noise_index = -1
                    if noise_index > -1:
                        y.append(float(experiments[noise_index + 1]))
                    else:
                        y.append(0.0)
                    z.append(float(cols[1]))
    ##Plot experiment data
    graph.scatter(x, y, z, marker = 'o', linewidth = 0.25)
    plt.show()
    
    