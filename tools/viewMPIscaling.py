'''

MPI STRONG AND WEAK SCALING
Tables and visualization

Tests strong and weak scaling for a given mesh or set of meshes.

'''

import argparse
import numpy as np
import sys
import matplotlib
from copy import deepcopy
import subprocess
from matplotlib.ticker import MultipleLocator, FuncFormatter

matplotlib.use("Agg")

from matplotlib import pyplot as plt
import os

def get_node_count(infile):
    '''
    Returns the node count of an AVS-UCD mesh.
    '''
    assert infile.split('.')[-1].lower() == 'inp','Must be an AVS-UCD mesh'

    with open(infile,'r') as f:
        header = f.readline().strip().split()

    return int(header[0])

def get_runtime(exe,runtime_file='._runtime'):
    '''
    Returns the runtime of an executable.
    '''
    with open(os.path.join(os.path.dirname(exe),runtime_file),'r') as rtf:
        elapsed_time = float(rtf.read())
    return elapsed_time

def strong_scaling(infile,mpi_nodes=[1,2,4,8,16,32,64],outtypes=['fehm','pflotran','tough2','hdf5'],imgout=None,view=True,exe=None):
    '''
    Generates a strong scaling graph and STDOUT table.

    :param infile: Mesh to test against
    :type infile: str
    :param mpi_nodes: A list containing the proc scaling
    :type mpi_nodes: list<int>
    :param outtypes: A list of the output formats to test on
    :type outtypes: list<str>
    :param imgout: Path to save scaling image to
    :type imgout: str
    :param view: View the image on screen or not
    :type view: bool
    :param exe: path to VORONOI executable
    :type exe: str
    '''

    # Set up initial parameters
    ext = infile.split('.')[-1].lower()

    #if otype is None: otype = 'fehm'
    if exe is None: exe = '../voronoi'
    outfile = '._temp'

    # Verify that the correct filetype is being passed in
    # Change to exclusively avs?
    assert ext in ['inp','lgi'],'Unsupported file format!'
    ext = 'avs' if ext == 'inp' else 'lg'

    # Construct the Matplotlib figure
    fig, ax = plt.subplots()

    total_runtimes = []

    colors = ['black','red','green','yellow']

    mesh_nodes = get_node_count(infile)
    for (k,otype) in enumerate(outtypes):
        runtimes = []
        cmds = ['mpiexec -np %d %s -%s %s -o %s -type %s' % (proc,exe,ext,infile,outfile,otype) for proc in mpi_nodes]
        for (ll,cmd) in enumerate(cmds):
            print(cmd)
            try:
                out = subprocess.check_output(cmd.split(),stderr=subprocess.STDOUT)
                if 'error' in out.lower():
                    print('ERROR')
                    raise Exception
            except subprocess.CalledProcessError as e:
                print('sdf')
                raise 
            
            runtimes.append(get_runtime(exe))

        # Calculate the scaling efficiency as a function of:
        # e(N) = t(1) / (N * t(N)) * 100%
        scaling_efficiency = [100.*runtimes[0]/((i+1)*runtimes[i]) for i in range(len(runtimes))]
        total_runtimes.append(runtimes)

        ax.plot(mpi_nodes, scaling_efficiency, color=colors[k], linestyle='solid', marker='o', markerfacecolor=colors[k], markersize=3, label=otype.upper())

    plt.xticks(mpi_nodes,[("%.0f" % x)  for x in mpi_nodes])
    ax.set_xlabel("Number of processing elements")
    ax.set_ylabel("Parallel Efficiency (% of linear scaling)")
    ax.set_title("Parallel Scaling: Strong Scaling on $10^6$ Vertex Mesh")
    legend = ax.legend(loc='upper right', shadow=False, fontsize='medium')

    ax.grid(True)

    # Save the figure?
    if imgout is not None:
        plt.savefig(imgout,format='eps', dpi=1000)

    for j in range(len(mpi_nodes)):
        a += '  %2d     ' % mpi_nodes[j] + '    '.join(['%4.4f' % a[i][j] for i in range(len(a))])

    a = '''
   strong scaling efficiency   
         meshinput.inp
===============================
   n         t(n)       E_fehm(n)
-------------------------------
   1        0.763s     1.0000
   2        0.763s     0.7540
   4        0.763s     0.7540
   8        0.763s     0.7540
   16       0.763s     0.7540
   32       0.763s     0.7540
   64       0.763s     0.7540'''

    # Show the figure?
    if view:
        plt.show()

strong_scaling("/scratch/sft/livingston/voronoi_lagrit/test/2D/mpi_test_1_delete/mesh.inp",imgout="TESTME.eps")


'''
if (__name__ == "__main__"):

    # Check this out for 'true' SM vis: https://github.com/nschloe/betterspy

    #applet_description = 

    parser = argparse.ArgumentParser(description=applet_description,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('filename',help='input filename (.stor; .lgi,.avs.)',type=str)
    parser.add_argument('-c','--cmap',help='Name of a PyPlot colormap (default: \'binary\')',type=str,default='binary')
    parser.add_argument('-s','--save',help='Filepath to save image (default: visualize only)',type=str)
    parser.add_argument('-v','--view',help='Flag to display image (default: true)',type=str,default='true')
    parser.add_argument('-b','--binary',help='Flag to view sparse matrix as binary or gradient (default: true)',type=str,default='true')
    args = parser.parse_args()

    if args.filename is None:
        parser.print_help()
        sys.exit()

    view = (args.view.lower()[0] == 't' or args.view.lower()[0] == 'y' or args.view.lower()[0] == '1')
    binary = (args.binary.lower()[0] == 't' or args.view.lower()[0] == 'y' or args.view.lower()[0] == '1')

    m = get_sparse_matrix(args.filename)
    #m = np.random.rand(100,100)*100.
    plot_sparse_matrix(m,outfile=args.save,binary=binary,cmap=args.cmap,viewfig=view)
'''