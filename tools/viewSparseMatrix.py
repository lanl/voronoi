import argparse
import numpy as np
import sys
from matplotlib import pyplot as plt

def get_header(f, skip_two):

    if skip_two:
        f.readline() # fehmstor ascir8i4 LaGriT Sparse Matrix Voronoi Coefficients
        f.readline() # Mon Oct 23 17:19:23 2017

    header = f.readline().strip().split()

    NUM_WRITTEN_COEFS = int(header[0])
    NEQ = int(header[1])
    NCOEF_NEQ_1 = int(header[2])
    NUM_AREA_COEF = int(header[3])
    try:
        NCON_MAX = int(header[4])
    except:
        NCON_MAX = None
    
    NCOEF = NCOEF_NEQ_1 - NEQ - 1
    return NUM_WRITTEN_COEFS,NEQ,NCOEF,NUM_AREA_COEF,NCON_MAX

def get_matrix(infile):
    
    f = open(infile)

    NUM_WRITTEN_COEFS,NEQ,NCOEF,NUM_AREA_COEF,NCON_MAX = get_header(f,True)
    
    matrix = np.zeros((NEQ,NEQ),dtype=np.double)
    vols = []
    row_count = []
    values = ""#np.array((0,0))#[]
    row_entries = []
    pointers = []
    zeros = []
    ptrs_diag = []

    for line in f:
        #values = values + line.strip().split()
        #values = np.append(values,line.strip().split())
        values = values + line
    f.close()

    values = values.strip().split()
    
    #------------------------------------------------#

    # There are seven sections of a STOR file after the header
    # Record all of them, and leave the remainder as area coeffs.
    for i in range(0,NEQ):
        vols.append(float(values.pop(0)))

    for i in range(0,NEQ+1):
        row_count.append(int(values.pop(0)))

    for i in range(0,NCOEF):
        row_entries.append(int(values.pop(0)))
    
    for i in range(0,NCOEF):
        pointers.append(int(values.pop(0)))

    for i in range(0,NEQ+1):
        zeros.append(int(values.pop(0)))

    for i in range(0,NEQ):
        ptrs_diag.append(int(values.pop(0)) - 5)

    #------------------------------------------------#

    tmp = row_count[:]
    for i in range(1,len(row_count)):
        row_count[i] = row_count[i] - tmp[i-1]
    row_count.pop(0)

    # Fill the matrix with area coeffs.
    a = 0
    bb = 0

    for i in range(0,len(row_count)):
        b = row_count[i]
        aa = row_entries[a:b+a]
        a = a + row_count[i]
        
        for j in range(0,len(aa)):
            matrix[i][aa[j]-1] = float(values[pointers[j+bb]-1])
        
        bb = bb + len(aa)

    # Finally, fill the diagonal with volume coeffs.
    for i in range(0,len(vols)):
        matrix[i][i] = vols[i]

    return matrix

def get_sparse_matrix(infile):

    ext = infile.split('.')[-1].lower()

    if ext == 'stor':
        return get_matrix(infile)
    elif ext == 'lgi':
        pass
    elif ext == 'avs':
        pass
    else:
        raise ValueError('File extension \'%s\' not recognized' % ext)

def plot_sparse_matrix(m,outfile=None,viewfig=True,binary=True,cmap='binary'):

    if isinstance(cmap,str):
        cmap = plt.get_cmap(cmap)

    fig = plt.figure()
    extent = [0,np.shape(m)[1],0,np.shape(m)[0]]

    if not binary:
        im = plt.imshow(m,cmap=cmap,interpolation='none',extent=extent)
        plt.colorbar(im,ticks=np.linspace(0,np.max(m),num=5))
    else:
        plt.imshow(m > 2,cmap=cmap,interpolation='none',extent=extent)

    # Set labels & title
    plt.title('Sparse Matrix Coefficients')
    plt.xlabel("node i")
    plt.ylabel("node j")

    plt.axis(extent)

    # Save the figure?
    if outfile is not None:
        fig.savefig(outfile)
        #plt.imsave("/tmp/foo.png", data, format="png", cmap="hot")

    # View the figure?
    if viewfig:
        plt.show()

if (__name__ == "__main__"):

    # Check this out for 'true' SM vis: https://github.com/nschloe/betterspy

    applet_description = '''
    SPARSE MATRIX VISUALIZATION: 
    This script visualizes a sparse matrix for a given mesh - 
    this is a handy way to quickly view the i,j connections between nodes.

    Voronoi cell volumes are represented on the matrix diagonal, while the 
    off-diagonals represent: (Voronoi area_ij/interface_ij length). 

    In the current implementation, an NxN sparse matrix is converted to an
    NxN image - be mindful of very large meshes as to not run into memory
    issues.

    References:

    [1] https://ieeexplore.ieee.org/document/112361/
    [2] https://en.wikipedia.org/wiki/Sparse_matrix
    '''

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

