"""

The VORONOI program only accepts meshes in a select few formats.
See [1] for a list of supported filetypes.

This script is intended to generate an AVS-UCD mesh for use in VORONOI,
from a collection of nodes and connectivity.

[1]: http://lagrit.lanl.gov/docs/commands/READ.html

"""

import argparse
import numpy as np
import sys

def read_nodes(infile,delimiter=' '):
    return np.loadtxt(infile,delimiter=delimiter,dtype=np.double)

def read_elements(infile,delimiter=' '):
    return np.loadtxt(infile,delimiter=delimiter,dtype=np.int)

def read_attributes(infile,delimiter=' '):
    pass

def write_avs(nodes,outfile,elements=None,ttype=None):
    
    if ttype is None:
        ttype = 'tri' if len(mo.elements[0]) == 3 else 'tet'
    else:
        assert ttype in ['tri','tet','quad','hex'],'Unknown element type'

    count_nodes = len(nodes)
    count_elems = len(elements) if elements is not None else 0
    count_natts = len(mo.natts[0]) if natts is not None else 0
    count_catts = len(mo.catts[0]) if catts is not None else 0

    assert count_nodes > 0,'nodes array cannot be empty'

    with open(outfile,'w') as f:
        # Write out header
        f.write("{:11d}{:11d}{:11d}{:11d}{:11d}\n".format(count_nodes,count_elems,count_natts,count_catts,0))
        
        # Write out nodes
        for i in range(0,count_nodes): f.write("{:03d}  {:010E}  {:010E}  {:010E}\n".format(i+1,mo.nodes[i][0],mo.nodes[i][1],mo.nodes[i][2]))
        
        # Write out elements
        if count_elems != 0:

            #f.write('\n'.join(["{:03d}        1   {}  {}".format(i+1,ttype,'  '.join(mo.elements[i])) for i in range(count_elems)]))

            if ttype == 'tet':
                for i in range(0,count_elems): f.write("{:03d}        1   tet  {}  {}  {}  {}\n".format(i+1,mo.elements[i][0],mo.elements[i][1],mo.elements[i][2],mo.elements[i][3]))
            elif ttype == 'tri':
                for i in range(0,count_elems): f.write("{:03d}        1   tri  {}  {}  {}\n".format(i+1,mo.elements[i][0],mo.elements[i][1],mo.elements[i][2]))
        
        # Write out node attributes
        if count_natts != 0:
            count = len(mo.natts[0])

            # Write out attribute header (data types needs to be fixed)
            f.write("{:05d}".format(count) + "  1"*count + "\n")
            names = ["imt{}, {}\n".format(i+1, "real") for i in range(0,count)]
            for i in range(0,len(names)): f.write(names[i])

            # Write out att data
            for i in range(0,len(mo.natts)):
                formatted_atts = '  '.join(["{0:010E}".format(mo.natts[i][j]) for j in range(0,count)])
                f.write('{0:010d}  '.format(i+1) + formatted_atts + "\n")

        # Write out element attributes
        if count_catts != 0:
            count = len(mo.catts[0])

            # Write out attribute header (data types needs to be fixed)
            f.write("{:05d}".format(count) + "  1"*count + "\n")
            names = ["isn{}, {}\n".format(i+1, "real") for i in range(0,len(count))]
            for i in range(0,len(names)): f.write(names[i])

            # Write out att data
            for i in range(0,len(mo.catts)):
                formatted_atts = '  '.join(["{0:010E}".format(mo.catts[i][j]) for j in range(0,count)])
                f.write('{0:010d}  '.format(i+1) + formatted_atts + "\n")

if (__name__ == "__main__"):

    applet_description = '''
    This script converts a file of node (x,y,z) values and element (i,j,k,...) indices into
    an AVS-UCD mesh.

    Node input file should be in the form:

        x0 y0 z0
        x1 y1 z1
        ...
        xN yN zN

    Element input file should be in the form:

        i1 j1 k1 ...
        i2 j2 k2 ...
        ...
        iN jN kN ...

    A delimiter between entries may be specified with the -d argument. Defaults to space (` `).

    It is recommended that the element file reference the node list with 1-based indexing; that is, 
    the element input file should reference the node `x0 y0 z0` as `1`.
    If you use zero- or N-based indexing, use the --index flag to indicate this.

    This script will automatically assume that an element file with 3 integers on a line refers to
    a triangle element, and 4 integers refers to tetrahedrons. If you wish to manually specify what the 
    element type is (i.e., quad) then use --type ['tet','tri','quad','hex']. Note that only one element
    type per file is allowed - mixing of element types in a single mesh is not supported in Voronoi.

    If you only have a list of nodes, this script will still write out a file - but no element connectivity
    will be defined. You can import the mesh into LaGriT and triangulate it, or use the SciPy Delaunay function.

    If you wish to view the created mesh, it is recommended that you use ParaView: https://www.paraview.org
    '''

    parser = argparse.ArgumentParser(description=applet_description,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o','--outfile',help='Filepath to save mesh to (defaults to \'mesh.avs\')',type=str,default='mesh.avs')
    parser.add_argument('-n','--nodes',help='Filepath to node input file',type=str)
    parser.add_argument('-e','--elements',help='Filepath to elements input file',type=str)
    parser.add_argument('-t','--type',help='Type of elements (tri,tet,quad,...)',type=str)
    parser.add_argument('-d','--delimiter',help='Delimiter between entries in node and elements files. Defaults to space.',type=str,default=' ')
    parser.add_argument('-na','--nodeatts',help='Filepath to node attributes',type=str)
    parser.add_argument('-ca','--cellatts',help='Filepath to cell attributes',type=str)
    parser.add_argument('-i','--index',help='First node reference integer in elements file',type=int,default=1)
    args = parser.parse_args()

    if args.nodes is None:
        parser.print_help()
        print("\nERROR: Node input file is required")
        sys.exit()

    if args.elements is None:
        print("WARNING: No element file provided. Writing only nodes...")

    nodes = read_nodes(args.nodes,delimiter=delimiter)
    elements = read_elements(args.elements,delimiter=delimiter) if args.elements is not None else None
    natts = read_node_attributes(args.nodeatts,delimiter=delimiter) if args.nodeatts is not None else None
    catts = read_node_attributes(args.cellatts,delimiter=delimiter) if args.cellatts is not None else None

    write_avs(nodes,args.outfile,elements=elements,ttype=args.type)

    print('Mesh successfully written to \'%s\'' % args.outfile)
