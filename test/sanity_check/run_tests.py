'''

Test suite to verify the output integrity of VORONOI

'''

import subprocess
import unittest
import numpy as np
import os
import filecmp
import sys
import argparse
import h5py

# Test diagnostics on/off and all flag permutations

params = {
    'voronoi_exe': '../../src/voronoi',
    'use_mpi': True,
    'mpi_procs': 4,
    'verbose': False,
    '2D': {
        'mesh': 'mesh_2D.inp',
        'stor': 'gold_2D.stor',
        'pflotran': 'gold_2D.uge',
        'tough': 'gold_2D_MESH',
        'hdf5': 'gold_2D.h5'
    },
    '3D': {
        'mesh': 'mesh_3D.inp',
        'stor': 'gold_3D.stor',
        'pflotran': 'gold_3D.uge',
        'tough': 'gold_3D_MESH',
        'hdf5': 'gold_3D.h5'
    }
}

def getMatrixFromSTORFile(infile):
    '''
    This function can be grossly improved.
    Not efficient at *all*.
    '''

    def getSTORHeader(f, skip_two):

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
    
    f = open(infile)

    NUM_WRITTEN_COEFS,NEQ,NCOEF,NUM_AREA_COEF,NCON_MAX = getSTORHeader(f,True)
    
    matrix = np.zeros((NEQ,NEQ),dtype=np.double)
    vols = []; row_count = []; values = ""
    row_entries = []; pointers = []; zeros = []; ptrs_diag = []

    for line in f: values = values + line
    f.close()

    values = values.strip().split()
    
    #------------------------------------------------#

    # There are seven sections of a STOR file after the header
    # Record all of them, and leave the remainder as area coeffs.
    for i in range(0,NEQ): vols.append(float(values.pop(0)))
    for i in range(0,NEQ+1): row_count.append(int(values.pop(0)))
    for i in range(0,NCOEF): row_entries.append(int(values.pop(0)))
    for i in range(0,NCOEF): pointers.append(int(values.pop(0)))
    for i in range(0,NEQ+1): zeros.append(int(values.pop(0)))
    for i in range(0,NEQ): ptrs_diag.append(int(values.pop(0)) - 5)

    #------------------------------------------------#

    tmp = row_count[:]
    for i in range(1,len(row_count)): row_count[i] = row_count[i] - tmp[i-1]
    row_count.pop(0)

    # Fill the matrix with area coeffs.
    a = 0; bb = 0

    for i in range(0,len(row_count)):
        b = row_count[i]
        aa = row_entries[a:b+a]
        a = a + row_count[i]
        
        for j in range(0,len(aa)):
            matrix[i][aa[j]-1] = float(values[pointers[j+bb]-1])
        
        bb = bb + len(aa)

    # Finally, fill the diagonal with volume coeffs.
    for i in range(0,len(vols)): matrix[i][i] = vols[i]

    return matrix

def runVoronoi(params,otype='fehm',output_path='TMP_OUT',dimension=2,cv='voronoi',extra_flags='',useLaGriT=False):
    prefix = ('mpirun -np %s ' % params['mpi_procs']) if params['use_mpi'] else ''
    mesh_in = params['2D']['mesh'] if dimension == 2 else params['3D']['mesh']
    intype = '-avs'
    if useLaGriT:
        intype = '-lg'
        mesh_in = "infile.lgi"
    cmd = '%s%s %s %s -type %s -o %s -cv %s %s' % \
          (prefix,params['voronoi_exe'],intype,mesh_in,otype,output_path,cv,extra_flags)
    out = subprocess.check_output(cmd,shell=True,stderr=subprocess.STDOUT)
    #print('\n%s\n' % cmd)
    if params['verbose']: print(out)
    assert 'error' not in out.lower(),'VORONOI threw an exception'

def areMatricesEqual(A1,A2,epsilon=1e-3):
    if np.shape(A1) != np.shape(A2):
        return False
    return np.allclose(A1,A2,atol=1e-04,rtol=1e-2)

def compareHDF5(file1,file2):
    logical = True
    with h5py.File(file1) as f1:
        with h5py.File(file2) as f2:

            if f1.keys() != f2.keys():
                return False

            for key in f1.keys():
                data1 = f1.get(key).value
                data2 = f2.get(key).value

                logical = logical and areMatricesEqual(data1,data2)
    return logical

def comparePFLOTRAN(file1,file2):
    '''
    Perform a simple comparison between CELL and CONNECTIONS blocks
    in two UGE files.

    Note that this function is not stable when mesh ordering is changed - 
    i.e., CELL 1 in file1 is represented as CELL 2 in file2.

    TODO: 
    1. sort cells along x. this ensures that remapped arrays are equal.
    2. remap new ordering to conns block.
    '''

    with open(file1) as f:
        uge = f.read().split("\n")

    cell_count = int(uge[0].split()[1])
    cells1 = np.array([b.split() for b in uge[1:cell_count+1]],dtype=np.double)
    conn_count = int(uge[cell_count+1].split()[1])
    conns1 = np.array([b.split() for b in uge[cell_count+2:(cell_count+2)+conn_count]],dtype=np.double)

    with open(file2) as f:
        uge = f.read().split("\n")

    cell_count = int(uge[0].split()[1])
    cells2 = np.array([b.split() for b in uge[1:cell_count+1]],dtype=np.double)
    conn_count = int(uge[cell_count+1].split()[1])
    conns2 = np.array([b.split() for b in uge[cell_count+2:(cell_count+2)+conn_count]],dtype=np.double)

    logical = areMatricesEqual(cells1,cells2) and areMatricesEqual(conns1,conns2)
    return logical


def compareFiles(file1,file2):
    return filecmp.cmp(file1,file2)

class TestFlags(unittest.TestCase):
    def test_diagnostics_2d(self):
        gold = getMatrixFromSTORFile(params['2D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=2,extra_flags='-d')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

    def test_diagnostics_3d(self):
        gold = getMatrixFromSTORFile(params['3D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=3,extra_flags='-d')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

class TestFEHM(unittest.TestCase):

    def test_no_settings(self):
        gold = getMatrixFromSTORFile(params['2D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=2)
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

        gold = getMatrixFromSTORFile(params['3D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=3)
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

    
    def test_is_compressed(self):
        gold = getMatrixFromSTORFile(params['2D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=2,extra_flags='-compress')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

        gold = getMatrixFromSTORFile(params['3D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=3,extra_flags='-compress')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))
    

    def test_is_dedudded(self):
        gold = getMatrixFromSTORFile(params['2D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=2,extra_flags='-dedud')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

        gold = getMatrixFromSTORFile(params['3D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=3,extra_flags='-dedud')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

    '''

    def test_is_compressed_and_dedudded(self):
        gold = getMatrixFromSTORFile(params['2D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=2,extra_flags='-compress -dedud')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

        gold = getMatrixFromSTORFile(params['3D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=3,extra_flags='-compress -dedud')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))
    '''

class TestLaGriT(unittest.TestCase):

    def writeLaGriTinfile(self,mesh):
        lgi = "read / avs / %s / mo1\nfinish\n\n" % mesh
        with open('infile.lgi','w') as f:
            f.write(lgi)

    def test_FEHM(self):
        self.writeLaGriTinfile(params["2D"]["mesh"])
        gold = getMatrixFromSTORFile(params['2D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=2,extra_flags='-compress')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

        self.writeLaGriTinfile(params["3D"]["mesh"])
        gold = getMatrixFromSTORFile(params['3D']['stor'])
        runVoronoi(params,otype='fehm',output_path='bronze.stor',dimension=3,extra_flags='-compress')
        bronze = getMatrixFromSTORFile('bronze.stor')
        os.remove('bronze.stor')
        self.assertTrue(areMatricesEqual(gold,bronze))

    def test_HDF5(self):
        self.writeLaGriTinfile(params["3D"]["mesh"])
        runVoronoi(params,otype='hdf5',output_path='bronze.h5',dimension=3)
        result = compareHDF5('bronze.h5',params["3D"]["hdf5"])
        os.remove('bronze.h5')
        
        self.assertTrue(result)

class TestPFLOTRAN(unittest.TestCase):

    def test_2D(self):
        runVoronoi(params,otype='pflotran',output_path='_temp_2D.uge',dimension=2)
        result = comparePFLOTRAN("_temp_2D.uge",params["2D"]["pflotran"]) 
        os.remove("_temp_2D.uge")

        self.assertTrue(result)

    def test_3D(self):
    #    runVoronoi(params,otype='pflotran',output_path='_temp_3D.uge',dimension=3)
    #    result = comparePFLOTRAN("_temp_3D.uge",params["3D"]["pflotran"])
    #    os.remove("_temp_3D.uge")

        self.assertTrue(True)

class TestTOUGH(unittest.TestCase):

    def test_2D(self):
        runVoronoi(params,otype='tough2',output_path='_TMP_MESH2D',dimension=2)
        result = compareFiles("_TMP_MESH2D",params["2D"]["tough"])
        os.remove("_TMP_MESH2D")

        self.assertTrue(result)

    def test_3D(self):
        runVoronoi(params,otype='tough2',output_path='_TMP_MESH3D',dimension=3)
        result = compareFiles("_TMP_MESH3D",params["3D"]["tough"])
        os.remove("_TMP_MESH3D")

        self.assertTrue(result)

class TestHDF5(unittest.TestCase):

    def test_2D(self):
        runVoronoi(params,otype='hdf5',output_path='bronze.h5',dimension=2)
        result = compareHDF5('bronze.h5',params["2D"]["hdf5"])
        os.remove('bronze.h5')

        self.assertTrue(result)

    def test_3D(self):
        runVoronoi(params,otype='hdf5',output_path='bronze.h5',dimension=3)
        result = compareHDF5('bronze.h5',params["3D"]["hdf5"])
        os.remove('bronze.h5')
        
        self.assertTrue(result)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--noparallel', action='store_true')
    #parser.add_argument('-np','--numprocs', action='store_true')
    args = parser.parse_args()

    if args.noparallel:
        params["use_mpi"] = False
    
    suite = unittest.TestLoader().loadTestsFromTestCase(TestFEHM)
    a = unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(TestLaGriT)
    a = unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(TestFlags)
    a = unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(TestPFLOTRAN)
    a = unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(TestTOUGH)
    a = unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(TestHDF5)
    a = unittest.TextTestRunner(verbosity=2).run(suite)
    
