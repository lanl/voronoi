import os
import sys
import json
import numpy as np
import subprocess

CONFIG_FILE = "config.json"


def init_parameters():
    with open(CONFIG_FILE, "r") as fp:
        params = json.load(fp)

    return params["voronoi_exe"], params["mpiexec_exe"]


VORONOI_EXE, MPIEXEC_EXE = init_parameters()
DATA_DIR = "./data"


class SampleMeshes:
    class TRI:
        def SIMPLE():
            return os.path.join(DATA_DIR, 'mesh_2D.inp')

    class TET:
        def SIMPLE():
            return os.path.join(DATA_DIR, 'mesh_3D.inp')


class Voronoi:
    def __init__(
        self,
        nprocs: int,
        infile: str,
        output: str = None,
        voronoi_exe: str = VORONOI_EXE,
        mpiexec_exe: str = MPIEXEC_EXE,
        infile_is_lgi: bool = None,
        output_type: str = None,
        compute_median: bool = False,
        show_diagnostics: bool = False,
        aperture_file: str = None,
    ):
        """
        If `infile_is_lgi` is set to `None`, then it will parse the file
        extension, or default to `avs` if that fails.
        """

        if output is None:
            output = "test_output.stor"
            output_type = "fehm"

        if infile_is_lgi is None:
            ext = infile.split(".")[-1].strip().lower()
            if not ext:
                infile_is_lgi = False
            elif ext in ["avs", "inp"]:
                infile_is_lgi = False
            elif ext in ["lgi", "lg"]:
                infile_is_lgi = True
            else:
                infile_is_lgi = False

        if output_type is None:
            ext = infile.split(".")[-1].strip().lower()
            if not ext:
                output_type = "fehm"
            elif ext in ["stor"]:
                output_type = "fehm"
            elif ext in ["uge"]:
                output_type = "pflotran"
            elif ext in ["h5", "hdf5"]:
                output_type = "hdf5"
            else:
                output_type = "tough2"

        self.exe = voronoi_exe
        self.mpiexec = mpiexec_exe

        self.nprocs = nprocs
        self.infile = infile
        self.output = output
        self.infile_is_lgi = infile_is_lgi
        self.output_type = output_type
        self.compute_median = compute_median
        self.show_diagnostics = show_diagnostics
        self.aperture_file = aperture_file

    def run(self, silent: bool = False, extra_flags: str = None):
        """
        Runs VORONOI with the configured parameters.
        """

        mpi_prefix = f"{self.mpiexec} -np {self.nprocs}" if self.nprocs > 1 else ""
        control_volume = "median" if self.compute_median else "voronoi"

        if self.infile_is_lgi:
            infile = f"-lg {self.infile}"
        else:
            infile = f"-avs {self.infile}"

        diagnostics = "-d" if self.show_diagnostics else ""

        if extra_flags is None:
            extra_flags = ""

        if self.aperture_file is not None:
            set_aperture_file = f"-aperture_file {self.aperture_file}"
        else:
            set_aperture_file = ""

        cmd = " ".join([
            mpi_prefix,
            self.exe,
            f"{infile}",
            f"-type {self.output_type}",
            f"-cv {control_volume}",
            f"-o {self.output}",
            f"{set_aperture_file}",
            f"{diagnostics}" f"{extra_flags}",
        ])

        out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

        try:
            out = out.decode("ascii")
        except AttributeError:
            pass
        assert "error" not in out.lower(), "VORONOI threw an exception: \n%s\n" % out

    def sparse_matrix(self):
        '''
        Returns the sparse matrix.
        '''
        if self.output_type != "fehm":
            raise AssertionError("Must be run with `-type fehm`")

        skip_two = True

        with open(self.output, 'r') as f:
    
            if skip_two:
                f.readline() # fehmstor ascir8i4 LaGriT Sparse Matrix Voronoi Coefficients
                f.readline() # Mon Oct 23 17:19:23 2017
    
            header = [int(x) for x in f.readline().strip().split()]
            NUM_WRITTEN_COEFS = header[0]
            NEQ = header[1]
            NCOEF_NEQ_1 = header[2]
            NUM_AREA_COEF = header[3]
    
            try:
                NCON_MAX = int(header[4])
            except IndexError:
                NCON_MAX = None
            
            NCOEF = NCOEF_NEQ_1 - NEQ - 1
    
            matrix = np.zeros((NEQ,NEQ),dtype=np.double)
            vols = []
            row_count = []
            values = ""#np.array((0,0))#[]
            row_entries = []
            pointers = []
            zeros = []
            ptrs_diag = []
    
            for line in f:
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

    def __del__(self):
        os.remove(self.output)

def compare_arrays(bronze, gold, rtol=1e-07, atol=0):
    np.testing.assert_allclose(bronze, gold, rtol=rtol, atol=atol)

