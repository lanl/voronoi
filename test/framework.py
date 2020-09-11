import os
import sys
import json
import numpy as np

CONFIG_FILE = "config.json"


def init_parameters():
    with open(CONFIG_FILE, "r") as fp:
        params = json.load(fp)

    return params["voronoi_exe"]


VORONOI_EXE = init_parameters()


class SampleMeshes:
    class TRI:
        def SIMPLE():
            raise ValueError('test')

    class TET:
        SIMPLE = None


class Voronoi:
    def __init__(
        self,
        nprocs: int,
        infile: str,
        output: str = None,
        infile_is_lgi: bool = None,
        output_type: str = None,
        compute_median: bool = False,
        show_diagnostics: bool = False,
        extra_flags: str = "",
        voronoi_exe: str = VORONOI_EXE,
    ):
        """
        If `infile_is_lgi` is set to `None`, then it will parse the file
        extension, or default to `avs` if that fails.
        """

        if output is None:
            output = "test_output"

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
        self.mpiexec = None

        self.nprocs = nprocs
        self.infile = infile
        self.output = output
        self.infile_is_lgi = infile_is_lgi
        self.output_type = output_type
        self.compute_median = compute_median
        self.show_diagnostics = show_diagnostics
        self.extra_flags = extra_flags

    def run(self, silent: bool = False):
        """
        Runs VORONOI with the configured parameters.
        """

        mpi_prefix = f"{mpiexec} -np {self.nprocs}" if self.nprocs > 1 else ""
        control_volume = "median" if self.compute_median else "voronoi"

        if self.infile_is_lgi:
            infile = f"-lg {self.infile}"
        else:
            infile = f"-avs {self.infile}"

        diagnostics = "-d" if self.diagnostics else ""

        cmd = " ".join(
            mpi_prefix,
            self.exe,
            f"{infile}" f"-type {self.output_type}",
            f"-cv {control_volume}",
            f"-o {self.output}",
            f"{diagnostics}" f"{extra_flags}",
        )

        out = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)

        try:
            out = out.decode("ascii")
        except AttributeError:
            pass
        assert "error" not in out.lower(), "VORONOI threw an exception: \n%s\n" % out

    def __del__(self):
        os.remove(self.output)
