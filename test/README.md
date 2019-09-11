## VORONOI Test Suite ##

Within this folder are two subfolders, containing validation tests for a VORONOI build.

### `sanity_check` ###

Within this folder is a Python script named `run_test.py`. This script:

1. Runs VORONOI against multiple permutations of options
2. Tests the written files against known 'gold standards'
3. Returns a non-zero exit code if the written coefficients do not match the 'gold' file (within an epsilon)

That is, this test verifies the integrity of written FEHM, PFLOTRAN, TOUGH, etc. files under variable settings.

### `edge_cases` ###

Within this folder is a c-shell script named `run_tests.scr`. This script:

1. Builds a number of very large 2D and 3D meshes with LaGriT
2. Reads the mesh in with VORONOI and writes to each output type
3. Returns a non-zero exit code if the test failed

The purpose of this test is to ensure that VORONOI can handle large and unconventional mesh geometeries.

### Contributing Tests ###

To contribute new tests, please submit a pull request with your test enclosed in a subfolder. In the pull request, explain the purpose of your test and instructions on running. To keep a compact repo, please limit the file size of included files.