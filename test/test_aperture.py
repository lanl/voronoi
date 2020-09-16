from framework import SampleMeshes, Voronoi, compare_arrays


def test_aperture():
    v = Voronoi(
        nprocs=4, 
        infile=SampleMeshes.TRI.SIMPLE(),
    )
    v.run()
    gold = v.sparse_matrix()

    v = Voronoi(
        nprocs=4, 
        infile=SampleMeshes.TRI.SIMPLE(),
        aperture_file=
    )
    v.run()
    bronze = v.sparse_matrix()

    compare_arrays(bronze, gold, rtol=1e-7, atol=0)