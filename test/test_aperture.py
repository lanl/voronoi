from framework import SampleMeshes, Voronoi


def test_aperture():
    v = Voronoi(
        nprocs=4, 
        infile=SampleMeshes.TRI.SIMPLE(),
        )
    #v.run()
    #gold = v.sparse_matrix()

    #v = Voronoi()
    #v.run()
    #bronze = v.sparse_matrix()

    #assert gold == bronze
