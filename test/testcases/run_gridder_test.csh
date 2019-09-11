#! /bin/csh
#
# Quick
foreach file ( \
input_2d_uniform_10.gridder \
input_2d_uniform_100.gridder \
input_2d_loglog_100.gridder \
)
#
# Not Quick
# foreach file ( input_2d*.gridder)
date
#
#   Run gridder to create 2D quad mesh
#
	echo INPUT FILE $file
	gridder < $file
	if(-e $file:r.inp) rm $file:r.inp
	mv grid.inp $file:r.inp
	if(-e input.grid) rm input.grid
	if(-e input.tmp) rm input.tmp
#
#   Run LaGriT to create 2D triangle mesh and geometric coefficient file
#
	if(-e mesh.inp) rm mesh.inp
	ln -s $file:r.inp mesh.inp
	lagrit < mesh_2d_to_geom_coef.lgi
	if(-e $file:r.stor) rm $file:r.stor
	mv mesh.stor $file:r.stor
	if(-e $file:r.tri.inp) rm $file:r.tri.inp
	mv mesh_tri.inp $file:r.tri.inp
	if(-e $file:r.outx3dgen) rm $file:r.outx3dgen
	mv outx3dgen $file:r.outx3dgen
	if(-e mesh.inp) rm mesh.inp
end

