testPath = @__DIR__;
aggPath = dirname(testPath);

include(aggPath * "/src/AgglomerationMultigrid2D.jl")
import .AgglomerationMultigrid2D as aggmg

mesh = aggmg.import_obj( testPath * "/std_plane.obj" );
mesh2 = aggmg.refine_quad( mesh );