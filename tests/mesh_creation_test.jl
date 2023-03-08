testPath = @__DIR__;
aggPath = dirname(testPath);

include(aggPath * "/src/AgglomerationMultigrid2D.jl")
import .AgglomerationMultigrid2D as aggmg

include(testPath * "/plot_matlab.jl")

mesh = aggmg.import_obj( testPath * "/std_plane.obj" );
mesh2 = aggmg.refine_quad( mesh );

s = get_default_msession();
plot_mesh( s, mesh; vertIndex=true, edgeIndex=true, faceIndex=true );
plot_mesh( s, mesh2; vertIndex=true, edgeIndex=true, faceIndex=true );