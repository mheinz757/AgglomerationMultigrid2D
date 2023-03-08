import MATLAB as ml

############################################################################################
#	----------------------------------------------------------------------------------------
#	HOW TO USE
#	----------------------------------------------------------------------------------------
#		1) Create a new session
#			e.g. s = create_session()
#		2) For standard MATLAB plotting functions, call them as per usual with s as added 
#		   1st argument
#			e.g. clf( s ), hold_on( s ), xlim( s, [0,1] )
#		3) For other more specific commands look at the "Plot Julia constructs section"
#		4) Close the session once you are done
#			e.g. close_session( s )
############################################################################################

############################################################################################
# Start/stop session
############################################################################################
function create_session()
	return ml.MSession( 0 );
end

function get_default_msession()
	return ml.get_default_msession();
end

function close_session( s )
	ml.close( s );
end

############################################################################################
# Basic tools
############################################################################################
function clf( s )
	ml.eval_string( s, "clf" );
end

function xlim( s, arr )
	ml.eval_string( s, string( "xlim(", arr, ")" ) );
end

function ylim( s, arr )
	ml.eval_string( s, string( "ylim(", arr, ")" ) );
end

function xlabel( s, label )
	ml.eval_string( s, string( "xlabel(\" ", label, "\" )" ) );
end

function ylabel( s, label )
	ml.eval_string( s, string( "ylabel(\" ", label, "\" )" ) );
end

function axis_equal( s )
	ml.eval_string( s, "axis equal" );
end

function hold_on( s )
	ml.eval_string( s, "hold on" );
end

function hold_off( s )
	ml.eval_string( s, "hold off" );
end

function colorbar( s )
	ml.eval_string( s, "colorbar()" );
end

function plot( s, x, y; tags="o" )
	ml.put_variable( s, :xT, x );
	ml.put_variable( s, :yT, y );
	ml.eval_string( s, string( "plot( xT,yT, \" ", tags, " \" )" ) );
end

function plot3( s, x, y, z; tags="o" )
	ml.put_variable( s, :xT, x );
	ml.put_variable( s, :yT, y );
	ml.put_variable( s, :zT, z );
	ml.eval_string( s, string( "plot3( xT,yT,zT, \" ", tags, " \" )" ) );
end

function scatter( s, x, y, z )
	ml.put_variable( s, :xT, x );
	ml.put_variable( s, :yT, y );
	ml.put_variable( s, :zT, z );
	ml.eval_string( s, "scatter( xT,yT,20,zT )" );
end

############################################################################################
# More complicated plotting
############################################################################################

# ###	Plot computational nodes of a mesh
# #	@s: MATLAB session
# #	@elemMesh: elemMesh object
# function plot_nodes( s, elemMesh::ElemMesh; nodeIndex=false )

# 	#Send to MATLAB and plot
# 	ml.put_variable( s, :xPts, elemMesh.mX[:,1] );
# 	ml.put_variable( s, :yPts, elemMesh.mX[:,2] );
# 	ml.eval_string( s, "plot( xPts, yPts, \"bo\" )" );  

# 	#Plot node indices
# 	if nodeIndex
# 		labels = Vector{String}( undef, elemMesh.nNodes2D );
# 		for ii = 1:elemMesh.nNodes2D
# 			labels[ii] = string( ii );
# 		end

# 		hold_on( s );
# 		ml.put_variable( s, :labels, labels ); 
# 		ml.eval_string( s, "text( xPts, yPts, labels )" ); 
# 		hold_off( s );
# 	end

# end

"""
	plot_velocity_field( s, x::Vector{Float64}, y::Vector{Float64}, u::Vector{Float64} )

Plot velocity field

# Arguments
- `s::MSession`: MATLAB session
- `x::Vector{Float64}`: x points
- `y::Vector{Float64}`: y points
- `u::Vector{Float64}`: velocity vector field
"""
function plot_velocity_field( s, x::Vector{Float64}, y::Vector{Float64}, 
	u::Vector{Float64} )

	#Send to MATLAB and plot
	ml.put_variable( s, :xPts, x );
	ml.put_variable( s, :yPts, y );

	ml.put_variable( s, :vX, u[1:length(x)] );
	ml.put_variable( s, :vY, u[length(x)+1:length(x)*2] );

	#Plot each of velocity
	ml.eval_string( s, "quiver( xPts, yPts, vX, vY, 1.0 )" );

end