############################################################################################
# Plot mesh construct
############################################################################################

"""
	plot_mesh( s, mesh::Mesh; vertIndex=false, edgeIndex=false, faceIndex=false )

Plot edges and vertices of an input mesh

# Arguments
- `s::MSession`: MATLAB session
- `mesh::Mesh`: mesh object
"""
function plot_mesh( s, mesh::Mesh; vertIndex=false, edgeIndex=false, faceIndex=false )

	#Get mesh lines
	xArr = zeros( 2, length(mesh.mEdges) );
	yArr = zeros( size(xArr) );
	for ( ii,edge ) in enumerate( mesh.mEdges )
		xArr[1,ii] = edge.mVertices[1].mX[1];
		xArr[2,ii] = edge.mVertices[2].mX[1];

		yArr[1,ii] = edge.mVertices[1].mX[2];
		yArr[2,ii] = edge.mVertices[2].mX[2];
	end

	#Send to MATLAB and plot
	ml.put_variable( s, :xArr, xArr );
	ml.put_variable( s, :yArr, yArr );
	ml.eval_string( s, "plot( xArr, yArr, \"ko-\" )" ); 

	#Show vertex indices
	if vertIndex
		ml.eval_string( s, "hold on" );
		for vert in mesh.mVertices
			ml.put_variable( s, :cX, vert.mX[1] ); 
			ml.put_variable( s, :cY, vert.mX[2] ); 
			ml.put_variable( s, :label, string(vert.mIndex) );
			ml.eval_string( s, "text(cX,cY,label)" );
		end
		ml.eval_string( s, "hold off" );
	end

	#Show edge indices
	if edgeIndex
		ml.eval_string( s, "hold on" );
		for edge in mesh.mEdges
			cX, cY = 0.5 .* ( edge.mVertices[1].mX + edge.mVertices[2].mX );
			ml.put_variable( s, :cX, cX ); 
			ml.put_variable( s, :cY, cY ); 
			ml.put_variable( s, :label, string(edge.mIndex) );
			ml.eval_string( s, "text(cX,cY,label)" );
		end
		ml.eval_string( s, "hold off" );
	end

	#Show face indices
	if faceIndex
		ml.eval_string( s, "hold on" );
		cX = 0.0; cY = 0.0;
		for face in mesh.mFaces
			cX = 0.0; cY = 0.0;
			for vert in face.mVertices
				cX += vert.mX[1]; cY += vert.mX[2];
			end
			cX /= length(face.mVertices); cY /= length(face.mVertices);
			ml.put_variable( s, :cX, cX ); 
			ml.put_variable( s, :cY, cY ); 
			ml.put_variable( s, :label, string(face.mIndex) );
			ml.eval_string( s, "text(cX,cY,label)" );
		end
		ml.eval_string( s, "hold off" );
	end

end

############################################################################################
# Plot mesh construct
############################################################################################

"""
	plot_ref_el( s, refEl::ReferenceElement; vertIndex=false )

Plot edges and vertices of a reference element

# Arguments
- `s::MSession`: MATLAB session
- `mesh::Mesh`: mesh object
"""
function plot_ref_el( s, refEl::ReferenceElement; vertIndex=false )

	# plot edge
	if refEl.mShape == :tri
		xEdge = [ 0.0, 1.0, 0.0, 0.0 ];
		yEdge = [ 0.0, 0.0, 1.0, 0.0 ];
	elseif refEl.mShape == :quad
		xEdge = [ -1.0, 1.0, 1.0, -1.0, -1.0 ];
		yEdge = [ -1.0, -1.0, 1.0, 1.0, -1.0 ];
	end

	ml.put_variable( s, :xEdge, xEdge );
	ml.put_variable( s, :yEdge, yEdge );
	ml.eval_string( s, "plot( xEdge, yEdge, \"k-\" )" ); 

	# plot nodes
	x = refEl.mNodesX[:,1];
	y = refEl.mNodesX[:,2];

	ml.eval_string( s, "hold on" );
	ml.put_variable( s, :x, x );
	ml.put_variable( s, :y, y );
	ml.eval_string( s, "plot( x, y, \"ko\" )" ); 

	# Show vertex indices
	if vertIndex
		ml.eval_string( s, "hold on" );
		for i in eachindex(x)
			ml.put_variable( s, :cX, x[i] ); 
			ml.put_variable( s, :cY, y[i] ); 
			ml.put_variable( s, :label, string(i) );
			ml.eval_string( s, "text(cX,cY,label)" );
		end
		ml.eval_string( s, "hold off" );
	end

end