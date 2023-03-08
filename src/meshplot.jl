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