############################################################################################
# 2D Refinement
############################################################################################

"""
    refine_tri( mesh::Mesh )

Refines a triangular mesh, returning a new triangular mesh.
"""
function refine_tri( mesh::Mesh )

	#Constructs of mesh to be refined
	oldFaces = mesh.mFaces; oldEdges = mesh.mEdges; oldVertices = mesh.mVertices;
	nVerts = length( oldVertices ); nEdges = length( oldEdges ); 
    nFaces = length( oldFaces );
	
	#Constructs for new mesh
	vIndex = 0; fIndex = 0;
	nVertices = Vector{Vertex}( undef, nVerts+nEdges ); 
    nFaces = Vector{Face}( undef, 4*nFaces );

	#Vertex points
	for vertex in oldVertices
		vIndex += 1; nVertices[vIndex] = Vertex( vIndex, vertex.mX );
	end

	#Edge points
	for edge in oldEdges
		vIndex += 1; 
        nVertices[vIndex] = Vertex( vIndex, 
            0.5 * (edge.mVertices[1].mX + edge.mVertices[2].mX) );
	end

	#Iterate over all old faces to form new faces
	vEdges = [0,0];
	for face in oldFaces

		#Faces involving vertex points
		for ii = 1:3
			fIndex += 1; 
            nF = Face( fIndex, 3 ); nFaces[fIndex] = nF; 
            foundEdges = 1; 
            vertex = face.mVertices[ii];

			for edge in face.mEdges
				if foundEdges < 3
					if edge.mVertices[1].mIndex == vertex.mIndex || edge.mVertices[2].mIndex == vertex.mIndex
						vEdges[foundEdges] = edge.mIndex; foundEdges += 1;
					end
				end
			end

			nF.mVertices[1] = nVertices[ vertex.mIndex ];
			push!( nF.mVertices[1].mFaces, fIndex );

			nF.mVertices[2] = nVertices[ nVerts+vEdges[1] ];
			push!( nF.mVertices[2].mFaces, fIndex );

			nF.mVertices[3] = nVertices[ nVerts+vEdges[2] ];
			push!( nF.mVertices[3].mFaces, fIndex );
		end

		#Face involving only edge points
		fIndex += 1; nF = Face( fIndex, 3 ); nFaces[fIndex] = nF; 

		nF.mVertices[1] = nVertices[ nVerts + face.mEdges[1].mIndex ];
		push!( nF.mVertices[1].mFaces, fIndex );

		nF.mVertices[2] = nVertices[ nVerts + face.mEdges[2].mIndex ];
		push!( nF.mVertices[2].mFaces, fIndex );

		nF.mVertices[3] = nVertices[ nVerts + face.mEdges[3].mIndex ];
		push!( nF.mVertices[3].mFaces, fIndex );

	end

	#Deal with face neighbor lists
	adjFaces = Dict{Int64,Int64}();
	for cFace in nFaces
		empty!( adjFaces );
		#Find all possible adjacencies
		for fVert in cFace.mVertices
			for face in fVert.mFaces
				if haskey( adjFaces, face )
					adjFaces[face] += 1;
				else
					push!( adjFaces, face => 1 );
				end
			end
		end

		#Figure out which faces are really adjacent
		nIndex = 1;
		for ( f, count ) in adjFaces
			if !(f == cFace.mIndex) && count == 2
				cFace.mNeighbors[nIndex] = f; nIndex += 1;
			end
		end
	end

	#Store faces in half-edge order and generate edges
	nEdges = resolve_edges!( nFaces );

	#Return index of the mesh
	return Mesh( nVertices, nEdges, nFaces );

end

"""
    refine_quad( mesh::Mesh )

Refines a quadrilateral mesh using "Catmull-Clark" like rules, returning a new 
quadrilateral mesh.
"""
function refine_quad( mesh::Mesh )
	
	#Constructs of mesh to be refined
	oldFaces = mesh.mFaces; oldEdges = mesh.mEdges; oldVertices = mesh.mVertices;
	nVerts = length( oldVertices ); nEdges = length( oldEdges ); 
    nFaces = length( oldFaces );
	
	#Constructs for new mesh
	vIndex = 0; eIndex = 0; fIndex = 0;
	nVertices = Vector{Vertex}( undef, nVerts+nEdges+nFaces); 
    nFaces = Vector{Face}( undef, 4*nFaces );

	#Vertex points
	for vertex in oldVertices
		vIndex += 1; nVertices[vIndex] = Vertex( vIndex, vertex.mX );
	end

	#Edge points
	for edge in oldEdges
		vIndex += 1; nVertices[vIndex] = Vertex( vIndex, 
            0.5 * (edge.mVertices[1].mX + edge.mVertices[2].mX) );
	end

	#Face points
	vEdges = zeros( Int64, 2 );
	for face in oldFaces
		vIndex += 1; nVertices[vIndex] = Vertex( vIndex, 
            0.25 * (face.mVertices[1].mX + face.mVertices[2].mX + face.mVertices[3].mX + 
            face.mVertices[4].mX) );

		#Create new faces
		for vertex in face.mVertices
			foundEdges = 1;
			for edge in face.mEdges
				if foundEdges < 3
					if edge.mVertices[1].mIndex == vertex.mIndex || edge.mVertices[2].mIndex == vertex.mIndex
						vEdges[foundEdges] = edge.mIndex; foundEdges += 1;
					end
				end
			end

			#Sanity check
			if foundEdges != 3
				error( "In refine_quad, each point of each face should have found 2 edges 
                    of face it is adjacent to. ")
			end

			fIndex += 1;  nF = Face( fIndex, 4 ); nFaces[fIndex] = nF;
			nF.mVertices[1] = nVertices[ vertex.mIndex ];
			push!( nF.mVertices[1].mFaces, fIndex );
			nF.mVertices[2] = nVertices[ nVerts+vEdges[1] ];
			push!( nF.mVertices[2].mFaces, fIndex );
			nF.mVertices[3] = nVertices[ vIndex ];
			push!( nF.mVertices[3].mFaces, fIndex );
			nF.mVertices[4] = nVertices[ nVerts+vEdges[2] ];
			push!( nF.mVertices[4].mFaces, fIndex );

		end
	end

	#Deal with face neighbor lists
	adjFaces = Dict{Int64,Int64}();
	for cFace in nFaces
		empty!( adjFaces );
		#Find all possible adjacencies
		for fVert in cFace.mVertices
			for face in fVert.mFaces
				if haskey( adjFaces, face )
					adjFaces[face] += 1;
				else
					push!( adjFaces, face => 1 );
				end
			end
		end

		#Figure out which faces are really adjacent
		nIndex = 1;
		for ( f, count ) in adjFaces
			if !(f == cFace.mIndex) && count == 2
				cFace.mNeighbors[nIndex] = f; nIndex += 1;
			end
		end
	end

	#Store faces in half-edge order and generate edges
	nEdges = resolve_edges!( nFaces );

	#Return index of the mesh
	return Mesh( nVertices, nEdges, nFaces );

end