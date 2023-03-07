############################################################################################
# Overwrite mesh terminal I/O
############################################################################################
function Base.show( io::IO, mesh::Mesh )
	print( "Mesh has: ", length(mesh.mVertices), " vertices, ", length(mesh.mEdges), 
        " edges, ", length(mesh.mFaces), " faces." );
end

############################################################################################
# OBJ File I/O
############################################################################################

"""
    import_obj( file::String )

Read in obj file and save as mesh
"""
function import_obj( file::String )

	# Keep track of constructed verts and faces
	vIndex = 0; tIndex = 0;

	# Find how many vertices and elements
	open( file, "r" ) do io
		for ln in eachline( io )
			if length( ln ) > 0 && ln[1] == 'v' && isspace( ln[2] )
				# Case for vertices
				vIndex += 1;
			elseif length( ln ) > 0 && ln[1] == 'f'
				# Case for faces
				tIndex += 1;
			end
		end
	end

	# Storage for all vertices of this mesh
	mVertices = Vector{Vertex}( undef, vIndex ); vIndex = 0;
	mFaces = Vector{Face}( undef, tIndex ); tIndex = 0;

	# I/O loop
	open( file, "r" ) do io
		for ln in eachline( io )
			if length( ln ) > 0 && ln[1] == 'v' && isspace( ln[2] )
				# Case for vertices
				tokens = split( ln, (r"\s+") ); vIndex += 1;
				x = parse( Float64, tokens[2] ); 
                y = parse( Float64, tokens[3] ); 
                z = parse( Float64, tokens[4] );
				nV = Vertex( vIndex, [x, y] ); mVertices[vIndex] = nV;

			elseif length( ln ) > 0 && ln[1] == 'f'
				# Case for faces
				tokens = split( ln, (r"\s+") ); tIndex += 1; 
                nF = Face( tIndex, length(tokens)-1 ); mFaces[tIndex] = nF;

				for i = 2:length(tokens)
					temp = split( tokens[i], ("/") ); vert = parse( Int64, temp[1] );

					# Tell each face which vertices are part of it, and vice versa
					nF.mVertices[i-1] = mVertices[vert]; 
                    push!( mVertices[vert].mFaces, tIndex );
				end
			end
		end
	end

	# Deal with face neighbor lists
	adjFaces = Dict{Int64,Int64}();
	for cFace in mFaces
		empty!( adjFaces );
		# Find all possible adjacencies
		for fVert in cFace.mVertices
			for face in fVert.mFaces
				if haskey( adjFaces, face )
					adjFaces[face] += 1;
				else
					push!( adjFaces, face => 1 );
				end
			end
		end

		# Figure out which faces are really adjacent
		nIndex = 1;
		for ( f, count ) in adjFaces
			if !(f == cFace.mIndex) && count == 2
				cFace.mNeighbors[nIndex] = f; nIndex += 1;
			end
		end
	end

	# Store faces in half-edge order and generate edges
	mEdges = resolve_edges!( mFaces );

	# Return index of the mesh
	return Mesh( mVertices, mEdges, mFaces );
end

"""
    export_obj( mesh::Mesh, file::String )

Export mesh as obj file
"""
function export_obj( mesh::Mesh, file::String )

	meshVerts = mesh.mVertices; meshFaces = mesh.mFaces;

	# I/O loop
	open( file, "w" ) do io
		for vert in meshVerts
			write( io, string( "v ", vert.mX[1], "\t", vert.mX[2], "\t", 0.0, "\n" ) );
		end
		for face in meshFaces
			if length( face.mVertices ) > 0 && length( face.mEdges ) > 0
				i = 1;
				write( io, string("f ") );
				for v in face.mVertices
					if i == length( face.mVertices )
						write( io, string( v.mIndex ) );
					else
						write( io, string( v.mIndex, "\t" ) );
					end
					i += 1;
				end
				write( io, string("\n") );
			end
		end
	end

end