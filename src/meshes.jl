"""
    Vertex

Stores information about a vertex on a 2-dimensional grid

# Fields
- `mIndex::Int64`
- `mX::Vector{Float64}`: X,Y position
- `mEdgeIDs::Vector{Vector{Int64}}`: IDs of incident edges
- `mEdges::Vector{Int64}`: vector of indices of incident edges
- `mFaces::Vector{Int64}`: vector of indices of faces adjacent to vertex
"""
struct Vertex

    mIndex::Int64					#Index
	mX::Vector{Float64}				#X,Y position
	mEdgeIDs::Vector{Vector{Int64}}	#IDs of incident edges
	mEdges::Vector{Int64}			#Indices of incidient edges
	mFaces::Vector{Int64}			#Indices of faces adjacent to vertex

	Vertex( mIndex, mX, mEdgeIDs=Vector{Vector{Int64}}(), mEdges=Vector{Int64}(), 
        mFaces=Vector{Int64}() ) = new( mIndex, mX, mEdgeIDs, mEdges, mFaces );
end

"""
    Edge

Stores information about an edge on a 2-dimensional grid

# Fields
- `mIndex::Int64`
- `mVertices::Vector{Float64}`: pointer to vertices at ends of edge
- `mFaces::Vector{Int64}`: vector of indices of faces adjacent to vertex
"""

struct Edge
	mIndex::Int64
	mVertices::Vector{Vertex}
	mFaces::Vector{Int64}

	Edge( mIndex, mVertices, mFaces=zeros(Int64,2) ) = new( mIndex, mVertices, mFaces );
end

"""
    Face

Stores information about a face on a 2-dimensional grid

# Fields
- `mIndex::Int64`
- `mVertices::Vector{Vertex}`: vector of pointers to vertices of the face
- `mEdges::Vector{Edge}`: vector of pointers to edges of the face
- `mNeighbors::Vector{Int64}`: vector of indices of neighbor faces
"""
struct Face
    mIndex::Int64
    mVertices::Vector{Vertex}
    mEdges::Vector{Edge}
    mNeighbors::Vector{Int64}

    Face( mIndex, mVertices, mEdges, mNeighbors ) = new( mIndex, mVertices, mEdges, 
        mNeighbors );
    Face( mIndex, nV::Int64 ) = new( mIndex, Vector{Vertex}( undef, nV ), 
        Vector{Edge}( undef, nV ), zeros(Int64, nV) );
end

"""
    Mesh

Stores information about a 2-dimensional mesh

# Fields
- `mVertices::Vector{Vertex}`: vector of pointers to vertices of the mesh
- `mEdges::Vector{Edge}`: vector of pointers to edges of the mesh
- `mFaces::Vector{Face}`: vector of pointers to faces of the mesh
"""
struct Mesh
    mVertices::Vector{Vertex}
    mEdges::Vector{Edge}
	mFaces::Vector{Face}
end

############################################################################################
# Mesh Creation
############################################################################################

"""
    resolve_edges!( faces::Vector{Face} )

Stores face vertices in half-edge order and creates all edges of the mesh. Faces will have 
their mVertices vector filled in the correct order. Also, the vertices belonging to the 
faces will have their mEdgeIDS and mEdges vectors filled.

# Outputs
- `createdEdges::Vector{Edge}`: vector of pointers to created edges
"""
function resolve_edges!( faces::Vector{Face} )
	# List of faces ready to be processed
	readyList = zeros( Int64, ( 3, length(faces) ) ); 
    readyInt = 1; readyUpdate = 2; readyList[1,1] = 1;

	# Dictionary keeping track of edges already been created
	createdEdges = Vector{Edge}(); 
    share = [ 0, 0 ]; shareFind = 0; v1Arr = [ 0, 0 ]; v2Arr = [ 0, 0 ];
	eIndex = 1;

	# Iterate through all faces to resolve
	while readyInt <= length(faces)
		# Get current face
		currFace = readyList[ 1, readyInt ]; readyInt += 1;
		face = faces[ currFace ]; fVerts = face.mVertices;

		if readyUpdate == 2
			# Case for the first face

			for (i, v) in enumerate( fVerts )
				# Find next vertex and create new edge
				vN = fVerts[ i%length(fVerts) + 1 ];
				check = [ min(v.mIndex, vN.mIndex), max(v.mIndex, vN.mIndex) ];
				nE = Edge( eIndex, [ v, vN ] );

				# Set face edges
				face.mEdges[i] = nE;

				# Set neighboring face
				if nE.mFaces[1] == 0
					nE.mFaces[1] = face.mIndex;
				elseif nE.mFaces[2] == 0
					nE.mFaces[2] = face.mIndex;
				else
					error( "Edge can only neighbor two faces." )
				end

				push!( createdEdges, nE );

				# Update vertex edge vectors
				push!( v.mEdgeIDs, check ); push!( v.mEdges, eIndex );
				push!( vN.mEdgeIDs, check ); push!( vN.mEdges, eIndex );

				eIndex += 1;
			end
		else
			# Case for all other faces, must follow orientation of first face

			# Find a neighboring face which has been resolved and follow that
			follow = face; followIndex = 0;
			for neighbor in face.mNeighbors
				follow = faces[ neighbor ]; followIndex = neighbor;
				if readyList[ 2, neighbor ] > 0
					break
				end
			end
			
			# Find the two vertices the faces share
			while shareFind < 2
				for ( vi, v ) in enumerate( follow.mVertices )
					if shareFind < 2
						for ( wi, w ) in enumerate( face.mVertices )
							if v.mIndex == w.mIndex
								shareFind += 1;
								share[ shareFind ] = v.mIndex; 
                                v1Arr[ shareFind ] = vi; 
                                v2Arr[ shareFind ] = wi;
								break
							end
						end
					end
				end
			end
			shareFind = 0; 
            length1 = length( follow.mVertices ); 
            length2 = length( face.mVertices );

			# Check whether face follows the same orientation
			v1a = v1Arr[1]; v1b = v1Arr[2]; 
			orderA = ( v1a > v1b && !( v1b == 1 && v1a == length1 ) ) || 
                ( v1a == 1 && v1b == length1 );

			v2a = v2Arr[1]; v2b = v2Arr[2];
			orderB = ( v2a > v2b && !( v2b == 1 && v2a == length2 ) ) || 
                ( v2a == 1 && v2b == length2 );

			# Enforce correct orientation
			if !( orderA ⊻ orderB )
				reverse!( face.mVertices );
			end

			# Create new edges
			for (i, v) in enumerate( fVerts ) 
				vN = fVerts[i%length(fVerts) + 1];

				# Check to see whether edge has already been created
				check = [ min(v.mIndex, vN.mIndex), max(v.mIndex, vN.mIndex) ]; vEIndex = 0;
				for (ii, id) in enumerate( v.mEdgeIDs )
					if check == id
						vEIndex = ii;
						break
					end
				end

				# Create edges
				if vEIndex == 0
					nE = Edge( eIndex, [ v, vN ] );
					push!( createdEdges, nE );

					# Set face edges
					face.mEdges[i] = nE;

					# Set neighboring face
					if nE.mFaces[1] == 0
						nE.mFaces[1] = face.mIndex;
					elseif nE.mFaces[2] == 0
						nE.mFaces[2] = face.mIndex;
					else
						error( "Edge can only neighbor two faces." )
					end

					# Update vertex edge vectors
					push!( v.mEdgeIDs, check ); push!( v.mEdges, eIndex );
					push!( vN.mEdgeIDs, check ); push!( vN.mEdges, eIndex );

					eIndex += 1;
				else					
					nE = createdEdges[ v.mEdges[vEIndex] ];

					# Set face edges
					face.mEdges[i] = nE;

					# Set neighboring face
					if nE.mFaces[1] == 0
						nE.mFaces[1] = face.mIndex;
					elseif nE.mFaces[2] == 0
						nE.mFaces[2] = face.mIndex;
					else
						error( "Edge can only neighbor two faces." )
					end
				end
			end
		end

		# Update list of stuff to still be done/was already done
		readyList[ 3, face.mIndex ] = 1;
		for f in face.mNeighbors
			if f > 0 && readyList[ 3, f ] == 0
				readyList[ 1, readyUpdate ] = f; readyUpdate += 1; readyList[ 3, f ] = 1;
			end
		end
		readyList[ 2, face.mIndex ] = 1;
	end

	# Return vector of edges
	return createdEdges
end

############################################################################################
# Mesh Utility
############################################################################################

"""
    is_boundary( vertex::Vertex )

Determines if a vertex is a boundary vertex of the 2-dimensional mesh
"""
function is_boundary( vertex::Vertex )
    return length( vertex.mEdges ) != length( vertex.mFaces );
end

"""
    is_boundary( edge::Edge )

Determines if an edge is a boundary edge of the 2-dimensional mesh
"""
function is_boundary( edge::Vertex )
    return edge.mFaces[2] < 1;
end

"""
    is_boundary( face::Face )

Determines if a face is a boundary face of the 2-dimensional mesh
"""
function is_boundary( face::Face )
	return face.mNeighbors[end] == 0;
end

"""
    is_extraordinary( vertex::Vertex )

Determines if a vertex is an extraordinary vertex of the 2-dimensional mesh
"""
function is_extraordinary( vertex::Vertex )
	return !is_boundary( vertex ) && length( vertex.mEdges ) != 4
end

"""
    is_extraordinary( face::Face )

Determines if a face is an extraordinary face of the 2-dimensional mesh
"""
function is_extraordinary( face::Face )
	return length( face.mEdges ) != 4
end

"""
    get_from_half_edges( mesh::Mesh, vertex::Vertex, face::Face, args::Vector{String} )

Uses underlying half-edge structure to traverse the mesh to get mesh information.

# Inputs
- `mesh`: underlying mesh
- `vertex`: vertex at base of half-edge
- `face`: face which half-edge belongs to
- `args`: array of strings of half-edge commands in order (next/prev/opp) from left to right

# Outputs
- `(currBase, currFace)`: a tuple which contains the vertex at the base of the current 
    half-edge and the current face
"""
function get_from_half_edges( mesh::Mesh, vertex::Vertex, face::Face, args::Vector{String} )

	# Get mesh constructs
	faces = mesh.mFaces;

	# Check that the vertex is part of the face
	currFace = face;
	currBase = vertex;
	currEnd = vertex;
	faceLen = length( currFace.mVertices );
	checkSum = 0;
	for v in currFace.mVertices
		if v.mIndex == vertex.mIndex
			checkSum = 1;
		end
	end

	if checkSum == 0
		throw( ArgumentError( "Must specify a vertex that lies on the specified face." ) )
	end

	# Loop through all arguments to get desired construct
	for arg in args
		# Find vertex on face
		baseIndex = 1;
		i = 1;
		for v in currFace.mVertices
			if v.mIndex == currBase.mIndex    
				baseIndex = i;
				break
			end
			i += 1;
		end

		# Deal with each case of the halfEdge construct
		if arg == "next"
			currBase = currFace.mVertices[ (baseIndex % faceLen) + 1 ];
		elseif arg == "prev"
			if baseIndex == 1
				currBase = currFace.mVertices[ faceLen ];
			else
				currBase = currFace.mVertices[ baseIndex - 1 ];
			end
		elseif arg == "opp"
			for f in currFace.mNeighbors
				fT = faces[f];
				currEnd = currFace.mVertices[ (baseIndex % faceLen) + 1 ];
				fCheck = 0;
				for v in fT.mVertices
					if v.mIndex == currEnd.mIndex || v.mIndex == currBase.mIndex
						fCheck += 1;
					end
				end
				if fCheck == 2
					currFace = fT;
					currBase = currEnd;
					faceLen = length( currFace.mVertices );
					break
				end
			end
		else
			print( string("In get_from_half_edges, arguments must be next/prev/opp. 
                Ignoring: ", arg) );
		end
	end

	return ( currBase, currFace );
end