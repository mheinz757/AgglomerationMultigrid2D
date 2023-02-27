module AgglomerationMultigrid2D

# import needed packages
import LinearAlgebra as la
import LinearAlgebra: ldiv!, \
import SparseArrays as sp
import SuiteSparse
# import BlockDiagonals as bd

abstract type AbstractElement end
abstract type AbstractAgglomeratedDgElement <: AbstractElement end
abstract type AbstractMesh end
abstract type AbstractSmoother end

include("meshes.jl")
# include("boundary_conditions.jl")
# include("legendre.jl")
# include("gaussquad.jl")
# include("reference_element.jl")

# include("cgMesh.jl")
# include("dgMesh.jl")
# include("agglomeratedDgMesh.jl")

# include("smoother.jl")
# include("interpolation.jl")

# include("meshHeirarchy.jl")
# include("solvers.jl")

end # end of module