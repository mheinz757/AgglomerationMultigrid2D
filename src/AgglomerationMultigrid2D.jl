module AgglomerationMultigrid2D

# import needed packages
import LinearAlgebra as la
import SparseArrays as sp
import SuiteSparse
import MATLAB as ml

import Base: *, \
import Base: Array, Matrix, show, similar, size
import LinearAlgebra: ldiv!, lu, mul!
import SparseArrays: sparse

abstract type AbstractElement end
abstract type AbstractAgglomeratedDgElement <: AbstractElement end
abstract type AbstractMesh end
abstract type AbstractSmoother end

include("meshes.jl")
include("meshio.jl")
include("refinement.jl")

include("legendre.jl")
include("bivariate_monomial.jl")
include("gauss_quad.jl")
include("reference_element.jl")

# include("boundary_conditions.jl")
# include("cgMesh.jl")
# include("dgMesh.jl")
# include("agglomeratedDgMesh.jl")

# include("smoother.jl")
# include("interpolation.jl")

# include("meshHeirarchy.jl")
# include("solvers.jl")

include("meshplot.jl")

end # end of module