testPath = @__DIR__;
aggPath = dirname(testPath);

include(aggPath * "/src/AgglomerationMultigrid2D.jl")
include(testPath * "/plot_matlab.jl")

import LinearAlgebra as la
import .AgglomerationMultigrid2D as aggmg


p = 2;
refEl = aggmg.ReferenceElement( p, :tri );
@show refEl.mNodesX
@show refEl.mGaussQuadNodes
@show refEl.mGaussQuadWeights
@show refEl.mBasisFunCoeff

s = get_default_msession();
aggmg.plot_ref_el( s, refEl; vertIndex = true );

funVal, gradVal = aggmg.evaluate_nodal_basis_fun_and_grad( refEl.mShape, 
    refEl.mBasisFunCoeff, refEl.mNodesX );
funVal2 = aggmg.evaluate_nodal_basis_fun( refEl.mShape, refEl.mBasisFunCoeff, 
    refEl.mNodesX );

println( " " );
println( la.norm( funVal - la.I( size(funVal, 1) ) ) );
println( la.norm( funVal - funVal2 ) );