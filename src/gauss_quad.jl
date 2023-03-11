include("gauss_quad_2d_dat.jl")

############################################################################################
# 1d Gaussian Quadrature
############################################################################################

"""
    gauss_quad_1d( p::Integer )

Gaussian quadrature on [-1,1] for given degree of precision `p`
"""
function gauss_quad_1d( p::Integer )
    p = max(p,1);
    n = ceil((p+1)/2);
    b = 1:n-1;
    b = @. b / sqrt(4*b^2 - 1);
    eval, evec = la.eigen(la.diagm(1 => b, -1 => b));

    nodes = eval;
    weights = 2*evec[1,:].^2;

    return nodes, weights;
end

############################################################################################
# 2d Gaussian Quadrature
############################################################################################

"""
    gauss_quad_2d( p::Integer, shape::Symbol )

Gaussian quadrature on reference quadrilateral [-1, 1]^2 or reference triangle with 
vertices (0,0), (1,0), and (0,1) for given degree of precision `p`.
"""
function gauss_quad_2d( p::Integer, shape::Symbol )
    p = max(p,1);
    if shape == :tri
        if p <= length( gaussQuadNodesTri )
            nodes = gaussQuadNodesTri[ p ];
            weights = gaussQuadWeightsTri[ p ] / 2.0;
        else
            x1, w1 = gauss_quad_1d( p ); # x1 = 2*x1 .- 1.0; w1 = 2.0 * w1;
            x2 = [ x1[i] for i in eachindex(x1), j in eachindex(x1) ];
            y2 = [ x1[j] for i in eachindex(x1), j in eachindex(x1) ];
            x0 = [ x2[:] y2[:] ];
            x = ( 1 .+ x0[:,1] - x0[:,2] - x0[:,1] .* x0[:,2] ) / 4.0;
            y = ( 1 .+ x0[:,2] ) / 2.0;
            nodes = [ x y ];
            w0 = w1 * w1';
            weights = w0[:] .* ( 1 .- x0[:,2] ) / 8.0;
        end   
    elseif shape == :quad
        x1, w1 = gauss_quad_1d( p );
        x2 = [ x1[i] for i in eachindex(x1), j in eachindex(x1) ];
        y2 = [ x1[j] for i in eachindex(x1), j in eachindex(x1) ];
        nodes = [ x2[:] y2[:] ];
        weights = ( w1 * w1' )[:];
    else
        throw( ArgumentError( "Shape must be triangle (:tri) or quadrilateral (:quad)." ) );
    end

    return nodes, weights;
end