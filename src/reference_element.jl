struct ReferenceElement
    mShape::Symbol
    mP::Int64
    mNodesX::Matrix{Float64}

    mGaussQuadNodes::Matrix{Float64}
    mGaussQuadWeights::Vector{Float64}

    mBasisFunCoeff::Union{ Matrix{Float64}, Array{Float64, 3} }
    mBasisGQFunVal::Matrix{Float64}
    mBasisGQGradVal::Array{Float64, 3}

    mMassMatrix::Matrix{Float64}
end

############################################################################################
# ReferenceElement creation 
############################################################################################

function ReferenceElement( mP, mShape )

    if mShape == :tri
        numNodes = Int( (mP+1)*(mP+2)/2 );

        mNodesX = zeros( numNodes, 2 );
        if mP >= 1
            x1d = range(0, 1, mP+1);

            # Initialize corners first
            mNodesX[1:3, :] = [ 0.0 0.0; 1.0 0.0; 0.0 1.0 ];

            if mP >= 2
                # Then go along edges
                i = 1;
                mNodesX[ (3 + (i-1)*(mP-1) + 1):(3 + i*(mP-1)), 1 ] = x1d[ 2:mP ];
                mNodesX[ (3 + (i-1)*(mP-1) + 1):(3 + i*(mP-1)), 2 ] .= 0.0;

                i = 2;
                mNodesX[ (3 + (i-1)*(mP-1) + 1):(3 + i*(mP-1)), 1 ] = x1d[ mP:-1:2 ];
                mNodesX[ (3 + (i-1)*(mP-1) + 1):(3 + i*(mP-1)), 2 ] = x1d[ 2:mP ];

                i = 3;
                mNodesX[ (3 + (i-1)*(mP-1) + 1):(3 + i*(mP-1)), 1 ] .= 0.0;
                mNodesX[ (3 + (i-1)*(mP-1) + 1):(3 + i*(mP-1)), 2 ] = x1d[ mP:-1:2 ];

                if mP >= 3
                    # Then do the interior
                    ind = 3 + 3*(mP-1);
                    for j = 2:(mP-1), i = 2:(mP+1-j)
                        ind += 1;
                        mNodesX[ ind, : ] = [ x1d[i] x1d[j] ]';
                    end
                end
            end
        else
            mNodesX[1,:] = [ 1.0/3.0, 1.0/3.0 ];
            mBasisFunCoeff[1,1] = 1.0;
        end

        V = monomial_val( mNodesX, mP );
        mBasisFunCoeff = la.inv(V);
    elseif mShape == :quad
        numNodes = (mP + 1)^2;

        x1d = zeros( mP + 1 );
        if mP >= 1
            x1d[ 1:2 ] = [ -1.0, 1.0 ];
            x1d[ 3:(mP+1) ] = cos.( pi * ( (mP-1):-1:1 ) / mP );
        else
            x1d[1] = 0.0;
        end

        V = legendre_val( x1d, mP );
        basisFunCoeff1d = la.inv(V);

        mNodesX = zeros( numNodes, 2 );
        mBasisFunCoeff = zeros( mP + 1, numNodes, 2 );
        if mP >= 1
            # Initialize corners first
            mNodesX[1:4,:] = [ -1.0 -1.0; 1.0 -1.0; 1.0 1.0; -1.0 1.0 ];

            mBasisFunCoeff[ :, 1, 1 ] = basisFunCoeff1d[:,1];
            mBasisFunCoeff[ :, 1, 2 ] = basisFunCoeff1d[:,1];

            mBasisFunCoeff[ :, 2, 1 ] = basisFunCoeff1d[:,2];
            mBasisFunCoeff[ :, 2, 2 ] = basisFunCoeff1d[:,1];

            mBasisFunCoeff[ :, 3, 1 ] = basisFunCoeff1d[:,2];
            mBasisFunCoeff[ :, 3, 2 ] = basisFunCoeff1d[:,2];

            mBasisFunCoeff[ :, 4, 1 ] = basisFunCoeff1d[:,1];
            mBasisFunCoeff[ :, 4, 2 ] = basisFunCoeff1d[:,2];

            if mP >= 2
                # Then go along edges
                i = 1;
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 1 ] = x1d[ 3:(mP+1) ];
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 2 ] .= -1.0;
                for j in 1:(mP-1)
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 1 ] = basisFunCoeff1d[:,j+2];
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 2 ] = basisFunCoeff1d[:,1];
                end

                i = 2;
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 1 ] .= 1.0;
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 2 ] = x1d[ 3:(mP+1) ];
                for j in 1:(mP-1)
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 1 ] = basisFunCoeff1d[:,2];
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 2 ] = basisFunCoeff1d[:,j+2];
                end

                i = 3;
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 1 ] = x1d[ (mP+1):-1:3 ];
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 2 ] .= 1.0;
                for j in 1:(mP-1)
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 1 ] = 
                        basisFunCoeff1d[:,mP+1-j+1];
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 2 ] = basisFunCoeff1d[:,2];
                end

                i = 4;
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 1 ] .= -1.0;
                mNodesX[ (4 + (i-1)*(mP-1) + 1):(4 + i*(mP-1)), 2 ] = x1d[ (mP+1):-1:3 ];
                for j in 1:(mP-1)
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 1 ] = basisFunCoeff1d[:,1];
                    mBasisFunCoeff[ :, 4 + (i-1)*(mP-1) + j, 2 ] = 
                        basisFunCoeff1d[:,mP+1-j+1];
                end

                # Then do the interior
                mNodesX[ (4 + 4*(mP-1) + 1):(mP + 1)^2, 1 ] = 
                    [ x1d[i] for j in 3:(mP+1) for i in 3:(mP+1) ];
                mNodesX[ (4 + 4*(mP-1) + 1):(mP + 1)^2, 2 ] = 
                    [ x1d[j] for j in 3:(mP+1) for i in 3:(mP+1) ];
                for j in 1:(mP-1), i in 1:(mP-1)
                    mBasisFunCoeff[ :, 4 + (4+j-1)*(mP-1) + i, 1 ] = basisFunCoeff1d[:,i+2];
                    mBasisFunCoeff[ :, 4 + (4+j-1)*(mP-1) + i, 2 ] = basisFunCoeff1d[:,j+2];
                end
            end
        else 
            mNodesX[1,:] = [ 0.0, 0.0 ];
            mBasisFunCoeff[1,1,:] = [ 1.0, 1.0 ];
        end
    else 
        throw( ArgumentError( "Shape must be triangle (:tri) or quadrilateral (:quad)." ) );
    end
    

    mGaussQuadNodes, mGaussQuadWeights = gauss_quad_2d( 2*mP, mShape );
    mBasisGQFunVal, mBasisGQGradVal = evaluate_nodal_basis_fun_and_grad( mShape, 
        mBasisFunCoeff, mGaussQuadNodes );

    mMassMatrix = zeros( numNodes, numNodes );

    for j = 1:numNodes
        for i = 1:j
            for (l, w) in enumerate( mGaussQuadWeights )
                mMassMatrix[i,j] += w * mBasisGQFunVal[l,i] * mBasisGQFunVal[l,j];
            end
        end
    end

    for j = 1:numNodes
        for i = (j+1):numNodes
            mMassMatrix[i,j] = mMassMatrix[j,i];
        end
    end

    return ReferenceElement( mShape, mP, mNodesX, mGaussQuadNodes, mGaussQuadWeights, 
        mBasisFunCoeff, mBasisGQFunVal, mBasisGQGradVal, mMassMatrix );
end

############################################################################################
# ReferenceElement utility
############################################################################################

function evaluate_nodal_basis_fun( shape::Symbol, basisFunCoeff, nodes )

    if shape == :tri
        basisFunVal = zeros( size( nodes, 1 ), size( basisFunCoeff, 1 ) );
        p = Int( ( sqrt( 1 + 8 * size(basisFunCoeff , 1) ) - 3.0 ) / 2.0 );

        monomialVal = monomial_val( nodes, p );
        basisFunVal[:,:] = monomialVal * basisFunCoeff;

    elseif shape == :quad
        basisFunVal = zeros( size( nodes, 1 ), size( basisFunCoeff, 2 ) );
        p = size( basisFunCoeff, 1 ) - 1;
    
        legendreValX = legendre_val( nodes[:,1], p );
        legendreValY = legendre_val( nodes[:,2], p );
        basisFunVal[:,:] = ( legendreValX * basisFunCoeff[:,:,1] ) .* 
            ( legendreValY * basisFunCoeff[:,:,2] );
    else 
        throw( ArgumentError( "Shape must be triangle (:tri) or quadrilateral (:quad)." ) );
    end

    return basisFunVal;
end

function evaluate_nodal_basis_fun_and_grad( shape::Symbol, basisFunCoeff, nodes )

    if shape == :tri
        basisFunVal = zeros( size( nodes, 1 ), size( basisFunCoeff, 1 ) );
        basisGradVal = zeros( size( nodes, 1 ), size( basisFunCoeff, 1 ), 2 );
        p = Int( ( sqrt( 1 + 8 * size(basisFunCoeff , 1) ) - 3.0 ) / 2.0 );

        monomialVal, monomialGrad = monomial_val_and_grad( nodes, p );
        basisFunVal[:,:] = monomialVal * basisFunCoeff;
        basisGradVal[:,:,1] = monomialGrad[:,:,1] * basisFunCoeff;
        basisGradVal[:,:,2] = monomialGrad[:,:,2] * basisFunCoeff;
    elseif shape == :quad
        basisFunVal = zeros( size( nodes, 1 ), size( basisFunCoeff, 2 ) );
        basisGradVal = zeros( size( nodes, 1 ), size( basisFunCoeff, 2 ), 2 );
        p = size( basisFunCoeff, 1 ) - 1;
    

        legendreValX, legendreDerivX = legendre_val_and_deriv( nodes[:,1], p );
        legendreValY, legendreDerivY = legendre_val_and_deriv( nodes[:,2], p );
        basisFunVal[:,:] = ( legendreValX * basisFunCoeff[:,:,1] ) .* 
            ( legendreValY * basisFunCoeff[:,:,2] );
        basisGradVal[:,:,1] = ( legendreDerivX * basisFunCoeff[:,:,1] ) .* 
            ( legendreValY * basisFunCoeff[:,:,2] );
        basisGradVal[:,:,2] = ( legendreValX * basisFunCoeff[:,:,1] ) .* 
            ( legendreDerivY * basisFunCoeff[:,:,2] );
    else
        throw( ArgumentError( "Shape must be triangle (:tri) or quadrilateral (:quad)." ) );
    end

    return basisFunVal, basisGradVal;
end
