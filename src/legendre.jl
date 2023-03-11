"""
    legendre_val( x::Union{Real, Vector{Real}}, n::Integer )

Evaluates all the Legendre polynomials up to degree `n` at points `x`

# Arguments
- `x`: points to evaluate at
- `n`: highest degree of Legendre polynomial to evaluate

# Outputs
- `funVal`: values of all the Legendre polynomials up to degree `n` at points `x`
"""
function legendre_val( x::Union{Real, Vector{<:Real}}, n::Integer )
    funVal = zeros( length(x), n+1 );
    funVal[:,1] .= 1.0;
    if n >= 1
        funVal[:,2] .= x;
        for i = 2:n
            funVal[:,i+1] .= ( (2*i-1) * x .* funVal[:,i] - (i-1) * funVal[:,i-1] ) / i;
        end
    end

    return funVal
end


"""
    legendre_val_and_deriv( x::Union{Real, Vector{Real}}, n::Integer )

Evaluates values of all the Legendre polynomials up to degree `n` and their derivatives at 
points `x`

# Arguments
- `x`: points to evaluate at
- `n`: highest degree of Legendre polynomial to evaluate

# Outputs
- `funVal`: values of all the Legendre polynomials up to degree `n` at points `x`
- `derivVal`: values of the derivative of all the Legendre polynomials up to degree `n` at 
    points `x`
"""
function legendre_val_and_deriv( x::Union{Real, Vector{<:Real}}, n::Integer )
    funVal = zeros( length(x), n+1 );
    derivVal = zeros( length(x), n+1 );
    funVal[:,1] .= 1.0;
    derivVal[:,1] .= 0.0;
    if n >= 1
        funVal[:,2] .= x;
        derivVal[:,2] .= 1.0;
        for i = 2:n
            funVal[:,i+1] .= ( (2*i-1) * x .* funVal[:,i] - (i-1) * funVal[:,i-1] ) / i;
            derivVal[:,i+1] .= (2*i-1) * funVal[:,i] + derivVal[:,i-1];
        end
    end

    return funVal, derivVal
end