"""
    monomial_val( x::Union{Vector{<:Real}, Matrix{<:Real}}, n::Integer )

Evaluates all bivariate monomials up to degree `n` at points `x`

# Arguments
- `x`: points to evaluate at
- `n`: highest degree of bivariate monomials to evaluate

# Outputs
- `funVal`: values of all bivariate monomials up to degree `n` at points `x`
"""
function monomial_val( x::Union{Vector{<:Real}, Matrix{<:Real}}, n::Integer )
    if typeof(x) <: Vector{<:Real}
        x = x';
    end

    numFun = Int( (n+1)*(n+2)/2 );
    funVal = zeros( size( x, 1 ), numFun );

    col = 0;
    for k = 0:n, j = 0:k
        col += 1;
        funVal[:,col] .= x[:,1].^(k-j) .* x[:,2].^j;
    end

    return funVal
end


"""
    monomial_val_and_grad( x::Union{Vector{<:Real}, Matrix{<:Real}}, n::Integer )

Evaluates values of all bivariate monomials up to degree `n` and their gradients at 
points `x`

# Arguments
- `x`: points to evaluate at
- `n`: highest degree of bivariate monomials to evaluate

# Outputs
- `funVal`: values of all bivariate monomials up to degree `n` at points `x`
- `gradVal`: values of the gradient of all bivariate monomials up to degree `n` at 
    points `x`
"""
function monomial_val_and_grad( x::Union{Vector{<:Real}, Matrix{<:Real}}, n::Integer )
    if typeof(x) <: Vector{<:Real}
        x = x';
    end

    numFun = Int( (n+1)*(n+2)/2 );
    funVal = zeros( size( x, 1 ), numFun );
    gradVal = zeros( size( x, 1 ), numFun, 2 );

    col = 0;
    for k = 0:n, j = 0:k
        col += 1;
        funVal[:,col] .= x[:,1].^(k-j) .* x[:,2].^j;
        
        if j == k
            gradVal[:,col,1] .= 0;
        else
            gradVal[:,col,1] .= (k-j) * x[:,1].^(k-j-1) .* x[:,2].^j;
        end

        if j == 0
            gradVal[:,col,2] .= 0;
        else
            gradVal[:,col,2] .= j * x[:,1].^(k-j) .* x[:,2].^(j-1);
        end
    end

    return funVal, gradVal
end