module Quadrature

using LinearAlgebra

struct quadrule
    label::String
    n::Int
    dim::Int
    ξ::Array
    w::Array
    iterator::Any
end


function gauss_legendre_1d(n)
    β = @. 0.5 / sqrt(1 - (2 * (1:(n - 1)))^(-2))
    T = diagm(-1 => β, 1 => β)
    λ, V = eigen(T)
    p = sortperm(λ)

    # nodes
    ξ = λ[p]

    # weights
    w = @. 2V[1, p]^2

    return quadrule("1D GL", n, 1, ξ, w,
                    zip(eachrow(ξ), w))
end




# not used - will do in script itself
function gauss_legendre_2d(n)
    # n is number of nodes in one dimension
    # getting 2d nodes and weights via 2x 1d nodes and weights
    

    # Get 1D nodes and weights
    x_nodes, x_weights = gauss_legendre_1d(n)
    y_nodes, y_weights = gauss_legendre_1d(n)


    return 



end

# remove? might be nonfunctional...
function gauss_legendre_2d_old(m,n)
    β = @. 0.5 / sqrt(1 - (2 * (1:(m - 1)))^(-2))
    γ = @. 0.5 / sqrt(1 - (2 * (1:(n - 1)))^(-2))

    T = diagm(-1 => β, 1 => β)
    U = diagm(-1 => γ, 1 => γ)

    λ1, V1 = eigen(T)
    λ2, V2 = eigen(U)

    p1 = sortperm(λ1)
    p2 = sortperm(λ1)


    # nodes
    ξ = λ1[p1]
    η = λ2[p2]

    # weights
    w1 = @. 2V[1, p1]^2
    w2 = @. 2V[1, p2]^2


    return quadrule2("2D GL", m, n, ξ, η, w,
                    zip(eachrow(ξ), w1), zip(eachrow(η), w2))

end

end