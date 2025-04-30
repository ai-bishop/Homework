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

struct quadrule2
    label::String
    m::Int
    n::Int
    dim::Int
    ξ::Array
    η::Array
    w1::Array
    # w2::Array
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






function gauss_legendre_2d(n)
    quad1 = gauss_legendre_1d(n)

    nint = n^2
    w = zeros(n,n)
    r = zeros(n,n)
    s = zeros(n,n)

    for (i,ξ) in enumerate(quad1.ξ)
        for (j,η) in enumerate(quad1.ξ)
            r[i,j] = ξ
            s[i,j] = η
            w[i,j] = quad1.w[i]*quad1.w[j]
        end
    end
    X = hcat( r[:], s[:] )
    W = w[:]

    return quadrule("2D GL", nint, 2, X, W,
                    zip(eachrow(X), W))

end

function gauss_legendre_2dim(m,n)
    quad1 = gauss_legendre_1d(m)
    quad2 = gauss_legendre_1d(n)

    nint = n^2
    w = zeros(n,n)
    r = zeros(n,n)
    s = zeros(n,n)

    for (i,ξ) in enumerate(quad1.ξ)
        for (j,η) in enumerate(quad2.ξ)
            r[i,j] = ξ
            s[i,j] = η
            w[i,j] = quad1.w[i]*quad2.w[j]
        end
    end
    X = hcat( r[:], s[:] )
    W = w[:]

    return quadrule("2D GL", nint, 2, X, W,
                    zip(eachrow(X), W))
end



function gauss_legendre_3d(n)
    quad1 = gauss_legendre_1d(n)

    nint = n^3
    w = zeros(n,n,n)
    r = zeros(n,n,n)
    s = zeros(n,n,n)
    t = zeros(n,n,n)

    for (i,ξ) in enumerate(quad1.ξ)
        for (j,η) in enumerate(quad1.ξ)
            for (k, ρ) in enumerate(quad1.ξ)
                r[i,j] = ξ
                s[i,j] = η
                t[i,j] = ρ
                w[i,j] = quad1.w[i]*quad1.w[j]*quad1.w[k]
            end
        end
    end
    X = hcat( r[:], s[:], t[:] )
    W = w[:]

    return quadrule("3D GL", nint, 3, X, W,
                    zip(eachrow(X), W))

end




end