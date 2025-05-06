module Plate2D

using LinearAlgebra
include("Shapefunctions.jl")
using .Shapefunctions
include("Preprocess.jl")
using .Preprocess

# copied from Bar1D
# update to confirm works in 2dims

function assemble_stiffness(mesh, properties, quad_rules)
    ned = mesh.ned
    totaldofs = ned * mesh.nnp
    K = zeros(totaldofs, totaldofs)

    ## loop over each element from each part of the mesh
    for (element_type, ien) in mesh.IEN
        nel = size(ien, 1)
        nen = std_element_defs[element_type].nen
        nee = ned * nen
        LM = mesh.LM[element_type]
        element_quad_rule = quad_rules[element_type]
        N = shapefunc(element_type)

        


        for e in 1:nel
            A = ien[e, 1:nen]
            xe = mesh.x[A]
            ye = mesh.y[A]

            ke = element_stiffness(xe, ye, N, properties, element_quad_rule)

            # assemble element stiffness into the Global stiffness Matrix
            for loop1 in 1:nee
                i = LM[loop1, e]
                for loop2 in 1:nee
                    j = LM[loop2, e]
                    K[i, j] += ke[loop1, loop2]
                end
            end
        end
    end
    return K
end

function element_stiffness(xe, ye, N, properties, quad_rules) # looks like element_forcing
    
    # check probs
    E = properties.E
    T = properties.T

    ned = 1
    nen = length(xe)
    nee = ned * nen
    ke = zeros(nee, nee)

    # integration loop
    for (X, W) in quad_rules.iterator

        # Evaluate the shape function
        Ne, Nξ, Nη = N(X) # shape function (not used) and shape function derivative

        # build Jacobian
        Jac = [(Nξ' * xe) (Nη' * xe); (Nξ' * ye) (Nη' * ye)]

        B = inv(Jac) * [Nξ'; Nη']
        
        ke += B' * E * T * B * det(Jac) * W

        # ke += Nξ * Nξ' * E * T * Nη * Nη' * W / det(Jac)

    end

    return ke
end

function assemble_rhs(mesh, external_forcing, quad_rules)
    ned = mesh.ned
    totaldofs = ned * mesh.nnp
    F = zeros(totaldofs)

    ## loop over each element from each part of the mesh
    for (element_type, ien) in mesh.IEN
        nel = size(ien, 1)
        N = shapefunc(element_type)
        element_quad_rule = quad_rules[element_type]
        nen = std_element_defs[element_type].nen
        nee = ned * nen
        LM = mesh.LM[element_type]

        for e in 1:nel
            A = ien[e, 1:nen]
            xe = mesh.x[A]
            ye = mesh.y[A]

            fe = element_forcing(xe, ye, N, element_quad_rule, external_forcing)

            # assemble element stiffness into the Global stiffness Matrix
            for loop1 in 1:nee
                i = LM[loop1, e]
                F[i] += fe[loop1]
            end
        end
    end
    return F
end

function element_forcing(xe, ye, N, element_quad_rule, external_forcing)
    ned = 1
    nen = length(xe)
    nee = ned * nen
    fe = zeros(nee)

    # integration loop
    for (ξ, w) in element_quad_rule.iterator

        # Evaluate the shape function
        Ne, Nξ, Nη = N(ξ)

        # evaluate the external loading at x(ξ)
        x = dot(Ne, xe)
        y = dot(Ne, ye)

        fext = external_forcing(x, y)

        # build Jacobian
        detJ = dot(Nξ, xe) - dot(Nη, ye)

        # dV0
        dV0 = detJ * w

        # integrate ke:
        fe += Ne * fext * dV0
    end

    return fe
end

end