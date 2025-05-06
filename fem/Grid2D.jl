# modify original code of Bar1D.jl for 2dim
module Grid2D

using LinearAlgebra
include("Shapefunctions.jl")
using .Shapefunctions
include("Preprocess.jl")
using .Preprocess





function assemble_resistance(mesh, BC_g_list, nnp)

    if nnp == 25
        # is 5x5 mesh
        # write temperatures as functions of adjacent nodes for interior
        locns = mesh.x

        # TempMatrix(i) = (locns(i-1) + locns(i+1) + locns(i-5) + locns(i+5)) / 4
        InteriorMatrix = zeros(3,3);

        InteriorMatrix[1] = (locns[6] + locns[8] + locns[2] + locns[12]) / 4 # 7
        InteriorMatrix[2] = (locns[7] + locns[9] + locns[3] + locns[13]) / 4 # 8
        InteriorMatrix[3] = (locns[8] + locns[10] + locns[4] + locns[14]) / 4 # 9

        InteriorMatrix[4] = (locns[11] + locns[13] + locns[7] + locns[17]) / 4 # 12
        InteriorMatrix[5] = (locns[12] + locns[14] + locns[8] + locns[18]) / 4 # 13
        InteriorMatrix[6] = (locns[13] + locns[15] + locns[9] + locns[19]) / 4 # 14

        InteriorMatrix[7] = (locns[16] + locns[18] + locns[12] + locns[22]) / 4 # 17
        InteriorMatrix[8] = (locns[17] + locns[19] + locns[13] + locns[23]) / 4 # 18
        InteriorMatrix[9] = (locns[18] + locns[20] + locns[14] + locns[24]) / 4 # 19

        return InteriorMatrix

    else

        error("size unaccounted for. please expand code.")

    end

end

function assemble_stiffness(mesh, quad_rules)
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

            hf = element_stiffness(xe, ye, N, element_quad_rule)

            # assemble element stiffness into the Global stiffness Matrix
            for loop1 in 1:nee
                i = LM[loop1, e]
                for loop2 in 1:nee
                    j = LM[loop2, e]
                    K[i, j] += hf[loop1, loop2]
                end
            end
        end
    end
    return K
end

function element_stiffness(xe, ye, N, quad_rules) # looks like element_forcing

    ned = 1
    nen = length(xe)
    nee = ned * nen
    ke = zeros(nee, nee)

    # integration loop
    for (ξ, w) in quad_rules.iterator

        # Evaluate the shape function
        Ne, Nξ, Nη = N(ξ) # shape function and shape function derivative
        
        Jac = [(Nξ' * xe) (Nη' * xe); (Nξ' * ye) (Nη' * ye)]

        B = inv(Jac) * [Nξ'; Nη']

        K = 1
        
        ke += B' * K * B * det(Jac) * w




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
    ned = length(ye)
    nen = length(xe)
    nee = ned * nen
    fe = zeros(nee)

    # integration loop
    for (ξ, η, w) in element_quad_rule.iterator

        # Evaluate the shape function
        Ne, Nξ, Nη = N(ξ, η)

        # evaluate the external loading at x(ξ)
        x = dot(Ne, xe)
        fext = external_forcing(x)

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