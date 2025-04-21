using Pkg
Pkg.activate(".")

include("Preprocess.jl")
include("Quadrature.jl")
include("Grid2D.jl")
include("SimpleVisualization.jl")
using GLMakie

quad_rules = Dict("quad2d" => Quadrature.gauss_legendre_1d(3)) # use for both directions/dimensions



# is a 5x5 mesh

## mesh connectivity
# IEN(e,a)
nnp = 25 # number of nodes
nez = 5 # sqrt of nnp - # nodes in xdim or ydim
nel = ((nnp^0.5 - 1) |> Int )^2 # number of elements; updated formula
nee = 4 # number of equations per element
IEN = Dict("quad2d" => zeros(Int, nel, nee)) 
# standard ordering
#
IEN["quad2d"][1,:] = [1, 2, 6, 7]
IEN["quad2d"][2,:] = [2, 3, 7, 8]
IEN["quad2d"][3,:] = [3, 4, 8, 9]
IEN["quad2d"][4,:] = [4, 5, 9, 10]
IEN["quad2d"][5,:] = [6, 7, 11, 12]
IEN["quad2d"][6,:] = [7, 8, 12, 13]
IEN["quad2d"][7,:] = [8, 9, 13, 14]
IEN["quad2d"][8,:] = [9, 10, 14, 15]
IEN["quad2d"][9,:] = [11, 12, 16, 17]
IEN["quad2d"][10,:] = [12, 13, 17, 18]
IEN["quad2d"][11,:] = [13, 14, 18, 19]
IEN["quad2d"][12,:] = [14, 15, 19, 20]
IEN["quad2d"][13,:] = [16, 17, 21, 22]
IEN["quad2d"][14,:] = [17, 18, 22, 23]
IEN["quad2d"][15,:] = [18, 19, 23, 24]
IEN["quad2d"][16,:] = [19, 20, 24, 25]

# make the list of node locations - is linearized to comply w IEN
nodes = LinRange(0, 1, nnp) |> collect
# y = LinRange(0, 1, nez) |> collect

## boundary Conditions

# Essential Boundary Conditions: [i,A]
BC_fix_list = zeros(Bool, 1, nnp)

# bc values
BC_g_list = zeros(1, nnp)
# have to write them all individually :(
# unlike MATLAB, cannot write to multiple array locations at once

# reads right then update
# so 1 is bottom left, 5 is bottom right, 21 is top left, 25 is top right

# corners - 0
BC_g_list[1] = 0.0 # first corner
BC_g_list[5] = 0.0 # second corner
BC_g_list[21] = 0.0 # third corner
BC_g_list[25] = 0.0 # fourth corner

# top - 100
BC_g_list[21] = 100.0
BC_g_list[22] = 100.0
BC_g_list[23] = 100.0

# left - 75
BC_g_list[6] = 75.0
BC_g_list[11] = 75.0 
BC_g_list[16] = 75.0

# right - 50
BC_g_list[10] = 50.0
BC_g_list[15] = 50.0
BC_g_list[20] = 50.0

# bottom - 0
BC_g_list[2] = 0.0
BC_g_list[3] = 0.0
BC_g_list[4] = 0.0



mesh = Preprocess.build_mesh(x, [], [], IEN, 1, BC_fix_list, BC_g_list)



## Assemble the global matrices


## Solve the system


# apply essential BCs in q[r2]



# solve




