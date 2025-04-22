using Pkg
Pkg.activate(".")

include("Preprocess.jl")
include("Quadrature.jl")
include("Grid2D.jl")
include("SimpleVisualization.jl")
using GLMakie

quad_rules = Dict("quad2d" => Quadrature.gauss_legendre_2d(3,3)) # use for both directions/dimensions



# is a 5x5 mesh

## mesh connectivity
# IEN(e,a)
nnp = 25 # number of nodes
nez = 5 # sqrt of nnp - # nodes in xdim or ydim
nel = ((nnp^0.5 - 2) |> Int )^2 # number of elements; updated formula
nee = 9 # number of equations per element
IEN = Dict("quad2d" => zeros(Int, nel, nee)) 
# standard ordering
#
IEN["quad2d"][1,:] = [1, 3, 13, 11, 2, 8, 12, 6, 7]
IEN["quad2d"][2,:] = [2, 4, 14, 12, 3, 9, 13, 7, 8]
IEN["quad2d"][3,:] = [3, 5, 15, 13, 4, 10, 14, 8, 9]

IEN["quad2d"][4,:] = [6, 8, 18, 16, 7, 13, 17, 11, 12]
IEN["quad2d"][5,:] = [7, 9, 19, 17, 8, 14, 18, 12, 13]
IEN["quad2d"][6,:] = [8, 10, 20, 18, 9, 15, 19, 13, 14]

IEN["quad2d"][7,:] = [11, 13, 23, 21, 12, 18, 22, 16, 17]
IEN["quad2d"][8,:] = [12, 14, 24, 22, 13, 19, 23, 17, 18]
IEN["quad2d"][9,:] = [13, 15, 25, 23, 14, 20, 24, 18, 19]




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


## Build mesh
mesh = Preprocess.build_mesh(nodes, [], [], IEN, 1, BC_fix_list, BC_g_list)

## Assemble the global matrices
K = Grid2D.assemble_stiffness(m, prop, quad_rules)
F = Grid2D.assemble_rhs(m, ___, quad_rules)









## Solve the system


# apply essential BCs in q[r2]



# solve




