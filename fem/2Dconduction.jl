using Pkg
Pkg.activate(".")

include("Preprocess.jl")
include("Quadrature.jl")
include("Grid2D.jl")
include("SimpleVisualization.jl")
using GLMakie

quad_rules = Dict("quad" => Quadrature.gauss_legendre_2d(3,3)) # use for both directions/dimensions



# is a 5x5 mesh

## mesh connectivity
# IEN(e,a)
nnp = 25 # number of nodes
nez = 5 # sqrt of nnp - # nodes in xdim or ydim
nel = ((nnp^0.5 - 1) |> Int )^2 # number of elements; updated formula
nee = 4 # number of equations per element
IEN = Dict("quad" => zeros(Int, nel, nee)) 
# standard ordering
#
IEN["quad"][1,:] = [1, 2, 7, 6]
IEN["quad"][2,:] = [2, 3, 8, 7]
IEN["quad"][3,:] = [3, 4, 9, 8]
IEN["quad"][4,:] = [4, 5, 10, 9]


IEN["quad"][5,:] = [6, 7, 12, 11]
IEN["quad"][6,:] = [7, 8, 13, 12]
IEN["quad"][7,:] = [8, 9, 14, 13]
IEN["quad"][8,:] = [9, 10, 15, 14]


IEN["quad"][9,:] = [11, 12, 17, 16]
IEN["quad"][10,:] = [12, 13, 18, 17]
IEN["quad"][11,:] = [13, 14, 19, 18]
IEN["quad"][12,:] = [14, 15, 20, 19]

IEN["quad"][13,:] = [16, 17, 22, 21]
IEN["quad"][14,:] = [17, 18, 23, 22]
IEN["quad"][15,:] = [18, 19, 24, 23]
IEN["quad"][16,:] = [19, 20, 25, 24]



# make the list of node locations - is linearized to comply w IEN
nodes = LinRange(0, 1, nnp) |> collect
xnodes = LinRange(0, 1, nez) |> collect
ynodes = LinRange(0, 1, nez) |> collect

## boundary Conditions

# Essential Boundary Conditions: [i,A]
BC_fix_list = zeros(Bool, 1, nnp)
BC_fix_list2 = zeros(Bool, nez, nez)
BC_fix_list2[1,2] = true
BC_fix_list2[1,3] = true
BC_fix_list2[1,4] = true
BC_fix_list2[2,1] = true
BC_fix_list2[3,1] = true
BC_fix_list2[4,1] = true
BC_fix_list2[5,2] = true
BC_fix_list2[5,3] = true
BC_fix_list2[5,4] = true
BC_fix_list2[2,5] = true
BC_fix_list2[3,5] = true
BC_fix_list2[4,5] = true

# bc values
BC_g_list2 = zeros(nez, nez)

# have to write them all individually :(
# unlike MATLAB, cannot write to multiple array locations at once


# corners, bottom - already set

# top - 100
BC_g_list2[1,2] = 100.0
BC_g_list2[1,3] = 100.0
BC_g_list2[1,4] = 100.0

# left - 75
BC_g_list2[2,1] = 75.0
BC_g_list2[3,1] = 75.0
BC_g_list2[4,1] = 75.0

# right - 50
BC_g_list2[2,5] = 50.0
BC_g_list2[3,5] = 50.0
BC_g_list2[4,5] = 50.0


# reads right then up
# so 1 is bottom left, 5 is bottom right, 21 is top left, 25 is top right

BC_g_list = zeros(1, nnp)

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
mesh = Preprocess.build_mesh(xnodes, ynodes, [], IEN, 1, BC_fix_list2, BC_g_list2)

## Assemble the global matrix
T_matrix = Grid2D.assemble_stiffness(mesh, quad_rules)

## Solve the system
q = zeros(m.nnp * m.ned)

 

# apply essential BCs in q[r2]
 




# K = stiffness

# solve
r1 = mesh.free_range
r2 = mesh.freefix_range
q[r1] = T_matrix[r1, r1] \ (T_matrix[r1, r2] * q[r2])



