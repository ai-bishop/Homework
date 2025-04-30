using Pkg
Pkg.activate(".")

include("Preprocess.jl")
include("Quadrature.jl")
include("Grid2D.jl")
include("SimpleVisualization.jl")
using GLMakie
using Printf

quad_rules = Dict("quad" => Quadrature.gauss_legendre_2d(3)) # use for both directions/dimensions



# is a 5x5 mesh

## mesh connectivity
# IEN(e,a)
nnp = 25 # number of nodes`
nez = 5 # sqrt of nnp - 1 - # nodes in xdim or ydim
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



# make the list of node locations
x = LinRange(0, 1, nez) |> collect
append!(x, x)
append!(x, x)
append!(x, LinRange(0, 1, nez) |> collect)
y = zeros(nez) |> collect
append!(y, 1/4 .+ y) 
append!(y, 2/4 .+ y)
append!(y, ones(nez)) 

## boundary Conditions

# Essential Boundary Conditions: [i,A]
BC_fix_list = zeros(Bool, 2, nnp)
BC_fix_list[1,2] = true
BC_fix_list[1,3] = true
BC_fix_list[1,4] = true
BC_fix_list[1,6] = true
BC_fix_list[1,11] = true
BC_fix_list[1,16] = true
BC_fix_list[1,10] = true
BC_fix_list[1,15] = true
BC_fix_list[1,20] = true
BC_fix_list[1,22] = true
BC_fix_list[1,23] = true
BC_fix_list[1,24] = true
BC_fix_list[1,1] = true
BC_fix_list[1,5] = true
BC_fix_list[1,21] = true
BC_fix_list[1,25] = true

# bc values

# have to write them all individually :(
# unlike MATLAB, cannot write to multiple array locations at once

BC_g_list = zeros(1, nnp)

# corners - 0
BC_g_list[1,1] = 37.5 # first corner
BC_g_list[1,5] = 25.0 # second corner
BC_g_list[1,21] = 87.5 # third corner
BC_g_list[1,25] = 75.0 # fourth corner

# top - 100
BC_g_list[1,21] = 100.0
BC_g_list[1,22] = 100.0
BC_g_list[1,23] = 100.0

# left - 75
BC_g_list[1,6] = 75.0
BC_g_list[1,11] = 75.0 
BC_g_list[1,16] = 75.0

# right - 50
BC_g_list[1,10] = 50.0
BC_g_list[1,15] = 50.0
BC_g_list[1,20] = 50.0

# bottom - 0
BC_g_list[1,2] = 0.0
BC_g_list[1,3] = 0.0
BC_g_list[1,4] = 0.0


## Build mesh
mesh = Preprocess.build_mesh(x, y, [], IEN, 1, BC_fix_list, BC_g_list)

## Assemble the global matrix
T_matrix = Grid2D.assemble_stiffness(mesh, quad_rules)

## Solve the system
q = zeros(mesh.nnp * mesh.ned)

# apply essential BCs in q[r2]
idx = findall(BC_fix_list)
q[ mesh.ID[idx] ] = BC_g_list[idx]

# solve
r1 = mesh.free_range
r2 = mesh.freefix_range
q[r1] = T_matrix[r1, r1] \ (  - T_matrix[r1, r2] * q[r2])

# order q
qt = q[r1]
T_interior = [  qt[7] qt[8] qt[9];
                qt[4] qt[5] qt[6];
                qt[1] qt[2] qt[3]]

T_interior_matlab = [   78.5714   76.1161   69.6429;
                        63.1696   56.2500   52.4554;
                        42.8571   33.2589   33.9286]

@printf "Temperature calculated via FEM Code"
display(T_interior)

@printf "" # display empty LinRange

@printf "Temperature calculated via MATLAB code in HW3"
display(T_interior_matlab)