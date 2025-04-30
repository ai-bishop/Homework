## simple plate

using Pkg
Pkg.activate(".")

include("Preprocess.jl")
include("Quadrature.jl")
include("Plate2D.jl")
include("SimpleVisualization.jl")
using GLMakie

# initialize dimensions - # nodes each dimension
# dimensions *must* be odd
m = 7 # xdim
if iseven(m)
    m = m-1
    printf("x dimension reduced by one for compatibility with meshing code")
end
n = 9 # ydim
if iseven(n)
    n = n-1
    printf("y dimension reduced by one for compatibility with meshing code")
end

quad_rules = Dict("quad2d" => Quadrature.gauss_legendre_2dim(m,n)) # changed input

# define properties
E = 1.0 # young's modulus, [F/m2]
T = 0.001 # thickness of plate [m]
L = 3.0 # length of plate [m]
W = 2.0 # width of plate [m]
ρ = 1, # density of plate [kg/m3]
grav = 9.81 # m/s^2

# insert parts into prop
prop = (
    E = E, # young's modulus, [F/m2]
    T = T, # thickness of plate [m]
    L = L, # length of plate [m]
    W = W, # width of plate [m]
    ρ = ρ, # density of plate [kg/m3]
    f0 = ρ * grav * T # loading force = density * volume * grav

)

# VERIFY NEL

## mesh connectivity
# IEN(e,a)
nnp = m * n # number of nodes
nel = ((m-1)/2 |> Int) * ((n-1)/2 |> Int) # number of elements; updated formula
nee = 9 # number of equations per element
IEN = Dict("quad2d" => zeros(Int, nel, nee)) 

# create IEN
# example for quad:
# IEN["quad"][1,:] = [1, 2, 7, 6]


# indexing correct
# need node creation
# node creation correct?
# if 3x3
# 4 7 3 
# 8 9 6
# 1 5 2
for ey in 1:((n-1)/2 |> Int) # ydim
    for ex in 1:((m-1)/2 |> Int) # xdim

# 2*ex - 1 + 2*ey, 2*ex + 2*ey, 2*ex+1 + 2*ey
# 2*ex - 1 + ey, 2*ex + ey, 2*ex + 1 + ey
# 2*ex - 1, 2*ex, 2*ex + 1
        IEN["quad2d"][(ex + ex * (ey - 1)),:] = [  
        2*ex - 1,       2*ex + 1,       2*ex+1 + 2*ey,
        2*ex - 1,       2*ex,           2*ex + 1 + ey,
        2*ex + 2*ey,    2*ex - 1 + ey,  2*ex + ey
        ]

    end
end

# make the list of node locations
x = LinRange(0, 1, m) |> collect
for zed in 2:n
    append!(x, LinRange(0, 1, m) |> collect)
end

y = zeros(n)
for zed in 2:m
    append!(y, (zed-1)/(m-1) .+ zeros(n))

end

# make Plate2D ref both dims

## boundary Conditions
# Essential Boundary Conditions: [i,A]
BC_fix_list = zeros(Bool, 2, nnp)

# ex of setting condition true:
# BC_fix_list[1,2] = true
# BC_fix_list[1,3] = true
# BC_fix_list[1,4] = true

BC_fix_list[1, 1] = true # bottom left corner only corner with fixed displacement



BC_g_list = zeros(1, nnp)
BC_g_list[1,1] = 0 # fixed deflection of 0

## Loading force
f(x) = prop.f0 * x * y

## Build mesh
mesh = Preprocess.build_mesh(x, y, [], IEN, 1, BC_fix_list, BC_g_list)

## Assemble the global matrix
K_global = Plate2D.assemble_stiffness(mesh, quad_rules)
F_global = Plate1D.assemble_rhs(mesh, f, quad_rules)
















