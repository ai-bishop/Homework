## simple bar

##
using Pkg
Pkg.activate(".")

include("Preprocess.jl")
include("Quadrature.jl")
include("Bar1D.jl")
include("SimpleVisualization.jl")
using GLMakie

quad_rules = Dict("parb" => Quadrature.gauss_legendre_1d(3)) # changed input

A = 1.0 # cross-sectional area [L^2]
L = 1.0 # length [L]
ρ = 1.0 # mass per length [kg/L^3]
Ω = 1.0 # angular velocity [rad/s]

## Properties of the uniform body
prop = (E=1.0, # Young's Modulus [F/m^2]
    A=A, # cross-sectional area [L^2]
    L=L, # length [L]
    ρ = ρ, # mass per length [kg/L^3]
    Ω = Ω, # angular velocity [rad/s]
    f0 = ρ * A * Ω^2 # loading [F/L] - change to f = ρAΩ^2 x
    # syntax? how define as a changing variable?
) 

## mesh connectivity
# IEN(e,a)
nnp = 5 # number of nodes
nel = (nnp - 1) / 2 |> Int # number of elements
nee = 3 # number of equations per element
IEN = Dict("parb" => zeros(Int, nel, nee)) 
# same ordering as from class 
IEN["parb"][1,:] = [1, 3, 2] 
# IEN["parb"][2,:] = [3, 4, 5] 
IEN["parb"][2,:] = [3, 5, 4] 

# make the list of node locations
x = LinRange(0, prop.L, nnp) |> collect

## Essential Boundary Conditions: [i,A]
BC_fix_list = zeros(Bool, 1, nnp)
BC_fix_list[1, 1] = true

# BC values
BC_g_list = zeros(1, nnp)
BC_g_list[1, 1] = 0.0

## Loading force
f(x) = prop.f0 * x

## Build mesh
m = Preprocess.build_mesh(x, [], [], IEN, 1, BC_fix_list, BC_g_list)

## Assemble the global matrices
K = Bar1D.assemble_stiffness(m, prop, quad_rules)
F = Bar1D.assemble_rhs(m, f, quad_rules)

## Solve the system
q = zeros(m.nnp * m.ned)

# apply essential BCs in q[r2]
idx = findall(BC_fix_list)
q[ m.ID[idx] ] = BC_g_list[idx]  

# solve
r1 = m.free_range
r2 = m.freefix_range
q[r1] = K[r1, r1] \ (F[r1] - K[r1, r2] * q[r2])

## Plot the result
fig, ax = SimpleVisualization.init_plot()
SimpleVisualization.draw_element(ax, m, q, linewidth=3, linestyle=:dash, color=:red)
u_true(x) = prop.f0^2 * x^3 / (6 * A)
lines!(ax, LinRange(0, prop.L, 20), u_true)
fig

## Plot the stress
σ_true(x) = prop.f0 * x /prop.A
fig, ax = SimpleVisualization.init_plot(ylabel="Stress")
SimpleVisualization.draw_element_stress(ax, m, q, prop)
lines!(ax, LinRange(0, prop.L, 20), σ_true)
fig