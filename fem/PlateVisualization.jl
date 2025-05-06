module PlateVisualization

include("Shapefunctions.jl")
using .Shapefunctions

using GLMakie
using LinearAlgebra

function init_plot(; xlabel="Position x", ylabel="Displacement u")
    fig = Figure()
    ax = Axis(fig[1, 1], autolimitaspect=1,
        xlabel=xlabel,
        ylabel=ylabel)

    return fig, ax
end

function draw_element(ax, mesh, q; linestyle=nothing, linewidth=2,
    color=GLMakie.Makie.wong_colors()[1])

    # plot cells
    for (element_type, ien) in mesh.IEN
        if element_type == "quad"
            nel = size(ien, 1)
            N = Shapefunctions.shapefunc(element_type)
            ξ = LinRange(-1, 1, 5)
            η = LinRange(-1, 1, 5)


            for e in 1:nel
                A = ien[e, :]
                xe = @view mesh.x[A]
                ye = @view mesh.y[A]
                idx = mesh.ID[1, A]
                qe = @view q[idx]
                vertices = Point2f[(x,y) for (x,y) in zip(xe, ye)]


                X = zeros(length(ξ))
                Y = zeros(length(ξ))
                u = zeros(length(ξ))
                for (i, ξ, η) in enumerate(ξ)
                    NN, Nξ, Nη = N(ξ, η)
                    X[i] = dot(NN, xe)
                    Y[i] = dot(NN, ye)
                    u[i] =  dot(NN, qe)
                end
            end

                contour(ax, X, Y, u, color=color,
                    linestyle=linestyle,
                    linewidth=linewidth)

            


        elseif element_type == "quad2d"


            nel = size(ien, 1)
            N = Shapefunctions.shapefunc(element_type)
            ξ = LinRange(-1, 1, 9)
            η = LinRange(-1, 1, 9)


            for e in 1:nel
                A = ien[e, :]
                xe = @view mesh.x[A]
                ye = @view mesh.y[A]
                idx = mesh.ID[1, A]
                qe = @view q[idx]
                vertices = Point2f[(x,y) for (x,y) in zip(xe, ye)]


                X = zeros(length(ξ), length(η))
                Y = zeros(length(ξ), length(η))
                u = zeros(length(ξ), length(η))
                for (i, ξ) in enumerate(ξ)
                    for (j,η) in enumerate(η)
                    NN, Nξ, Nη = N(ξ, η)
                    X[i,j] = dot(NN, xe)
                    Y[i,j] = dot(NN, ye)
                    u[i,j] =  dot(NN, qe)
                    end
                end

                contour(X, Y, u, color=color,
                    linestyle=linestyle,
                    linewidth=linewidth)
            end
        end
    end
end



function draw_element_stress(ax, mesh, q, properties;
    linestyle=nothing, linewidth=2,
    color=GLMakie.Makie.wong_colors()[1])

    # plot cells
    for (element_type, ien) in mesh.IEN
        if element_type == "line" || element_type == "parb"
            nel = size(ien, 1)
            N = Shapefunctions.shapefunc(element_type)
            ξ = LinRange(-1, 1, 2)
            for e in 1:nel
                A = ien[e, :]
                xe = @view mesh.x[A]
                idx = mesh.ID[1, A]
                qe = @view q[idx]
                X = zeros(length(ξ))
                σ = zeros(length(ξ))
                for (i, ξ) in enumerate(ξ)
                    NN, Nξ = N(ξ)
                    X[i] = dot(NN, xe)
                    σ[i] = properties.E * dot(Nξ, qe) / dot(Nξ, xe)
                end

                lines!(ax, X, σ, color=color,
                    linestyle=linestyle,
                    linewidth=linewidth)

            end
            
        elseif element_type == "quad2d"
            nel = size(ien, 1)
            N = Shapefunctions.shapefunc(element_type)
            ξ = LinRange(-1, 1, 5)
            η = LinRange(-1, 1, 5)
            for e in 1:nel
                A = ien[e, :]
                xe = @view mesh.x[A]
                ye = @view mesh.y[A]
                idx = mesh.ID[1, A]
                qe = @view q[idx]
                X = zeros(length(ξ))
                Y = zeros(length(η))
                σ = zeros(length(ξ), length(η))
                for (i, ξ) in enumerate(ξ)
                    for (j, η) in enumerate(η)
                        NN, Nξ, Nη = N(ξ,η)
                    X[i] = dot(NN, xe)
                    Y[i] = dot(NN, ye)
                    σ[i] = properties.E * (dot(Nξ, qe) / dot(Nξ, xe)) * (dot(Nη, qe) / dot(Nη, ye))
                    end
                end

                # check syntax - may not be correct
                # check abt contour plot?
                lines!(ax, X, Y, σ, color=color,
                    linestyle=linestyle,
                    linewidth=linewidth) 

            end            
        end
    end
end

end
