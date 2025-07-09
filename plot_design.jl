using Base.Iterators
function plot_design(ρ_new, cell_coords, iter, L, W)
    # Range of limit of colors
    range_limits = (0, 1)
    # Polygon for FE MESH
    shape_polygon = [Point2f[(cell_coords[e][1][1], cell_coords[e][1][2]), # Node 1
        (cell_coords[e][2][1], cell_coords[e][2][2]), # Node 2
        (cell_coords[e][4][1], cell_coords[e][4][2]),
        (cell_coords[e][3][1], cell_coords[e][3][2])] for e in 1:no_ele]
    fontsize = 16
    # ==========================
    fig1 = Figure(fontsize=16)
    ax1 = Axis(fig1[1, 1], title="TOPOLOGY OPTIMIZATION; Iteration = " *string(iter), xlabelsize=fontsize, ylabelsize=fontsize, xlabel="X", ylabel="Y", aspect=DataAspect())

    Makie.poly!(ax1, shape_polygon, color=ρ_new, colormap=:jet)

    Makie.Colorbar(fig1[1, 2], colormap=:jet, colorrange=range_limits)
    Label(fig1[1, 2, Top()], L"$\rho$", padding=(0, 0, 2, 0), fontsize=18)
    xlims!(ax1, 0, L)
    ylims!(ax1, 0, W)

    colsize!(fig1.layout, 1, Aspect(1, L / W))
    Makie.resize_to_layout!(fig1)
    return fig1
end