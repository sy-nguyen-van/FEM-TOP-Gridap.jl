include("fd_check_cost.jl")
# include("fd_check_constraint.jl")
function run_finite_difference_check(type_FEA, dv, fd_step_size, H, p, k, q_in, s, Ω, dΩ, dΓ_flux, Uₕ, Vₕ)
    # This function performs a finite difference check of the analytical
    path_output = "FD_check/"
    # sensitivities of the cost and/or constraint functions by invoking the
    # corresponding routines.
    grad_cost_0, grad_cost_i = fd_check_cost(dv, fd_step_size, H, p, k, q_in, s, Ω, dΩ, dΓ_flux, Uₕ, Vₕ)
    # # Vector of index of dv
    n_dv = num_cells(Ω) # Number of elements
    x = LinRange(1, n_dv, n_dv)
    if type_FEA == "Thermal"
        Name_FEA = "Thermal Comliance"
    else
        Name_FEA = "Comliance"
    end
    # # Name of plots
    name_obj = "FD of cost function: " * Name_FEA
    # # --
    fig_FD1 = Figure()
    ax_FD1 = Axis(fig_FD1[1, 1], xlabel="Design variable: v", ylabel="dz/dv", title=name_obj)
    Makie.lines!(ax_FD1, x, vec(grad_cost_0), color=:red, label="Analytical")
    Makie.scatter!(ax_FD1, x, vec(grad_cost_i), color=:blue, label="FD")
    axislegend(position=:rb)
    figs = fig_FD1
    # ------
    Makie.save(path_output * "FD_cost_" * Name_FEA * ".png", fig_FD1)
    return figs
end