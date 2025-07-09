include("obj.jl")
function fd_check_cost(dv, fd_step_size, H, p, k, q_in, s, Ω, dΩ, dΓ_flux, Uₕ, Vₕ)
    #
    # This function performs a finite difference check of the sensitivities of
    # the COST function with respect to the bar design variables.
    #
    # ===============================
    # FINITE DIFFERENCE SENSITIVITIES
    # ===============================
    num_ele = num_cells(Ω) # Number of elements
    n_dv = num_ele
    grad_theta_i = zeros(n_dv, 1)

    fd_step = fd_step_size

    max_error = 0.0
    max_rel_error = 0.0
    max_error_dv = 0
    max_rel_error_dv = 0
    dv_0 = copy(dv)
    dv_i = copy(dv)

    theta_0, grad_theta_0 = obj(dv, H, p, k, q_in, s, Ω, dΩ, dΓ_flux, Uₕ, Vₕ)
    # Finite differences
    println("Computing finite difference sensitivities of cost...")
    # Do this for all design variables or only a few
    up_to_dv = n_dv

    for i in 1:up_to_dv
        dv_i[i] = dv_0[i] + fd_step
        theta_i, _,  = obj(dv_i, H, p, k, q_in, s, Ω, dΩ, dΓ_flux, Uₕ, Vₕ)
        grad_theta_i[i] = (theta_i - theta_0) / fd_step
        error = grad_theta_0[i] - grad_theta_i[i]

        if abs(error) > abs(max_error)
            max_error = error
            max_error_dv = i
        end
        rel_error = error / theta_0
        if abs(rel_error) > abs(max_rel_error)
            max_rel_error = rel_error
            max_rel_error_dv = i
        end
        dv_i = copy(dv_0)

    end
    dv = dv_0
    @printf("Max. ABSOLUTE error is:  %.5e\n", max_error)
    @printf("It occurs at variable of  %d\n", max_error_dv)

    @printf("Max. RELATIVE error is:  %.5e\n", max_rel_error)
    @printf("It occurs at variable of  %d\n", max_rel_error_dv)

    #------------------------------------
    return grad_theta_0, grad_theta_i
end