include("nonlcon.jl")
function fd_check_constraint(OPT::OPT_struct, FE::FE_struct,SGs::SGs_struct)
    #
    # This function performs a finite difference check of the sensitivities of
    # the COST function with respect to the bar design variables.
    #
    # ===============================
    # FINITE DIFFERENCE SENSITIVITIES
    # ===============================
    n_dv = OPT.n_dv
    grad_theta_i = zeros(n_dv, 1)

    fd_step = OPT.fd_step_size

    max_error = 0.0
    max_rel_error = 0.0
    max_error_dv = 0
    max_rel_error_dv = 0
    dv_0 = copy(OPT.dv)
    dv_i = copy(OPT.dv)

    theta_0, grad_theta_0, _, _ = nonlcon(OPT, FE, SGs, dv_0)
    # Finite differences
    println("Computing finite difference sensitivities of constraint...")
    # Do this for all design variables or only a few
    up_to_dv = n_dv

    for i in 1:up_to_dv
        dv_i[i] = dv_0[i] + fd_step
        theta_i, _, _, _ = nonlcon(OPT, FE, SGs, dv_i)
        grad_theta_i[i] = (theta_i- theta_0)/ fd_step
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
    OPT.dv = dv_0
    @printf("Max. ABSOLUTE error is:  %.5e\n", max_error)
    @printf("It occurs at variable of  %d\n", max_error_dv)

    @printf("Max. RELATIVE error is:  %.5e\n", max_rel_error)
    @printf("It occurs at variable of  %d\n", max_rel_error_dv)
    #------------------------------------
    return grad_theta_0, grad_theta_i
end