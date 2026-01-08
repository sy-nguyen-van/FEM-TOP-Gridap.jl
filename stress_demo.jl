
# --- 3. Optimization Loop (SIMP Iteration) ---

# Initialization
ρ_e = fill(0.5, nelems) # Elemental densities
ρ_tilde = copy(ρ_e)     # Filtered densities
ρ_prev = copy(ρ_e)
λ_oc = L_oc             # Lagrange multiplier for the volume constraint (used in compliance)
μ_oc = 0.0              # Lagrange multiplier for the stress constraint

println("Starting SIMP Stress-Constrained Optimization...")

for iter in 1:Max_Iter

    # 3.1. FEA: Update and Solve for Displacement (u)
    a(u, v) = ∫(ε(v) ⊙ (C_simp ∘ ρ_tilde) ⊙ ε(u)) * dΩ
    op = AffineFEOperator(a, (v) -> 0, U, V0)
    K, b = get_matrix(op), get_vector(op)
    uh = solve(op, K, f_vec)
    u_dof = get_free_dof_values(uh)

    # 3.2. Stress Aggregation and Constraint (g)
    ε_uh = ε(uh)
    σ_vM_field = σ_vM ∘ ε_uh ∘ ρ_tilde
    σ_vM_vals = get_cell_dof_values(σ_vM_field, Ω)

    # Compute P-Norm Aggregated Stress (σ_P)
    # The sum is over normalized stress to the power P
    p_norm_sum = sum(ρ_tilde[e] * (σ_vM_vals[e] / σ_yield)^P_norm for e in 1:nelems)
    σ_P = p_norm_sum^(1 / P_norm)

    # Constraint value (g <= 0)
    g_val = σ_P - σ_yield

    # 3.3. Sensitivity Analysis (Volume and Stress)

    # Objective Gradient (Volume): d(Vol)/d(rho_e)
    dVol_drho = V_e # Simple elemental volumes

    # Constraint Gradient (Stress): d(g)/d(rho_e) - Requires Adjoint Method

    # i) Adjoint Load Vector (f_adj) - Complex term to derive in Gridap
    # This involves the derivative of the P-norm stress function w.r.t nodal displacements
    # Gridap's automatic differentiation isn't feasible for this high-level structure, 
    # so we define a symbolic function for the adjoint load form (l_adj(v)).

    # Derivative of Von Mises stress w.r.t stress vector [σxx, σyy, τxy]
    function d_von_mises_d_sigma(σ_vec)
        σₓ, σᵧ, τₓᵧ = σ_vec[1], σ_vec[2], σ_vec[3]
        σ_vm = sqrt(σₓ^2 + σᵧ^2 - σₓ * σᵧ + 3 * τₓᵧ^2)
        if σ_vm < 1e-12
            return VectorValue(0.0, 0.0, 0.0)
        end
        # Derivative: [ (2σx - σy)/2σvm, (2σy - σx)/2σvm, 3τxy/σvm ]
        return VectorValue((2 * σₓ - σᵧ) / (2 * σ_vm), (2 * σᵧ - σₓ) / (2 * σ_vm), (3 * τₓᵧ) / σ_vm)
    end

    function d_sigma_P_d_u(ε_v, ρ_e, ε_uh_val)
        D_e = C_simp(ρ_e)
        σ_vec = D_e ⋅ ε_uh_val # Stress vector from current displacement

        # Chain rule terms:
        # 1. d(σ_P)/d(Sum) = (1/P) * (Sum)^(1/P - 1)
        # 2. d(Sum)/d(σ_vm_e) = ρ_e * P * (σ_vm_e / σ_yield)^(P-1) * (1/σ_yield)
        # 3. d(σ_vm_e)/d(σ_vec) = d_von_mises_d_sigma(σ_vec)
        # 4. d(σ_vec)/d(ε_v) = D_e

        # Combine 1 & 2:
        # coeff = (Sum)^(1/P - 1) * ρ_e * (σ_vm_e / σ_yield)^(P-1) * (1/σ_yield)
        # Note: p_norm_sum is the Sum. σ_P = Sum^(1/P).
        # So coeff = σ_P^(1-P) * ρ_e * ...

        σ_vm_e = σ_vM(ε_uh_val, ρ_e)

        term1 = (p_norm_sum)^(1 / P_norm - 1)
        term2 = ρ_e * P_norm * (σ_vm_e / σ_yield)^(P_norm - 1) * (1.0 / σ_yield)

        d_vm_d_sigma = d_von_mises_d_sigma(σ_vec)

        # Adjoint load integrand: coeff * (d_vm_d_sigma ⋅ (D_e ⋅ ε_v))
        return (term1 * term2) * (d_vm_d_sigma ⋅ (D_e ⋅ ε_v))
    end

    # Apply the adjoint load as a weak form (linear form l_adj(v))
    l_adj(v) = ∫(v ⋅ (d_sigma_P_d_u ∘ (ε(v), ρ_tilde, ε(uh)))) * dΩ

    # ii) Solve Adjoint Equation (K * u_adj = f_adj)
    # The stiffness matrix K is the same as for the primal problem
    op_adj = AffineFEOperator(a, l_adj, U, V0)
    K_adj = get_matrix(op_adj)
    f_adj = get_vector(op_adj) # This vector now contains the adjoint load from l_adj(v)

    # Solve for adjoint displacement
    u_adj = solve(op_adj, K_adj, f_adj)
    u_adj_dof = get_free_dof_values(u_adj)

    # iii) Compute the Final Stress Sensitivity d(g)/d(rho_e)
    dG_drho = zeros(nelems)

    # Function to compute elemental stiffness derivative (dK_e/drho_e)
    function d_C_d_rho(ρ_e)
        dE_drho = penal * ρ_e^(penal - 1) * (E₀ - ρ_min^penal * E₀)
        # Derivative of D matrix w.r.t E_eff
        D_prime = zeros(3, 3)
        D_prime[1, 1] = D_prime[2, 2] = 1 / (1 - ν^2)
        D_prime[1, 2] = D_prime[2, 1] = ν / (1 - ν^2)
        D_prime[3, 3] = 1 / (2 * (1 + ν))
        return dE_drho * D_prime
    end

    # Compute the implicit term integral (u_adj^T * dK/drho_e * u)
    # This must be done element by element and requires elemental interpolation
    dK_drho_field = d_C_d_rho ∘ ρ_tilde

    for e in 1:nelems
        # Simplified integration of the implicit term: ∫ (ε(u_adj) : dD/drho : ε(u)) dΩ
        implicit_term = sum(∫(ε(u_adj) ⊙ dK_drho_field ⊙ ε(uh)) * dΩ(CellMeasure([e])))

        # Explicit term (related to the derivative of the P-norm stress w.r.t. rho_e)
        # This is a small term related to the density filter in the P-norm sum.
        explicit_term = (1 / P_norm) * (σ_P / p_norm_sum) * (σ_vM_vals[e] / σ_yield)^P_norm


        # The full sensitivity: dG_drho = Explicit - Implicit
        dG_drho[e] = explicit_term - implicit_term
    end

    # 3.4. Optimization Update (Modified OC for Inequality Constraints)
    # Use a modified OC scheme to handle the two objectives (min Vol) and constraint (max Stress).

    # Only update the Lagrange multiplier μ_oc for the stress constraint.
    # We fix λ_oc (the volume multiplier) or use a basic search.

    # Update μ_oc based on the KKT condition for the stress constraint
    if g_val > 0 || μ_oc > 1e-6 # If constraint is violated or active
        μ_oc = μ_oc + 0.1 * g_val # Simple gradient ascent update for Lagrange multiplier
        μ_oc = max(0.0, μ_oc)
    end

    # The combined sensitivity (Objective + Constraint)
    # L_total = Vol + μ * g
    # dL/drho = dVol/drho + μ * dG/drho

    # --- Optimality Criteria Factor (B) ---
    # We solve for the KKT stationarity condition dL/drho = 0, leading to:
    # dVol/drho + μ * dG/drho = λ * rho^(p-1) * sensitivity_factor (for compliance, NOT here)
    # The simplest update is a direct search on a fictitious multiplier.

    l1, l2 = 0.0, 1e9 # Fictitious multiplier range

    function update_rho_stress_based(l1, l2, dVol_drho, dG_drho, μ_oc, ρ_old)
        while (l2 - l1) / (l1 + l2) > 1e-4
            λ = 0.5 * (l1 + l2) # Fictitious multiplier (to satisfy volume *change* limit)

            ρ_new = zeros(nelems)
            for e in 1:nelems
                # The sensitivity ratio (B) balances dVol, dG, and the fictitious multiplier λ
                B = (-dG_drho[e] * μ_oc) / (dVol_drho[e] * λ)

                # Update formula (similar to compliance, but using dG/drho)
                ρ_e_new = max(ρ_old[e] * (1 - move_limit), min(1.0, ρ_old[e] * (1 + move_limit), B^(1 / penal)))
                ρ_new[e] = max(ρ_min, ρ_e_new)
            end

            # Check if the volume reduction is "too much" (i.e., fictitious constraint is active)
            # Since this is volume minimization, we just check total volume.
            current_V = sum(ρ_new .* V_e)
            if current_V > sum(ρ_old .* V_e) # If we are increasing volume (bad)
                l2 = λ
            else
                l1 = λ
            end
        end
        return ρ_new
    end

    ρ_new = update_rho_stress_based(l1, l2, dVol_drho, dG_drho, μ_oc, ρ_e)

    # 3.5. Update and Check Convergence
    change = norm(ρ_new - ρ_prev, Inf) / norm(ρ_new, Inf)
    ρ_prev = ρ_e
    ρ_e = ρ_new
    ρ_tilde = ρ_e # No filter applied here

    Current_Vol = sum(ρ_e .* V_e)

    println("Iter: $iter, Volume: $(Current_Vol), Max Stress: $(σ_P), Constraint (g): $(g_val), Change: $(change), Mu: $(μ_oc)")

    if iter % 10 == 0
        writevtk(Ω, "cantilever_stress_simp_iter_$iter", cellfields=["density" => ρ_e, "displacement" => uh, "stress_vm" => σ_vM_vals])
    end

    if change < 0.001
        println("\nStress-constrained SIMP converged in $iter iterations.")
        break
    end
end