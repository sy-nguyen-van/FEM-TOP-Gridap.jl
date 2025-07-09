
function obj(dv, H, p, k, q_in, s, Ω, dΩ, dΓ_flux, Uₕ, Vₕ)

    dv_fil = H * dv # Filter
    # SIMP
    ρ_SIMP, dρ_dp = SIMP_Inter(dv_fil, p, Ω)
    # ===FE-ANALYSIS
    a(u,v) = ∫(ρ_SIMP * k * ∇(v) ⋅ ∇(u) ) * dΩ
    l(v) = ∫( -1*q_in * v ) * dΓ_flux + ∫( s * v ) * dΩ
    # Solution of the FE problem
    op = AffineFEOperator(a, l, Uₕ, Vₕ)
    Th = solve(op)
    C = sum(∫( ρ_SIMP *k * ∇(Th) ⋅ ∇(Th) ) * dΩ) # Compliance
    dC_fil = get_contribution(∫( -dρ_dp*k * ∇(Th) ⋅ ∇(Th) ) * dΩ, Ω) # Sensitivity of compliance w.r.t rho
    dC = transpose(H) * (dC_fil[:])
    return C, dC
end