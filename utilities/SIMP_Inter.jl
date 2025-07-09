# This function is for The SIMP (Solid Isotropic Material with Penalization)
# Reference: Topology optimization: theory, methods, and applications, MP Bendsoe, O Sigmund - 2013
function SIMP_Inter(ρ, p, Ω)
    ρ_SIMP = CellField(ρ .^ p, Ω) # SIMP: 
    dρ_dp = p * CellField(ρ .^ (p - 1), Ω) # SIMP: 
    return ρ_SIMP, dρ_dp
end