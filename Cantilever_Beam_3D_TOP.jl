using Gridap
using Gridap
using Gridap.CellData
using Gridap.Visualization
function OC(x, volfrac, dc, nUm_elements)
    x_min = 0.001
    l1 = 0.0
    l2 = 100000.0
    move = 0.2
    xnew = copy(x)
    while (l2 - l1 > 1e-4)
        lmid = 0.5 * (l2 + l1)
        x_candidate = x .* sqrt.(-dc ./ lmid)
        xnew = max.(x_min, max.(x .- move, min.(1.0, min.(x .+ move, x_candidate))))
        if sUm(xnew) - volfrac * nUm_elements > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end
    return xnew
end
function SIMP_Inter(ρ₀,p,Ω)
    ρ_SIMP= CellField(ρ₀.^p, Ω) # SIMP: 
    dρ_dp= p*CellField(ρ₀.^(p-1), Ω) # SIMP: 
    return ρ_SIMP, dρ_dp
end
# ========
L = 120 # Length
W = 24 # Width and Height
domain = (0, L, 0, W, 0, W) # sizes along X, Y, Z
partition = (120, 24, 24) # mesh sizes
model = CartesianDiscreteModel(domain, partition)
# Tag the left face (x = 0); right face (x = L)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"left", [1,3,5,7,13,15,17,19,25]) # left face: corners (1,3,5,7); edges (13,15,17,19); others (25) 
add_tag_from_tags!(labels, "right", [2,4,6,8,14,16,18,20,26]) # Right face: corners (2,4,6,8); edges (14,16,18,20); others (26) 
# Integrating on the domain Ω
degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
# neUmann boUndaries
Γ = BoundaryTriangulation(model,tags = "right")
dΓ = Measure(Γ,degree)
# Vector-valUed FE space
order = 1
reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)
V = TestFESpace(Ω, reffe; conformity=:H1, dirichlet_tags = "left")
U = TrialFESpace(V)
# ConstitUtive law for LINEAR ELASTICITY 
const E = 119e3 # YoUng's modUlUs
const ν = 0.3   # Poisson's ratio
λ = (E*ν)/((1+ν)*(1-2*ν)) # Lamé parameters
μ = E/(2*(1+ν)) # Lamé parameters
ε₀(U) = 0.5 * ( ∇(U) + transpose(∇(U))) # Strain; In Gridap this can be aUtomatically defined in Gridap.ε
σ(ε₀) = (λ * tr(ε₀) * one(ε₀) + 2μ * ε₀) # Stress
# The weak form
p = 3
volfrac = 0.3
num_elements = num_cells(Ω)
ρ₀ =rand(num_elements)
ρ_SIMP= CellField(ρ₀.^p, Ω) # SIMP: 
# SIMP: Bendsøe, Martin P., and Ole SigmUnd. "Material interpolation schemes in topology optimization." Archive of applied mechanics 69 (1999): 635-654.
# ====================================
F(x) = VectorValue(0.0, -1, 0.0) # I want to apply a Unit load along Y axis to each node on the right face
a(U,v) = ∫( ρ_SIMP * ( (σ∘ε₀(U)) ⊙ ε₀(v) ) ) * dΩ # ; (∘) Composite fUnctions
l(v) = ∫(F ⋅ v) * dΓ # Right-hand size
# SolUtion of the FE problem
op = AffineFEOperator(a,l,U,V)
uh = solve(op)
# Export the displacement, strain, stress
writevtk(Ω,"Cantilever_Beam_3D_Uh", cellfields=["rho"=>ρ_SIMP,"Uh"=>uh,"epsi"=>ε₀(uh),"sigma"=>σ∘ε₀(uh)])

C = ρ_SIMP * ( (σ∘ε₀(uh)) ⊙ ε₀(uh) )
get_free_dof_values(uh)

# ρ_dC = CellField(ρ₀.^(p-1), Ω) # SIMP: 
# dC =  -p*ρ_dC * ( (σ∘ε₀(Uh)) ⊙ ε₀(Uh) ) 


