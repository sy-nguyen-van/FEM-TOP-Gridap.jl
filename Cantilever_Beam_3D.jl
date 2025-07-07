using Gridap
L = 100 # Length
W = 24 # Width
domain = (0, L, 0, W, 0, W) # sizes along X, Y, Z
partition = (100, 24, 24) # mesh sizes
model = CartesianDiscreteModel(domain, partition)
# Tag the left face (x = 0); right face (x = L)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"left", [1,3,5,7,13,15,17,19,25]) # left face: corners (1,3,5,7); edges (13,15,17,19); others (25) 
add_tag_from_tags!(labels, "right", [2,4,6,8,14,16,18,20,26]) # Right face: corners (2,4,6,8); edges (14,16,18,20); others (26) 
# Integrating on the domain Ω
degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω,degree)
# neumann boundaries
Γ = BoundaryTriangulation(model,tags = "right")
dΓ = Measure(Γ,degree)
# Vector-valued FE space
order = 1
reffe = ReferenceFE(lagrangian, VectorValue{3,Float64}, order)
Vₕ = TestFESpace(Ω, reffe; conformity=:H1, dirichlet_tags = "left")
Uₕ = TrialFESpace(Vₕ) 
# Constitutive law for LINEAR ELASTICITY 
const E = 119e3 # Young's modulus
const ν = 0.3   # Poisson's ratio
const λ = (E*ν)/((1+ν)*(1-2*ν)) # Lamé parameters
const μ = E/(2*(1+ν)) # Lamé parameters
ε₀(u) = 0.5 * ( ∇(u) + transpose(∇(u))) # Strain; In Gridap this can be automatically defined in Gridap.ε
σ(ε₀) = λ * tr(ε₀) * one(ε₀) + 2μ * ε₀ # Stress
# The weak form
F(x) = VectorValue(0.0, -1, 0.0) # I want to apply a unit load along Y axis to each node on the right face
a(u,v) = ∫((σ∘ε₀(u)) ⊙ ε₀(v)) *dΩ # Right-hand size; (∘) Composite functions
l(v) = ∫(F ⋅ v) * dΓ # Left-hand size
# Solution of the FE problem
op = AffineFEOperator(a,l,Uₕ,Vₕ)
uₕ = solve(op)
# Export the displacement
writevtk(Ω,"3D_Problem_uh",cellfields=["u"=>uₕ])




