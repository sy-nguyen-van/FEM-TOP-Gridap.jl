using Gridap, Gridap.CellData, SparseArrays
using NearestNeighbors, Printf  # make sure this is at the top of your file
using Makie,GLMakie, LaTeXStrings, Makie.GeometryBasics
include("plot_design.jl")
include("OC.jl")
include("SIMP_Inter.jl")
# ==============================================
L = 60 # Length
W = 20 # Height
domain = (0, L, 0, W) # sizes along X, Y
partition = (L*2, W*2) # mesh sizes
model = CartesianDiscreteModel(domain, partition)
writevtk(model, "Cantilever_Beam_2D")
# Tag the left face (x = 0); right face (x = L)
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "left", [1, 3, 7]) # 
add_tag_from_tags!(labels, "right", [2, 4, 8]) # 
# Integrating on the domain Ω
degree = 2
Ω = Triangulation(model)
dΩ = Measure(Ω, degree)
# neumann boundaries
Γ = BoundaryTriangulation(model, tags="right")
dΓ = Measure(Γ, degree)
# Vector-valued FE space
order = 1
reffe = ReferenceFE(lagrangian, VectorValue{2,Float64}, order)
Vₕ = TestFESpace(Ω, reffe; conformity=:H1, dirichlet_tags="left")
Uₕ = TrialFESpace(Vₕ)
# Constitutive law for LINEAR ELASTICITY 
const E = 119e3 # Young's modulus
const ν = 0.3   # Poisson's ratio
const λ = (E * ν) / ((1 + ν) * (1 - 2 * ν)) # Lamé parameters
const μ = E / (2 * (1 + ν)) # Lamé parameters
ε₀(u) = 0.5 * (∇(u) + transpose(∇(u))) # Strain; In Gridap this can be automatically defined in Gridap.ε
σ(ε₀) = λ * tr(ε₀) * one(ε₀) + 2μ * ε₀ # Stress

num_ele = num_cells(Ω) # Number of elements
volfrac = 0.4
p = 3
max_iter = 50
# Compute centroids of elements and store in a matrix
cell_coords = get_cell_coordinates(Ω) # Get coordinates of vertices for each cell
no_ele = length(cell_coords) # Number of elements
centroid_matrix = zeros(2, no_ele) # Matrix to store centroids (no_ele × 2)

for (cell_id, coords) in enumerate(cell_coords)
    # Compute centroid by averaging vertex coordinates
    centroid = sum(coords) / length(coords)
    centroid_matrix[:, cell_id] = [centroid[1], centroid[2]] # Store x, y coordinates
end

# Number of neighbors for each element
filter_radius = 2.5
tol = 1e-3
search_dist = filter_radius - tol
balltree = BallTree(centroid_matrix)
OPT_H = zeros(num_ele, num_ele)
for iel in 1:num_ele
    near_ele = inrange(balltree, centroid_matrix[:, iel], search_dist)
    dist = sqrt.(sum((centroid_matrix[:, near_ele] .- centroid_matrix[:, iel]) .^ 2, dims=1))
    num = 1 .- dist ./ filter_radius
    den = sum(num)
    OPT_H[iel, near_ele] = num / den
end
H = sparse(OPT_H)

for iter = 1:max_iter
    if iter == 1
        ρ_new = 0.5 * ones(num_ele)
    end
    ρ_old = copy(ρ_new)
    ρ_new_fil = H * ρ_new # Filter
    # SIMP
    ρ_SIMP, dρ_dp = SIMP_Inter(ρ_new_fil, p, Ω)
    # ===FE-ANALYSIS
    F(x) = VectorValue(0.0, -1) # I want to apply a unit load along Y axis to each node on the right face
    l(v) = ∫(F ⋅ v) * dΓ # Right-hand size
    a(u, v) = ∫(ρ_SIMP * (σ ∘ ε₀(u)) ⊙ ε₀(v)) * dΩ # Left-hand size; (∘) Composite functions
    # Solution of the FE problem
    op = AffineFEOperator(a, l, Uₕ, Vₕ)
    uh = solve(op)

    C = sum(∫(ρ_SIMP * (σ ∘ ε₀(uh)) ⊙ ε₀(uh)) * dΩ) # Compliance

    dC_fil = get_contribution(∫(-dρ_dp * (σ ∘ ε₀(uh)) ⊙ ε₀(uh)) * dΩ, Ω) # Sensitivity of compliance w.r.t rho
    dC = transpose(H) * (dC_fil[:])

    ρ_new = OC(ρ_new, volfrac, dC, num_ele)
    change_obj = maximum(abs.(ρ_new .- ρ_old))
    figs = plot_design(ρ_new, cell_coords, L, W)

    display(figs)

    println("It. ", iter, ", Obj = ", @sprintf("%.5e", C), ", Vol.: ", @sprintf("%.5e", sum(ρ_new) / num_ele), ", ch.: ", @sprintf("%.5e", change_obj))
    if iter == max_iter
        writevtk(Ω, "Cantilever_Beam_2D_TOP", cellfields=["rho" => ρ_SIMP])
    end
end



