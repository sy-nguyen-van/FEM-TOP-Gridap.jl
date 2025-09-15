using Gridap
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces: TestBasis, TrialBasis
using Gridap.FESpaces: SingleFieldFEBasis
using FillArrays

# Mesh
domain = (0, 2, 0, 2)
cells = (2, 2)
model = CartesianDiscreteModel(domain, cells)

Ω = Triangulation(model)
Qₕ = CellQuadrature(Ω, 2)
Qₕ_cell_data = get_data(Qₕ)

Qₕ_cell_point = get_cell_points(Qₕ)
qₖ = get_data(Qₕ_cell_point)
q = Qₕ_cell_data[1]
p = get_coordinates(q)
w = get_weights(q)

# Reference shape functions
k = 1
reffe = ReferenceFE(QUAD, lagrangian, VectorValue{2,Float64}, k)
Vₕ = FESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
Uₕ = TrialFESpace(Vₕ)

dv = get_fe_basis(Vₕ)
du = get_trial_fe_basis(Uₕ)

shape_test = evaluate(dv[1], Qₕ_cell_point)
shape_trial = evaluate(du[1], Qₕ_cell_point)

M = zeros(8, 8)
for i = 1:4
  detJK_wi = w[i]
  for j = 1:8
    for k = 1:8
      M[j, k] += shape_test[i, j] * shape_trial[i, k] * detJK_wi

    end
  end
end


M


# Integration measure

dv_at_Qₕ[1]
using Gridap.TensorValues
# Define material constants
E = 119e3  # Young's modulus (Pa)
ν = 0.3    # Poisson's ratio
λ = E * ν / ((1 + ν) * (1 - 2 * ν))  # Lamé constant
μ = E / (2 * (1 + ν))               # Shear modulus
Ce = SymFourthOrderTensorValue(
  λ + 2 * μ, λ, 0.0,   # Row 1: C1111, C1122, C1112
  λ, λ + 2 * μ, 0.0,   # Row 2: C2211, C2222, C2212
  0.0, 0.0, μ      # Row 3: C1211, C1222, C1212
) # Plain
σ(ε₀) = λ * tr(ε₀) * one(ε₀) + 2μ * ε₀; # Stress



# ε₀(u) = 0.5 * ( ∇(u) + transpose(∇(u)))
# a1 = ∫((ε₀(u)⊙ C(λ,μ)) ⊙  ε₀(v))*dΩ
# Cell-wise stiffness matrix
