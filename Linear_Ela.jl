using Gridap, Gridap.CellData, SparseArrays, BenchmarkTools
using NearestNeighbors, Printf  # make sure this is at the top of your file
using Makie,GLMakie, LaTeXStrings, Makie.GeometryBasics
using Gridap
using GridapGmsh
# Discrete model
model = GmshDiscreteModel("Lbracket2d_Sy.msh")
writevtk(model,"model")

# Vector-valued FE space
order = 1
reffe = ReferenceFE(lagrangian,VectorValue{3,Float64},order)
V0 = TestFESpace(model,reffe;
  conformity=:H1,
  dirichlet_tags="BC_Surf")

U = TrialFESpace(V0)

# Constitutive law for LINEAR ELASTICITY 
const ν = 0.3   # Poisson's ratio
const E = 2100 # Young's modulus

# lame parameters
const λ = (E*ν)/((1+ν)*(1-2*ν))
const μ = E/(2*(1+ν))

# stress definition
σ(ε) = λ * tr(ε) * one(ε) + 2 * μ * ε; # Stress

# Construct weak form
degree = 2*order

Ω = Triangulation(model)
dΩ = Measure(Ω, degree)

# neumann boundaries
Γ = BoundaryTriangulation(model,tags = "Load_Surf")
dΓ = Measure(Γ,degree)

# force application
F(x) = VectorValue(0.0, -1, 0.0) # A unit load along Y axis applies to each node on the right face

# bi-linear form
a(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ

# linear form
l(v) = ∫(F ⋅ v) * dΓ # Right-hand size

# Solution of the FE problem
op = AffineFEOperator(a,l,U,V0)
uh = solve(op)

# write the resukt to a file
writevtk(Ω,"results",cellfields=["uh"=>uh,"epsi"=>ε(uh),"sigma"=>σ∘ε(uh)])