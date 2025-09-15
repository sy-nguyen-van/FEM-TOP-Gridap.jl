using GridapGmsh
using Gridap, GridapPETSc
using GridapDistributed
using PartitionedArrays

function main_ex4(nparts,distribute)
  parts  = distribute(LinearIndices((nparts,)))
  options = "-ksp_type cg -pc_type gamg -ksp_monitor"
  GridapPETSc.with(args=split(options)) do
    model = GmshDiscreteModel(parts,"demo.msh")
    order = 2
    u((x,y)) = (x+y)^order
    f(x) = -Δ(u,x)
    reffe = ReferenceFE(lagrangian,Float64,order)
    V = TestFESpace(model,reffe,dirichlet_tags=["boundary1","boundary2"])
    U = TrialFESpace(u,V)
    Ω = Triangulation(model)
    dΩ = Measure(Ω,2*order)
    a(u,v) = ∫( ∇(v)⋅∇(u) )dΩ
    l(v) = ∫( v*f )dΩ
    op = AffineFEOperator(a,l,U,V)
    solver = PETScLinearSolver()
    uh = solve(solver,op)
    writevtk(Ω,"results_ex4",cellfields=["uh"=>uh,"grad_uh"=>∇(uh)])
  end
end

nparts = 4
with_mpi() do distribute
  main_ex4(nparts,distribute)
end

