using SparseArrays, LinearAlgebra
using IterativeSolvers
using BenchmarkTools
using Test
using IterativeSolvers
using Test
using LinearAlgebra
using SparseArrays
using Random
# # Create a large sparse symmetric positive definite matrix A
# n = 10000
# A = sprand(n, n, 0.0005)  # very sparse
# A = A + A' + n * I        # make it symmetric and SPD

# b = randn(n)

# println("---- Conjugate Gradient (IterativeSolvers.jl) ----")
# @btime cg(A, b; maxiter=10000, reltol=1e-6);

# println("---- Direct Solver (factorize(Symmetric(A))) ----")

# @btime factorize(Symmetric(A)) \ b;

# # @test norm(x1) ≈ norm(x2)

T = Float64
  n = 1000
  A = sprand(n, n,0.1)
  A = A' * A + I
  b = rand(n)
  reltol = √eps(real(T))

  # x, ch = 
  @btime cg(A, b; reltol=reltol, maxiter=10000, log=true);
  @btime cg(A, b; maxiter=10000, reltol=1e-6);
  x1,ch=cg(A, b; reltol=reltol, maxiter=10000, log=true);
  x2=cg(A, b; maxiter=10000, reltol=1e-6);
   @test norm(x1) ≈ norm(x2)

