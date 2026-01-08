using Gridap
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using LinearAlgebra
using FillArrays

# Mesh
nodes = Point{2,Float64}[
(0,0),(0,1),(0,2),
(1,0),(1,1),(1,2),
(2,0),(2,1),(2,2)]
conn = [[1,2,4,5],[2,3,5,6],[4,5,7,8],[5,6,8,9]]

# Reference shape functions
filter(e,p) = true
ref_nodes = Point{2,Float64}[
(0,0),(1,0),(0,1),(1,1)]
m = MonomialBasis{2}(Float64,1,filter)
l = LagrangianDofBasis(Float64,ref_nodes)
change = inv(evaluate(l,m))
s = linear_combination(change,m)

# Reference quadrature rule
q = Point{2,Float64}[(.5,.5)]
w = [1.]

# Cell-wise shape functions
ncells = length(conn)
cell_s = Fill(s,ncells)

# Cell-wise quadrature
cell_q = Fill(q,ncells)
cell_w = Fill(w,ncells)

# Geometrical map and Jacobian
cell_nodes = lazy_map(
Broadcasting(Reindex(nodes)),conn)
cell_φ = lazy_map(
linear_combination,cell_nodes,cell_s)
cell_Jt = lazy_map(∇,cell_φ)

# Cell-wise shape function derivatives
cell_∇ref_s = lazy_map(Broadcasting(∇),cell_s)
cell_invJt = lazy_map(Operation(inv),cell_Jt)
cell_∇s = lazy_map(
Broadcasting(Operation(⋅)),cell_invJt,cell_∇ref_s)

# Cell-wise stifness matrix
cell_∇st = lazy_map(transpose,cell_∇s)
cell_∇s∇st = lazy_map(
Broadcasting(Operation(⋅)),cell_∇s,cell_∇st)

cell_mat = lazy_map(
integrate,cell_∇s∇st,cell_q,cell_w,cell_Jt)