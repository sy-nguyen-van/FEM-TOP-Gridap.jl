using Gridap
using Gridap.Geometry
using LinearAlgebra: norm, tr, dot
using WriteVTK

try
    model = CartesianDiscreteModel((0, 1), (10,))
    grid = get_grid(model)
    coords = Gridap.Geometry.get_node_coordinates(grid)
    println("Gridap.Geometry.get_node_coordinates(grid) works")
catch e
    println("Error: $e")
end
