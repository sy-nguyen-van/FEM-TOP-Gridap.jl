function filters(filter_radius, Ω)
    # Compute centroids of elements and store in a matrix
    cell_coords = get_cell_coordinates(Ω) # Get coordinates of vertices for each cell
    no_ele = length(cell_coords) # Number of elements
    centroid_matrix = zeros(3, no_ele) # Matrix to store centroids (no_ele × 2)

    for (cell_id, coords) in enumerate(cell_coords)
        # Compute centroid by averaging vertex coordinates
        centroid = sum(coords) / length(coords)
        centroid_matrix[:, cell_id] = [centroid[1], centroid[2], centroid[3]] # Store x, y coordinates
    end

    # Number of neighbors for each element
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
    return H
end