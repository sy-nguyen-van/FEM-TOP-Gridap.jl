# ============================
function Surface_Vec(X, Y, degree_x, degree_y)
    A = ones(size(X))
    dA = zeros(size(X))
    for i_x in 0:degree_x
        for i_y in 0:degree_y
            if i_x != 0 || i_y != 0
                A = vcat(A, (X .^ i_x) .* (Y .^ i_y))
               dA = vcat(dA, (i_x*X .^ (i_x-1)) .* (Y .^ i_y))
            end
        end
    end
    return vec(A), vec(dA)
end