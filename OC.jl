function OC(ρ, volfrac, dc, num_ele)
    ρ_min = 0.001
    l1 = 0.0
    l2 = 100000.0
    move = 0.2
    ρ_new = copy(ρ)
    while (l2 - l1 > 1e-4)
        lmid = 0.5 * (l2 + l1)
        ρ_candidate = ρ .* sqrt.(-dc ./ lmid)
        ρ_new = max.(ρ_min, max.(ρ .- move, min.(1.0, min.(ρ .+ move, ρ_candidate))))
        if sum(ρ_new) - volfrac * num_ele > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end
    return ρ_new
end
