using Gridap
for n in names(Gridap; all=true)
    if occursin("coord", string(n))
        println(n)
    end
end
