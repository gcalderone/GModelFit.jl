# ====================================================================
function free_params_indices(mevals::Vector{ModelEval})
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(mevals)
        nn = length(mevals[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


function free_params(mevals::Vector{ModelEval})
    out = Vector{Parameter}()
    for id in 1:length(mevals)
        append!(out, free_params(mevals[id]))
    end
    return out
end
nfree(mevals::Vector{ModelEval}) = sum(nfree.(mevals))


function set_pvalues!(mevals::Vector{ModelEval}, pvalues::Vector{Float64})
    for (id, i1, i2) in free_params_indices(mevals)
        set_pvalues!(mevals[id], pvalues[i1:i2])
    end
end

function update!(mevals::Vector{ModelEval})
    if length(mevals[1].pvmulti) == 0
        pvmulti = [mevals[i].pvalues for i in 1:length(mevals)]
        for i in 1:length(mevals)
            append!(mevals[i].pvmulti, pvmulti)
        end
    end
    return update!.(mevals)
end

function set_bestfit!(mevals::Vector{ModelEval}, pvalues::Vector{Float64}, uncerts::Vector{Float64})
    for (id, i1, i2) in free_params_indices(mevals)
        set_bestfit!(mevals[id], pvalues[i1:i2], uncerts[i1:i2])
    end
end
