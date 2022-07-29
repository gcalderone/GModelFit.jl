# ====================================================================
# MultiModel
#
struct MultiModel <: AbstractMultiModel
    models::Vector{Model}
    pvalues::Vector{HashHashVector{Float64}}
    function MultiModel(v::Vararg{Model})
        multi = new([v...], Vector{HashHashVector{Float64}}())
        for m in v
            m.parent = multi
        end
        evaluate(multi)
        return multi
    end
end

Base.getindex(m::MultiModel, id::Int) = m.models[id]
Base.length(m::MultiModel) = length(m.models)

function push!(multi::MultiModel, model::Model)
    push!(multi.models, model)
    model.parent = multi
    evaluate(multi)
end


function evaluate(multi::MultiModel)
    eval_step0(multi)
    # eval_step1(multi)
    eval_step2(multi)
    eval_step3(multi)
    eval_step4(multi)
    return multi
end


function free_params_indices(multi::MultiModel)
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(multi)
        nn = length(multi[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


function eval_step0(multi::MultiModel)
    for id in 1:length(multi)
        eval_step0(multi.models[id])
    end
    empty!(multi.pvalues)
    for id in 1:length(multi)
        push!(multi.pvalues, multi.models[id].pvalues)
    end
end

function eval_step1(multi::MultiModel, pvalues::Vector{Float64})
    for (id, i1, i2) in free_params_indices(multi)
        eval_step1(multi[id], pvalues[i1:i2])
    end
end

function eval_step2(multi::MultiModel)
    for id in 1:length(multi)
        eval_step2(multi.models[id])
    end
end

function eval_step3(multi::MultiModel)
    for id in 1:length(multi)
        eval_step3(multi.models[id])
    end
end

function eval_step4(multi::MultiModel, uncerts=Vector{Float64}[])
    for (id, i1, i2) in free_params_indices(multi)
        if length(uncerts) > 0
            eval_step4(multi[id], pvalues[i1:i2])
        else
            eval_step4(multi[id])
        end
    end
end

function free_params(multi::MultiModel)
    out = Vector{Parameter}()
    for id in 1:length(multi)
        append!(out, free_params(multi[id]))
    end
    return out
end
