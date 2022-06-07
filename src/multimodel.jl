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
    eval_step1(multi)
    eval_step2(multi)
    eval_step3(multi)
    eval_step4(multi)
    return multi
end

function eval_step1(multi::MultiModel)
    for id in 1:length(multi.models)
        eval_step1(multi.models[id])
    end
    empty!(multi.pvalues)
    for id in 1:length(multi.models)
        push!(multi.pvalues, multi.models[id].pvalues)
    end
end

function eval_step2(multi::MultiModel)
    for id in 1:length(multi.models)
        eval_step2(multi.models[id])
    end
end

function eval_step3(multi::MultiModel)
    for id in 1:length(multi.models)
        eval_step3(multi.models[id])
    end
end

function eval_step4(multi::MultiModel)
    for id in 1:length(multi.models)
        eval_step4(multi.models[id])
    end
end
