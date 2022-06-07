# ====================================================================
# MultiModel
#
mutable struct MultiModel <: AbstractMultiModel
    models::Vector{Model}
    patchfuncts::Vector{λFunct}
    patchcomps::Vector{OrderedDict{Symbol, HashVector{Float64}}}

    function MultiModel(v::Vararg{Model})
        multi = new([v...], Vector{λFunct}(),
                    Vector{OrderedDict{Symbol, HashVector{Float64}}}())
        for m in v
            m.parent = multi
        end
        evaluate(multi)
        return multi
    end
end

Base.getindex(m::MultiModel, id::Int) = m.models[id]
Base.length(m::MultiModel) = length(m.models)


function evaluate(multi::MultiModel)
    eval_step1(multi)
    eval_step2(multi)
    eval_step3(multi)
    eval_step4(multi)
end

function eval_step1(multi::MultiModel)
    for id in 1:length(multi.models)
        eval_step1(multi.models[id])
    end
    empty!(multi.patchcomps)
    for id in 1:length(multi.models)
        push!(multi.patchcomps, multi.models[id].meval.patchcomps)
    end
end

function eval_step2(multi::MultiModel)
    for id in 1:length(multi.models)
        eval_step2(multi.models[id], fromparent=true)
    end
    for pf in multi.patchfuncts
        pf.funct(multi.patchcomps)
    end
    nothing
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

function push!(multi::MultiModel, model::Model)
    push!(multi.models, model)
    model.parent = multi
    evaluate(multi)
end

function patch!(func::λFunct, multi::MultiModel)
    evaluate(multi)  # ensure sub models are evaluated before adding new patch functions
    push!(multi.patchfuncts, func)
    evaluate(multi)
    return multi
end
