# ====================================================================
# MultiModel
#
mutable struct MultiModel <: AbstractMultiModel
    models::Vector{Model}
    patchfuncts::Vector{ExprFunction}
    patchcomps::Vector{OrderedDict{Symbol, PatchComp}}

    function MultiModel(v::Vararg{Model})
        multi = new([v...], Vector{ExprFunction}(),
                    Vector{OrderedDict{Symbol, PatchComp}}())
        for m in v
            m.parent = multi
        end
        evaluate!(multi)
        return multi
    end
end

Base.getindex(m::MultiModel, id::Int) = m.models[id]
Base.length(m::MultiModel) = length(m.models)


function evaluate!(multi::MultiModel)
    for id in 1:length(multi.models)
        evaluate!(multi.models[id])
    end
    empty!(multi.patchcomps)
    for id in 1:length(multi.models)
        push!(multi.patchcomps, multi.models[id].peval.patchcomps)
    end
    patch_params(multi)
    quick_evaluate(multi)
end

function patch_params(multi::MultiModel)
    for id in 1:length(multi.models)
        model = multi.models[id]
        model.peval.patched .= model.peval.pvalues  # copy all values by default
        for pf in model.patchfuncts
            pf.funct(model.peval.patchcomps)
        end
    end
    for pf in multi.patchfuncts
        pf.funct(multi.patchcomps)
    end
    nothing
end

function quick_evaluate(multi::MultiModel)
    for id in 1:length(multi.models)
        quick_evaluate(multi.models[id])
    end
end

function push!(multi::MultiModel, model::Model)
    push!(multi.models, model)
    model.parent = multi
    evaluate!(multi)
end

function patch!(multi::MultiModel, exfunc::ExprFunction)
    push!(multi.patchfuncts, exfunc)
    evaluate!(multi)
    return multi
end
