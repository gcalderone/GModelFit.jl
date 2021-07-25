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
    eval1!(multi)
    eval2!(multi)
    eval3!(multi)
    eval4!(multi)
end

function eval1!(multi::MultiModel)
    for id in 1:length(multi.models)
        eval1!(multi.models[id])
    end
    empty!(multi.patchcomps)
    for id in 1:length(multi.models)
        push!(multi.patchcomps, multi.models[id].meval.patchcomps)
    end
end

function eval2!(multi::MultiModel)
    for id in 1:length(multi.models)
        eval2!(multi.models[id], fromparent=true)
    end
    for pf in multi.patchfuncts
        pf.funct(multi.patchcomps)
    end
    nothing
end

function eval3!(multi::MultiModel)
    for id in 1:length(multi.models)
        eval3!(multi.models[id])
    end
end

function eval4!(multi::MultiModel)
    for id in 1:length(multi.models)
        eval4!(multi.models[id])
    end
end

function push!(multi::MultiModel, model::Model)
    push!(multi.models, model)
    model.parent = multi
    evaluate!(multi)
end

function patch!(multi::MultiModel, exfunc::ExprFunction)
    evaluate!(multi)  # ensure sub models are evaluated before adding
                      # new patch functions
    push!(multi.patchfuncts, exfunc)
    evaluate!(multi)
    return multi
end
