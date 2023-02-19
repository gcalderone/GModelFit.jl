# ====================================================================
"""
    MultiModel

A structure containing a multi-model description, whose evaluation is suitable to be compared to a set of empirical datasets. A `MultiModel` is very similar to a vector of `Model` objects, with the ability to trasparently handle the patch constraints between one model and the other.

Constructor is: `MultiModel(model1, model2, ...)`.

You may access the individual `Model` objects the indexing syntax, as if it was a `Vector{Model}`.  Also, you may add new model to a `MultiModel` after it has been created using the `push!()` function. Finally, you may retrieve the length of the vector with `length()`.

Just like a `Model` object you may need to manually trigger a `MultiModel` evaluation using the `update!()` function.
"""
struct MultiModel
    models::Vector{Model}
    function MultiModel(v::Vararg{Model})
        multi = new([v...])
        update!(multi)
        return multi
    end
end

Base.getindex(m::MultiModel, id::Int) = m.models[id]

"""
    length(multi::MultiModel)

Returns how many `Model` objects are contained in a `MultiModel`.
"""
Base.length(multi::MultiModel) = length(multi.models)

"""
    push!(multi::MultiModel, model::Model)

Push a new `Model` object into a `MultiModel`.
"""
function push!(multi::MultiModel, model::Model)
    push!(multi.models, model)
    update!(multi)
end


"""
    update!(multi::MultiModel)

Update! a `MultiModel` and update internal structures.
"""
function update!(multi::MultiModel)
    update_step0(multi)
    # update_step1(multi)
    update_step2(multi)
    update_step3(multi)
    update_step4(multi)
    return multi
end


function free_params_indices(multi::MultiModel)
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(multi)
        nn = length(multi[id].pv.ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


function update_step0(multi::MultiModel)
    for id in 1:length(multi)
        update_step0(multi.models[id])
    end
    for i in 1:length(multi)
        for j in 1:length(multi)
            push!(multi.models[i].pv.mvalues, multi.models[j].pv.values)
        end
    end
end

function update_step1(multi::MultiModel, pvalues::Vector{Float64})
    for (id, i1, i2) in free_params_indices(multi)
        update_step1(multi[id], pvalues[i1:i2])
    end
end

function update_step2(multi::MultiModel)
    for id in 1:length(multi)
        update_step2(multi.models[id])
    end
end

function update_step3(multi::MultiModel)
    for id in 1:length(multi)
        update_step3(multi.models[id])
    end
end

function update_step4(multi::MultiModel, uncerts=Vector{Float64}[])
    for (id, i1, i2) in free_params_indices(multi)
        if length(uncerts) > 0
            update_step4(multi[id], pvalues[i1:i2])
        else
            update_step4(multi[id])
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
