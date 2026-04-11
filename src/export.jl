using DataFrames
export export_params, export_folded_eval, export_unfolded_eval

function export_params(obj::Union{ <: AbstractModel, <: AbstractGroup, <: AbstractCompSlot})
    replace(::Nothing) = ""
    replace(f::Symbol) = string(f)
    replace(f::FunctDesc) = f.display

    evaluate!(isa(obj, AbstractModel)  ?  obj  :  obj.rootmodel)
    df = DataFrame()
    for (pname, par) in getparams(obj)
        cols = (:pname, propertynames(par)...)
        vals = (string(pname), [getproperty(par, n) for n in propertynames(par)]...)
        push!(df, NamedTuple{cols}(vals), promote=true)
    end
    df[!, :patch] = replace.(df.patch)
    df[!, :reparam] = replace.(df.reparam)
    return df
end

export_domain(dom::CartesianDomain) = export_domain(Domain(dom))
function export_domain(dom::Domain{N}) where N
    df = DataFrame()
    for i in 1:N
        df[!, Symbol(:coord, i)] = coords(dom, i)
    end
    return df
end

export_data(data::GaussianData)  = DataFrame(:data_vals => reshape(values(data), :), :data_unc => reshape(uncerts(data), :))
export_data(data::PoissonCounts) = DataFrame(:data_vals => reshape(values(data)))

export_folded_eval(model::AbstractModelEval{Tuple{}}) = export_folded_eval(model[()])
function export_folded_eval(group::AbstractGroupEval)
    evaluate!(group.rootmodel)
    df = export_domain(domain(group))
    isa(group, AbstractGroupLH)  &&  (df = hcat(df, export_data(getdata(group))))
    df[!, :model] = reshape(group(), :)
    df[!, :resid] = reshape(getresiduals(group), :)
    return df
end

export_unfolded_eval(model::AbstractModelEval{Tuple{}}) = export_unfolded_eval(model[()])
function export_unfolded_eval(group::AbstractGroupEval)
    evaluate!(group.rootmodel)
    df = export_domain(unfolded_domain(getresp(group)))
    for (cname, comp) in getcomps(group)
        df[!, cname[end]] = reshape(comp(), :)
    end
    return df
end
