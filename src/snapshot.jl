"""
    ModelSnapshot

A structure containing a *snapshot* (i.e. a "*frozen*" state) of a `Model`.  A snapshot contains the same parameters and component evaluations of the original model, and provide the same user interface.  Moreover, a `ModelSnapshot` can be serialized to a file and de-serialized in another Julia session (see `GModelFit.serialize()`).

The best fit model and parameter values returned by the `fit()` function are provided as a `ModelSnapshot` object .
"""
struct ModelSnapshot
    domain::AbstractDomain
    params::PVModel{Parameter}
    buffers::OrderedDict{Symbol, Vector{Float64}}
    maincomp::Symbol
    comptypes::OrderedDict{Symbol, String}
    isfreezed::OrderedDict{Symbol, Bool}
    deps::OrderedDict{Symbol, Vector{Symbol}}
    evalcounters::OrderedDict{Symbol, Int}
end
function ModelSnapshot(meval::ModelEval, bestfit::PVModel{Parameter})
    deps = OrderedDict{Symbol, Vector{Symbol}}()
    for cname in keys(meval.cevals)
        deps[cname] = dependencies(meval.model, cname)
    end

    ModelSnapshot(deepcopy(meval.domain), deepcopy(bestfit),
                  OrderedDict([Pair(cname, ceval.tpar.buffer) for (cname, ceval) in meval.cevals]),
                  meval.seq[end],
                  comptypes(meval.model),
                  OrderedDict([Pair(cname, isfreezed(meval.model, cname)) for cname in keys(meval.cevals)]),
                  deps, evalcounters(meval))
end

domain(model::ModelSnapshot) = model.domain
Base.keys(model::ModelSnapshot) = collect(keys(model.buffers))
(model::ModelSnapshot)() = reshape(domain(model), model.buffers[model.maincomp])
(model::ModelSnapshot)(name::Symbol) = reshape(domain(model), model.buffers[name])
function deptree(model::ModelSnapshot)
    function deptree(model, cname::Symbol, level::Int, parent::Union{Nothing, Symbol})
        out = DependencyNode(cname, level, parent)
        for d in dependencies(model, cname)
            push!(out.childs, deptree(model, d, level+1, cname))
        end
        return out
    end
    return deptree(model, model.maincomp, 1, nothing)
end
isfreezed(model::ModelSnapshot, cname::Symbol) = model.isfreezed[cname]
dependencies(model::ModelSnapshot, cname::Symbol) = model.deps[cname]
evalcounter(model::ModelSnapshot, cname::Symbol) = model.evalcounters[cname]
comptype(model::ModelSnapshot, cname::Symbol) = model.comptypes[cname]
comptypes(model::ModelSnapshot) = model.comptypes
Base.haskey(m::ModelSnapshot, name::Symbol) = haskey(m.buffers, name)
function Base.getindex(model::ModelSnapshot, name::Symbol)
    @assert name in keys(model.buffers)
    if name in keys(model.params)
        return model.params[name]
    end
    return PVComp(model.params)
end

Base.length(model::ModelSnapshot) = length(model.buffers)

function iterate(model::ModelSnapshot, i=1)
    k = collect(keys(model))
    (i > length(k))  &&  return nothing
    return (k[i] => model[k[i]], i+1)
end

function getparams(comp::GModelFit.PV.PVComp{GModelFit.Parameter})
    out = OrderedDict{Symbol, Parameter}()
    for pname in propertynames(comp)
        out[pname] = getproperty(comp, pname)
    end
    return out
end
