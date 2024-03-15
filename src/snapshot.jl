"""
    ModelSnapshot

A structure containing a *snapshot* (i.e. a "*frozen*" state) of a `Model`.  A snapshot contains the same parameters and component evaluations of the original model, and provide the same user interface.  Moreover, a `ModelSnapshot` can be serialized to a file and de-serialized in another Julia session (see `GModelFit.serialize()`).

The best fit model and best fit parameter values returned as a `ModelSnapshot` object by the `fit()` function.
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
function ModelSnapshot(meval::ModelEval)
    deps = OrderedDict{Symbol, Vector{Symbol}}()
    for cname in keys(meval.cevals)
        deps[cname] = dependencies(meval.model, cname)
    end
    ModelSnapshot(deepcopy(meval.domain), deepcopy(meval.pv.params),
                  OrderedDict([Pair(cname, ceval.buffer) for (cname, ceval) in meval.cevals]),
                  meval.maincomp,
                  comptypes(meval.model),
                  OrderedDict([Pair(cname, isfreezed(meval.model, cname)) for cname in keys(meval.cevals)]),
                  deps, evalcounters(meval))
end

domain(model::ModelSnapshot) = model.domain
Base.keys(model::ModelSnapshot) = collect(keys(model.buffers))
(model::ModelSnapshot)() = reshape(domain(model), model.buffers[model.maincomp])
(model::ModelSnapshot)(name::Symbol) = reshape(domain(model), model.buffers[name])
find_maincomp(model::ModelSnapshot) = model.maincomp
isfreezed(model::ModelSnapshot, cname::Symbol) = model.isfreezed[cname]
dependencies(model::ModelSnapshot, cname::Symbol) = model.deps[cname]
evalcounter(model::ModelSnapshot, cname::Symbol) = model.evalcounters[cname]
comptype(model::ModelSnapshot, cname::Symbol) = model.comptypes[cname]
comptypes(model::ModelSnapshot) = model.comptypes
Base.haskey(m::ModelSnapshot, name::Symbol) = haskey(m.params, name)
function Base.getindex(model::ModelSnapshot, name::Symbol)
    if name in keys(model.params)
        return model.params[name]
    end
    error("Name $name not defined")
end

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
