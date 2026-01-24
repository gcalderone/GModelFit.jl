struct ComponentSnapshot
    comptype::String
    isfreezed::Bool
    deps::Vector{Symbol}
    evalcounter::Int
    params::OrderedDict{Symbol, Parameter}
    buffer::Vector{Float64}
end
getproperty(comp::ComponentSnapshot, name::Symbol) = getfield(comp, :params)[name]
getindex(comp::ComponentSnapshot, name::Symbol)    = getfield(comp, :params)[name]
propertynames(comp::ComponentSnapshot) = collect(keys(getfield(comp, :params)))

getparams(comp::ComponentSnapshot) = getfield(comp, :params)

"""
    ModelSnapshot

A structure containing a *snapshot* (i.e. a "*frozen*" state) of a `Model`.  A snapshot contains the same parameters and component evaluations of the original model, and provide the same user interface.  Moreover, a `ModelSnapshot` can be serialized to a file and de-serialized in another Julia session (see `GModelFit.serialize()`).

The best fit model and parameter values returned by the `fit()` function are provided as a `ModelSnapshot` object .
"""
struct ModelSnapshot
    domain::AbstractDomain
    comps::OrderedDict{Symbol, ComponentSnapshot}
    maincomp::Symbol
    folded_domain::AbstractDomain
    folded::Vector{Float64}
end
function ModelSnapshot(meval::ModelEval, bestfit::PVModel{Parameter})
    comps = OrderedDict{Symbol, ComponentSnapshot}()
    for cname in keys(meval.cevals)
        params = OrderedDict{Symbol, Parameter}()
        for (pname, par) in bestfit[cname]
            params[pname] = par
        end
        comps[cname] = ComponentSnapshot(comptype(meval.model, cname),
                                         isfreezed(meval.model, cname),
                                         dependencies(meval.model, cname),
                                         evalcounter(meval, cname),
                                         params,
                                         meval.cevals[cname].buffer)
    end

    ModelSnapshot(deepcopy(meval.domain), comps,
                  meval.seq[end], meval.folded_domain, meval.folded)
end

domain(model::ModelSnapshot) = model.domain
Base.keys(model::ModelSnapshot) = collect(keys(model.comps))
(model::ModelSnapshot)() = model(model.maincomp)

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

Base.haskey(m::ModelSnapshot, name::Symbol) = haskey(m.comps, name)
function Base.getindex(model::ModelSnapshot, name::Symbol)
    @assert name in keys(model.comps) "$name is not a component in the model"
    return model.comps[name]
end

Base.length(model::ModelSnapshot) = length(model.comps)

function iterate(model::ModelSnapshot, i=1)
    k = collect(keys(model))
    (i > length(k))  &&  return nothing
    return (k[i] => model[k[i]], i+1)
end

(model::ModelSnapshot)(cname::Symbol) = reshape(domain(model), getfield(model.comps[cname], :buffer))
isfreezed(model::ModelSnapshot, cname::Symbol) =               getfield(model.comps[cname], :isfreezed)
dependencies(model::ModelSnapshot, cname::Symbol) =            getfield(model.comps[cname], :deps)
evalcounter(model::ModelSnapshot, cname::Symbol) =             getfield(model.comps[cname], :evalcounter)
comptype(model::ModelSnapshot, cname::Symbol) =                getfield(model.comps[cname], :comptype)


folded_domain(model::ModelSnapshot) = model.folded_domain
folded(model::ModelSnapshot) = model.folded
