struct ComponentSnapshot
    comptype::String
    params::OrderedDict{Symbol, Parameter}
    isfrozen::Bool
    deps::Vector{Symbol}
    evalcounter::Int
    buffer::Array{Float64}
end

getparams(comp::ComponentSnapshot) = comp.params

# ====================================================================
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
    folded::Array{Float64}
end
function ModelSnapshot(meval::ModelEval, params::OrderedDict{NTuple{2, Symbol}, Parameter})
    comps = OrderedDict{Symbol, ComponentSnapshot}()
    for cname in keys(meval.model.comps)
        # TODO is this needed ? (cname in keys(meval.cevals))  ||  continue
        p = OrderedDict([k[2] => p for (k,p) in
                             filter(p -> p.first[1] == cname, params)])
        comps[cname] = ComponentSnapshot(comptype(meval.model, cname), p,
                                         isfrozen(meval.model, cname),
                                         dependencies(meval.model, cname),
                                         evalcounter(meval, cname),
                                         deepcopy(meval.cevals[cname].buffer))
    end
    ModelSnapshot(deepcopy(meval.domain), comps,
                  meval.seq[end], meval.folded_domain, meval.folded)
end

domain(model::ModelSnapshot) = model.domain
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

function Base.getindex(model::ModelSnapshot, cname::Symbol)
    @assert cname in keys(model.comps) "$cname is not a component in the model"
    return model.comps[cname]
end

(model::ModelSnapshot)(cname::Symbol) =             getfield(model.comps[cname], :buffer)
isfrozen(model::ModelSnapshot, cname::Symbol) =     getfield(model.comps[cname], :isfrozen)
dependencies(model::ModelSnapshot, cname::Symbol) = getfield(model.comps[cname], :deps)
evalcounter(model::ModelSnapshot, cname::Symbol) =  getfield(model.comps[cname], :evalcounter)
comptype(model::ModelSnapshot, cname::Symbol) =     getfield(model.comps[cname], :comptype)


folded_domain(model::ModelSnapshot) = model.folded_domain
folded(model::ModelSnapshot) = model.folded


# ====================================================================
struct ModelSetSnapshot
    dict::OrderedDict{Symbol, ModelSnapshot}
end

function ModelSetSnapshot(mseval::ModelSetEval)
    out = ModelSetSnapshot(OrderedDict{Symbol, ModelSnapshot}())
    i = 0
    for (mname, model) in mseval.ms
        i += 1
        params = OrderedDict([(k[2], k[3]) => p for (k,p) in
                                  filter(p -> p.first[1] .== mname, getparams(mseval))])
        out.dict[mname] = ModelSnapshot(mseval.vec[i], params)
    end
    return out
end
