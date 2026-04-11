struct RespSnapshot <: AbstractResponse
    type::String
    unfolded_domain::AbstractDomain
    folded_domain::AbstractDomain
end
RespSnapshot(   resp::AbstractResponse) = RespSnapshot(stype(resp), folded_domain(resp), unfolded_domain(resp))
folded_domain(  resp::RespSnapshot) = resp.folded_domain
unfolded_domain(resp::RespSnapshot) = resp.unfolded_domain
fold_model!(    resp::RespSnapshot, unfolded::Array, folded::Array) = error("Can't fold model using a RespSnapshot object")

"""
    ModelSnapshot

A structure containing a *snapshot* (i.e. a "*frozen*" state) of a `Model`.  A snapshot contains the same parameters and component evaluations of the original model, and provide the same user interface.  Moreover, a `ModelSnapshot` can be serialized to a file and de-serialized in another Julia session (see `GModelFit.serialize()`).

The best fit model and parameter values returned by the `fit()` function are provided as a `ModelSnapshot` object .
"""
struct ModelSnapshot{GT, CT, PT} <: AbstractModelEval{GT, CT, PT}
    dt::DictTree
    ifree::Vector{Int}

    function ModelSnapshot{GT, CT, PT}() where {GT, CT, PT}
        dt = DictTree()
        add_layer!(dt, SDTree{GT, GroupSnapshot}(on_delete=forbid_delete), label=:groups)
        add_layer!(dt, SDTree{CT, CompSnapshot}( on_delete=forbid_delete), label=:comps)
        add_layer!(dt, SDTree{PT, ParameterEval}(on_delete=forbid_delete), label=:params)
        return new{GT, CT, PT}(dt, Vector{Int}())
    end
end

function ModelSnapshot(dt::DictTree)   # deserialization method
    T = fill(Tuple{}, 3)
    for l in keys(getlabels(dt))
        if l == :groups
            T[1] = StaticDictTrees.getKT(get_layer(dt, l))
        elseif l == :comps
            T[2] = StaticDictTrees.getKT(get_layer(dt, l))
        else
            @assert l == :params
            T[3] = StaticDictTrees.getKT(get_layer(dt, l))
        end
    end
    out = ModelSnapshot{T...}()
    i = 1
    for (k, v) in dt
        if isa(v, GroupSnapshot)
            out.dt[k] = GroupSnapshot(out, k, [getfield(v, f) for f in fieldnames(typeof(v))[3:end]]...)
        elseif isa(v, CompSnapshot)
            out.dt[k] =  CompSnapshot(out, k, [getfield(v, f) for f in fieldnames(typeof(v))[3:end]]...)
        else
            @assert isa(v, ParameterEval)
            out.dt[k] = v
            !v.fixed  &&  !v.actually_fixed  &&  push!(out.ifree, i)
            i += 1
        end
    end
    return out
end

function snapshot(model::ModelEval{T, GT, CT, PT}) where {T, GT, CT, PT}
    out = ModelSnapshot{GT, CT, PT}()
    for (gkey, group) in getgroups(model)
        out.dt[gkey] = GroupSnapshot(out, gkey, group)
    end
    for (ckey, comp) in getcomps(model)
        out.dt[ckey] = CompSnapshot(out, ckey, comp)
    end
    for (pkey, par) in getparams(model)
        out.dt[pkey] = deepcopy(par)
    end
    return out
end

struct LikelihoodSnapshot{GT, CT, PT} <: AbstractModelEval{GT, CT, PT}
    dt::DictTree
    ifree::Vector{Int}
    evalbuf::Array{Float64}
    loglh::Float64
    gofstat::Float64

    function LikelihoodSnapshot{GT, CT, PT}(args...) where {GT, CT, PT}
        dt = DictTree()
        add_layer!(dt, SDTree{GT, GroupLHSnapshot}(on_delete=forbid_delete), label=:groups)
        add_layer!(dt, SDTree{CT, CompSnapshot}(   on_delete=forbid_delete), label=:comps)
        add_layer!(dt, SDTree{PT, ParameterEval}(  on_delete=forbid_delete), label=:params)
        return new{GT, CT, PT}(dt, Vector{Int}(), args...)
    end
end

function LikelihoodSnapshot(dt::DictTree, args...)   # deserialization method
    T = fill(Tuple{}, 3)
    for l in keys(getlabels(dt))
        if l == :groups
            T[1] = StaticDictTrees.getKT(get_layer(dt, l))
        elseif l == :comps
            T[2] = StaticDictTrees.getKT(get_layer(dt, l))
        else
            @assert l == :params
            T[3] = StaticDictTrees.getKT(get_layer(dt, l))
        end
    end
    out = LikelihoodSnapshot{T...}(args...)
    i = 1
    for (k, v) in dt
        if isa(v, GroupLHSnapshot)
            out.dt[k] = GroupLHSnapshot(out, k, [getfield(v, f) for f in fieldnames(typeof(v))[3:end]]...)
        elseif isa(v, CompSnapshot)
            out.dt[k] =    CompSnapshot(out, k, [getfield(v, f) for f in fieldnames(typeof(v))[3:end]]...)
        else
            @assert isa(v, ParameterEval)
            out.dt[k] = v
            !v.fixed  &&  !v.actually_fixed  &&  push!(out.ifree, i)
            i += 1
        end
    end
    return out
end

function snapshot(lh::Likelihood{Float64, GT, CT, PT}) where {GT, CT, PT}
    out = LikelihoodSnapshot{GT, CT, PT}(lh.evalbuf, loglh(lh), gofstat(lh))
    for (gkey, group) in getgroups(lh)
        out.dt[gkey] = GroupLHSnapshot(out, gkey, group)
    end
    for (ckey, comp) in getcomps(lh)
        out.dt[ckey] = CompSnapshot(out, ckey, comp)
    end
    for (pkey, par) in getparams(lh)
        out.dt[pkey] = deepcopy(par)
    end
    return out
end

struct CompSnapshot <: AbstractCompSlotEval
    rootmodel::Union{Nothing, ModelSnapshot, LikelihoodSnapshot}
    key::Tuple
    stype::String
    frozen::Bool
    domain::AbstractDomain
    counter::Int
    deps::Vector{Symbol}
    evalbuf::Array{Float64}
end

CompSnapshot(root::Union{ModelSnapshot, LikelihoodSnapshot}, key::Tuple, comp::CompSlotEval{Float64}) =
    CompSnapshot(root, key, stype(getcomp(comp)), comp.frozen, deepcopy(comp.domain),
                 comp.counter, deepcopy(dependencies(getcomp(comp))), deepcopy(comp.evalbuf))

struct GroupSnapshot <: AbstractGroupEval
    rootmodel::Union{Nothing, ModelSnapshot}
    key::Tuple
    domain::AbstractDomain
    seq::Vector{Symbol}
    evalbuf::Array{Float64}
end

GroupSnapshot(root::ModelSnapshot, key::Tuple, group::GroupEval) =
    GroupSnapshot(root, key, deepcopy(folded_domain(getresp(group))),
                  deepcopy(group.seq), deepcopy(group.evalbuf))

struct GroupLHSnapshot <: AbstractGroupLH
    rootmodel::Union{Nothing, LikelihoodSnapshot}
    key::Tuple
    domain::AbstractDomain
    seq::Vector{Symbol}
    evalbuf::Array{Float64}
    data::AbstractData
    residuals::Vector{Float64}
    loglh::Float64
    gofstat::Float64
    gofstat_name::String
end

GroupLHSnapshot(root::LikelihoodSnapshot, key::Tuple, group::GroupLH{Float64}) =
    GroupLHSnapshot(root, key, deepcopy(folded_domain(getresp(group))),
                    deepcopy(group.seq), deepcopy(group.evalbuf),
                    deepcopy(group.data), deepcopy(group.residuals),
                    loglh(group), gofstat(group), gofstat_name(group))
