
module PV

using DataStructures

import Base.getindex, Base.setindex!, Base.push!,
Base.getproperty, Base.setproperty!, Base.propertynames,
Base.length, Base.keys, Base.empty!, Base.iterate

export  PVComp, PVModel, PVMulti, items, set_items!

struct PVComp{T}
    params::OrderedDict{Symbol, Int}
    data::Vector{T}
end

struct PVModel{T}
    comps::OrderedDict{Symbol, PVComp}
    data::Vector{T}
end


PVComp(parent::PVModel{T}) where T = PVComp{T}(OrderedDict{Symbol, Int}(), parent.data)
internal_vector(pv::PVComp)  = getfield(pv, :data)
internal_dict(  pv::PVComp)  = getfield(pv, :params)
propertynames(pv::PVComp) = collect(keys(internal_dict(pv)))

getindex(pv::PVComp, key::Symbol) =
    internal_vector(pv)[internal_dict(pv)[key]]
getproperty(pv::PVComp, key::Symbol) = getindex(pv, key)

function setindex!(pv::PVComp, value, key::Symbol)
    if haskey(internal_dict(pv), key)
        setindex!(internal_vector(pv), value, internal_dict(pv)[key])
    else
        i = length(internal_vector(pv)) + 1
        push!(internal_vector(pv), value)
        internal_dict(pv)[key] = i
    end
    nothing
end
setproperty!(pv::PVComp, key::Symbol, value) =
    setindex!(pv, value, key)

items(pv::PVComp) =
    view(internal_vector(pv), collect(values(internal_dict(pv))))

function iterate(pv::PVComp, state...)
    out = iterate(internal_dict(pv), state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => getproperty(pv, out[1][1]), out[2])
end





PVModel{T}() where T = PVModel{T}(OrderedDict{Symbol, PVComp{T}}(), Vector{T}())

internal_vector(pv::PVModel) = pv.data

function empty!(pv::PVModel)
    empty!(pv.comps)
    empty!(pv.data)
end

keys(pv::PVModel) = collect(keys(pv.comps))

function getindex(pv::PVModel{T}, key::Symbol) where T
    if !haskey(pv.comps, key)
        pv.comps[key] = PVComp(pv)
    end
    return pv.comps[key]
end

function items(pv::PVModel)
    ids = Vector{Int}()
    for (cname, comp) in pv
        append!(ids, collect(values(internal_dict(comp))))
    end
    return view(internal_vector(pv), ids)
end

iterate(pv::PVModel, state...) = iterate(pv.comps, state...)



# Multi
struct PVMulti{T}
    models::Vector{PVModel{T}}
    PVMulti{T}() where T = new(Vector{PVModel{T}}())
end

empty!(pv::PVMulti) = empty!(pv.models)
push!(pv::PVMulti, m::PVModel) = push!(pv.models, m)

getindex(pv::PVMulti, i::Int) = pv.models[i]
length(pv::PVMulti) = length(pv.models)

# Note: the following returns a copy of the items, rather than a view
function items(pv::PVMulti{T}) where T
    out = Vector{T}()
    for i in 1:length(pv)
        append!(out, items(pv[i]))
    end
    return out
end

function set_items!(pv::PVMulti, values::Vector)
    c = 1
    for i in 1:length(pv)
        v = items(pv[i])
        v .= values[c:(c+length(v)-1)]
        c += length(v)
    end
    nothing
end

end
