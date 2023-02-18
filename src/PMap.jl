
module PMap

using DataStructures

import Base.getindex, Base.setindex!,
Base.getproperty, Base.setproperty!, Base.propertynames,
Base.length, Base.keys, Base.empty!, Base.iterate

export  PMapComponent, PMapModel, PMapMultiModel, items, set_items!

struct PMapComponent{T}
    params::OrderedDict{Symbol, Int}
    data::Vector{T}
end

struct PMapModel{T}
    comps::OrderedDict{Symbol, PMapComponent}
    data::Vector{T}
end


PMapComponent(parent::PMapModel{T}) where T = PMapComponent{T}(OrderedDict{Symbol, Int}(), parent.data)
internal_vector(pmap::PMapComponent)  = getfield(pmap, :data)
internal_dict(  pmap::PMapComponent)  = getfield(pmap, :params)
propertynames(pmap::PMapComponent) = collect(keys(internal_dict(pmap)))

getindex(pmap::PMapComponent, key::Symbol) =
    internal_vector(pmap)[internal_dict(pmap)[key]]
getproperty(pmap::PMapComponent, key::Symbol) = getindex(pmap, key)

function setindex!(pmap::PMapComponent, value, key::Symbol)
    if haskey(internal_dict(pmap), key)
        setindex!(internal_vector(pmap), value, internal_dict(pmap)[key])
    else
        i = length(internal_vector(pmap)) + 1
        push!(internal_vector(pmap), value)
        internal_dict(pmap)[key] = i
    end
    nothing
end
setproperty!(pmap::PMapComponent, key::Symbol, value) =
    setindex!(pmap, value, key)

items(pmap::PMapComponent) =
    view(internal_vector(pmap), collect(values(internal_dict(pmap))))

function iterate(pmap::PMapComponent, state...)
    out = iterate(internal_dict(pmap), state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => getproperty(pmap, out[1][1]), out[2])
end





PMapModel{T}() where T = PMapModel{T}(OrderedDict{Symbol, PMapComponent{T}}(), Vector{T}())

internal_vector(pmap::PMapModel) = pmap.data

function empty!(pmap::PMapModel)
    empty!(pmap.comps)
    empty!(pmap.data)
end

keys(pmap::PMapModel) = collect(keys(pmap.comps))

function getindex(pmap::PMapModel{T}, key::Symbol) where T
    if !haskey(pmap.comps, key)
        pmap.comps[key] = PMapComponent(pmap)
    end
    return pmap.comps[key]
end

function items(pmap::PMapModel)
    ids = Vector{Int}()
    for (cname, comp) in pmap
        append!(ids, collect(values(internal_dict(comp))))
    end
    return view(internal_vector(pmap), ids)
end

iterate(pmap::PMapModel, state...) = iterate(pmap.comps, state...)



# MultiModel
struct PMapMultiModel{T}
    models::Vector{PMapModel{T}}
end

getindex(pmap::PMapMultiModel, i::Int) = pmap.models[i]
length(pmap::PMapMultiModel) = length(pmap.models)

# Note: the following returns a copy of the items, rather than a view
function items(pmap::PMapMultiModel{T}) where T
    out = Vector{T}()
    for i in 1:length(pmap)
        append!(out, items(pmap[i]))
    end
    return out
end

function set_items!(pmap::PMapMultiModel, values::Vector)
    c = 1
    for i in 1:length(pmap)
        v = items(pmap[i])
        v .= values[c:(c+length(v)-1)]
        c += length(v)
    end
    nothing
end

end
