import Base.length
import Base.keys
import Base.propertynames
import Base.getindex
import Base.setindex!
import Base.getproperty
import Base.setproperty!
import Base.empty!
import Base.iterate
import Base.values
using DataStructures

# A wrapper for a vector allowing access to elements via a symbol name
struct HashVector{T}
    dict::OrderedDict{Symbol, Int}
    data::Vector{T}

    HashVector{T}() where T =
        new{T}(OrderedDict{Symbol, Int}(), Vector{T}())

    HashVector{T}(data::Vector{T}) where T =
        new{T}(OrderedDict{Symbol, Int}(), data)
end

internal_data(hv::HashVector) = getfield(hv, :data)

function empty!(hv::HashVector)
    dict = getfield(hv, :dict)
    data = getfield(hv, :data)
    deleteat!(data, collect(values(dict)))
    empty!(dict)
    nothing
end

length(hv::HashVector) = length(getfield(hv, :dict))
values(hv::HashVector) = getfield(hv, :data)[collect(values(getfield(hv, :dict)))]
keys(hv::HashVector) = propertynames(hv)
propertynames(hv::HashVector) = collect(keys(getfield(hv, :dict)))

getindex(hv::HashVector, key::Symbol) = getproperty(hv, key)
function getproperty(hv::HashVector, key::Symbol)
    id = getfield(hv, :dict)[key]
    return getfield(hv, :data)[id]
end

setindex!(hv::HashVector, value, key::Symbol) = setproperty!(hv, key, value)
function setproperty!(hv::HashVector{T}, key::Symbol, value::T) where T
    dict = getfield(hv, :dict)
    data = getfield(hv, :data)

    id = get(dict, key, nothing)
    if isnothing(id)
        id = length(data) + 1
        dict[key] = id
        push!(data, value)
    else
        data[id] = value
    end
    return value
end

function iterate(hv::HashVector, state...)
    dict = getfield(hv, :dict)
    data = getfield(hv, :data)
    out = iterate(dict, state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => data[out[1][2]], out[2])
end


# --------------------------------------------------------------------
# A Dict{HashVector{V}} with automatic generation of empty HashVector values

struct HashHashVector{V}
    dict::OrderedDict{Symbol, HashVector{V}}
    data::Vector{V}

    HashHashVector{V}() where V =
        new{V}(OrderedDict{Symbol, HashVector{V}}(), Vector{V}())
end

internal_data(hhv::HashHashVector) = hhv.data
keys(hhv::HashHashVector) = collect(keys(getfield(hhv, :dict)))

function empty!(hhv::HashHashVector)
    empty!(hhv.dict)
    empty!(hhv.data)
end

function getindex(hhv::HashHashVector{V}, key::Symbol) where V
    if !haskey(hhv.dict, key)
        hhv.dict[key] = HashVector{V}(hhv.data)
    end
    return hhv.dict[key]
end

iterate(hhv::HashHashVector, state...) =
    iterate(hhv.dict, state...)
