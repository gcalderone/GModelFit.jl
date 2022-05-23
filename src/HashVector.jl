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

struct HashVector{V}
    dict::OrderedDict{Symbol, Int}
    data::Vector{V}

    HashVector{V}() where V =
        new{V}(OrderedDict{Symbol, Int}(), Vector{V}())

    HashVector{V}(data::Vector{V}) where V =
        new{V}(OrderedDict{Symbol, Int}(), data)
end

internal_dict(hv::HashVector) = getfield(hv, :dict)
internal_data(hv::HashVector) = getfield(hv, :data)
indices(hv::HashVector) = collect(values(internal_dict(hv)))

length(hv::HashVector) = length(internal_dict(hv))
keys(hv::HashVector{V}) where V = keys(internal_dict(hv))
values(hv::HashVector{V}) where V = internal_data(hv)[indices(hv)]  # view(internal_data(hv), indices(hv))
propertynames(hv::HashVector) = keys(internal_dict(hv))


function getindex(hv::HashVector, key::Symbol)
    id = internal_dict(hv)[key]
    return internal_data(hv)[id]
end

function setindex!(hv::HashVector{V}, value::V, key::Symbol) where V
    if haskey(internal_dict(hv), key)
        id = internal_dict(hv)[key]
        internal_data(hv)[id] = value
    else
        id = length(internal_data(hv)) + 1
        internal_dict(hv)[key] = id
        push!(internal_data(hv), value)
    end
    return value
end

getproperty(hv::HashVector, key::Symbol) = getindex(hv, key)
setproperty!(hv::HashVector{V}, key::Symbol, value::V) where V = setindex!(hv, value, key)

function iterate(hv::HashVector, state...)
    out = iterate(internal_dict(hv), state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => internal_data(hv)[out[1][2]], out[2])
end


function empty!(hv::HashVector)
    deleteat!(internal_data(hv), indices(hv))
    empty!(internal_dict(hv))
    nothing
end



# ====================================================================
struct HashHashVector{V}
    dict::OrderedDict{Symbol, HashVector{V}}
    data::Vector{V}

    HashHashVector{V}() where V =
        new{V}(OrderedDict{Symbol, HashVector{V}}(), Vector{V}())

    HashHashVector{V}(data::Vector{V}) where V =
        new{V}(OrderedDict{Symbol, HashVector{V}}(), data)
end


internal_data(hhv::HashHashVector) = hhv.data


function getindex(hhv::HashHashVector{V}, key::Symbol) where V
    if !haskey(hhv.dict, key)
        hhv.dict[key] = HashVector{V}(hhv.data)
    end
    return hhv.dict[key]
end


iterate(hhv::HashHashVector, state...) = 
    iterate(hhv.dict, state...)


function empty!(hhv::HashHashVector)
    for (key, hv) in hhv.dict
        empty!(internal_dict(hv))
    end
    empty!(hhv.data)
    nothing
end
