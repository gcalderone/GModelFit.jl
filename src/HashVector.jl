import Base.length
import Base.propertynames
import Base.getproperty
import Base.setproperty!
import Base.push!
import Base.iterate
using DataStructures

struct HashVector{V}
    dict::OrderedDict{Symbol, Int}
    data::Vector{V}

    HashVector{V}() where V =
        new{V}(OrderedDict{Symbol, Int}(), Vector{V}())

    function HashVector{V}(data::Vector{V}) where V
        @assert length(data) == 0
        new{V}(OrderedDict{Symbol, Int}(), data)
    end
end

length(hv::HashVector) = length(getfield(hv, :dict))
propertynames(hv::HashVector) = keys(getfield(hv, :dict))

function getproperty(hv::HashVector{V}, key::Symbol) where V
    id = getfield(hv, :dict)[key]
    return getfield(hv, :data)[id]
end

function setproperty!(hv::HashVector{V}, key::Symbol, value::V) where V
    if haskey(getfield(hv, :dict), key)
        id = getfield(hv, :dict)[key]
        getfield(hv, :data)[id] = value
    else
        id = length(getfield(hv, :data)) + 1
        getfield(hv, :dict)[key] = id
        push!(getfield(hv, :data), value)
    end
    return value
end

function iterate(hv::HashVector{V}, state...) where V
    out = iterate(getfield(hv, :dict), state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => getproperty(hv, out[1][1]), out[2])
end

