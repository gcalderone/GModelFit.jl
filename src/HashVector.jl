import Base.length
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

internal_dict(hv::HashVector{V}) where V = getfield(hv, :dict)
internal_data(hv::HashVector{V}) where V = getfield(hv, :data)


length(hv::HashVector) = length(internal_dict(hv))
propertynames(hv::HashVector) = keys(internal_dict(hv))

values(hv::HashVector{V}) where V = internal_data(hv)[collect(values(internal_dict(hv)))]

getindex(hv::HashVector{V}, key::Symbol) where V = getproperty(hv, key)
function getproperty(hv::HashVector{V}, key::Symbol) where V
    id = internal_dict(hv)[key]
    return internal_data(hv)[id]
end

setindex!(hv::HashVector{V}, value::V, key::Symbol) where V = setproperty!(hv, key, value)
function setproperty!(hv::HashVector{V}, key::Symbol, value::V) where V
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

function iterate(hv::HashVector{V}, state...) where V
    out = iterate(internal_dict(hv), state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => getproperty(hv, out[1][1]), out[2])
end


function empty!(hv::HashVector{V}) where V
    empty!(internal_dict(hv))
    empty!(internal_data(hv))
    nothing
end
