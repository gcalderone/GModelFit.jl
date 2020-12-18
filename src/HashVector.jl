import Base.length
import Base.propertynames
import Base.getproperty
import Base.setproperty!
import Base.iterate
using DataStructures

struct HashVector{K,V}
    dict::OrderedDict{K, Union{Int, Vector{Int}}}
    data::Vector{V}
end

function HashVector(in_dict::AbstractDict{K, T}) where {K,T}
    Eltypes = unique(eltype.(values(in_dict)))
    @assert length(Eltypes) == 1
    V = Eltypes[1]
    dict = OrderedDict{K, Union{Int, Vector{Int}}}()
    data = Vector{V}()
    for (key, values) in in_dict
        i0 = length(data) + 1
        i1 = i0 + length(values) - 1
        append!(data, values)
        if isa(values, AbstractArray)
            dict[key] = collect(i0:i1)
        else
            dict[key] = i0
        end
    end
    return HashVector{K,V}(dict, data)
end

length(hv::HashVector) = length(getfield(hv, :dict))
propertynames(hv::HashVector) = keys(getfield(hv, :dict))

function _getproperty(hv::HashVector{K,V}, key::K) where {K,V}
    id = getfield(hv, :dict)[key]
    if isa(id, Array)
        return view(getfield(hv, :data), id)
    end
    return getfield(hv, :data)[id]
end

function _setproperty!(hv::HashVector{K,V}, key::K, value::V) where {K,V}
    id = get(getfield(hv, :dict), key, nothing)

    # Avoid raising an error here to simplify patch functions
    isnothing(id)  &&  (return value)

    @assert !isa(id, Array) "Can't set a scalar value to a Vector{$V}: try with dot (broadcast) notation."
    getfield(hv, :data)[id] = value
end

function iterate(hv::HashVector{K,V}, state...) where {K,V}
    out = iterate(getfield(hv, :dict), state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => getproperty(hv, out[1][1]), out[2])
end


#=
dd = Dict(:a => 1, :b => [2,3]);
hv = HashVector(dd)

@assert length(hv) == 2
@assert all(propertynames(hv) .== [:a, :b])

@assert hv.a == 1
@assert hv.b == [2,3]

hv.a = 10
@assert hv.a == 10

hv.b .= 10
@assert hv.b == [10,10]

for (key, val) in hv
    println(key, "  ", val)
end
=#
