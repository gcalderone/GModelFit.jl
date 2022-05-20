import Base.length
import Base.propertynames
import Base.getproperty
import Base.setproperty!
import Base.push!
import Base.iterate
using DataStructures

struct HashVector{V}
    dict::OrderedDict{Symbol, Union{Int, Vector{Int}}}
    data::Vector{V}

    HashVector{V}() where V =
        new{V}(OrderedDict{Symbol, Union{Int, Vector{Int}}}(), Vector{V}())

    HashVector{V}(data::Vector{V}) where V =
        new{V}(OrderedDict{Symbol, Union{Int, Vector{Int}}}(), data)
end

length(hv::HashVector) = length(getfield(hv, :dict))
propertynames(hv::HashVector) = keys(getfield(hv, :dict))

function getproperty(hv::HashVector{V}, key::Symbol) where V
    id = getfield(hv, :dict)[key]
    if isa(id, Array)
        return view(getfield(hv, :data), id)
    end
    return getfield(hv, :data)[id]
end

function setproperty!(hv::HashVector{V}, key::Symbol, value::V) where V
    id = get(getfield(hv, :dict), key, nothing)

    # Avoid raising an error here to simplify patch functions
    isnothing(id)  &&  (return value)

    @assert !isa(id, Array) "Can't set a scalar value to a Vector{$V}: try with dot (broadcast) notation."
    getfield(hv, :data)[id] = value
end

function push!(hv::HashVector{V}, k::Symbol, v::V) where V
    dict = getfield(hv, :dict)
    data = getfield(hv, :data)
    
    i1 = length(data) + 1
    if k in keys(dict)
        if isa(dict[k], Int)
            dict[k] = [dict[k], i1]
        else
            push!(dict[k], i1)
        end
    else
        dict[k] = i1
    end
    push!(getfield(hv, :data), v)
end


function iterate(hv::HashVector{V}, state...) where V
    out = iterate(getfield(hv, :dict), state...)
    isnothing(out)  &&  (return nothing)
    return (out[1][1] => getproperty(hv, out[1][1]), out[2])
end


#=
hv = GFit.HashVector{Int}()
push!(hv, :a, 1)
push!(hv, :b, 2)
push!(hv, :b, 3)

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



#=
function HashVector(in_dict::AbstractDict{Symbol, T}) where T
    Eltypes = unique(eltype.(values(in_dict)))
    @assert length(Eltypes) == 1
    V = Eltypes[1]
    dict = OrderedDict{Symbol, Union{Int, Vector{Int}}}()
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
    return HashVector{V}(dict, data)
end


function force_vector(hv::HashVector{V}, k::Symbol) where V
    dict = getfield(hv, :dict)
    @assert k in keys(dict)
    if isa(dict[k], Int)
        dict[k] = [dict[k]]
    end
    return hv
end

=#
