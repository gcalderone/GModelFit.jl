abstract type AbstractDomain{N} end

struct Domain{N} <: AbstractDomain{N}
    axis::NTuple{N, Vector{Float64}}

    function Domain(axis::Vararg{AbstractVector{T},N}) where {T <: Real, N}
        if N > 1
            @assert all(length(axis[1]) .== [length.(axis)...])
        end
        return new{N}(deepcopy(axis))
    end
end

ndims(dom::Domain{N}) where N = N
length(dom::Domain) = length(dom.axis[1])
getindex(dom::Domain, dim) = dom.axis[dim]


struct CartesianDomain{N} <: AbstractDomain{N}
    axis::NTuple{N, Vector{Float64}}
    size::NTuple{N, Int}
    index::Vector{Int}

    function CartesianDomain(axis::Vararg{AbstractVector{T},N}; index=nothing) where {T <: Real, N}
        @assert N >= 2
        isnothing(index)  &&  (index = collect(1:prod(length.(axis))))
        return new{N}(deepcopy(axis), tuple(length.(axis)...), index)
    end
end

ndims(dom::CartesianDomain{N}) where N = N
length(dom::CartesianDomain) = length(dom.index)
size(dom::CartesianDomain) = dom.size
getindex(dom::CartesianDomain, dim) = dom.axis[dim]


# ====================================================================
# Measures and Counts types
#
abstract type AbstractData{T,N} <: AbstractArray{T,N} end

struct Measures{N} <: AbstractData{Float64,N}
    meta::MDict
    val::Array{Float64,N}
    unc::Array{Float64,N}

    Measures(val::AbstractArray{T, N}, unc::AbstractArray{T, N}) where {T <: Real, N} =
        new{N}(MDict(),
            deepcopy(convert(Array{Float64, N}, val)),
            deepcopy(convert(Array{Float64, N}, val)))
end

Measures(val::AbstractArray{T, N}, unc::T) where {T <: Real, N} =
    Measures(val, fill(unc, size(val)))

Measures(val::AbstractArray{T, N}) where {T <: Real, N} =
    Measures(val, one(Float64))


struct Counts{N} <: AbstractData{Int,N}
    meta::MDict
    val::Array{Int,N}

    Counts(val::AbstractArray{T, N}) where {T <: Integer, N} =
        new{N}(MDict(),
            deepcopy(convert(Array{Int, N}, val)))
end


ndims(d::Measures{N}) where N = N
size(d::Measures) = size(d.val)
length(d::Measures) = length(d.val)
getindex(d::Measures, args...) where N = getindex(d.val, args...)
iterate(d::Measures, args...) = iterate(d.val, args...)
uncerts(d::Measures) = d.uncerts

ndims(d::Counts{N}) where N = N
size(d::Counts) = size(d.val)
length(d::Counts) = length(d.val)
getindex(d::Counts, args...) where N = getindex(d.val, args...)
iterate(d::Counts, args...) = iterate(d.val, args...)


# ====================================================================
# Methods to "flatten" a multidimensional object into a 1D one
#
flatten(dom::Domain) = dom

function flatten(dom::CartesianDomain{N}) where N
    out = Matrix{Float64}(undef, N, length(dom.index))
    ci = Tuple.(CartesianIndices(dom.size))[dom.index]
    for i = 1:N
        out[i, :] .= dom.axis[i][getindex.(ci, i)]
    end
    return Domain(getindex.(Ref(out), 1:N, :)...)
end

function flatten(data::Measures{N}, dom::Domain{N}) where N
    @assert length(domain) == length(data)
    (N == 1)  &&  return data
    return Measures(data.val[:], data.unc[:])
end

function flatten(data::Counts  {N}, dom::Domain{N}) where N
    @assert length(domain) == length(data)
    (N == 1)  &&  return data
    return Counts(data.val)
end


#=
"""
# reshape

Reshape an array according to the size of a CartesianDomain object
"""
function reshape(array::AbstractArray, dom::AbstractCartesianDomain)
    @assert length(array) == length(dom) "Domain and array must have the same length"
    out = Array{Float64}(undef, size(dom)...)
    out .= NaN
    out[dom.index] .= array
    return out
end
=#
