abstract type AbstractDomain{N} end

struct Domain{N} <: AbstractDomain{N}
    coords::NTuple{N, Vector{Float64}}

    function Domain(coords::Vararg{AbstractVector{T},N}) where {T <: Real, N}
        (N > 1)  &&  (@assert all(length(coords[1]) .== [length.(coords)...]))
        return new{N}(deepcopy(coords))
    end
    function Domain(length::Integer)
        @assert length >= 1
        return new{1}((collect(1.:length),))
    end
end

ndims(d::Domain{N}) where N = N
length(d::Domain) = length(d.coords[1])
getindex(d::Domain, dim) = d.coords[dim]
function iterate(d::Domain{N}, dim=1) where N
    (dim > N)  &&  (return nothing)
    return (d[dim], dim+1)
end

struct CartesianDomain{N} <: AbstractDomain{N}
    axis::NTuple{N, Vector{Float64}}
    roi::Vector{Int}
    ldomain::Domain{N}

    function CartesianDomain(axis::Vararg{AbstractVector{T},N}; roi=nothing) where {T <: Real, N}
        @assert N >= 2
        len = prod(length.(axis))
        ss = tuple(length.(axis)...)
        isnothing(roi)  &&  (roi = collect(1:len))

        # Prepare corresponding linear domain
        matrix = Matrix{Float64}(undef, N, length(roi))
        ci = Tuple.(CartesianIndices(ss))[roi]
        for i = 1:N
            matrix[i, :] .= axis[i][getindex.(ci, i)]
        end
        ldomain = Domain(getindex.(Ref(matrix), 1:N, :)...)
        return new{N}(deepcopy(axis), roi, ldomain)
    end

    function CartesianDomain(lengths::Vararg{T,N}; kw...) where {T <: Integer, N}
        @assert N >= 2
        @assert all(lengths .>= 1)
        axis = [collect(1.:lengths[i]) for i in 1:N]
        return CartesianDomain(axis...; kw...)
    end
end

# Forward methods to ldomain field
ndims(d::CartesianDomain{N}) where N = N
length(d::CartesianDomain) = length(d.ldomain)
getindex(d::CartesianDomain, dim) = getindex(d.ldomain, dim)
iterate(d::CartesianDomain, args...) = iterate(d.ldomain, args...)

# Cartesian-only methods
size(d::CartesianDomain) = tuple(length.(d.axis)...)
axis(d::CartesianDomain, dim) = d.axis[dim]
roi(d::CartesianDomain) = d.roi


# ====================================================================
# Measures and Counts types
#
abstract type AbstractData{T,N} <: AbstractArray{T,N} end

struct Measures{N} <: AbstractData{Float64,N}
    meta::MDict
    val::Array{Float64,N}
    unc::Array{Float64,N}

    Measures(val::AbstractArray{T, N}, unc::AbstractArray{T, N}) where {T <: Real, N} =
        new{N}(MDict(), deepcopy(val), deepcopy(unc))
end

Measures(val::AbstractArray{T, N}, unc::T) where {T <: Real, N} =
    Measures(val, fill(unc, size(val)))

Measures(val::AbstractArray{T, N}) where {T <: Real, N} =
    Measures(val, one(Float64))


struct Counts{N} <: AbstractData{Int,N}
    meta::MDict
    val::Array{Int,N}

    Counts(val::AbstractArray{T, N}) where {T <: Integer, N} =
        new{N}(MDict(), deepcopy(val))
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
# TODO: check these
function flatten(data::Measures{N}, dom::Domain{N}) where N
    @assert length(dom) == length(data)
    (N == 1)  &&  return data
    return Measures(data.val[:], data.unc[:])
end

function flatten(data::Counts{N}, dom::Domain{N}) where N
    @assert length(dom) == length(data)
    (N == 1)  &&  return data
    return Counts(data.val)
end


function flatten(data::Measures{N}, dom::CartesianDomain{N}) where N
    @assert length(dom) == length(data)
    (N == 1)  &&  return data
    return Measures(data.val[roi(dom)], data.unc[roi(dom)])
end

function flatten(data::Counts{N}, dom::CartesianDomain{N}) where N
    @assert length(dom) == length(data)
    (N == 1)  &&  return data
    return Counts(data.val[roi(dom)])
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
