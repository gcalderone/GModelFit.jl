abstract type AbstractDomain{N} <: AbstractVector{Union{Float64, Vector{Float64}}} end

struct Domain{N} <: AbstractDomain{N}
    meta::MDict
    coords::AbstractMatrix{Float64}

    Domain(coords::Matrix{T}) where {T <: Real} =
        new{size(coords)[2]}(MDict(), coords)

    function Domain(coords::Vararg{AbstractVector{T},N}) where {T <: Real, N}
        (N > 1)  &&  (@assert all(length(coords[1]) .== [length.(coords)...]))
        mat = deepcopy(coords[1])
        if N > 1
            for i in 2:N
                mat = hcat(mat, coords[i])
            end
        else
            mat = reshape(mat, length(mat), 1)
        end
        return new{N}(MDict(), mat)
    end

    function Domain(length::Integer)
        @assert length >= 1
        return new{1}(MDict(), reshape(collect(1.:length), length, 1))
    end
end

ndims(d::Domain{N}) where N = N
length(d::Domain) = size(d.coords)[1]
size(d::Domain{N}) where N = (length(d),)

getindex(d::Domain{1}, ::Colon) = d.coords[:, 1]
getindex(d::Domain{1}, index::Integer) = d.coords[index, 1]
getindex(d::Domain, dim::Integer) = d.coords[:, dim]

function iterate(d::Domain{1}, ii=1)
    (ii > size(d.coords)[1])  &&  (return nothing)
    return (d.coords[ii, 1], ii+1)
end
function iterate(d::Domain, ii=1)
    (ii > size(d.coords)[1])  &&  (return nothing)
    return (d.coords[ii, :], ii+1)
end

coords(d::Domain{N}) where N = [d.coords[:, i] for i in 1:N]


struct CartesianDomain{N} <: AbstractDomain{N}
    meta::MDict
    axis::NTuple{N, Vector{Float64}}
    roi::Vector{Int}
    ldomain::Domain{N}

    function CartesianDomain(axis::Vararg{AbstractVector{T},N}; roi=nothing) where {T <: Real, N}
        @assert N >= 2 "A cartesian domain requires at least 2 dimensions"
        len = prod(length.(axis))
        ss = tuple(length.(axis)...)
        isnothing(roi)  &&  (roi = collect(1:len))

        # Prepare corresponding linear domain
        mat = Matrix{Float64}(undef, length(roi), N)
        ci = Tuple.(CartesianIndices(ss))[roi]
        for i = 1:N
            mat[:, i] .= axis[i][getindex.(ci, i)]
        end
        ldomain = Domain(mat)
        return new{N}(MDict(), deepcopy(axis), roi, ldomain)
    end

    function CartesianDomain(lengths::Vararg{T,N}; kw...) where {T <: Integer, N}
        @assert all(lengths .>= 1)
        axis = [collect(1.:lengths[i]) for i in 1:N]
        return CartesianDomain(MDict(), axis...; kw...)
    end
end

# Forward methods to ldomain field
ndims(d::CartesianDomain{N}) where N = N
length(d::CartesianDomain) = length(d.ldomain)
size(d::CartesianDomain) = size(d.ldomain)
getindex(d::CartesianDomain, ii) = getindex(d.ldomain, ii)
iterate(d::CartesianDomain, args...) = iterate(d.ldomain, args...)
coords(d::CartesianDomain) = coords(d.ldomain)

# Cartesian-only methods
orig_size(d::CartesianDomain) = tuple(length.(d.axis)...)
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
#
function flatten(data::Measures{N}, dom::Domain{N}) where N
    @assert length(dom) == length(data) "Domain and dataset have incompatible lengths"
    (N == 1)  &&  return data
    return Measures(data.val[:], data.unc[:])
end

function flatten(data::Counts{N}, dom::Domain{N}) where N
    @assert length(dom) == length(data) "Domain and dataset have incompatible lengths"
    (N == 1)  &&  return data
    return Counts(data.val)
end

function flatten(data::Measures{N}, dom::CartesianDomain{N}) where N
    @assert (length(dom) <= length(data)  &&  prod(orig_size(dom)) == length(data)) "Domain and dataset have incompatible lengths"
    return Measures(data.val[roi(dom)], data.unc[roi(dom)])
end

function flatten(data::Counts{N}, dom::CartesianDomain{N}) where N
    @assert (length(dom) <= length(data)  &&  prod(orig_size(dom)) == length(data)) "Domain and dataset have incompatible lengths"
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
