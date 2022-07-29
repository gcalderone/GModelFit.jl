abstract type AbstractDomain{N} end

# A non-cartesian domain has the same number of points in all dimensions.
struct Domain{N} <: AbstractDomain{N}
    axis::NTuple{N, Vector{Float64}}

    function Domain(coords::Vararg{AbstractVector{T},N}) where {T <: Real, N}
        @assert N >= 1
        @assert all(length(coords[1]) .== [length.(coords)...])
        return new{N}(convert(NTuple{N, Vector{Float64}}, coords))
    end

    function Domain(length::Integer)
        @assert length >= 1
        return Domain(collect(1.:length))
    end
end


# A cartesian domain has the coordinates specified independently for each axis, and can always be trasformed into a non-cartesian one.  Cartesian domains also supports region-of-interest (ROI),
struct CartesianDomain{N} <: AbstractDomain{N}
    axis::NTuple{N, Vector{Float64}}
    roi::Vector{Int}
    ldomain::Domain{N}

    function CartesianDomain(axis::Vararg{AbstractVector{T},N}; roi=nothing) where {T <: Real, N}
        @assert N >= 2 "A cartesian domain requires at least 2 dimensions"
        isnothing(roi)  &&  (roi = collect(1:prod(length.(axis))))

        # Prepare corresponding linear domain
        ss = tuple(length.(axis)...)
        ci = Tuple.(CartesianIndices(ss))[roi]
        vv = Vector{Vector{Float64}}()
        for i = 1:N
            push!(vv, axis[i][getindex.(ci, i)])
        end
        ldomain = Domain(vv...)
        return new{N}(deepcopy(axis), roi, ldomain)
    end

    function CartesianDomain(lengths::Vararg{T,N}; kw...) where {T <: Integer, N}
        @assert all(lengths .>= 1)
        axis = [collect(1.:lengths[i]) for i in 1:N]
        return CartesianDomain(axis...; kw...)
    end
end


ndims(d::Union{Domain{N}, CartesianDomain{N}}) where N = N

length(d::Domain) = length(d.axis[1])
length(d::CartesianDomain) = length(linear_domain(d))
size(d::CartesianDomain) = [length(v) for v in d.axis]

getindex(d::Domain{1}, ::Colon) = d.axis[1]
getindex(d::Union{Domain, CartesianDomain}, dim::Integer) = d.axis[dim]

function iterate(d::Union{Domain, CartesianDomain}, ii=1)
    (ii > ndims(d))  &&  (return nothing)
    return (d[ii], ii+1)
end

linear_domain(d::CartesianDomain) = d.ldomain


# ====================================================================
# Measures and Counts types
#
abstract type AbstractData{T,N} end

struct Measures{N} <: AbstractData{Float64,N}
    domain::AbstractDomain{N}
    val::Array{Float64,N}
    unc::Array{Float64,N}

    function Measures(domain::Domain{N}, val::AbstractArray{T, N}, unc::AbstractArray{T, N}) where {T <: Real, N}
        @assert length(domain) == length(val) "Domain and dataset have incompatible lengths"
        @assert length(unc) == length(val) "Dataset values and uncertainties have incompatible lengths"
        return new{N}(deepcopy(domain), deepcopy(val), deepcopy(unc))
    end

    function Measures(domain::CartesianDomain{N}, val::AbstractArray{T, N}, unc::AbstractArray{T, N}) where {T <: Real, N}
        @assert prod(orig_size(domain)) == length(val) "Domain and dataset have incompatible lengths"
        return new{N}(deepcopy(domain), deepcopy(val), deepcopy(unc))
    end
end

Measures(domain::AbstractDomain{N}, val::AbstractArray{T, N}, unc::T) where {T <: Real, N} =
    Measures(domain, val, fill(unc, size(val)))


struct Counts{N} <: AbstractData{Int,N}
    domain::AbstractDomain{N}
    val::Array{Int,N}

    function Counts(domain::Domain{N}, val::AbstractArray{T, N}) where {T <: Integer, N}
        @assert length(domain) == length(val) "Domain and dataset have incompatible lengths"
        new{N}(deepcopy(domain), deepcopy(val))
    end

    function Counts(domain::CartesianDomain{N}, val::AbstractArray{T, N}) where {T <: Real, N}
        @assert prod(orig_size(domain)) == length(val) "Domain and dataset have incompatible lengths"
        return new{N}(deepcopy(domain), deepcopy(val))
    end
end


ndims(d::Measures{N}) where N = N
size(d::Measures) = size(d.val)
length(d::Measures) = length(d.val)
uncerts(d::Measures) = d.uncerts

ndims(d::Counts{N}) where N = N
size(d::Counts) = size(d.val)
iterate(d::Counts, args...) = iterate(d.val, args...)


# ====================================================================
# Methods to "flatten" a multidimensional object into a 1D one
#
function flatten(data::Measures{N}) where N
    (N == 1)  &&  return data
    if isa(data.domain, CartesianDomain)
        return Measures(data.val[roi(data.domain)], data.unc[roi(data.domain)])
    end
    return Measures(data.val[:], data.unc[:])
end

function flatten(data::Counts{N}) where N
    (N == 1)  &&  return data
    if isa(data.domain, CartesianDomain)
        return Counts(data.val[roi(data.domain)])
    end
    return Counts(data.val[:])
end
