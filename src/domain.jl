abstract type AbstractDomain{N} end

# A non-cartesian domain has the same number of points in all
# dimensions.
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


# A cartesian domain has the coordinates specified independently for
# each axis, and can always be trasformed into a non-cartesian one.
# Cartesian domains also support region-of-interest (ROI).
struct CartesianDomain{N} <: AbstractDomain{N}
    axis::NTuple{N, Vector{Float64}}
    roi::Vector{Int}
    ldomain::Domain{N}

    function CartesianDomain(axis::Vararg{AbstractVector{T},N}; roi=nothing) where {T <: Real, N}
        @assert N >= 2 "A cartesian domain requires at least 2 dimensions"
        isnothing(roi)  &&  (roi = collect(1:prod(length.(axis))))

        # Pre-compute corresponding linear domain
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
length(d::CartesianDomain) = length(coords(d, 1))

# Return coordinates of all points along a given dimension
coords(d::Domain{1}) = d.axis[1]
coords(d::Domain, dim::Integer) = d.axis[dim]
coords(d::CartesianDomain, dim::Integer) = coords(flatten(d), dim)
getindex(d::Union{Domain, CartesianDomain}, dim::Integer) = coords(d, dim)

# Iterate through domain dimensions returning coordinates
function iterate(d::Union{Domain, CartesianDomain}, ii=1)
    (ii > ndims(d))  &&  (return nothing)
    return (coords(d, ii), ii+1)
end

# Cartesian-only methods
flatten(d::CartesianDomain) = d.ldomain
size(d::CartesianDomain) = tuple([length(v) for v in d.axis]...)
axis(d::CartesianDomain, dim::Integer) = d.axis[dim]



# ====================================================================
# Measures and Counts types
#
abstract type AbstractData{T,N} end

struct Measures{N} <: AbstractData{Float64,N}
    domain::AbstractDomain{N}
    val::Vector{Float64}
    unc::Vector{Float64}

    function Measures(domain::Domain{N}, val::AbstractArray{T, N}, unc::AbstractArray{T, N}) where {T <: Real, N}
        @assert length(domain) == length(val) == length(unc) "Domain and dataset have incompatible lengths"
        return new{N}(deepcopy(domain), deepcopy(val[:]), deepcopy(unc[:]))
    end

    function Measures(domain::CartesianDomain{N}, val::AbstractArray{T, N}, unc::AbstractArray{T, N}) where {T <: Real, N}
        @assert size(domain) == size(val) == size(unc) "Domain and dataset have incompatible size"
        return new{N}(deepcopy(domain), deepcopy(val[domain.roi]), deepcopy(unc[domain.roi]))
    end
end

Measures(domain::AbstractDomain{N}, val::AbstractArray{T, N}, unc::T) where {T <: Real, N} =
    Measures(domain, val, fill(unc, size(val)))


struct Counts{N} <: AbstractData{Int,N}
    domain::AbstractDomain{N}
    val::Vector{Int}

    function Counts(domain::Domain{N}, val::AbstractArray{T, N}) where {T <: Real, N}
        @assert length(domain) == length(val) "Domain and dataset have incompatible lengths"
        return new{N}(deepcopy(domain), deepcopy(val[:]))
    end

    function Counts(domain::CartesianDomain{N}, val::AbstractArray{T, N}) where {T <: Real, N}
        @assert size(domain) == size(val) == size(unc) "Domain and dataset have incompatible size"
        return new{N}(deepcopy(domain), deepcopy(val[domain.roi]))
    end
end

domain(d::AbstractData) = d.domain
values(d::AbstractData) = d.val
uncerts(d::Measures) = d.unc
ndims(d::AbstractData) where N = ndims(domain(d))
length(d::AbstractData) = length(domain(d))
size(d::AbstractData) = size(domain(d))
