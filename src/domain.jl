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
        @assert N >= 2 "A cartesian domain requires at least 2 dimensions"
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

abstract type AbstractMeasures{N} end

domain(d::AbstractMeasures) = d.domain
ndims(d::AbstractMeasures) = ndims(domain(d))
length(d::AbstractMeasures) = length(domain(d))
size(d::AbstractMeasures) = size(domain(d))

values(d::AbstractMeasures, index=1) = d.values[index]
function original_shape(d::AbstractMeasures{N}, index=1) where N
    if isa(domain(d), CartesianDomain)
        out = fill(NaN, size(d))
        out[domain(d).roi] .= values(d, index)
        return out
    end
    out = values(d, index)
    return out
end


struct Measures{N} <: AbstractMeasures{N}
    domain::AbstractDomain{N}
    values::NTuple{2, Vector{Float64}}
    labels::NTuple{2, String}

    # Measures with linear domain are built using 1D vector(s).
    function Measures(domain::Domain{N}, values::AbstractVector{T}, uncerts::AbstractVector{T}) where {T <: AbstractFloat, N}
        @assert length(domain) == length(values) == length(uncerts) "Domain and dataset have incompatible length"
        return new{N}(deepcopy(domain), tuple(deepcopy(values), deepcopy(uncerts)), ("values", "uncerts"))
    end

    # Measures with cartesian domain are built using N-dim arrays(s).
    function Measures(domain::CartesianDomain{N}, values::AbstractArray{T, N}, uncerts::AbstractArray{T, N}) where {T <: AbstractFloat, N}
        @assert size(domain) == size(values) == size(uncerts) "Domain and dataset have incompatible size"
        return new{N}(deepcopy(domain), tuple(deepcopy(values[domain.roi]), deepcopy(uncerts[domain.roi])), ("values", "uncerts"))
    end
end
uncerts(d::Measures) = d.values[2]
Measures(dom::AbstractDomain, values::AbstractArray, uncert::Real) = Measures(dom, values, fill(uncert, size(values)))


struct PoissonCounts{N} <: AbstractMeasures{N}
    domain::AbstractDomain{N}
    values::NTuple{1, Vector{Int}}
    labels::NTuple{1, String}

    # Measures with linear domain are built using 1D vector(s).
    function PoissonCounts(domain::Domain{N}, values::AbstractVector{T}) where {T <: Integer, N}
        @assert length(domain) == length(values) "Domain and dataset have incompatible length"
        return new{N}(deepcopy(domain), tuple(deepcopy(values)), ("counts", ))
    end

    # Measures with cartesian domain are built using N-dim arrays(s).
    function PoissonCounts(domain::CartesianDomain{N}, values::AbstractArray{T, N}) where {T <: Integer, N}
        @assert size(domain) == size(values) "Domain and dataset have incompatible size"
        return new{N}(deepcopy(domain), tuple(deepcopy(values[domain.roi])), ("counts", ))
    end
end
