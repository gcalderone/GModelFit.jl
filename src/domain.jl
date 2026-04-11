abstract type AbstractDomain{N} end

# A non-cartesian domain has the same number of points in all
# dimensions.

"""
    Domain{N}

An object representing a N-dimensional linear domain.

Available constructors:
- `Domain(coords...)`: each argument is a vector, one for each dimension (all arguments must have same lengths);
- `Domain(length)`: returns a 1-dim `Domain` object of the given length.

Coordinates for all points along a given axis can be obtained with the `coords` function.
"""
struct Domain{N} <: AbstractDomain{N}
    coords::NTuple{N, Vector{Float64}}

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
"""
    CartesianDomain{N}

An object representing a model, or a dataset, N-dimensional cartesian domain (i.e. a grid).

Available constructors:
- `CartesianDomain(axis...)`: each argument is a vector containing the coordinates on a given axis (arguments may have different lengths);
- `CartesianDomain(lengths...)`: returns a N-dim `CartesianDomain` object whose axis lengths are specified in the arguments.

Note that a `CartesianDomain` requires at least 2 dimensions.

Coordinates for all points along a given axis can be obtained with the `coords()` function, while the coordinates of the grid can be obtained with `gridcoords()`.
"""
struct CartesianDomain{N} <: AbstractDomain{N}
    coords::NTuple{N, Vector{Float64}}

    function CartesianDomain(coords::Vararg{AbstractVector{T},N}) where {T <: Real, N}
        @assert N >= 2 "A cartesian domain requires at least 2 dimensions"
        return new{N}(deepcopy(coords))
    end

    function CartesianDomain(lengths::Vararg{T,N}) where {T <: Integer, N}
        @assert N >= 2 "A cartesian domain requires at least 2 dimensions"
        @assert all(lengths .>= 1)
        coords = [collect(1.:lengths[i]) for i in 1:N]
        return CartesianDomain(coords...)
    end
end

ndims(d::Union{Domain{N}, CartesianDomain{N}}) where N = N

length(d::Domain) = length(d.coords[1])
length(d::CartesianDomain) = prod(length.(d.coords))

"""
    coords(d::Domain{1})
    coords(d::Domain, dim::Integer)
    coords(d::CartesianDomain, dim::Integer)

Returns coordinates of all points along a given dimension as a `Vector{Float64}`.
"""
coords(d::Domain{1}) = d.coords[1]
coords(d::Domain, dim::Integer) = d.coords[dim]

# Iterate through domain dimensions returning coordinates
function iterate(d::Domain, ii=1)
    (ii > ndims(d))  &&  (return nothing)
    return (coords(d, ii), ii+1)
end

function iterate(d::CartesianDomain, ii=1)
    (ii > ndims(d))  &&  (return nothing)
    return (gridcoords(d, ii), ii+1)
end

# Cartesian-only methods
size(d::CartesianDomain) = tuple([length(v) for v in d.coords]...)

function Domain(d::CartesianDomain{N}) where N
    v = reshape(collect(Iterators.product(d.coords...)), :)
    return Domain([getindex.(v, i) for i in 1:N]...)
end

"""
    gridcoords(d::CartesianDomain, dim::Integer)

Returns the coordinates of the grid along a given dimension as a `Vector{Float64}`.
"""
gridcoords(d::CartesianDomain, dim::Integer) = d.coords[dim]

domain2buffer(::Type{T}, domain::Domain)             where {T <: Real}    = Array{T, 1}(undef, length(domain))
domain2buffer(::Type{T}, domain::CartesianDomain{N}) where {T <: Real, N} = Array{T, N}(undef,   size(domain))

# ====================================================================
abstract type AbstractResponse end

struct IdentityResp <: AbstractResponse
    domain::AbstractDomain
end
folded_domain(  resp::IdentityResp) = resp.domain
unfolded_domain(resp::IdentityResp) = resp.domain
fold_model!(    resp::IdentityResp, unfolded::Array, folded::Array) = folded .= unfolded

# ====================================================================
# Data types, used to identify the data statistics and the response
#
abstract type AbstractData{N, R <: AbstractResponse} end

"""
    getresp(d::AbstractData)

Return the instrument response associated to an AbstractData object.
"""
getresp(d::AbstractData) = d.resp

"""
    domain(d::AbstractData)

Return the domain associated to an AbstractData object.
"""
domain(d::AbstractData) = folded_domain(getresp(d))
ndims(d::AbstractData) = ndims(domain(d))
length(d::AbstractData) = length(domain(d))
size(d::AbstractData) = size(domain(d))

"""
    GaussianData{N}

An object representing a set of empirical measurements (with Gaussian uncertainties) as measured on a specific domain.

Available constructors:
- `GaussianData(domain::Domain{N},
            values::AbstractVector{T},
            uncerts::AbstractVector{T}) where {T <: AbstractFloat, N}`
- `GaussianData(domain::CartesianDomain{N},
            values::AbstractArray{T, N},
            uncerts::AbstractArray{T, N}) where {T <: AbstractFloat, N}`
- `GaussianData(values::AbstractVector, uncerts)`

In the above constructor methods the last argument may also be a scalar value, to set the same uncertainty for all the measurements. The method accepting a `CartesianDomain` requires arrays with at least 2 dimensions.  In the last constructor the `Domain` object is automatically built depending on the length of the `values` vector.

The domain, values and uncertainties for a `GaussianData` object can be retrieved using the `domain`, `values` and `uncerts` functions respectively.
"""
struct GaussianData{N, R <: AbstractResponse} <: AbstractData{N, R}
    resp::R
    values::Array{Float64, N}
    uncerts::Array{Float64, N}

    function GaussianData(resp::R, values::AbstractArray{T}, uncerts::AbstractArray{T}) where {R <: AbstractResponse, T <: AbstractFloat}
        dom = folded_domain(resp)
        @assert ndims(dom) == ndims(values) == ndims(uncerts)
        if isa(dom, Domain)
            @assert length(dom) == length(values) == length(uncerts) "Domain and dataset have incompatible length"
        else
            @assert size(dom) == size(values) == size(uncerts) "Domain and dataset have incompatible size"
        end
        return new{ndims(values), R}(resp, values, uncerts)
    end
end

GaussianData(values::AbstractVector, uncerts) = GaussianData(IdentityResp(Domain(length(values))), values, uncerts)
GaussianData(domain::Domain{N}, values::AbstractVector{T}, uncerts::AbstractVector{T}) where {T <: AbstractFloat, N} = GaussianData(IdentityResp(domain), values, uncerts)
GaussianData(domain::CartesianDomain{N}, values::AbstractArray{T, N}, uncerts::AbstractArray{T, N}) where {T <: AbstractFloat, N} = GaussianData(IdentityResp(domain), values, uncerts)
GaussianData(domain::AbstractDomain{N}, values::AbstractArray{T, N}, uncert::T) where {T <: AbstractFloat, N} = GaussianData(domain, values, fill(uncert, size(values)))
GaussianData(resp::AbstractResponse, values::AbstractArray, uncert::Real) = GaussianData(resp, values, fill(uncert, size(values)))
GaussianData(dom::AbstractDomain, values::AbstractArray, uncert::Real) = GaussianData(dom, values, fill(uncert, size(values)))

"""
    values(d::GaussianData)

Returns the measurement values as a `Vector{Float64}`.
"""
values(d::GaussianData) = d.values

"""
    uncerts(d::GaussianData)

Returns the measurement uncertainties as a `Vector{Float64}`.
"""
uncerts(d::GaussianData) = d.uncerts
getresp(d::GaussianData) = d.resp

const Measures = GaussianData

# --------------------------------------------------------------------
struct PoissonCounts{N, R <: AbstractResponse} <: AbstractData{N, R}
    resp::R
    values::Array{Float64}

    function PoissonCounts(resp::R, values::AbstractArray{T}) where {R <: AbstractResponse, T <: Integer}
        dom = folded_domain(resp)
        @assert ndims(dom) == ndims(values)
        if isa(dom, Domain)
            @assert length(dom) == length(values) "Domain and dataset have incompatible length"
        else
            @assert size(dom) == size(values) "Domain and dataset have incompatible size"
        end
        return new{ndims(values), R}(resp, convert(Array{Float64}, values))
    end
end

PoissonCounts(values::AbstractVector) = PoissonCounts(IdentityResp(Domain(length(values))), values)
PoissonCounts(domain::Domain{N}, values::AbstractVector{T}) where {T <: Integer, N} = PoissonCounts(IdentityResp(domain), values)
PoissonCounts(domain::CartesianDomain{N}, values::AbstractArray{T, N}) where {T <: Integer, N} = PoissonCounts(IdentityResp(domain), values)

values(d::PoissonCounts) = d.values
getresp(d::PoissonCounts) = d.resp
