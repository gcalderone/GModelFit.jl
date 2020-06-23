# ====================================================================
# Abstract types
#
abstract type AbstractDomain end
abstract type AbstractLinearDomain    <: AbstractDomain end
abstract type AbstractCartesianDomain <: AbstractDomain end

abstract type AbstractData end
abstract type AbstractMeasures <: AbstractData end
abstract type AbstractCounts   <: AbstractData end


# ====================================================================
# Macro to define domain, data and associated methods for any
# dimensionality
#
macro define_ndim(ndim::Int)
    @assert ndim >= 1 "Number of dimensions must be >= 1"
    out = Expr(:block)

    # Structure def.
    if ndim > 1
        name = Symbol(:CartesianDomain_, ndim, :D)
        push!(out.args, :(
            struct $name <: AbstractCartesianDomain
                axis::Vector{Vector{Float64}}
                size::NTuple{$(ndim), Int}
                index::Vector{Int}
            end;
            Base.ndims(dom::$name) = $ndim;
            Base.size(dom::$name) = dom.size;
            Base.size(dom::$name, dim::Int) = dom.size[dim];
            Base.getindex(dom::$name, dim::Int) = dom.axis[dim];
        ))

        # Constructors
        a = ["d$(i)::AbstractVector{Float64}" for i in 1:ndim]
        b = ["deepcopy(d$(i))" for i in 1:ndim]
        c = ["length(d$(i))" for i in 1:ndim]
        s = "CartesianDomain(" * join(a, ", ") * "; index=1:" * join(c, " * ") *
            ") = CartesianDomain_$(ndim)D([" * join(b, ", ") * "], (" * join(c, ", ") * ",), index)"
        push!(out.args, Meta.parse(s))

        a = ["n$(i)::Int" for i in 1:ndim]
        b = ["one(Float64):one(Float64):n$(i)" for i in 1:ndim]
        c = ["n$(i)" for i in 1:ndim]
        s = "CartesianDomain(" * join(a, ", ") * "; index=1:" * join(c, " * ") *
            ") = CartesianDomain(" * join(b, ", ") * ", index=index)"
        push!(out.args, Meta.parse(s))
    end

    # Structure def.
    name = Symbol(:Domain_, ndim, :D)
    tmp = (ndim > 1  ?  2  :  1)
    push!(out.args, :(
        struct $name <: AbstractLinearDomain
            axis::Array{Float64, $tmp}
            vmin::Vector{Float64}
            vmax::Vector{Float64}
            length::Int
        end;
        Base.ndims(dom::$name) = $ndim;
    ))

    # getindex
    if ndim == 1
        push!(out.args, :(Base.getindex(dom::$name, dim::Int) = dom.axis))
    else
        push!(out.args, :(Base.getindex(dom::$name, dim::Int) = dom.axis[dim,:]))
    end

    # Constructors
    a = ["d$(i)::AbstractVector{Float64}" for i in 1:ndim]
    b = ["deepcopy(d$(i))" for i in 1:ndim]
    c = ["d$(i)" for i in 1:ndim]
    s = "function Domain(" * join(a, ", ") * ")\n"
    if ndim > 1
        s *= "  @assert " * join("length(" .* c .* ")", " == ") * " \"Arrays must have same length\" \n"
    end
    s *= "  return Domain_$(ndim)D(" *
        (ndim == 1  ?  "d1"  :  "[" * join(b, " ") * "]'") * ", " *
        "[" * join("minimum(" .* c .* ")", ", ") * "], " *
        "[" * join("maximum(" .* c .* ")", ", ") * "], " *
        "length(d1))\n"
    s *= "end"
    push!(out.args, Meta.parse(s))

    a = ["n$(i)::Int" for i in 1:ndim]
    b = ["one(Float64):one(Float64):n$(i)" for i in 1:ndim]
    s = "Domain(" * join(a, ", ") * ") = Domain(" * join(b, ", ") * ")"
    push!(out.args, Meta.parse(s))

    # flatten
    if ndim > 1
        s = "function flatten(dom::CartesianDomain_$(ndim)D)::Domain_$(ndim)D\n"
        s *= "  out = Matrix{Float64}(undef, $ndim, length(dom))\n"
        s *= "  for i in 1:length(dom)\n"
        a = ["d$(i)" for i in 1:ndim]
        s *= "  (" * join(a, ", ") * ") = Tuple(CartesianIndices(size(dom))[dom.index[i]])\n"
        for i in 1:ndim
            s *= "  out[$i, i] = dom.axis[$i][d$(i)]\n"
        end
        s *= "  end\n"
        a = ["out[$(i), :]" for i in 1:ndim]
        s *= "  return Domain(" * join(a, ", ") * ")\n"
        s *= "end\n"
        push!(out.args, Meta.parse(s))
    end

    name = Symbol(:Measures_, ndim, :D)
    push!(out.args, :(
        struct $name <: AbstractMeasures
            val::Array{Float64, $(ndim)}
            unc::Array{Float64, $(ndim)}
        end;
        Base.ndims(dom::$name) = $ndim;
    ))

    push!(out.args, :(
        function Measures(measure::Array{Float64, $ndim}, uncert::Array{Float64, $(ndim)})
            @assert length(measure) == length(uncert) "Measure and uncertainty arrays must have same size"
            $(name)(measure, uncert)
        end;
    ))

    push!(out.args, :(
        function Measures(measure::Array{Float64, $ndim}, uncert::Float64=one(Float64))
            u = similar(measure)
            fill!(u, uncert)
            $(name)(measure, u)
        end;
    ))

    name = Symbol(:Counts_, ndim, :D)
    push!(out.args, :(
        struct $name <: AbstractCounts
            val::Array{Int, $(ndim)}
        end;
        Base.ndims(dom::$name) = $ndim;
    ))

    push!(out.args, :(
        function Counts(measure::Array{Int, $ndim})
            $(name)(measure)
        end
    ))

    return esc(out)
end

@define_ndim 1
@define_ndim 2
@define_ndim 3

getaxismin(dom::AbstractLinearDomain, dim::Int) = dom.vmin[dim]
getaxismax(dom::AbstractLinearDomain, dim::Int) = dom.vmax[dim]
getaxisextrema(dom::AbstractLinearDomain, dim::Int) = (dom.vmin[dim], dom.vmax[dim])
length(dom::AbstractLinearDomain) = dom.length
length(dom::AbstractCartesianDomain) = length(dom.index)
length(data::AbstractData) = length(data.val)
size(data::AbstractData) = size(data.val)
size(data::AbstractData, dim::Int) = size(data.val)[dim]

# Methods to "flatten" a multidimensional object <: AbstractData into a 1D one
flatten(dom::AbstractLinearDomain) = dom
flatten(data::Measures_1D, dom::AbstractLinearDomain)::Measures_1D = data
flatten(data::Counts_1D, dom::AbstractLinearDomain)::Counts_1D = data
flatten(data::AbstractMeasures, dom::AbstractLinearDomain)::Measures_1D    = Measures_1D(data.val[:], data.unc[:])
flatten(data::AbstractMeasures, dom::AbstractCartesianDomain)::Measures_1D = Measures_1D(data.val[dom.index], data.unc[dom.index])
flatten(data::AbstractCounts  , dom::AbstractLinearDomain)::Counts_1D      = Counts_1D(data.val)
flatten(data::AbstractCounts  , dom::AbstractCartesianDomain)::Counts_1D   = Counts_1D(data.val[dom.index])

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
