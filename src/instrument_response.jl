# ====================================================================
# Instrument Response
#
abstract type AbstractInstrumentResponse end

"""
    prepare!(IR::AbstractInstrumentResponse, folded_domain::AbstractDomain)

Invoked to precompute instrument response specific quantities.

This method is invoked when the instrument response is first evaluated hence it is the perfect place to pre-compute relevant quantities associated to the evaluation on a specific folded domain.

Extension of this method for user defined instrument response is optional, and the default implementation `nothing`.
"""
prepare!(IR::AbstractInstrumentResponse, folded_domain::AbstractDomain) = nothing


# ====================================================================
struct IREval{T <: AbstractInstrumentResponse, TDomain <: AbstractDomain}
    IR::T
    unfolded_domain::AbstractDomain
    folded_domain::TDomain
    folded::Vector{Float64}
    folded_ad::Vector{Union{Dual, Float64}}

    function IREval(IR::T, folded_domain::TDomain) where {T <: AbstractInstrumentResponse, TDomain <: AbstractDomain}
        prepare!(IR, folded_domain)
        d = unfolded_domain(IR, folded_domain)
        return new{T, TDomain}(IR, d, folded_domain,
                               Vector{Float64}(undef, length(folded_domain)),
                               Vector{Union{Dual, Float64}}(undef, length(folded_domain)))
    end
end

apply_ir!(ireval::IREval, unfolded::Vector{Float64}) =
    apply_ir!(ireval.IR, ireval.folded_domain, ireval.folded   , ireval.unfolded_domain, unfolded)

apply_ir!(ireval::IREval, unfolded::Vector) =
    apply_ir!(ireval.IR, ireval.folded_domain, ireval.folded_ad, ireval.unfolded_domain, unfolded)


unfolded_domain(ireval::IREval) = ireval.unfolded_domain
folded_domain(ireval::IREval) = ireval.folded_domain

function fold_model(ireval::IREval, unfolded::Vector{Float64})
    folded = deepcopy(ireval.folded)
    apply_ir!(ireval.IR, ireval.folded_domain, folded, ireval.unfolded_domain, unfolded)
    return folded
end

last_eval_folded(ireval::IREval) = ireval.folded


# ====================================================================
"""
    IdealInstrument

An instrument response representing an ideal instrument, for which unfolded and folded models are identical.

This structure has no fields.
"""
struct IdealInstrument <: AbstractInstrumentResponse
end

unfolded_domain(IR::IdealInstrument, folded_domain::AbstractDomain) = folded_domain
apply_ir!(IR::IdealInstrument,
          folded_domain::AbstractDomain  , folded::Vector,
          unfolded_domain::AbstractDomain, unfolded::Vector) = folded .= unfolded
