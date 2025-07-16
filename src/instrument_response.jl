# ====================================================================
# Instrument Response
#
abstract type AbstractInstrumentResponse end


struct IREval{T <: AbstractInstrumentResponse, TDomain <: AbstractDomain}
    IR::T
    unfolded_domain::AbstractDomain
    folded_domain::TDomain
    folded::Vector{Float64}
    folded_ad::Vector{Union{Dual, Float64}}

    function IREval(IR::T, folded_domain::TDomain) where {T <: AbstractInstrumentResponse, TDomain <: AbstractDomain}
        d = model_domain(IR, folded_domain)
        return new{T, TDomain}(IR, d, folded_domain,
                               Vector{Float64}(undef, length(folded_domain)),
                               Vector{Union{Dual, Float64}}(undef, length(folded_domain)))
    end
end

apply_ir!(ireval::IREval, unfolded::Vector{Float64}) =
    apply_ir!(ireval.IR, ireval.folded_domain, ireval.folded   , ireval.unfolded_domain, unfolded)

apply_ir!(ireval::IREval, unfolded::Vector) =
    apply_ir!(ireval.IR, ireval.folded_domain, ireval.folded_ad, ireval.unfolded_domain, unfolded)


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

model_domain(IR::IdealInstrument, folded_domain::AbstractDomain) = folded_domain
apply_ir!(IR::IdealInstrument,
          folded_domain::AbstractDomain  , folded::Vector,
          unfolded_domain::AbstractDomain, unfolded::Vector) = folded .= unfolded
