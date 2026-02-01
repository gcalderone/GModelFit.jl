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
"""
    IdealInstrument

An instrument response representing an ideal instrument, for which unfolded and folded models are identical.

This structure has no fields.
"""
struct IdealInstrument <: AbstractInstrumentResponse
end

unfolded_domain(IR::IdealInstrument, folded_domain::AbstractDomain) = folded_domain
apply_ir!(IR::IdealInstrument,
          folded_domain::AbstractDomain  , folded::Array,
          unfolded_domain::AbstractDomain, unfolded::Array) = folded .= unfolded
