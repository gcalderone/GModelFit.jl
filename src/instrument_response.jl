import ForwardDiff: Dual

# ====================================================================
# Instrument Response
#
abstract type AbstractInstrumentResponse end


struct IREval{T <: AbstractInstrumentResponse, TDomain <: AbstractDomain}
    IR::T
    data_domain::TDomain
    model_domain::AbstractDomain
    folded::Vector{Float64}
    folded_ad::Vector{Union{Dual, Float64}}

    function IREval(IR::T, data_domain::TDomain) where {T <: AbstractInstrumentResponse, TDomain <: AbstractDomain}
        mdomain = model_domain(IR, data_domain)
        return new{T, TDomain}(IR, data_domain, mdomain, 
                               Vector{Float64}(undef, length(data_domain)),
                               Vector{Union{Dual, Float64}}(undef, length(data_domain)))
    end
end


# ====================================================================
struct IdealInstrument <: AbstractInstrumentResponse
end

model_domain(IR::IdealInstrument, data_domain::AbstractDomain) = data_domain
apply_ir!(IR::IdealInstrument,
          data_domain::AbstractDomain, folded::Vector,
          model_domain::AbstractDomain, unfolded::Vector) = folded .= unfolded
