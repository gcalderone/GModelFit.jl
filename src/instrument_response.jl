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

apply_ir!(ireval::IREval, unfolded::Vector{Float64}) =
    apply_ir!(ireval.IR, ireval.data_domain, ireval.folded   , ireval.model_domain, unfolded)

apply_ir!(ireval::IREval, unfolded::Vector) =
    apply_ir!(ireval.IR, ireval.data_domain, ireval.folded_ad, ireval.model_domain, unfolded)


function fold_model(ireval::IREval, unfolded::Vector{Float64})
    folded = deepcopy(ireval.folded)
    apply_ir!(ireval.IR, ireval.data_domain, folded, ireval.model_domain, unfolded)
    return folded
end

last_eval_folded(ireval::IREval) = ireval.folded


# ====================================================================
struct IdealInstrument <: AbstractInstrumentResponse
end

model_domain(IR::IdealInstrument, data_domain::AbstractDomain) = data_domain
apply_ir!(IR::IdealInstrument,
          data_domain::AbstractDomain, folded::Vector,
          model_domain::AbstractDomain, unfolded::Vector) = folded .= unfolded
