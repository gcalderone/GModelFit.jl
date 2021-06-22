
# ====================================================================
# Fit results
#
struct BestFitParam
    val::Float64
    unc::Float64
    fixed::Bool
    patched::Float64  # value after transformation
end

const BestFitComp = HashVector{BestFitParam}


struct BestFitResult
    timestamp::DateTime
    elapsed::Float64
    mzer::AbstractMinimizerStatus
    comps::OrderedDict{Symbol, BestFitComp}
    mdc::MDComparison
end


Base.getindex(res::BestFitResult, cname::Symbol) = res.comps[cname]
Base.keys(res::BestFitResult) = collect(keys(res.preds))
Base.haskey(res::BestFitResult, name::Symbol) = haskey(res.comps, name)


struct BestFitMultiResult
    timestamp::DateTime
    elapsed::Float64
    mzer::AbstractMinimizerStatus
    models::Vector{OrderedDict{Symbol, BestFitComp}}
    mdc::MDMultiComparison
end

Base.getindex(res::BestFitMultiResult, id::Int) = res.models[id]
Base.length(res::BestFitMultiResult) = length(res.models)
