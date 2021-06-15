
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
    comps::OrderedDict{Symbol, BestFitComp}
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:OK, :Warn, :Error
    log10testprob::Float64
    timestamp::DateTime
    elapsed::Float64
end


Base.getindex(res::BestFitResult, cname::Symbol) = res.comps[cname]
Base.keys(res::BestFitResult) = collect(keys(res.preds))
Base.haskey(res::BestFitResult, name::Symbol) = haskey(res.comps, name)


struct BestFitMultiResult
    models::Vector{OrderedDict{Symbol, BestFitComp}}
    ndata::Int
    dof::Int
    cost::Float64
    status::Symbol      #:OK, :Warn, :Error
    log10testprob::Float64
    timestamp::DateTime
    elapsed::Float64
end

Base.getindex(res::BestFitMultiResult, id::Int) = res.models[id]
Base.length(res::BestFitMultiResult) = length(res.models)
