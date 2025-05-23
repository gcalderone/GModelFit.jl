import ForwardDiff: Dual

mutable struct CompEvalT{T <: Real}
    counter::Int
    lastparvalues::Vector{Union{T, Float64}}
    lastdepscounter::Vector{Int}
    deps::Vector{CompEvalT}
    buffer::Vector{Union{T, Float64}}

    CompEvalT{T}(npar::Int, nres::Int) where T <: Real =
        new{T}(0,
               Vector{Union{T, Float64}}( undef, npar),
               Vector{Int}(),
               Vector{CompEvalT}(),
               Vector{Union{T, Float64}}(undef, nres))

    # The following is used only for domain coordinates in _update_eval()
    CompEvalT{T}(buffer::Vector{T}) where T <: Real =
        new{T}(1,
               Vector{Union{T, Float64}}(),
               Vector{Int}(),
               Vector{CompEvalT}(),
               buffer)
end

mutable struct CompEval{TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    tpar::CompEvalT{Float64}
    tparad::CompEvalT{Dual}

    function CompEval(comp::AbstractComponent, domain::AbstractDomain)
        prepare!(comp, domain)
        npar = length(getparams(comp))
        nres = result_length(comp, domain)
        return new{typeof(comp), typeof(domain)}(
            comp, domain,
            CompEvalT{Float64}(npar, nres),
            CompEvalT{Dual}(   npar, nres))
    end
end


"""
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, param1, param2....
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, deps::AbstractVector, param1, param2....

Evaluate component `comp` on the given `domain` using `deps` dependencies and `param1`, `param2, ... parameters.  Output should be stored in the `output` vector.

If the component has no dependencies the `deps` argument should not be present.

The `evaluate!` function is called with the `output`, `deps` and parameter arguments containing either `Float64` values (to evaluate the component) or `ForwardDiff.Dual` values (to evaluate the component derivatives).
"""
evaluate!(::AbstractComponent, ::AbstractDomain, args...) =
    error("No evaluate method implemented for $(TComp), $(TDomain)")


update_eval!(ceval::CompEval, pvalues::AbstractVector{Float64}) = _update_eval!(ceval, ceval.tpar  , pvalues)
update_eval!(ceval::CompEval, pvalues::AbstractVector)          = _update_eval!(ceval, ceval.tparad, pvalues)

function _update_eval!(ceval::CompEval, tpar::CompEvalT, pvalues::AbstractVector)
    doeval = (tpar.counter == 0)                                  ||
        any(tpar.lastparvalues   .!= pvalues)                     ||
        length(tpar.lastdepscounter) != length(tpar.deps)         ||
        any(tpar.lastdepscounter .!= getfield.(tpar.deps, :counter))
    if doeval
        if length(tpar.deps) > 0
            evaluate!(ceval.comp, ceval.domain, tpar.buffer, getfield.(tpar.deps, :buffer), pvalues...)
            if length(  tpar.lastdepscounter) != length(tpar.deps)
                empty!( tpar.lastdepscounter)
                append!(tpar.lastdepscounter, getfield.(tpar.deps, :counter))
            else
               tpar.lastdepscounter .= getfield.(tpar.deps, :counter)
            end
        else
            evaluate!(ceval.comp, ceval.domain, tpar.buffer, pvalues...)
        end
        tpar.lastparvalues .= pvalues
        tpar.counter += 1
    end
    return tpar.buffer
end


# Evaluate component on the given domain.  Parameter values are the
# ones stored in the component unless a custom value is provided via a
# keyword.
function (comp::AbstractComponent)(domain::AbstractDomain; kws...)
    @assert length(dependencies(comp)) == 0 "Can't evaluate a component with dependencies as a stand-alone one."
    ceval = CompEval(comp, domain)
    pvalues = OrderedDict([(pname, par.val) for (pname, par) in getparams(comp)])
    for (pname, pval) in kws
        if pname in keys(pvalues)
            pvalues[pname] = pval
        else
            @warn "$pname is not a parameter name for $(typeof(comp)). Valid names are: " * join(string.(keys(pvalues)), ", ")
        end
    end
    return update_eval!(ceval, collect(values(pvalues)))
end


# ====================================================================
# Built-in components
include("components/FComp.jl")
include("components/OffsetSlope.jl")
include("components/Polynomial.jl")
include("components/Gaussian.jl")
include("components/Lorentzian.jl")
include("components/SumReducer.jl")


# ====================================================================
struct ModelEvalT{T <: Real}
    pvalues::PVModel{Union{T, Float64}}
    pactual::PVModel{Union{T, Float64}}
    pvmulti::Vector{PVModel{Union{T, Float64}}}

    ModelEvalT{T}() where T  <: Real =
        new(PVModel{Union{T, Float64}}(),
            PVModel{Union{T, Float64}}(),
            Vector{PVModel{Union{T, Float64}}}())
end

function empty!(v::ModelEvalT)
    empty!(v.pvalues)
    empty!(v.pactual)
    empty!(v.pvmulti)
end


struct ModelEval
    model::Model
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    ifree::Vector{Int}
    patched::Vector{NTuple{2, Symbol}}
    tpar::ModelEvalT{Float64}
    tparad::ModelEvalT{Dual}
    seq::Vector{Symbol}

    function ModelEval(model::Model, domain::AbstractDomain)
        meval = new(model, domain, OrderedDict{Symbol, CompEval}(),
                    Vector{Int}(), Vector{NTuple{2, Symbol}}(),
                    ModelEvalT{Float64}(), ModelEvalT{Dual}(),
                    Vector{CompEval}())
        scan_model!(meval, evaluate=false)
        return meval
    end
end


function scan_model!(meval::ModelEval; evaluate=true)
    function isParamFixed(par::Parameter)
        if !isnothing(par.patch)
            @assert isnothing(par.mpatch) "Parameter [$cname].$pname has both patch and mpatch fields set, while only one is allowed"
            if isa(par.patch, Symbol)  # use same param. value from a different component
                return true
            else                       # invoke a patch function
                @assert length(par.patch.args) in [1,2]
                return  length(par.patch.args) == 1
            end
        elseif !isnothing(par.mpatch)
            @assert length(par.mpatch.args) in [1,2]
            return  length(par.mpatch.args) == 1
        end
        return par.fixed
    end

    empty!(meval.ifree)
    empty!(meval.patched)
    empty!(meval.tpar)
    empty!(meval.tparad)

    isfixed = Vector{Bool}()
    for (cname, comp) in meval.model.comps
        for (pname, par) in getparams(comp)
            # Some solvers do not handle parameter limits, ensure values are in the allowed range
            if par.val < par.low
                s = string(par) * "[$(cname)].$(pname) value outside allowed range, using lower limit"
                @warn s
                par.val = par.low
            elseif par.val > par.high
                s = string(par) * "[$(cname)].$(pname) value outside allowed range, using upper limit"
                par.val = par.high
            end
            if isnan(par.low)  ||  isnan(par.high)  ||  isnan(par.val)
                s = "NaN value detected for param [$(cname)].$(pname):\n" * string(par)
                error(s)
            end

            push!(isfixed, isParamFixed(par)  ||  meval.model.fixed[cname])
            if !isnothing(par.patch)  ||  !isnothing(par.mpatch)
                push!(meval.patched, (cname, pname))
            end

            push!(meval.tpar.pvalues  , cname, pname, par.val)
            push!(meval.tpar.pactual  , cname, pname, par.val)
            push!(meval.tparad.pvalues, cname, pname, par.val)
            push!(meval.tparad.pactual, cname, pname, par.val)
        end
        if !(cname in keys(meval.cevals))
            ceval = CompEval(comp, meval.domain)
            meval.cevals[cname] = ceval
        end
    end
    append!(meval.ifree, findall(.! isfixed))

    for (cname, ceval) in meval.cevals
        if length(ceval.tpar.deps) == 0
            i = 1
            for d in dependencies(meval.model, cname, select_domain=true)
                push!(ceval.tpar.deps  , CompEvalT{Float64}(coords(meval.domain, i)))
                push!(ceval.tparad.deps, CompEvalT{Float64}(coords(meval.domain, i)))
                i += 1
            end
            for d in dependencies(meval.model, cname, select_domain=false)
                push!(ceval.tpar.deps  , meval.cevals[d].tpar)
                push!(ceval.tparad.deps, meval.cevals[d].tparad)
            end
        else
            d        = dependencies(meval.model, cname, select_domain=true)
            append!(d, dependencies(meval.model, cname, select_domain=false))
            deleteat!(d, 1:length(ceval.tpar.deps))
            for d in d
                push!(ceval.tpar.deps  , meval.cevals[d].tpar)
                push!(ceval.tparad.deps, meval.cevals[d].tparad)
            end
        end
    end

    function compeval_sequence!(meval::ModelEval, cname::Symbol)
        for d in dependencies(meval.model, cname)
            compeval_sequence!(meval, d)
        end
        push!(meval.seq, cname)
    end
    function compeval_sequence!(meval::ModelEval)
        empty!(meval.seq)
        compeval_sequence!(meval, find_maincomp(meval.model))
    end
    compeval_sequence!(meval)

    # Avoid early evaluation since ModelEval ina multi-model situation may not yet been created
    evaluate  &&  update_eval!(meval)
    nothing
end


function free_params(meval::ModelEval)
    out = Vector{Parameter}()
    for (cname, comp) in meval.model.comps
        for (pname, par) in getparams(comp)
            push!(out, par)
        end
    end
    return out[meval.ifree]
end
free_params_val(meval::ModelEval) = convert(Vector{Float64}, getfield.(free_params(meval), :val))  # convert() is useful when there are no free parameters
nfree(meval::ModelEval) = length(meval.ifree)


# Set new model parameters
function set_pvalues!(meval::ModelEval, pvalues::AbstractVector{Float64})
    items(meval.tpar.pvalues)[meval.ifree] .= pvalues
    items(meval.tpar.pactual)[meval.ifree] .= pvalues
end

function set_pvalues!(meval::ModelEval, pvalues::AbstractVector)
    items(meval.tparad.pvalues)[meval.ifree] .= pvalues
    items(meval.tparad.pactual)[meval.ifree] .= pvalues
end


function run_patch_functs!(meval::ModelEval, tpar::ModelEvalT)
    for (cname, pname) in meval.patched
        par = getfield(meval.model[cname], pname)
        if !isnothing(par.patch)
            @assert isnothing(par.mpatch) "Parameter [:$(cname)].$pname has both patch and mpatch fields set, while only one is allowed"
            if isa(par.patch, Symbol)  # use same param. value from a different component
                tpar.pactual[cname][pname] = tpar.pvalues[par.patch][pname]
            else                       # invoke a patch function
                if length(par.patch.args) == 1
                    tpar.pactual[cname][pname] = par.patch(tpar.pvalues)
                else
                    tpar.pactual[cname][pname] = par.patch(tpar.pvalues, tpar.pvalues[cname][pname])
                end
            end
        elseif !isnothing(par.mpatch)
            @assert length(tpar.pvmulti) > 0 "Parameter [:$(cname)].$pname has the mpatch field set but no other Model is being considered"
            if length(par.mpatch.args) == 1
                tpar.pactual[cname][pname] = par.mpatch(tpar.pvmulti)
            else
                tpar.pactual[cname][pname] = par.mpatch(tpar.pvmulti, tpar.pvalues[cname][pname])
            end
        end
    end
end


update_eval!(meval::ModelEval) = update_eval!(meval, items(meval.tpar.pvalues)[meval.ifree])
function update_eval!(meval::ModelEval, pvalues::AbstractVector{Float64})
    set_pvalues!(meval, pvalues)
    run_patch_functs!(meval, meval.tpar)
    for cname in meval.seq
        update_eval!(meval.cevals[cname], items(meval.tpar.pactual[cname]))
    end
    return meval.cevals[meval.seq[end]].tpar.buffer
end

function update_eval!(meval::ModelEval, pvalues::AbstractVector)
    set_pvalues!(meval, pvalues)
    run_patch_functs!(meval, meval.tparad)
    for cname in meval.seq
        update_eval!(meval.cevals[cname], items(meval.tparad.pactual[cname]))
    end
    return meval.cevals[meval.seq[end]].tparad.buffer
end


#=
evalcounter(meval::ModelEval, cname::Symbol)

Return the number of times the component with name `cname` has been evaluated.
=#
evalcounter(meval::ModelEval, cname::Symbol) = meval.cevals[cname].tpar.counter + meval.cevals[cname].tparad.counter
evalcounter(model::Model, cname::Symbol) = "???"


#=
evalcounters(meval::ModelEval)

Return a `OrderedDict{Symbol, Int}` with the number of times each component has been evaluated.
=#
evalcounters(meval::ModelEval) = OrderedDict([cname => evalcounter(meval, cname) for cname in keys(meval.cevals)])


#=
last_eval(meval::ModelEval)
last_eval(meval::ModelEval, name::Symbol)

Return last evaluation of a component whose name is `cname` in a `ModelEval` object.  If `cname` is not provided the evaluation of the main component is returned.
=#
last_eval(meval::ModelEval) = last_eval(meval, meval.seq[end])
function last_eval(meval::ModelEval, cname::Symbol)
    if meval.cevals[cname].tpar.counter == 0  # Ensure model is evaluated
        update_eval!(meval)
    end
    meval.cevals[cname].tpar.buffer
end

# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    meval = ModelEval(model, domain)
    update_eval!(meval)
    isnothing(cname)  &&  (cname = meval.seq[end])
    return last_eval(meval, cname)
end
