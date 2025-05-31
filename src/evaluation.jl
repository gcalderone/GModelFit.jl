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
evaluate!(::TComp, ::TDomain, args...) where {TComp <: AbstractComponent, TDomain <: AbstractDomain} =
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
struct SingleMEvalT{T <: Real}
    pvalues::PVModel{Union{T, Float64}}
    pactual::PVModel{Union{T, Float64}}
    pvmulti::Vector{PVModel{Union{T, Float64}}}

    SingleMEvalT{T}() where T  <: Real =
        new(PVModel{Union{T, Float64}}(),
            PVModel{Union{T, Float64}}(),
            Vector{PVModel{Union{T, Float64}}}())
end

function empty!(v::SingleMEvalT)
    empty!(v.pvalues)
    empty!(v.pactual)
    empty!(v.pvmulti)
end


struct SingleMEval
    model::Model
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    ifree::Vector{Int}
    patched::Vector{NTuple{2, Symbol}}
    tpar::SingleMEvalT{Float64}
    tparad::SingleMEvalT{Dual}
    seq::Vector{Symbol}

    function SingleMEval(model::Model, domain::AbstractDomain)
        meval = new(model, domain, OrderedDict{Symbol, CompEval}(),
                    Vector{Int}(), Vector{NTuple{2, Symbol}}(),
                    SingleMEvalT{Float64}(), SingleMEvalT{Dual}(),
                    Vector{CompEval}())
        return meval
    end
end


function scan_model!(meval::SingleMEval)
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

    function compeval_sequence!(meval::SingleMEval, cname::Symbol)
        for d in dependencies(meval.model, cname)
            compeval_sequence!(meval, d)
        end
        push!(meval.seq, cname)
    end
    function compeval_sequence!(meval::SingleMEval)
        empty!(meval.seq)
        compeval_sequence!(meval, find_maincomp(meval.model))
    end
    compeval_sequence!(meval)
    nothing
end


function free_params(meval::SingleMEval)
    out = Vector{Parameter}()
    for (cname, comp) in meval.model.comps
        for (pname, par) in getparams(comp)
            push!(out, par)
        end
    end
    return out[meval.ifree]
end
nfree(meval::SingleMEval) = length(meval.ifree)


# Set new model parameters
function set_pvalues!(meval::SingleMEval, pvalues::AbstractVector{Float64})
    items(meval.tpar.pvalues)[meval.ifree] .= pvalues
    items(meval.tpar.pactual)[meval.ifree] .= pvalues
end

function set_pvalues!(meval::SingleMEval, pvalues::AbstractVector)
    items(meval.tparad.pvalues)[meval.ifree] .= pvalues
    items(meval.tparad.pactual)[meval.ifree] .= pvalues
end


function run_patch_functs!(meval::SingleMEval, tpar::SingleMEvalT)
    for (cname, pname) in meval.patched
        par = getproperty(meval.model[cname], pname)
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


function update_eval!(meval::SingleMEval)
    run_patch_functs!(meval, meval.tpar)
    for cname in meval.seq
        update_eval!(meval.cevals[cname], items(meval.tpar.pactual[cname]))
    end
    return meval.cevals[meval.seq[end]].tpar.buffer
end

function update_eval_ad!(meval::SingleMEval)
    run_patch_functs!(meval, meval.tparad)
    for cname in meval.seq
        update_eval!(meval.cevals[cname], items(meval.tparad.pactual[cname]))
    end
    return meval.cevals[meval.seq[end]].tparad.buffer
end


#=
evalcounter(meval::SingleMEval, cname::Symbol)

Return the number of times the component with name `cname` has been evaluated.
=#
evalcounter(meval::SingleMEval, cname::Symbol) = meval.cevals[cname].tpar.counter + meval.cevals[cname].tparad.counter
evalcounter(model::Model, cname::Symbol) = "???"


#=
evalcounters(meval::SingleMEval)

Return a `OrderedDict{Symbol, Int}` with the number of times each component has been evaluated.
=#
evalcounters(meval::SingleMEval) = OrderedDict([cname => evalcounter(meval, cname) for cname in keys(meval.cevals)])


#=
last_eval(meval::SingleMEval)
last_eval(meval::SingleMEval, name::Symbol)

Return last evaluation of a component whose name is `cname` in a `SingleMEval` object.  If `cname` is not provided the evaluation of the main component is returned.
=#
last_eval(meval::SingleMEval) = last_eval(meval, meval.seq[end])
last_eval(meval::SingleMEval, cname::Symbol) = meval.cevals[cname].tpar.buffer


# ====================================================================
struct MEval{N}
    v::Vector{SingleMEval}

    MEval(model::Model, domain::AbstractDomain) = MEval([model], [domain])
    function MEval(models::Vector{Model}, domains::Vector{<: AbstractDomain})
        @assert length(models) == length(domains)
        out = new{length(models)}([SingleMEval(models[i], domains[i]) for i in 1:length(models)])
        scan_model!(out)
        return out
    end
end
length(meval::MEval) = length(meval.v)


function free_params_indices(meval::MEval)
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(meval)
        nn = length(meval.v[id].ifree)
        if nn > 0
            i2 = i1 + nn - 1
            push!(out, (id, i1, i2))
            i1 += nn
        end
    end
    return out
end


scan_model!(meval::MEval{1}) = scan_model!(meval.v[1])
function scan_model!(meval::MEval)
    for i in 1:length(meval)
        scan_model!(meval.v[i])
    end
    pv1 = [meval.v[i].tpar.pvalues for i in 1:length(meval)]
    for i in 1:length(meval)
        empty!( meval.v[i].tpar.pvmulti)
        append!(meval.v[i].tpar.pvmulti, pv1)
    end
    pv2 = [meval.v[i].tparad.pvalues for i in 1:length(meval)]
    for i in 1:length(meval)
        empty!( meval.v[i].tparad.pvmulti)
        append!(meval.v[i].tparad.pvmulti, pv2)
    end
end


function free_params(meval::MEval)
    out = Vector{Parameter}()
    for id in 1:length(meval)
        append!(out, free_params(meval.v[id]))
    end
    return out
end
free_params_val(meval::MEval) = getfield.(free_params(meval), :val)
nfree(meval::MEval) = sum(nfree.(meval.v))


function set_pvalues!(meval::MEval{N}, pvalues::AbstractVector) where N
    if N == 1
        set_pvalues!(meval.v[1], pvalues)
    else
        for (id, i1, i2) in free_params_indices(meval)
            set_pvalues!(meval.v[id], pvalues[i1:i2])
        end
    end
end

function update_eval!(meval::MEval, pvalues::AbstractVector{Float64})
    set_pvalues!(meval, pvalues)
    return update_eval!.(meval.v)
end

function update_eval!(meval::MEval{N}, pvalues::AbstractVector) where N
    set_pvalues!(meval, pvalues)
    return update_eval_ad!.(meval.v)
end

update_eval!(meval::MEval) = update_eval!.(meval.v)


last_eval(meval::MEval{1}) = last_eval(meval, 1)
last_eval(meval::MEval{1}, cname::Symbol) = last_eval(meval, 1, cname)
last_eval(meval::MEval, id::Int) = last_eval(meval.v[id])
last_eval(meval::MEval, id::Int, cname::Symbol) = last_eval(meval.v[id], cname)


# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    meval = MEval(model, domain)
    update_eval!(meval)
    isnothing(cname)  &&  (return last_eval(meval))
    return last_eval(meval, cname)
end
