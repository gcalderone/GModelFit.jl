# ====================================================================
"""
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, param1, param2....
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, deps::AbstractVector, param1, param2....

Evaluate component `comp` on the given `domain` using `deps` dependencies and `param1`, `param2, ... parameters.  Output should be stored in the `output` vector.

If the component has no dependencies the `deps` argument should not be present.

The `evaluate!` function is called with the `output`, `deps` and parameter arguments containing either `Float64` values (to evaluate the component) or `ForwardDiff.Dual` values (to evaluate the component derivatives).
"""
evaluate!(::TComp, ::TDomain, args...) where {TComp <: AbstractComponent, TDomain <: AbstractDomain} =
    error("No evaluate method implemented for $(TComp), $(TDomain)")


domain2buffer(T, domain::Domain) = Array{T, 1}(undef, length(domain))
domain2buffer(T, domain::CartesianDomain{N}) where N = Array{T, N}(undef, size(domain))


mutable struct CompEval{T <: Real, TComp <: AbstractComponent, TDomain <: AbstractDomain}
    comp::TComp
    domain::TDomain
    counter::Int
    lastparvalues::Vector{T}
    lastdepscounter::Vector{Int}
    deps::Vector{CompEval}
    buffer::Array{T}

    function CompEval{T}(comp::AbstractComponent, domain::AbstractDomain) where {T <: Real}
        prepare!(comp, domain)
        return new{T, typeof(comp), typeof(domain)}(
            comp, domain,
            0,
            Vector{T}( undef, length(getparams(comp))),
            Vector{Int}(),
            Vector{CompEval}(),
            domain2buffer(T, domain))
    end
end


function update_eval!(ceval::CompEval{T, TComp, TDomain}, pvalues::AbstractVector{T}) where {TComp, TDomain, T}
    doeval = (ceval.counter == 0)                                   ||
        any(ceval.lastparvalues   .!= pvalues)                      ||
        length(ceval.lastdepscounter) != length(ceval.deps)         ||
        any(ceval.lastdepscounter .!= getfield.(ceval.deps, :counter))
    if doeval
        if length(ceval.deps) > 0
            evaluate!(ceval.comp, ceval.domain, ceval.buffer, getfield.(ceval.deps, :buffer), pvalues...)
            if length(  ceval.lastdepscounter) != length(ceval.deps)
                empty!( ceval.lastdepscounter)
                append!(ceval.lastdepscounter, getfield.(ceval.deps, :counter))
            else
                ceval.lastdepscounter .= getfield.(ceval.deps, :counter)
            end
        else
            evaluate!(ceval.comp, ceval.domain, ceval.buffer, pvalues...)
        end
        ceval.lastparvalues .= pvalues
        ceval.counter += 1
    end
    return ceval.buffer
end


# Evaluate component on the given domain.  Parameter values are the
# ones stored in the component unless a custom value is provided via a
# keyword.
function (comp::AbstractComponent)(domain::AbstractDomain; kws...)
    @assert length(dependencies(comp)) == 0 "Can't evaluate a component with dependencies as a stand-alone one."
    ceval = CompEval{Float64}(comp, domain)
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
struct ModelEval{T <: Real}
    model::Model
    pvalues::PVModel{T}
    pactual::PVModel{T}
    patched::OrderedDict{NTuple{2, Symbol}, Parameter}
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    seq::Vector{Symbol}
    folded_domain::AbstractDomain
    folded::Array{T}

    function ModelEval{T}(model::Model, folded_domain::AbstractDomain,
                          pvalues::PVModel{T}, pactual::PVModel{T}) where {T <: Real}
        prepare!(model.IR, folded_domain)
        meval = new{T}(model, pvalues, pactual,
                       OrderedDict{NTuple{2, Symbol}, Parameter}(),
                       unfolded_domain(model.IR, folded_domain),
                       OrderedDict{Symbol, CompEval}(),
                       Vector{Symbol}(),
                       folded_domain,
                       domain2buffer(T, folded_domain))

        for (key, par) in getparams(model)
            if !isnothing(par.patch)  ||  !isnothing(par.reparam)
                meval.patched[(key[1], key[2])] = par
            end
        end

        # Update evaluation sequence
        ftree = flatten(deptree(model))
        append!(meval.seq, reverse(getfield.(ftree, :cname)))

        # Walk dependency tree and create relevant CompEval structures (if
        # not yet existing)
        for d in ftree
            cname = d.cname
            if !(cname in keys(meval.cevals))
                meval.cevals[cname] = CompEval{T}(model.comps[cname], meval.domain)  # all components share the same domain
            end
        end

        # Update references to dependencies
        for (cname, ceval) in meval.cevals
            deps = dependencies(ceval.comp)
            if isa(ceval.comp, FComp)
                _domdeps = Vector{Symbol}()
                _compdeps = Vector{Symbol}()
                for d in deps
                    if haskey(meval.cevals, d)  # dependency with known name
                        push!(_compdeps, d)
                    else # dependency with unknown name is intended as a domain dimension
                        @assert length(_compdeps) == 0 "Domain dependencies for component $cname must be listed before other component name(s)"
                        push!(_domdeps, d)
                    end
                end
                @assert length(_domdeps) == ndims(meval.domain) "Domain has $(ndims(meval.domain)) dimensions but component $cname accepts $(length(_domdeps)) domain arguments: " * join(string.(_domdeps), ", ")
                deps = _compdeps
            end
            for d in deps
                @assert haskey(meval.cevals, d) "No component has name $d"
                push!(ceval.deps, meval.cevals[d])
            end
        end
        return meval
    end
end


struct ModelSetEval{T <: Real}
    ms::ModelSet
    dict::OrderedDict{Symbol, ModelEval{T}}
    vec::Vector{ModelEval{T}}
    params::OrderedDict{NTuple{3, Symbol}, Parameter}
    pvalues::PVSet{T}
    pactual::PVSet{T}
    ifree::Vector{Int}

    function ModelSetEval{T}(ms::ModelSet, domains::Vector{<: AbstractDomain}) where {T <: Real}
        N = length(ms)
        @assert N == length(domains)
        out = new{T}(ms, OrderedDict{Symbol, ModelEval{T}}(),
                     Vector{ModelEval{T}}(),
                     OrderedDict{NTuple{3, Symbol}, Parameter}(),
                     PVSet{T}(),
                     PVSet{T}(allow_overwrite=true),
                     Vector{Int}())
        for (key, par) in getparams(ms)
            out.params[ key]    = deepcopy(par) # local copy, these will contain best fit nd uncertainties
            out.pvalues[key...] = par.val
            out.pactual[key...] = par.val
            if !par.actually_fixed
                push!(out.ifree, length(out.params))
            end
        end
        for (mname, model) in ms
            out.dict[mname] = ModelEval{T}(model, domains[length(out.vec) + 1],
                                        PVModel(mname, out.pvalues),
                                        PVModel(mname, out.pactual))
            push!(out.vec, out.dict[mname])
        end
        return out
    end
end


function update_eval!(meval::ModelEval)
    for (key, par) in meval.patched
        if !isnothing(par.patch)
            if isa(par.patch, Symbol)  # use same param. value from a different component
                meval.pactual[key...] = meval.pvalues[par.patch, key[2]]
            else                       # invoke a patch function
                meval.pactual[key...] = par.patch(meval.pvalues)
            end
        else
            meval.pactual[key...] = par.reparam(meval.pvalues, meval.pvalues[key...])
        end
    end
    for cname in meval.seq
        update_eval!(meval.cevals[cname], meval.pactual[cname])
    end
    unfolded = meval.cevals[meval.seq[end]].buffer
    apply_ir!(meval.model.IR, meval.folded_domain, meval.folded, meval.domain, unfolded)
end

function update_eval!(mseval::ModelSetEval, pvalues::Vector)
    i = mseval.ifree
    mseval.pvalues.vec[i] .= pvalues
    mseval.pactual.vec[i] .= pvalues
    update_eval!.(mseval.vec)
end

function update_eval!(mseval::ModelSetEval)
    mseval.pactual.vec .= mseval.pvalues.vec
    update_eval!.(mseval.vec)
end

getparams(mseval::ModelSetEval) = mseval.params
nfree(mseval::ModelSetEval) = length(mseval.ifree)

function set_bestfit!(mseval::ModelSetEval, pvalues::Vector{Float64}, puncerts::Vector{Float64})
    update_eval!(mseval, pvalues)
    i = 0
    for (key, par) in mseval.params
        par.val    = mseval.pvalues[key...]
        par.actual = mseval.pactual[key...]
        if !par.actually_fixed
            i += 1
            @assert par.val == pvalues[i]
            par.unc = puncerts[i]
        end
    end
end


# Return the number of times the component with name `cname` has been evaluated.
evalcounter(meval::ModelEval, cname::Symbol) = meval.cevals[cname].counter


# Return last evaluation of a component whose name is `cname` in a `ModelEval` object.  If `cname` is not provided the evaluation of the main component is returned.
last_eval(meval::ModelEval) = last_eval(meval, meval.seq[end])
last_eval(meval::ModelEval, cname::Symbol) = meval.cevals[cname].buffer
last_eval_folded(meval::ModelEval) = meval.folded
function fold_model(meval::ModelEval{T}, cname::Symbol) where T
    output = domain2buffer(T, meval.folded_domain)
    apply_ir!(meval.model.IR, meval.folded_domain, output, meval.domain, last_eval(meval, cname))
    return output
end

last_eval(mseval::ModelSetEval, mname::Symbol) = last_eval(mseval.dict[mname])
last_eval(mseval::ModelSetEval, mname::Symbol, cname::Symbol) = last_eval(mseval.dict[mname], cname)

#last_eval_folded(multi::ModelSetEval{1}) = last_eval_folded(multi, 1)
last_eval_folded(mseval::ModelSetEval, mname::Symbol) = last_eval_folded(mseval.dict[mname])

#fold_model(multi::ModelSetEval{1}, cname::Symbol) = fold_model(multi, 1, cname)
fold_model(mseval::ModelSetEval, mname::Symbol, cname::Symbol) = fold_model(mseval.dict[mname], cname)


# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    mseval = ModelSetEval{Float64}(ModelSet(:_ => model), [domain])
    update_eval!(mseval)
    isnothing(cname)  &&  (return last_eval_folded(mseval, :_))
    return fold_model(mseval, :_, cname)
end
