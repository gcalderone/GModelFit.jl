
"""
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, param1, param2....
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, deps::AbstractVector, param1, param2....

Evaluate component `comp` on the given `domain` using `deps` dependencies and `param1`, `param2, ... parameters.  Output should be stored in the `output` vector.

If the component has no dependencies the `deps` argument should not be present.

The `evaluate!` function is called with the `output`, `deps` and parameter arguments containing either `Float64` values (to evaluate the component) or `ForwardDiff.Dual` values (to evaluate the component derivatives).
"""
evaluate!(::TComp, ::TDomain, args...) where {TComp <: AbstractComponent, TDomain <: AbstractDomain} =
    error("No evaluate method implemented for $(TComp), $(TDomain)")


domain2buffer(T::DataType, domain::Domain) = Array{T, 1}(undef, length(domain))
domain2buffer(T::DataType, domain::CartesianDomain{N}) where N = Array{T, N}(undef, size(domain))


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
    domain::AbstractDomain
    cevals::OrderedDict{Symbol, CompEval}
    ifree::Vector{Int}
    patched::Vector{NTuple{2, Symbol}}
    pvalues::PVModel{T}
    pactual::PVModel{T}
    pvmulti::Vector{PVModel{T}}
    seq::Vector{Symbol}
    folded_domain::AbstractDomain
    folded::Array{T}

    function ModelEval{T}(model::Model, folded_domain::AbstractDomain) where {T <: Real}
        prepare!(model.IR, folded_domain)
        meval = new{T}(model, unfolded_domain(model.IR, folded_domain),
                       OrderedDict{Symbol, CompEval}(),
                       Vector{Int}(), Vector{NTuple{2, Symbol}}(),
                       PVModel{T}(),
                       PVModel{T}(),
                       Vector{PVModel{T}}(),
                       Vector{Symbol}(),
                       folded_domain,
                       domain2buffer(T, folded_domain))
        return meval
    end
end


function scan_model!(meval::ModelEval{T}) where {T <: Real}
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
    empty!(meval.pvalues)
    empty!(meval.pactual)
    empty!(meval.pvmulti)

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

            push!(meval.pvalues, cname, pname, par.val)
            push!(meval.pactual, cname, pname, par.val)
        end
    end
    append!(meval.ifree, findall(.! isfixed))

    # Update evaluation sequence
    tree = deptree(meval.model)
    ftree = flatten(tree)
    empty!(meval.seq)
    append!(meval.seq, reverse(getfield.(ftree, :cname)))

    # Walk dependency tree and create relevant CompEval structures (if
    # not yet existing)
    for d in ftree
        cname = d.cname
        if !(cname in keys(meval.cevals))
            ceval = CompEval{T}(meval.model.comps[cname], meval.domain)  # all components share the same domain
            meval.cevals[cname] = ceval
        end
    end

    # Update references to dependencies
    for (cname, ceval) in meval.cevals
        empty!(ceval.deps)
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
            @assert haskey(meval.cevals, d)  "No component has name $d"
            push!(ceval.deps, meval.cevals[d])
        end
    end

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
nfree(meval::ModelEval) = length(meval.ifree)


function run_patch_functs!(meval::ModelEval)
    for (cname, pname) in meval.patched
        par = getproperty(meval.model[cname], pname)
        if !isnothing(par.patch)
            @assert isnothing(par.mpatch) "Parameter [:$(cname)].$pname has both patch and mpatch fields set, while only one is allowed"
            if isa(par.patch, Symbol)  # use same param. value from a different component
                meval.pactual[cname][pname] = meval.pvalues[par.patch][pname]
            else                       # invoke a patch function
                if length(par.patch.args) == 1
                    meval.pactual[cname][pname] = par.patch(meval.pvalues)
                else
                    meval.pactual[cname][pname] = par.patch(meval.pvalues, meval.pvalues[cname][pname])
                end
            end
        elseif !isnothing(par.mpatch)
            @assert length(meval.pvmulti) > 0 "Parameter [:$(cname)].$pname has the mpatch field set but no other Model is being considered"
            if length(par.mpatch.args) == 1
                meval.pactual[cname][pname] = par.mpatch(meval.pvmulti)
            else
                meval.pactual[cname][pname] = par.mpatch(meval.pvmulti, meval.pvalues[cname][pname])
            end
        end
    end
end


#=
evalcounter(meval::ModelEval, cname::Symbol)

Return the number of times the component with name `cname` has been evaluated.
=#
evalcounter(meval::ModelEval, cname::Symbol) = meval.cevals[cname].counter
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
last_eval(meval::ModelEval, cname::Symbol) = meval.cevals[cname].buffer
last_eval_folded(meval::ModelEval) = meval.folded
function fold_model(meval::ModelEval{T}, cname::Symbol) where T
    output = domain2buffer(T, meval.folded_domain)
    apply_ir!(meval.model.IR, meval.folded_domain, output, meval.domain, last_eval(meval, cname))
    return output
end


# ====================================================================
struct MultiEval{N, T <: Real}
    v::Vector{ModelEval{T}}

    MultiEval{T}(model::Model, domain::AbstractDomain) where {T <: Real} = MultiEval{T}([model], [domain])
    function MultiEval{T}(models::Vector{Model}, domains::Vector{<: AbstractDomain}) where {T <: Real}
        @assert length(models) == length(domains)
        out = new{length(models), T}([ModelEval{T}(models[i], domains[i]) for i in 1:length(models)])
        scan_model!(out)
        return out
    end
end
length(multi::MultiEval) = length(multi.v)


function free_params_indices(multi::MultiEval)
    out = Vector{NTuple{3, Int}}()
    i1 = 1
    for id in 1:length(multi)
        nn = length(multi.v[id].ifree)
        #(nn == 0)  && continue   This is a bug (e.g. if in multimodel one model has zero free params
        i2 = i1 + nn - 1
        push!(out, (id, i1, i2))
        i1 += nn
    end
    return out
end


scan_model!(multi::MultiEval{1,T}) where {T <: Real} = scan_model!(multi.v[1])
function scan_model!(multi::MultiEval)
    for i in 1:length(multi)
        scan_model!(multi.v[i])
    end
    pv = [multi.v[i].pvalues for i in 1:length(multi)]
    for i in 1:length(multi)
        empty!( multi.v[i].pvmulti)
        append!(multi.v[i].pvmulti, pv)
    end
end


function free_params(multi::MultiEval)
    out = Vector{Parameter}()
    for id in 1:length(multi)
        append!(out, free_params(multi.v[id]))
    end
    return out
end
nfree(multi::MultiEval) = sum(nfree.(multi.v))


function update_eval!(multi::MultiEval{N, T}, pvalues::Vector) where {N, T  <: Real}
    for (i, i1, i2) in free_params_indices(multi)
        meval = multi.v[i]
        items(meval.pvalues)[meval.ifree] .= pvalues[i1:i2]
        items(meval.pactual)[meval.ifree] .= pvalues[i1:i2]
    end

    for i in 1:N
        meval = multi.v[i]
        run_patch_functs!(meval)
        for cname in meval.seq
            update_eval!(meval.cevals[cname], items(meval.pactual[cname]))
        end
        unfolded = meval.cevals[meval.seq[end]].buffer
        apply_ir!(meval.model.IR, meval.folded_domain, meval.folded, meval.domain, unfolded)
    end
    nothing
end

update_eval!(multi::MultiEval) = update_eval!(multi, getfield.(free_params(multi), :val))


last_eval(multi::MultiEval{1}) = last_eval(multi, 1)
last_eval(multi::MultiEval{1}, cname::Symbol) = last_eval(multi, 1, cname)
last_eval(multi::MultiEval, id::Int) = last_eval(multi.v[id])
last_eval(multi::MultiEval, id::Int, cname::Symbol) = last_eval(multi.v[id], cname)

last_eval_folded(multi::MultiEval{1}) = last_eval_folded(multi, 1)
last_eval_folded(multi::MultiEval, id::Int) = last_eval_folded(multi.v[id])

fold_model(multi::MultiEval{1}, cname::Symbol) = fold_model(multi, 1, cname)
fold_model(multi::MultiEval, id::Int, cname::Symbol) = fold_model(multi.v[id], cname)


# ====================================================================
# Evaluate Model on the given domain
function (model::Model)(domain::AbstractDomain, cname::Union{Nothing, Symbol}=nothing)
    multi = MultiEval{Float64}(model, domain)
    update_eval!(multi)
    isnothing(cname)  &&  (return last_eval_folded(multi))
    return fold_model(multi, cname)
end
