stype(::T) where T = replace(string(T), "GModelFit." => "")

# --------------------------------------------------------------------
# Data statistic methods

loglh(lh::AbstractLikelihood) = (evaluate!(lh); sum([loglh(group) for group in values(getgroups(lh))]))
loglh(group::GroupLHSnapshot) = group.loglh
loglh(group::GroupLH{T, <: GaussianData})  where T = -0.5 * gofstat(group)
loglh(group::GroupLH{T, <: PoissonCounts}) where T = -0.5 * gofstat(group)

gofstat(lh::AbstractLikelihood) = (evaluate!(lh); sum([gofstat(group) for group in values(getgroups(lh))]))

gofstat_name(group::GroupLHSnapshot) = group.gofstat_name
gofstat(group::GroupLHSnapshot) = group.gofstat

gofstat_name(group::GroupLH{T, <: GaussianData}) where T = "Chi-squared"
gofstat(group::GroupLH{T, <: GaussianData}) where T = sum(abs2, group.residuals)

gofstat_name(group::GroupLH{T, <: PoissonCounts}) where T = "C-stat"
function gofstat(group::GroupLH{T, <: PoissonCounts}) where T
    D = values(group.data)
    M = group.evalbuf
    cstat = zero(T)
    @inbounds for i in eachindex(D, M)
        d = D[i]
        m = M[i]

        # Protect against domain errors if the solver explores negative model space
        m_safe = max(m, eps(Float64))
        if d > 0
            cstat += m_safe - d + d * log(d / m_safe)
        else
            cstat += m_safe
        end
    end
    return 2 * cstat
end

evaluate_resid!(group::GroupLH{T, <: GaussianData})  where T = group.residuals .= view((group.evalbuf .- values(group.data)) ./ uncerts(group.data), :)
function evaluate_resid!(group::GroupLH{T, <: PoissonCounts}) where T
    D = values(group.data)
    M = group.evalbuf
    R = group.residuals

    @inbounds for i in eachindex(D, M, R)
        d = D[i]
        m = max(M[i], eps(Float64))

        # Calculate the C-stat contribution for this specific bin
        c_i = 2 * (d > 0 ? (m - d + d * log(d / m)) : m)

        # Residuals are the signed square root of the C-stat contribution (Anscombe transform)
        R[i] = sign(d - m) * sqrt(c_i)
    end
    return R
end

# --------------------------------------------------------------------
# AbstractComponent methods
"""
    prepare!(comp::AbstractComponent, domain::AbstractDomain)

Invoked to precompute component-specific quantities

This method is invoked only once when the component is first evaluated hence it is the perfect place to pre-compute quantities associated to a component evaluation on a specific domain.
"""
prepare!(comp::AbstractComponent, domain::AbstractDomain) = nothing

"""
    dependencies(comp::AbstractComponent)
    dependencies(comp::CompSnapshot)

Returns a `Vector{Symbol}` containing the name of dependecies.

Default implementation returns `Symbol[]` (i.e. no dependency).
"""
dependencies(comp::AbstractComponent) = Symbol[]  # fall back method, no dependencies

"""
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, param1, param2....
    evaluate!(comp::AbstractComponent, domain::AbstractDomain, output::Abstractvector, deps::AbstractVector, param1, param2....

Evaluate component `comp` on the given `domain` using `deps` dependencies and `param1`, `param2, ... parameters.  Output should be stored in the `output` vector.

If the component has no dependencies the `deps` argument should not be present.

The `evaluate!` function is called with the `output`, `deps` and parameter arguments containing either `Float64` values (to evaluate the component) or `ForwardDiff.Dual` values (to evaluate the component derivatives).
"""
evaluate!(::TComp, ::TDomain, args...) where {TComp <: AbstractComponent, TDomain <: AbstractDomain} =
    error("No evaluate method implemented for $(TComp), $(TDomain)")

# Collect parameter from a component and return a OrderedDict{Symbol, Parameter}
function getparams(comp::T) where T <: AbstractComponent
    dict_fields = Vector{Symbol}()
    for name in fieldnames(T)
        if fieldtype(T, name) <: AbstractDict{Symbol, Parameter}
            push!(dict_fields, name)
        end
    end

    out = OrderedDict{Symbol, Parameter}()
    for name in fieldnames(T)
        if fieldtype(T, name) == Parameter
            out[name] = getfield(comp, name)
        end
    end

    if (length(dict_fields) > 1)  ||  ((length(dict_fields) > 0)  &&  (length(out) > 0))
        error("Ambiguous set of Parameters in structure $T.  You need to implement a custom GModelFit.getparams(::$(T)) method")
    end

    if length(dict_fields) == 1
        getfield(comp, dict_fields[1])
    end
    return out
end

# Quick evaluation: parameter values are the ones stored in the component unless a
# custom value is provided via a keyword.
function (comp::AbstractComponent)(domain::AbstractDomain; kws...)
    @assert length(dependencies(comp)) == 0 "Can't evaluate a component with dependencies as a stand-alone one."
    model = ModelEval(domain, :main => comp)
    pvalues = OrderedDict{Symbol, Float64}()
    for (pname, par) in getparams(comp)
        @assert isnothing(par.patch)
        @assert isnothing(par.reparam)
    end
    for (pname, pval) in kws
        model[(:main, pname)].val = pval
    end
    evaluate!(model)
    return model[:main]()
end

# --------------------------------------------------------------------
# Parameter level
# Upon construction a Parameter is not dirty, the corresponding component is already flagged as dirty.
Parameter(value::Real) = Parameter(float(value), -Inf, +Inf, false, nothing, nothing)
ParameterEval(p::Parameter) = ParameterEval([getproperty(p, s) for s in fieldnames(Parameter)]..., false, false, NaN, NaN)

isdirty(object::ParameterEval) = getfield(object, :_dirty)
function setdirty!(object::ParameterEval, dirty::Bool=true)
    setfield!(object, :_dirty, dirty)
    if dirty
        setfield!(object, :mval, NaN)
        setfield!(object, :unc , NaN)
    end
end

Base.propertynames(object::ParameterEval) = (:val, :low, :high, :fixed, :actually_fixed, :patch, :reparam, :unc, :mval)
function Base.setproperty!(object::ParameterEval, field::Symbol, value)
    (field in fieldnames(Parameter))  ||  error("Field $field is not supposed to be modified by the user")
    setdirty!(object, true)
    return setfield!(object, field, convert(fieldtype(ParameterEval, field), value))
end

# --------------------------------------------------------------------
# Component level
isdirty(object::CompSlotEval) = object._dirty
setdirty!(object::CompSlot, dirty::Bool=true) = dirty
setdirty!(object::CompSlotEval, dirty::Bool=true) = object._dirty = dirty

getcomp(comp::AbstractCompSlot) = getfield(comp, :comp)
getcomp(comp::CompSlotEval) = (setdirty!(comp); getfield(comp, :comp))

dependencies(comp::AbstractCompSlot) = dependencies(getcomp(comp))
dependencies(comp::CompSnapshot) = getfield(comp, :deps)

stype(comp::AbstractCompSlot) = stype(getcomp(comp))
stype(comp::CompSnapshot) = comp.stype

"""
    isfrozen(comp)

Check whether a component is *frozen*.
"""
isfrozen(comp::AbstractCompSlot) = comp.frozen

"""
    freeze!(comp)

Freeze a component in the model (i.e. treat all component parameters as fixed for fitting).
"""
freeze!(comp::AbstractCompSlot) = (setdirty!(comp); comp.frozen = true; return comp)

"""
    thaw!(comp)

Thaw a frozen component in the model (i.e. treat component parameters as fixed only if explicitly set in the corresponding `Parameter` structure).
"""
thaw!(comp::AbstractCompSlot) = (setdirty!(comp); comp.frozen = false; return comp)

domain(comp::AbstractCompSlotEval) = comp.domain

function (comp::CompSlotEval)()
    evaluate!(comp.rootmodel)
    return comp.evalbuf
end
(comp::CompSnapshot)() = comp.evalbuf

# --------------------------------------------------------------------
# Group level
domain(group::AbstractGroupEval) = folded_domain(getresp(group))

getresp(group::AbstractGroupEval) = group.resp
getdata(group::AbstractGroupLH) = group.data
getresiduals(group::AbstractGroupLH) = (evaluate!(group.rootmodel); group.residuals)

function (group::AbstractGroupEval)()
    evaluate!(group.rootmodel)
    return group.evalbuf
end
(group::GroupLHSnapshot)() = group.evalbuf

# --------------------------------------------------------------------
# Model level
NumType(autodiff::Bool) = autodiff  ?  Union{Dual, Float64}  :  Float64
NumType(::ModelEval{T}) where T = T
NumType(::Likelihood{T}) where T = T

Model(p::Pair{CT}, args::Vararg{Pair{CT}}) where {CT <: Tuple} = Model(OrderedDict(p, args...))
function Model(dict::AbstractDict{CT}) where {CT <: Tuple}
    @assert fieldcount(CT) >= 1
    GT = Tuple{fieldtypes(CT)[1:(fieldcount(CT)-1)]...}

    model = Model{GT}()
    for (name, item) in dict
        if isa(item, FunctDesc)
            setindex!(model, FComp(item), name)
        else
            @assert isa(item, AbstractComponent)
            setindex!(model, item, name)
        end
    end
    return model
end

Model(p::Pair{Symbol}, args::Vararg{Pair{Symbol}}) = Model(OrderedDict(p, args...))
Model(dict::AbstractDict{Symbol})                  = Model(OrderedDict([(k,) => v for (k, v) in dict]...))
Model(arg::Union{<: AbstractComponent, FunctDesc}) = Model(OrderedDict(:main => arg))


ModelEval(resp::AbstractDict{GT, <: AbstractResponse}, p::Pair{CT}, args::Vararg{Pair{CT}}; kws...) where {GT <: Tuple, CT <: Tuple} = ModelEval(resp, OrderedDict(p, args...); kws...)
function ModelEval(resp::AbstractDict{GT, <: AbstractResponse}, dict::AbstractDict{CT}; kws...) where {GT <: Tuple, CT <: Tuple}
    @assert fieldcount(CT) >= 1
    @assert GT == Tuple{fieldtypes(CT)[1:(fieldcount(CT)-1)]...}
    model = ModelEval(resp; kws...)
    for (name, item) in dict
        if isa(item, FunctDesc)
            setindex!(model, FComp(item), name)
        else
            @assert isa(item, AbstractComponent)
            setindex!(model, item, name)
        end
    end
    return model
end

ModelEval(resp::AbstractDict{Tuple{}, <: AbstractResponse}, p::Pair{Symbol}, args::Vararg{Pair{Symbol}}; kws...) = ModelEval(resp, OrderedDict(p, args...)                       ; kws...)
ModelEval(resp::AbstractDict{Tuple{}, <: AbstractResponse}, dict::AbstractDict{Symbol}                 ; kws...) = ModelEval(resp, OrderedDict([(k,) => v for (k, v) in dict]...); kws...)
ModelEval(resp::AbstractDict{Tuple{}, <: AbstractResponse}, arg::Union{<: AbstractComponent, FunctDesc}; kws...) = ModelEval(resp, OrderedDict(:main => arg)                     ; kws...)
ModelEval(domain::AbstractDomain                      , args...; kws...)          = ModelEval(OrderedDict(() => IdentityResp(domain))                     , args...; kws...)
ModelEval(domains::AbstractDict{GT, <: AbstractDomain}, args...; kws...) where GT = ModelEval(OrderedDict([k => IdentityResp(v) for (k, v) in domains]...), args...; kws...)
ModelEval(resp::AbstractResponse                   , args...; kws...)          = ModelEval(OrderedDict(() => resp)                                        , args...; kws...)


Likelihood(data::AbstractDict{GT, <: AbstractData}, p::Pair{CT}, args::Vararg{Pair{CT}}; kws...) where {GT <: Tuple, CT <: Tuple} = Likelihood(resp, OrderedDict(args...); kws...)
function Likelihood(data::AbstractDict{GT, <: AbstractData}, dict::AbstractDict{CT}; kws...) where {GT <: Tuple, CT <: Tuple}
    @assert fieldcount(CT) >= 1
    @assert GT == Tuple{fieldtypes(CT)[1:(fieldcount(CT)-1)]...}
    model = Likelihood(data; kws...)
    for (name, item) in dict
        if isa(item, FunctDesc)
            setindex!(model, FComp(item), name)
        else
            @assert isa(item, AbstractComponent)
            setindex!(model, item, name)
        end
    end
    return model
end

Likelihood(data::AbstractDict{Tuple{}, <: AbstractData}, p::Pair{Symbol}, args::Vararg{Pair{Symbol}}; kws...) = Likelihood(data, OrderedDict(p, args...)                       ; kws...)
Likelihood(data::AbstractDict{Tuple{}, <: AbstractData}, dict::AbstractDict{Symbol}                 ; kws...) = Likelihood(data, OrderedDict([(k,) => v for (k, v) in dict]...); kws...)
Likelihood(data::AbstractDict{Tuple{}, <: AbstractData}, arg::Union{<: AbstractComponent, FunctDesc}; kws...) = Likelihood(data, OrderedDict(:main => arg)                     ; kws...)
Likelihood(data::AbstractData, args...; kws...) = Likelihood(OrderedDict(() => data), args...; kws...)


function isdirty(model::AbstractModelEval)
    for (ckey, comp) in getcomps(model)
        comp._dirty  &&  (return true)
    end
    for (pkey, par) in getparams(model)
        par._dirty   &&  (return true)
    end
    return false
end

function setdirty!(model::AbstractModelEval, dirty::Bool=true)
    for (ckey, comp) in getcomps(model)
        setfield!(comp, :_dirty, dirty)
    end
    for (pkey, par) in getparams(model)
        setfield!(par , :_dirty, dirty)
    end
end

domain(model::AbstractModel{Tuple{}}) = domain(model[()])
(model::AbstractModelEval{Tuple{}})() = model[()]()
(model::AbstractModel{Tuple{}})(dom::AbstractDomain) = ModelEval(dom, model)[()]()

function getfreepars(model::AbstractModelEval)
    evaluate!(model);
    out = OrderedDict(getparams(model))
    filter!(p -> p[1] in collect(keys(out))[model.ifree], out)
    return out
end

nfree(model::AbstractModelEval) = (evaluate!(model); length(model.ifree))
ndata(model::AbstractLikelihood) = length(model.evalbuf)

function dof(lh::AbstractLikelihood)
    for (gkey, group) in getgroups(lh)
        isa(getdata(group), GaussianData)  ||  (return nothing)  # Can't compute DOF for a non-chisquared statistic
    end
    return ndata(lh) - nfree(lh)
end

# --------------------------------------------------------------------
# Dict related methods
group_key( ::AbstractModel{GT, CT, PT}, key::PT) where {GT, CT, PT} = key[1:fieldcount(GT)]
group_key( ::AbstractModel{GT, CT, PT}, key::CT) where {GT, CT, PT} = key[1:fieldcount(GT)]
comp_key(  ::AbstractModel{GT, CT, PT}, key::PT) where {GT, CT, PT} = key[1:fieldcount(CT)]
copa_key(  ::AbstractModel{GT, CT, PT}, key::PT) where {GT, CT, PT} = key[( fieldcount(GT)+1):(fieldcount(CT)+1)]
comp_name( ::AbstractModel{GT, CT, PT}, key::PT) where {GT, CT, PT} = key[  fieldcount(CT)]
comp_name( ::AbstractModel{GT, CT, PT}, key::CT) where {GT, CT, PT} = key[end]
param_name(::AbstractModel{GT, CT, PT}, key::PT) where {GT, CT, PT} = key[end]

_tree(    object::AbstractModel) = object.dt
_tree(    object::Union{AbstractCompSlot, AbstractGroup}) = view(_tree(object.rootmodel), object.key)

keys(     object::Union{AbstractCompSlot, AbstractGroup, AbstractModel})          = keys(    _tree(object))
length(   object::Union{AbstractCompSlot, AbstractGroup, AbstractModel})          = length(  _tree(object))
haskey(   object::Union{AbstractCompSlot, AbstractGroup, AbstractModel}, key)     = haskey(  _tree(object), key)
values(   object::Union{AbstractCompSlot, AbstractGroup, AbstractModel})          = values(  _tree(object))
iterate(  object::Union{AbstractCompSlot, AbstractGroup, AbstractModel}, args...) = iterate( _tree(object), args...)
getindex( object::Union{AbstractCompSlot, AbstractGroup, AbstractModel}, key)     = getindex(_tree(object), key)

function getparams(object::Union{AbstractCompSlot, AbstractGroup})
    dt = get_layer(object.rootmodel.dt, :params)
    haspath(dt, object.key)  &&  (return view(dt, object.key))
    if isa(object, CompSlot)  ||  isa(object, Group)
        return SDTree{StaticDictTrees.getKT(dt), Parameter}()
    end
    return SDTree{StaticDictTrees.getKT(dt), ParameterEval}()
end
getparams(object::AbstractModel)                          =      get_layer(_tree(object)      , :params)
getcomps( object::AbstractGroup)                          = view(get_layer(object.rootmodel.dt, :comps ), object.key)
getcomps( object::AbstractModel)                          =      get_layer(_tree(object)      , :comps )
getgroups(object::AbstractModel)                          =      get_layer(_tree(object)      , :groups)

setindex!(m::AbstractModel{GT, CT} , f::FunctDesc, ckey::CT) where {GT, CT}  = setindex!(m, FComp(f), ckey)
setindex!(m::AbstractModel{Tuple{}}, f::FunctDesc, cname::Symbol)            = setindex!(m, FComp(f), (cname,))
setindex!(m::AbstractModel{Tuple{}}, comp::AbstractComponent, cname::Symbol) = setindex!(m, comp, (cname,))

function setindex!(model::Model{GT, CT}, comp::AbstractComponent, ckey::CT) where {GT, CT}
    gkey = group_key(model, ckey)
    haskey(model, gkey)  ||  (model.dt[gkey] = Group(model, gkey))
    model.dt[ckey] = CompSlot(model, ckey, deepcopy(comp))
    for (pname, par) in getparams(comp)
        model.dt[(ckey..., pname)] = par
    end
    return model.dt[ckey]
end

function setindex!(model::AbstractModelEval{GT, CT}, comp::AbstractComponent, ckey::CT) where {GT, CT}
    gkey = group_key(model, ckey)
    @assert haskey(model, gkey)  # we do not have a new response here, hence group must already exists
    model.dt[ckey] = CompSlotEval(model, ckey, deepcopy(comp), unfolded_domain(model[gkey].resp))
    for (pname, par) in getparams(comp)
        model.dt[(ckey..., pname)] = ParameterEval(par)
    end
    return model.dt[ckey]
end

# --------------------------------------------------------------------
# Actual fit
include("evaluation.jl")
include("solvers.jl")
using .Solvers

function fit!(lh::Likelihood{T}; solver::AbstractSolver=Solvers.lsqfit()) where T
    (nfree(lh) > 0)  ||  error("No free parameter in the model")
    (T == Float64)  &&  need_autodiff(solver)  &&  error("Solver $(typeof(solver)) requires automatic differentiation (autodiff=true)")
    return maximize!(lh, solver)
end

"""
    fit(mdesc::Model{Tuple{}}, data::AbstractData, solver=lsqfit())
    fit(mdesc::Model{GT}, data::AbstractDict{GT, <: AbstractData}, solver=lsqfit()) where GT

Fit a model to empirical data set(s) using the specified solver (default: `lsqfit()`).
"""
function fit(args...; solver::AbstractSolver=Solvers.lsqfit(), kws...)
    lh = Likelihood(args...; autodiff=need_autodiff(solver), kws...)
    fsumm = fit!(lh; solver=solver)
    return lh, fsumm
end

# --------------------------------------------------------------------
# Mock
"""
    mock(::Type{GaussianData}, mdesc::Model; keywords...)
    mock(::Type{GaussianData}, multi::Vector{Model}; keywords...)

Generate mock dataset(s) using a ground truth `Model` or `Vector{Model}` object. The first version returns a single `GaussianData` object, while the second returns a `Vector{GaussianData}`.

The measurement random errors added to the data points are drawn from a Normal distribution centered on the data value itself, and a width given by the sum of three contributions:
- *proportional* part: error proportional to each data point value;
- *range* part: error proportional to the range spanned by all values in a single dataset;
- *absolute* part: absolute error value.

No systematic error is considered when generating mock dataset(s).

# Accepted keywords:
- `properr=0.01`: proportional error;
- `rangeerr=0.05`: range error;
- `abserr=0.`: absolute error;
- `seed=nothing`: seed for the `Random.MersenneTwister` generator.
"""
# --------------------------------------------------------------------
# Mock
mock(::Type{T}, args...; kws...) where {T <: AbstractData} = mock(T, ModelEval(args...); kws...)
function mock(::Type{GaussianData}, model::ModelEval{Float64, GT}; properr=0.01, rangeerr=0.05, abserr=0., seed=nothing) where {GT}
    evaluate!(model)
    out = OrderedDict{GT, GaussianData}()
    rng = isnothing(seed)  ?  Random.default_rng()  :  MersenneTwister(seed)
    for (gkey, group) in getgroups(model)
        values = deepcopy(group.evalbuf)
        ee = extrema(values)
        range = ee[2] - ee[1]
        @assert range > 0
        err = (properr .* abs.(values) .+ rangeerr .* range .+ abserr)
        values .+= err .* randn(rng, size(values))
        out[gkey] = GaussianData(group.resp, values, err)
    end
    (GT == Tuple{})  &&  (return out[()])
    return out
end

function mock(::Type{PoissonCounts}, model::ModelEval{Float64, GT}; seed=nothing) where {GT}
    evaluate!(model)
    out = OrderedDict{GT, PoissonCounts}()
    rng = isnothing(seed) ? Random.default_rng() : MersenneTwister(seed)

    for (gkey, group) in getgroups(model)
        M = group.evalbuf
        counts = Array{Int}(undef, size(M))

        @inbounds for i in eachindex(M, counts)
            λ = max(M[i], 0.0) # Ensure expected rate is non-negative

            # Dependency-free Poisson random number generator
            if λ < 30.0  # Knuth's algorithm for small λ
                L = exp(-λ)
                k = 0
                p = 1.0
                while p > L
                    k += 1
                    p *= rand(rng)
                end
                counts[i] = k - 1
            else
                # Normal approximation for large λ
                counts[i] = max(0, round(Int, λ + sqrt(λ) * randn(rng)))
            end
        end

        out[gkey] = PoissonCounts(group.resp, counts)
    end

    (GT == Tuple{}) && (return out[()])
    return out
end
