# DummyComp is used to create serializable snapshots of actual components
struct DummyComp <: AbstractComponent
    tname::String
    deps::Vector{Symbol}
    params::OrderedDict{Symbol, Parameter}
end

getproperty(comp::DummyComp, pname::Symbol) = getfield(comp, :params)[pname]
getparams(comp::DummyComp) = getfield(comp, :params)
dependencies(comp::DummyComp) = getfield(comp, :deps)
original_type(comp::DummyComp) = getfield(comp, :tname)

function evaluate!(buffer::Vector{Float64}, comp::DummyComp, domain::AbstractDomain, pars...)
    @warn "Can't evaluate a dummy component!"
    nothing
end

# Named function used to create serializable FunctDesc objects
function dummyfunct(args...)
    @warn "Can't evaluate a dummy function!"
    nothing
end


# The `serializable` methods are supposed to create copies of GFit
# objects replacing all non-serializable items (such as
# component instances and FunctDesc) with dummy ones.
serializable(f::FunctDesc) = FunctDesc(dummyfunct, f.display, deepcopy(f.args), deepcopy(f.optargs))


function serializable!(par::Parameter)
    if isa(par.patch, FunctDesc)
        par.patch  = serializable(par.patch)
    end
    if isa(par.mpatch, FunctDesc)
        par.mpatch = serializable(par.mpatch)
    end
end


function serializable(ceval::CompEval, domain::AbstractDomain)
    tname = collect(split(string(typeof(ceval)), ['{', ',']))[2]
    (tname == "DummyComp")  &&  (return tname)

    deps = deepcopy(dependencies(ceval.comp))
    params = deepcopy(getparams(ceval.comp))
    for (pname, par) in params
        serializable!(par)
    end
    out = CompEval(DummyComp(tname, deps, params), domain)
    for fname in fieldnames(typeof(ceval))
        (fname == :comp)  &&  continue
        setfield!(out, fname, getfield(ceval, fname))
    end
    return out
end


function serializable(source::Model)
    model = deepcopy(source)
    model.parent = nothing
    for (cname, ceval) in source.cevals
        model.cevals[cname] = serializable(ceval, model.domain)
    end
    model.maincomp = find_maincomp(source)
    return model
end


function serializable(source::MultiModel)
    model = deepcopy(source)
    for i in 1:length(model)
        model.models[i] = serializable(source[i])
    end
    return model
end


serializable(v::AbstractDomain) = v
serializable(v::AbstractMeasures) = v
function serializable(source::FitResult)
    res = deepcopy(source)
    if isa(res.bestfit, Vector)
        for i in 1:length(res.bestfit)
            bb = res.bestfit[i]
            for p in getfield(bb, :data)
                serializable!(p)
            end
        end
    else
        bb = res.bestfit
        for p in getfield(bb, :data)
            serializable!(p)
        end
    end
    return res
end


"""
    snapshot(filename::String, args...)

Save a binary snapshot of one (or more) GFit object(s) such as `Model`, `MultiModel`, `Domain`, `Measures`, etc using the standard `Serialization` package.  The snapshot can be restored in a later session, and the objects will be similar to the original ones, with the following notable differences:
- in `Model` objects, all components are casted into `GFit.DummyComp` ones.  The original type is availble (as a string) via the `original_type()` function, while the content of the original structure is lost;
- all `FunctDesc` objects retain their textual representation, but the original function is lost;
- `Model` and `MultiModel` objects, as well as all the components, retain their last evaluated values but they can no longer be evaluated (an attempt to invoke `evaluate()` will result in an error);

The reason to introduce such differences is to ensure that all data structures being serialized are defined within the `GFit` package with no further external depencency, and to overcome limitation of the `Serialization` package related to, e.g., anonymous functions.

## Example:
```julia-repl
# Create GFit objects
using GFit
dom  = Domain(1:5)
model = Model(dom, :linear => @λ (x, b=2, m=0.5) -> (b .+ x .* m))
data = Measures(dom, [4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
res = fit!(model, data)

# Save a snapshot
GFit.snapshot("my_snapshot.dat", (model, data, res))

# Restore snapshot (possibly in a different Julia session)
using Serialization, GFit
(model, data, res) = deserialize("my_snapshot.dat")
```

!!! note
    The GFit binary serialization facility is **experimental**,

See also `GFit.snapshot_json()`.
"""
function snapshot(filename::String, arg)
    serialize(filename, serializable(args))
    return filename
end

function snapshot(filename::String, args::Tuple)
    serialize(filename, serializable.(args))
    return filename
end


function todict(vv)
    tt = typeof(vv)
    @assert isstructtype(tt) "Unsupported type: $(string(tt))"
    # @assert parentmodule(tt) == GFit
    if tt <: AbstractComponent
        @assert tt == DummyComp "Can't serialize instances of $(string(tt))"
    end
    out = OrderedDict{Symbol, Any}()
    out[:_structtype] = string(tt)
    for fname in fieldnames(tt)
        if  (tt == GFit.FunctDesc)  && 
            (fname == :funct)
            out[fname] = nothing
        else
            out[fname] = todict(getfield(vv, fname))
        end
    end

    # Add :show key
    plain = GFit.showsettings.plain
    GFit.showsettings.plain = true
    try
        io = IOBuffer()
        show(io, vv)
        out[:show] = String(take!(io))
    catch err
        @warn "Exception caught while invoking show($(tt))"
        println(err)
    end
    GFit.showsettings.plain = plain
    return out
end

function todict(vv::AbstractDict{Symbol,T}) where {T <: Any}
    out = OrderedDict{Symbol, Any}()
    out[:_dicttype] = string(typeof(vv))
    for (key, val) in vv
        out[key] = todict(val)
    end
    return out
end

todict(::Nothing) = nothing
todict(v::DateTime) = string(v)
todict(v::String) = v
todict(v::Symbol) = string(v)
todict(v::Number) = v
todict(v::AbstractArray) = todict.(v)
todict(v::Tuple) = todict.(v)


function snapshot_json(filename::String, args...; compress=false)
    filename = ensure_file_extension(filename, "json")
    if compress
        filename = ensure_file_extension(filename, "gz")
        io = GZip.open(filename, "w")
    else
        io = open(filename, "w")  # io = IOBuffer()
    end
    JSON.print(io, todict(serializable.(args)))
    close(io)                     # String(take!(io))
    return filename
end

