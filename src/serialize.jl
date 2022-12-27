using Serialization, JSON

struct DummyComp <: AbstractComponent
    tname::String
    deps::Vector{Symbol}
    params::OrderedDict{Symbol, Parameter}
end

getproperty(comp::DummyComp, pname::Symbol) = getfield(comp, :params)[pname]
getparams(comp::DummyComp) = getfield(comp, :params)
dependencies(comp::DummyComp) = getfield(comp, :deps)

function evaluate!(buffer::Vector{Float64}, comp::DummyComp, domain::AbstractDomain, pars...)
    @warn "Can't evaluate a dummy component!"
    nothing
end

function dummyfunct(args...)
    @warn "Can't evaluate a dummy function!"
    nothing
end


# Create a copy of a Model object replacing all non-serializable
# objects (such as component instances and λFunct) with dummy ones.
function serializable(source::Model)
    model = deepcopy(source)
    model.parent = nothing
    cevals = OrderedDict{Symbol, CompEval}()
    for (cname, ceval) in model.cevals
        tname = collect(split(string(typeof(ceval)), ['{', ',']))[2]
        deps = deepcopy(dependencies(ceval.comp))
        params = deepcopy(getparams(ceval.comp))
        for (pname, par) in params
            f = par.patch
            if isa(f, λFunct)
                par.patch  = λFunct(dummyfunct, f.display, deepcopy(f.args), deepcopy(f.optargs))
            end
            f = par.mpatch
            if isa(f, λFunct)
                par.mpatch = λFunct(dummyfunct, f.display, deepcopy(f.args), deepcopy(f.optargs))
            end
        end
        tmp = CompEval(DummyComp(tname, deps, params), model.domain)
        for fname in fieldnames(typeof(ceval))
            (fname == :comp)  &&  continue
            setfield!(tmp, fname, getfield(ceval, fname))
        end
        cevals[cname] = tmp
    end
    for (cname, ceval) in cevals
        model.cevals[cname] = cevals[cname]
    end
    return model
end


function todict(vv)
    @assert isstructtype(typeof(vv)) "Unsupported type: $(string(typeof(vv)))"
    out = OrderedDict{Symbol, Any}()
    tt = typeof(vv)
    if tt <: AbstractComponent
        @assert tt == DummyComp "Can't serialize instances of $(string(tt))"
    end
    out[:_structtype] = string(tt)
    for fname in fieldnames(tt)
        if  (tt == GFit.λFunct)  && 
            (fname == :funct)
            out[fname] = nothing
        else
            out[fname] = todict(getfield(vv, fname))
        end
    end
    return out
end

function todict(vv::AbstractDict{K,T}) where {K<:Any, T <: Any}
    @assert K == Symbol
    out = OrderedDict{Symbol, Any}()
    out[:_dicttype] = string(typeof(vv))
    for (key, val) in vv
        out[key] = todict(val)
    end
    return out
end

todict(::Nothing) = nothing
todict(v::String) = v
todict(v::Symbol) = string(v)
todict(v::Number) = v
todict(v::AbstractArray) = todict.(v)
todict(v::Tuple) = todict.(v)


function snapshot(filename::String, model::Model)
    serialize(filename, serializable(model))
end

function snapshot_json(filename::String, model::Model)
    io = open(filename, "w")  # io = IOBuffer()
    JSON.print(io, todict(serializable(model)))
    close(io)                 # String(take!(io))
end

