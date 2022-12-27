using Serialization, JSON, GZip

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
serializable(f::λFunct) = λFunct(dummyfunct, f.display, deepcopy(f.args), deepcopy(f.optargs))

# Create a copies of GFit objects replacing all non-serializable items
# objects (such as component instances and λFunct) with dummy ones.
function serializable!(par::Parameter)
    if isa(par.patch, λFunct)
        par.patch  = serializable(par.patch)
    end
    if isa(par.mpatch, λFunct)
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

function todict(vv)
    tt = typeof(vv)
    @assert isstructtype(tt) "Unsupported type: $(string(tt))"
    @assert parentmodule(tt) == GFit
    if tt <: AbstractComponent
        @assert tt == DummyComp "Can't serialize instances of $(string(tt))"
    end
    out = OrderedDict{Symbol, Any}()
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

function todict(vv::AbstractDict{Symbol,T}) where {T <: Any}
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


function snapshot(filename::String, args...)
    serialize(filename, serializable.(args))
end

function snapshot_json(filename::String, args...; compress=false)
    if compress
        fname = deepcopy(filename)
        if fname[end-2:end] != ".gz"
            fname *= ".gz"
        end
        io = GZip.open(fname, "w")
    else
        io = open(filename, "w")  # io = IOBuffer()
    end
    JSON.print(io, todict(serializable.(args)))
    close(io)                 # String(take!(io))
end

