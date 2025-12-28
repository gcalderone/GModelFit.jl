JSONDict(meta::String, obj=nothing) = Dict("meta" => meta, "obj" => obj)

# Serialization methods: convert a Julia value into a type suitable to
# be serialized using JSON.
_serialize(::Nothing) = JSONDict("Nothing")
_serialize(::Missing) = JSONDict("Missing")
function _serialize(v::Number)
    if isnan(v)
        return JSONDict("NaN")
    elseif isinf(v)
        if v > 0
            return JSONDict("+Inf")
        else
            return JSONDict("-Inf")
        end
    end
    return v
end

_serialize(v::Expr) = return JSONDict("Expr", string(v))
_serialize(v::Date) = return JSONDict("Date", string(v))
_serialize(v::DateTime) = return JSONDict("DateTime", string(v))
_serialize(v::String) = v
_serialize(v::Symbol) = return JSONDict("Symbol", string(v))
_serialize(v::AbstractVector) = JSONDict("Vector", _serialize.(v))
_serialize(v::Tuple) = JSONDict("Tuple", [_serialize.(v)...])

function _serialize(dict::OrderedDict{Symbol,T}) where {T}
    out = OrderedDict{String, Any}()
    for (key, val) in dict
        out[string(key)] = _serialize(val)
    end
    return JSONDict("OrderedDict{Symbol}", out)
end

function _serialize_struct(str::T) where {T}
    dict = OrderedDict{String, Any}()
    for field in fieldnames(T)
        ff = getfield(str, field)
        dict[string(field)] = _serialize(ff)
    end
    return JSONDict(string(parentmodule(T)) * "." * string(nameof(T)), dict)
end

function _serialize(str::T) where {T}
    @assert isstructtype(T) "No _serialize() method defined for type $T"
    return _serialize_struct(str)
end

# Can't serialize solver's status, replace it with nothing
_serialize(s::FitSummary) =
    _serialize_struct(FitSummary(s.start, s.elapsed, s.ndata, s.nfree, s.fitstat, s.status, nothing))


"""
    GModelFit.serialize(filename::String, ::ModelSnapshot[, ::FitSummary[, ::Measures]]; compress=false)
    GModelFit.serialize(filename::String, ::Vector{ModelSnapshot}[, ::FitSummary[, ::Vector{Measures}]]; compress=false)

Serialize GModelFit object(s) using a JSON format. The serializable objects are:
- `ModelSnapshot` and `Vector{ModelSnapshot}` (mandatory argument);
- `FitSummary` (optional);
- `Measures` and `Vector{Measures}` (optional);

If `compress=true` the resulting JSON file will be compressed using GZip.
Objects can later be deserialized in a different Julia session with `GModelFit.deserialize`.

Note: The `GModelFit.serialize` function also accepts `Model` and `Vector{Model}` but they will be internally converted to `ModelSnapshot`(s).


## Example:
```julia-repl
# Create GModelFit objects
using GModelFit
model = Model(:linear => @fd (x, b=2, m=0.5) -> (b .+ x .* m))
data = Measures([4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
bestfit, stats = fit(model, data)

# Serialize objects and save in a file
GModelFit.serialize("my_snapshot.json", bestfit, stats, data)

# Restore objects (possibly in a different Julia session)
using GModelFit
(bestfit, stats, data) = GModelFit.deserialize("my_snapshot.json")
```
"""
function serialize(filename::String, args...; compress=false)
    data = _serialize([args...])
    filename = ensure_file_extension(filename, "json")
    if compress
        filename = ensure_file_extension(filename, "gz")
        io = GZip.open(filename, "w")
    else
        io = open(filename, "w")
    end
    JSON.json(io, data)
    close(io)
    return filename
end


# ====================================================================
# Deserialization methods: transform a JSON parsed value into its
# original Julia value
_deserialize(::Val{:Nothing}, obj) = nothing
_deserialize(::Val{:Missing}, obj) = missing
_deserialize(::Val{:NaN}, obj) = NaN
_deserialize(::Val{Symbol("+Inf")}, obj) = +Inf
_deserialize(::Val{Symbol("-Inf")}, obj) = -Inf
_deserialize(::Val{:Expr}, obj) = Meta.parse(obj)
_deserialize(::Val{:Date}, obj) = Date(obj)
_deserialize(::Val{:DateTime}, obj) = DateTime(obj)
_deserialize(::Val{:Symbol}, obj) = Symbol(obj)
_deserialize(v::Number) = v
_deserialize(v::String) = v

function _deserialize(::Val{:Vector}, v)
    out = _deserialize.(v)
    #@assert length(unique(typeof.(out))) == 1
    return out
end

_deserialize(::Val{:Tuple}, obj) = tuple(_deserialize.(obj)...)

function _deserialize(::Val{Symbol("OrderedDict{Symbol}")}, dict::OrderedDict{K,V}) where {K,V}
    out = OrderedDict{Symbol, Any}()
    for (key, obj) in dict
        out[Symbol(key)] = _deserialize(obj)
    end
    return out
end

function _deserialize(dict::OrderedDict)
    @assert "meta" in keys(dict)
    return _deserialize(Val(Symbol(dict["meta"])), dict["obj"])
end

_deserialize(::Val{Symbol("GModelFit.ModelSnapshot")}, dict::OrderedDict) = GModelFit.ModelSnapshot(_deserialize.(values(dict))...)
_deserialize(::Val{Symbol("GModelFit.Domain")}, dict::OrderedDict) = Domain(_deserialize(dict["axis"])...)
_deserialize(::Val{Symbol("GModelFit.CartesianDomain")}, dict::OrderedDict) = CartesianDomain(_deserialize(dict["axis"])..., roi=_deserialize(dict["roi"]))
_deserialize(::Val{Symbol("GModelFit.Parameter")}, dict::OrderedDict) = GModelFit.Parameter(_deserialize.(values(dict))...)
_deserialize(::Val{Symbol("GModelFit.ComponentSnapshot")}, dict::OrderedDict) = GModelFit.ComponentSnapshot(_deserialize.(values(dict))...)
_deserialize(::Val{Symbol("GModelFit.Solvers.FitSummary")}, dict::OrderedDict) = GModelFit.Solvers.FitSummary(_deserialize.(values(dict))...)

function deserialized_function(args...)
    @warn "Can't evaluate a deserialized function"
    nothing
end

_deserialize(::Val{Symbol("GModelFit.FunctDesc")},
             dict::OrderedDict) =
                 FunctDesc(deserialized_function,
                           _deserialize(dict["display"]),
                           _deserialize(dict["args"]),
                           _deserialize(dict["optargs"]))

_deserialize(::Val{Symbol("GModelFit.Solvers.SolverStatusOK")}, dict::OrderedDict) = GModelFit.Solvers.SolverStatusOK()
_deserialize(::Val{Symbol("GModelFit.Solvers.SolverStatusWarn")}, dict::OrderedDict) = GModelFit.Solvers.SolverStatusWarn(_deserialize.(values(dict))...)
_deserialize(::Val{Symbol("GModelFit.Solvers.SolverStatusError")}, dict::OrderedDict) = GModelFit.Solvers.SolverStatusError(_deserialize.(values(dict))...)

function _deserialize(::Val{Symbol("GModelFit.Measures")}, dict::OrderedDict)
    dom = _deserialize(dict["domain"])
    tmp = _deserialize(dict["values"])
    if isa(dom, CartesianDomain)
        return Measures(dom, reshape(tmp[1], size(dom)), reshape(tmp[2], size(dom)))
    else
        return Measures(dom, tmp[1], tmp[2])
    end
end

function deserialize(filename::String)
    if filename[end-2:end] == ".gz"
        io = GZip.open(filename)
    else
        io = open(filename)
    end
    j = JSON.parse(io, dicttype=OrderedDict)
    close(io)
    return _deserialize(j)
end


function test_serialization(args...)
    function comparedata(A::TA, B::TB) where {TA, TB}
        isa(B, Function)  &&  return
        @assert TA == TB
        if (TA <: AbstractVector)  ||  (TA <: Tuple)  ||  (TA <: AbstractDict)
            @assert length(A) == length(B)
            for i in eachindex(A)
                comparedata(getindex(A, i), getindex(B, i))
            end
        elseif isstructtype(TA)
            for i in 1:nfields(A)
                comparedata(getfield(A, i), getfield(B, i))
            end
        else
            @assert isequal(A, B)
        end
    end
    dd = deserialize(serialize(tempname(), args...))
    comparedata(dd, [args...])
end
