#=
Serialization methods: output must be suitable to be stored in a file using the JSON format.
=#
_serialize(::Nothing) = nothing
_serialize(::Function) = nothing
_serialize(v::Expr) = "_TE_" * string(v)
_serialize(v::Date) = "_TD_" * string(v)
_serialize(v::DateTime) = "_TDT" * string(v)
_serialize(v::String) = v
_serialize(v::Symbol) = "_TS_" * string(v)
_serialize(v::AbstractVector) = _serialize.(v)
_serialize(v::Tuple) = [_serialize.(v)...]
_serialize(v::Number) = (isnan(v)  ||  isinf(v)  ?  "_TN_" * string(v)  :  v)

function _serialize(vv::AbstractDict{Symbol,T}) where {T <: Any}
    out = OrderedDict{String, Any}()
    for (key, val) in vv
        out[String(key)] = _serialize(val)
    end
    return out
end

function _serialize_struct(vv; add_show=false)
    out = OrderedDict{String, Any}()
    out["_structtype"] = string(typeof(vv))
    for field in fieldnames(typeof(vv))
        ff = getfield(vv, field)
        out[String(field)] = _serialize(ff)
    end
    if add_show
        io = IOBuffer()
        if showsettings.plain
            show(io , vv)
        else
            ctx = IOContext(io, :color => true)
            show(ctx, vv)
        end
        out["show"] = String(take!(io))
    end
    return out
end

_serialize(vv::PV.PVComp) = _serialize_struct(vv)
_serialize(vv::PV.PVModel) = _serialize_struct(vv)
_serialize(vv::Parameter) = _serialize_struct(vv)
_serialize(vv::FunctDesc) = _serialize_struct(vv)
_serialize(vv::FitStats) = _serialize_struct(vv, add_show=true)
_serialize(vv::ModelSnapshot) = _serialize_struct(vv, add_show=true)
_serialize(vv::AbstractMinimizerStatus) = _serialize_struct(vv, add_show=true)
_serialize(vv::AbstractDomain) = _serialize_struct(vv, add_show=true)
_serialize(vv::AbstractMeasures) = _serialize_struct(vv, add_show=true)


_serialize(model::ModelSnapshot, fitstats::FitStats                         ) = _serialize([model, fitstats])
_serialize(model::ModelSnapshot, fitstats::FitStats, data::AbstractMeasures ) = _serialize([model, fitstats, data])
_serialize(multi::Vector{Model        }                                     ) = _serialize(ModelSnapshot.(multi))
_serialize(multi::Vector{ModelSnapshot}, fitstats::FitStats                 ) = _serialize([multi, fitstats])
function _serialize(multi::Vector{ModelSnapshot}, fitstats::FitStats, data::Vector{T}) where T <: AbstractMeasures
    @assert length(multi) == length(data)
    _serialize([multi, fitstats, data])
end


"""
    GModelFit.serialize(filename::String, ::ModelSnapshot[, ::FitStats[, ::Measures]]; compress=false)
    GModelFit.serialize(filename::String, ::Vector{ModelSnapshot}[, ::FitStats[, ::Vector{Measures}]]; compress=false)

Serialize GModelFit object(s) using a JSON format. The serializable objects are:
- `ModelSnapshot` and `Vector{ModelSnapshot}` (mandatory argument);
- `FitStats` (optional);
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
(bestit, stats, data) = GModelFit.deserialize("my_snapshot.json")
```
"""
function serialize(filename::String, args...; compress=false)
    data = _serialize(args...)
    filename = ensure_file_extension(filename, "json")
    if compress
        filename = ensure_file_extension(filename, "gz")
        io = GZip.open(filename, "w")
    else
        io = open(filename, "w")
    end
    JSON.print(io, data)
    close(io)

    return filename
end


# ====================================================================
# Deserialization methods
_deserialize(::Nothing) = nothing
_deserialize(v::Number) = v
function _deserialize(v::AbstractVector)
    tmp = _deserialize.(v)
    tt = unique(typeof.(tmp))
    if length(tt) == 1
        out = similar(tmp, tt[1])
        out .= tmp
    else
        out = tmp
    end
    return out
end

function _deserialize(v::String)
    if length(v) > 4
        magic = v[1:4]
        if magic == "_TE_"
            return Meta.parse(v)
        elseif magic == "_TD_"
            return Date(v[5:end])
        elseif magic == "_TDT"
            return DateTime(v[5:end])
        elseif magic == "_TS_"
            return Symbol(v[5:end])
        elseif v == "_TN_Inf"
            return +Inf
        elseif v == "_TN_-Inf"
            return -Inf
        elseif v == "_TN_NaN"
            return NaN
        end
    end
    return v
end

function _deserialize(dd::AbstractDict)
    function deserialized_function(args...)
        @warn "Can't evaluate a deserialized function"
        nothing
    end

    if "_structtype" in keys(dd)
        if !isnothing(findfirst("PVComp{GModelFit.Parameter}", dd["_structtype"]))
            return PVComp{Parameter}(_deserialize(dd["pnames"]), _deserialize(dd["indices"]), _deserialize(dd["data"]))
        elseif !isnothing(findfirst("PVModel{GModelFit.Parameter}", dd["_structtype"]))
            return PVModel{Parameter}(_deserialize(dd["comps"]), _deserialize(dd["indices"]), _deserialize(dd["data"]))
        elseif dd["_structtype"] == "GModelFit.FunctDesc"
            return FunctDesc(deserialized_function,
                             _deserialize(dd["display"]),
                             _deserialize(dd["args"]),
                             _deserialize(dd["optargs"]))
        elseif dd["_structtype"] == "GModelFit.Parameter"
            return Parameter(_deserialize(dd["val"]),
                             _deserialize(dd["low"]),
                             _deserialize(dd["high"]),
                             _deserialize(dd["fixed"]),
                             _deserialize(dd["patch"]),
                             _deserialize(dd["mpatch"]),
                             _deserialize(dd["actual"]),
                             _deserialize(dd["unc"]))
        elseif dd["_structtype"] == "GModelFit.ModelSnapshot"
            return ModelSnapshot(_deserialize(dd["domain"]),
                                 _deserialize(dd["params"]),
                                 _deserialize(dd["buffers"]),
                                 _deserialize(dd["maincomp"]),
                                 _deserialize(dd["comptypes"]),
                                 _deserialize(dd["isfreezed"]),
                                 _deserialize(dd["deps"]),
                                 _deserialize(dd["evalcounters"]))
        elseif dd["_structtype"] == "GModelFit.FitStats"
            return FitStats(_deserialize(dd["elapsed"]),
                            _deserialize(dd["ndata"]),
                            _deserialize(dd["nfree"]),
                            _deserialize(dd["dof"]),
                            _deserialize(dd["fitstat"]),
                            _deserialize(dd["status"]))
        elseif dd["_structtype"] == "GModelFit.MinimizerStatusOK"
            return MinimizerStatusOK()
        elseif dd["_structtype"] == "GModelFit.MinimizerStatusDry"
            return MinimizerStatusDry()
        elseif dd["_structtype"] == "GModelFit.MinimizerStatusWarn"
            return MinimizerStatusWarn(_deserialize(dd["message"]))
        elseif dd["_structtype"] == "GModelFit.MinimizerStatusError"
            return MinimizerStatusError(_deserialize(dd["message"]))
        elseif !isnothing(findfirst("CartesianDomain", dd["_structtype"]))
            axis = _deserialize(dd["axis"])
            roi =  _deserialize(dd["roi"])
            return CartesianDomain(axis..., roi=roi)
        elseif !isnothing(findfirst("Domain", dd["_structtype"]))
            axis = _deserialize(dd["axis"])
            return Domain(axis...)
        elseif !isnothing(findfirst("Measures", dd["_structtype"]))
            dom = _deserialize(dd["domain"])
            tmp = _deserialize(dd["values"])
            if isa(dom, CartesianDomain)
                return Measures(dom, reshape(tmp[1], size(dom)), reshape(tmp[2], size(dom)))
            else
                return Measures(dom, tmp[1], tmp[2])
            end
        else
            error("Unrecognized structure in serialized data: " * dd["_structtype"])
        end
    else
        out = OrderedDict{Symbol, Any}()
        for (kk, vv) in dd
            out[Symbol(kk)] = _deserialize(vv)
        end
        return out
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
