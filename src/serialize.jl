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
        if hasmethod(_serialize, (typeof(ff),))
            out[String(field)] = _serialize(ff)
        end
    end
    if add_show
        io = IOBuffer()
        show(io, vv)
        out["show"] = String(take!(io))
    end
    return out
end

_serialize(vv::Union{HashVector, HashHashVector, FunctDesc, Parameter}) = _serialize_struct(vv)
_serialize(vv::Model) = _serialize(ModelBuffers(vv))
_serialize(vv::ModelBuffers) = _serialize_struct(vv, add_show=true)
_serialize(vv::FitResult) = _serialize_struct(vv, add_show=true)
_serialize(vv::MinimizerStatus) = _serialize_struct(MinimizerStatus(vv.code, vv.message, nothing))
_serialize(vv::MinimizerStatusCode) = Int(vv)
_serialize(vv::AbstractDomain) = _serialize_struct(vv, add_show=true)
_serialize(vv::AbstractMeasures) = _serialize_struct(vv, add_show=true)



# Deserialization methods
_deserialize(::Nothing) = nothing
_deserialize(v::AbstractVector) = _deserialize.(v)
_deserialize(v::Number) = v

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
        if dd["_structtype"] == "GFit.HashVector{GFit.Parameter}"
            out = HashVector{Parameter}(_deserialize(dd["data"]))
            for (k, v) in dd["dict"]
                getfield(out, :dict)[Symbol(k)] = _deserialize(v)
            end
            return out
        elseif dd["_structtype"] == "GFit.HashHashVector{GFit.Parameter}"
            out = HashHashVector{Parameter}()
            for (k, v) in dd["dict"]
                getfield(out, :dict)[Symbol(k)] = _deserialize(v)
            end
            append!(getfield(out, :data), _deserialize(dd["data"]))
            return out
        elseif dd["_structtype"] == "GFit.FunctDesc"
            return FunctDesc(deserialized_function,
                             _deserialize(dd["display"]),
                             _deserialize(dd["args"]),
                             _deserialize(dd["optargs"]))
        elseif dd["_structtype"] == "GFit.Parameter"
            return Parameter(_deserialize(dd["val"]),
                             _deserialize(dd["low"]),
                             _deserialize(dd["high"]),
                             _deserialize(dd["fixed"]),
                             _deserialize(dd["patch"]),
                             _deserialize(dd["mpatch"]),
                             _deserialize(dd["actual"]),
                             _deserialize(dd["unc"]))
        elseif dd["_structtype"] == "GFit.ModelBuffers"
            return ModelBuffers(_deserialize(dd["domain"]),
                                _deserialize(dd["buffers"]),
                                _deserialize(dd["maincomp"]),
                                _deserialize(dd["show"]))
        elseif dd["_structtype"] == "GFit.FitResult"
            return FitResult(_deserialize(dd["timestamp"]),
                             _deserialize(dd["elapsed"]),
                             _deserialize(dd["ndata"]),
                             _deserialize(dd["nfree"]),
                             _deserialize(dd["dof"]),
                             _deserialize(dd["fitstat"]),
                             _deserialize(dd["status"]),
                             _deserialize(dd["bestfit"]))
        elseif dd["_structtype"] == "GFit.MinimizerStatus"
            return MinimizerStatus(MinimizerStatusCode(_deserialize(dd["code"])),
                                   _deserialize(dd["message"]),
                                   _deserialize(dd["internal"]))
        elseif !isnothing(findfirst("CartesianDomain", dd["_structtype"]))
            axis = _deserialize(dd["axis"])
            roi =  _deserialize(dd["roi"])
            return CartesianDomain(axis..., roi=roi)
        elseif !isnothing(findfirst("Domain", dd["_structtype"]))
            axis = _deserialize(dd["axis"])
            return Domain(axis...)
        elseif !isnothing(findfirst("Measures", dd["_structtype"]))
            tmp = _deserialize(dd["values"])
            return Measures(_deserialize(dd["domain"]), tmp[1], tmp[2])
        end
    else
        out = OrderedDict{Symbol, Any}()
        for (kk, vv) in dd
            out[Symbol(kk)] = _deserialize(vv)
        end
        return out
    end

    error("TODO")
end



"""
    GFit.serialize(filename::String, args...; compress=false)

Save a snapshot of one (or more) GFit object(s) using a JSON format.  The snapshot can be restored in a different Julia session with `GFit.deserialize`, and the resulting objects will be similar to the original ones, with the following notable differences:
- `Model` objects are casted into `ModelBuffers` ones containing just the latest component evaluations;
- `MultiModel` objects are casted into `Vector{ModelBuffers}`;
- `FunctDesc` objects retain their textual representation, but the original function is lost;

## Example:
```julia-repl
# Create GFit objects
using GFit
dom  = Domain(1:5)
model = Model(dom, :linear => @λ (x, b=2, m=0.5) -> (b .+ x .* m))
data = Measures(dom, [4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
res = fit!(model, data)

# Save a snapshot
GFit.serialize("my_snapshot.json", [model, data, res])

# Restore snapshot (possibly in a different Julia session)
using GFit
(model, data, res) = GFit.deserialize("my_snapshot.json")
```

!!! note
    The GFit binary serialization facility is **experimental**,
"""
function serialize(filename::String, arg; compress=false)
    filename = ensure_file_extension(filename, "json")
    if compress
        filename = ensure_file_extension(filename, "gz")
        io = GZip.open(filename, "w")
    else
        io = open(filename, "w")
    end
    JSON.print(io, _serialize(arg))
    close(io)
    return filename
end


function deserialize(filename::String)
    if filename[end-2:end] == "gz"
        io = GZip.open(filename)
    else
        io = open(filename)
    end
    j = JSON.parse(io, dicttype=OrderedDict)
    close(io)
    return _deserialize(j)
end
