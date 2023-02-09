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
_serialize(v::Tuple) = _serialize.(v)
_serialize(v::Number) = (isnan(v)  ||  isinf(v)  ?  "_TN_" * string(v)  :  v)

function _serialize(vv::AbstractDict{Symbol,T}) where {T <: Any}
    out = OrderedDict{Symbol, Any}()
    out[:_dicttype] = string(typeof(vv))
    for (key, val) in vv
        out[key] = _serialize(val)
    end
    return out
end

function _serialize_struct(vv; add_show=false)
    out = OrderedDict{Symbol, Any}()
    out[:_structtype] = string(typeof(vv))
    for field in fieldnames(typeof(vv))
        ff = getfield(vv, field)
        if hasmethod(_serialize, (typeof(ff),))
            out[field] = _serialize(ff)
        end
    end
    if add_show
        plain = GFit.showsettings.plain
        GFit.showsettings.plain = false
        io = IOBuffer()
        show(io, vv)
        out[:show] = String(take!(io))
        GFit.showsettings.plain = plain
    end
    return out
end

_serialize(vv::Union{HashVector, HashHashVector, FunctDesc, Parameter}) = _serialize_struct(vv)
_serialize(vv::Model) = _serialize(ModelBuffers(vv))
_serialize(vv::ModelBuffers) = _serialize_struct(vv, add_show=true)
_serialize(vv::FitResult) = _serialize_struct(vv, add_show=true)
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

function _deserialize(dict::AbstractDict)
    function deserialized_function(args...)
        @warn "Can't evaluate a deserialized function"
        nothing
    end

    dd = OrderedDict{Symbol, Any}()
    for (kk, vv) in dict
        dd[Symbol(kk)] = _deserialize(vv)
    end
    if :_structtype in keys(dd)
        println()
        @info dd[:_structtype]
        dump(dd)
        println()
        if dd[:_structtype] == "GFit.HashVector{GFit.Parameter}"
            tmp = HashVector{Parameter}(dd[:data])
            for (k, v) in dd[:dict]
                @info k
                @info v
                getfield(tmp, :dict)[k] = v
            end
            return tmp
        elseif dd[:_structtype] == "GFit.HashHashVector{GFit.Parameter}"
            tmp = HashHashVector{Parameter}()
            for (k, v) in dd[:dict]
                getfield(tmp, :dict)[k] = v
            end
            append!(getfield(tmp, :data), dd[:data])
            return tmp
        elseif dd[:_structtype] == "GFit.FunctDesc"
            return FunctDesc(deserialized_function, dd[:display], dd[:args], d[:optargs])
        elseif dd[:_structtype] == "GFit.Parameter"
            return Parameter(dd[:val], dd[:low], dd[:high], dd[:fixed],
                             dd[:patch], dd[:mpatch], dd[:actual], dd[:unc])
        elseif dd[:_structtype] == "GFit.ModelBuffers"
            return ModelBuffers(dd[:domain], dd[:buffers], dd[:maincomp])
        elseif dd[:_structtype] == "GFit.FitResult"
            return FitResult(dd)
        elseif dd[:_structtype] == "GFit.AbstractDomain"
            # TODO
        elseif dd[:_structtype] == "GFit.AbstractMeasures"
            # TODO
        end
    end
    return dd
end



"""
    serialize(filename::String, args...; compress=false)

Save a snapshot of one (or more) GFit object(s) such as `Model`, `MultiModel`, `Domain`, `Measures`, etc using the JSON format.  The snapshot can be restored in a later session with `deserialize`, and the objects will be similar to the original ones, with the following notable differences:
- in `Model` objects, all components are casted into `GFit.DummyComp` ones.  The original type is availble (as a string) via the `original_type()` function, while the content of the original structure is lost;
- all `FunctDesc` objects retain their textual representation, but the original function is lost;
- `Model` and `MultiModel` objects, as well as all the components, retain their last evaluated values but they can no longer be evaluated (an attempt to invoke `evaluate()` will result in an error);

The reason to introduce such differences is to ensure that all data structures can be serialized to a JSON format and can be safely eserialized.

## Example:
```julia-repl
# Create GFit objects
using GFit
dom  = Domain(1:5)
model = Model(dom, :linear => @λ (x, b=2, m=0.5) -> (b .+ x .* m))
data = Measures(dom, [4.01, 7.58, 12.13, 19.78, 29.04], 0.4)
res = fit!(model, data)

# Save a snapshot
GFit.serialize("my_snapshot.json", model, data, res)

# Restore snapshot (possibly in a different Julia session)
using Serialization, GFit
(model, data, res) = deserialize("my_snapshot.json")
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
    JSON.print(io, _serialize(args))
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
    return _deserialize(j)
end
