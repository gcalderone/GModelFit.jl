# --------------------------------------------------------------------
# Misc. utils
forbid_delete(key, args...) = error("Object can not be deleted ($key)")
forbid_update(key, args...) = error("Object can not be replaced ($key)")

flattenkeys(dict::AbstractDict, sep="_") = OrderedDict([Symbol(join(string.(k), sep)) => v for (k, v) in dict]...)

function _validateparams(model::Model)
    for (name, par) in getparams(model)
        @assert !isnan(par.low)
        @assert !isnan(par.high)
        @assert isfinite(par.val)
        @assert par.low < par.high
        @assert par.val >= par.low  "Value for $name is smaller than the minimum allowed value"
        @assert par.val <= par.high "Value for $name is larger than the maximum allowed value"
        @assert isnothing(par.patch)  ||  isnothing(par.reparam) "Parameter $name can either be patched to another value, or re-parametrized (reparam), but not both."
    end
end
# version() = Pkg.TOML.parsefile(joinpath(pkgdir(GModelFit), "Project.toml"))["version"]

function ensure_file_extension(_filename, _ext)
    filename = deepcopy(_filename)
    ext = "." * _ext
    nn = length(ext)
    if  (length(filename) <= nn)  ||
        (filename[(end-nn+1):end] != ext)
        filename *= ext
    end
    return filename
end
