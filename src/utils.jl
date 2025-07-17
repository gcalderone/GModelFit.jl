version() = Pkg.TOML.parsefile(joinpath(pkgdir(GModelFit), "Project.toml"))["version"]


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


"""
    mock(::Type{Measures}, model::Model; keywords...)
    mock(::Type{Measures}, multi::Vector{Model}; keywords...)

Generate mock dataset(s) using a ground truth `Model` or `Vector{Model}` object. The first version returns a single `Measures` object, while the second returns a `Vector{Measures}`.

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
function mock(::Type{Measures}, multi::MultiEval; properr=0.01, rangeerr=0.05, abserr=0., seed=nothing)
    update_eval!(multi)
    out = Vector{Measures}()
    for i in 1:length(multi)
        values = last_eval_folded(multi, i)
        ee = extrema(values)
        range = ee[2] - ee[1]
        @assert range > 0
        err = (properr .* abs.(values) .+ rangeerr .* range .+ abserr)
        values .+= err .* randn(MersenneTwister(seed), size(values))
        push!(out, Measures(folded_domain(multi.v[i].ireval), values, err))
    end
    return out
end

mock(::Type{T}, model::Model, domain::AbstractDomain; kws...) where T =
    mock(T, MultiEval(model, domain); kws...)[1]

mock(::Type{T}, models::Vector{Model}, domains::Vector{<: AbstractDomain}; kws...) where T =
    mock(T, MultiEval(models, domains); kws...)
