version() = Pkg.TOML.parsefile(joinpath(pkgdir(GFit), "Project.toml"))["version"]


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


# ====================================================================
#= TODO
function print_param_covariance(fitres::FitStats;
                                select=nothing, sort=false, threshold=0.)
    @assert isa(fitres.status.internal, CMPFit.Result) "Minimizer is not CMPFit"

    parnames = String[]
    if isa(fitres.bestfit, Vector{GFit.HashHashVector{GFit.Parameter}})
        for i in 1:length(fitres.bestfit)
            for (cname, hv) in fitres.bestfit[i]
                for (pname, par) in hv
                    par.fixed  &&  continue
                    push!(parnames, "[$(i)][$(cname)].$(pname)")
                end
            end
        end
    else
        for (cname, hv) in fitres.bestfit
            for (pname, _) in hv
                par.fixed  &&  continue
                push!(parnames, "[$(cname)].$(pname)")
            end
        end
    end
    @assert length(parnames)^2 == length(fitres.status.internal.covar)
    ii = Vector{Int}()
    jj = Vector{Int}()
    covar = Vector{Float64}()
    for i in 1:length(parnames)
        for j in i+1:length(parnames)
            push!(covar, fitres.status.internal.covar[i, j])
            push!(ii, i)
            push!(jj, j)
        end
    end
    if sort
        ii    = ii[   sortperm(abs.(covar))]
        jj    = jj[   sortperm(abs.(covar))]
        covar = covar[sortperm(abs.(covar))]
    end
    for i in 1:length(ii)
        if !isnothing(select)
            (parnames[ii[i]] in select)  ||  continue
        end
        (abs(covar[i]) < threshold)  &&  continue
        @printf "%-30s  %-30s  %10.4f\n" parnames[ii[i]] parnames[jj[i]] covar[i]
    end
end
=#

"""
    mock(::Type{Measures}, model::Model; keywords...)
    mock(::Type{Measures}, multi::MultiModel; keywords...)

Generate mock dataset(s) using a ground truth `Model` or `MultiModel` object. The first version returns a single `Measures` object, while the second returns a `Vector{Measures}`.

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
function mock(::Type{Measures}, model::Model; properr=0.01, rangeerr=0.05, abserr=0., seed=nothing)
    rng = MersenneTwister(seed);
    update!(model)
    values = model()
    ee = extrema(values)
    range = ee[2] - ee[1]
    @assert range > 0
    err = (properr .* abs.(values) .+ rangeerr .* range .+ abserr)
    values .+= err .* randn(rng, size(values))
    Measures(domain(model), values, err)
end

function mock(T, multi::Vector{Model}; kws...)
    update!(multi)
    return [mock(T, multi[i]; kws...) for i in 1:length(multi)]
end

