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


# ====================================================================
#= TODO
function print_param_covariance(fitres::FitSummary;
                                select=nothing, sort=false, threshold=0.)
    @assert isa(fitres.status.internal, CMPFit.Result) "Solver is not CMPFit"

    parnames = String[]
    if isa(fitres.bestfit, Vector{GModelFit.HashHashVector{GModelFit.Parameter}})
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
function mock(::Type{Measures}, meval::MEval; properr=0.01, rangeerr=0.05, abserr=0., seed=nothing)
    update_eval!(meval)
    out = Vector{Measures}()
    for i in 1:length(meval)
        values = last_eval(meval, i)
        ee = extrema(values)
        range = ee[2] - ee[1]
        @assert range > 0
        err = (properr .* abs.(values) .+ rangeerr .* range .+ abserr)
        values .+= err .* randn(MersenneTwister(seed), size(values))
        push!(out, Measures(meval.v[i].domain, values, err))
    end
    return out
end

mock(::Type{T}, model::Model, domain::AbstractDomain; kws...) where T =
    mock(T, MEval(model, domain); kws...)[1]

mock(::Type{T}, models::Vector{Model}, domains::Vector{<: AbstractDomain}; kws...) where T =
    mock(T, MEval(models, domains); kws...)
