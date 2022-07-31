# ====================================================================
function print_param_covariance(fitres::FitResult;
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


function mockdata(model::Model; propnoise=0.01, rangenoise=0.05, absnoise=0., seed=nothing)
    rng = MersenneTwister(seed);
    evaluate(model)
    values = model()
    ee = extrema(values)
    range = ee[2] - ee[1]
    @assert range > 0
    noise = (propnoise .* values .+ rangenoise .* range .+ absnoise)
    Measures(domain(model), values .+ noise .* randn(rng, length(values)), noise)
end
