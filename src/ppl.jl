using Turing, Distributions, Statistics, HypothesisTests
using Gnuplot

struct PPLData{T <: Likelihood}
    params::OrderedDict{Symbol, Parameter}
    priors::Vector{Distribution{Univariate, Continuous}}
    lh::Vector{T}

    PPLData(model::Model, data::AbstractMeasures) = PPLData([model], [data])
    function PPLData(model::Vector{Model}, data::Vector{T}) where T <: AbstractMeasures
        @assert length(model) == length(data)
        PPLData(Likelihood(data, model, use_AD=true))
    end

    function PPLData(lh::T) where T <: Likelihood
        params = OrderedDict{Symbol, Parameter}()
        priors = Vector{Distribution{Univariate, Continuous}}()
        for (t, par) in getparams(lh; free=true)
            @assert isfinite(par.val)
            @assert isfinite(par.unc)
            push!(priors, Uniform(max(par.val - 5 * par.unc, par.low),
                                  min(par.val + 5 * par.unc, par.high)))
            params[Symbol(join(string.(t), "_"))] = par
        end
        new{T}(params, priors, [deepcopy(lh) for i in 1:(Threads.nthreads() + 1)])
    end
end
getlh(ppldata::PPLData) = ppldata.lh[Threads.threadid()]

@model function pplmodel(ppldata::PPLData{Likelihood{T, ChiSquared}}) where T <: Measures
    pvalues ~ arraydist(ppldata.priors)
    lh = getlh(ppldata)
    update_eval!(lh, pvalues)
    Turing.@addlogprob! loglikelihood(lh)
end


renameparams(ppldata::PPLData, chains::Chains) = 
    replacenames(chains, Dict(Pair.(Symbol.(:pvalues, "[", 1:length(ppldata.priors), "]"), collect(keys(ppldata.params)))))


function cornerplot(ppldata::PPLData, chains::Chains, fsumm::Union{Solvers.FitSummary, Nothing}=nothing; chain=:all)
    function correlation(x, y)
        test = CorrelationTest(x, y)
        p_val = pvalue(test)
        m = cov(x, y) / var(x)
        c = mean(y) - m * mean(x)
        line_x = [minimum(x), maximum(x)]
        line_y = m .* line_x .+ c
        return p_val, line_x, line_y
    end
    pp = get_params(chains)
    nchains = size(pp.logjoint)[2]
    if chain == :all
        chain = 1:nchains
    else
        @assert (1 <= chain <= nchains)
    end
    npar = length(ppldata.params)
    pnames = collect(keys(ppldata.params))
    @gp    "set size ratio 1" Gnuplot.palette(:viridis) :-
    @gp :- "set tmargin 0" :-
    @gp :- "set bmargin 0" :-
    @gp :- "set lmargin 1" :-
    @gp :- "set rmargin 1" :-
    @gp :- "set multiplot layout $(npar+1), $npar margins 0.1, 0.99, 0.15, 0.95 spacing 0.05,0.01" :-
    @gp :- "unset colorbox" :-
    @gp :- "set xlabel noenhanced" "set ylabel noenhanced" "set x2label noenhanced" "set y2label noenhanced rotate by -90" :-
    @gp :- "unset xtics" xlab="" :-
    count = 0
    for iy in 2:npar
        for ix in 1:npar
            count += 1
            if ix < iy
                parx = ppldata.params[pnames[ix]]
                pary = ppldata.params[pnames[iy]]
                x = collect(reshape(getproperty(pp, pnames[ix])[:, chain], :))
                y = collect(reshape(getproperty(pp, pnames[iy])[:, chain], :))
                xr = [min(parx.val, minimum(x)), max(parx.val, maximum(x))]
                yr = [min(pary.val, minimum(y)), max(pary.val, maximum(y))]
                hh = hist(x, y)
                @gp :- count hist_bins(hh, 1) hist_bins(hh, 2) hist_weights(hh) "w image notit" xr=xr yr=yr ylog=false :-

                p_val, line_x, line_y = correlation(x, y);  @gp :- line_x line_y "w l notit dt 2 lc rgb 'white'" :-
                p_val, line_y, line_x = correlation(y, x);  @gp :- line_x line_y "w l notit dt 2 lc rgb 'white'" :-
                @gp :- count parx.val .* [1,1] pary.val .* [1,1] "w p notit pt 6 ps 3 lc rgb 'red'" :-
            else
                @gp :- count "set multiplot next" :-
            end
            if (ix < iy)  &&  (ix+1 == iy)
                @gp :- count "set x2label '$(pnames[ix])'" "set y2label '$(pnames[iy])'" :-
            else
                @gp :- count "set x2label ''" "set y2label ''" :-
            end
            @gp :- count ylab=(ix == 1  ?  string(pnames[iy])  :  "") :-
        end
    end

    for ix in 1:npar
        count += 1 
        x = collect(reshape(getproperty(pp, pnames[ix])[:, chain], :))
        hh = hist(x)
        yr = [extrema(hist_weights(hh))...]
        @gp :- count hist_bins(hh) hist_weights(hh) "w steps notit lc rgb 'blue'" xr=[extrema(x)...] yr=yr ylog=false :-
        @gp :- count quantile(x, (1 - 0.9) / 2) .* [1, 1] yr "w l notit dt 3 lc rgb 'black'" :-
        @gp :- count quantile(x, (1 + 0.9) / 2) .* [1, 1] yr "w l notit dt 3 lc rgb 'black'" :-
        x = range(minimum(x), maximum(x), length=100)
        parx = ppldata.params[pnames[ix]]
        if isfinite(parx.unc)
            @gp :- count x maximum(hist_weights(hh)) * exp.(-0.5 .* ((x .- parx.val) ./ parx.unc).^2) "w l notit lc rgb 'red'" :-
        end
        @gp :- count ylab=(ix == 1  ?  "Counts"  :  "") :-
    end

    for ix in 1:npar
        count += 1
        x = collect(reshape(getproperty(pp, pnames[ix])[    :, chain], :))
        y = collect(reshape(getproperty(pp, :loglikelihood)[:, chain], :) ./ (-0.5 * dof(getlh(ppldata))))
        hh = hist(x, y)
        yr = isnothing(fsumm)  ?  [extrema(y)...]  :  [min(fsumm.fitstat, minimum(y)), maximum(y)]
        @gp :- count x y "w dots notit lc rgb 'blue'" xr=[extrema(x)...] yr=yr ylog=true :-
        @gp :- count quantile(x, (1 - 0.9) / 2) .* [1, 1] yr "w l notit dt 3 lc rgb 'black'" :-
        @gp :- count quantile(x, (1 + 0.9) / 2) .* [1, 1] yr "w l notit dt 3 lc rgb 'black'" :-
        parx = ppldata.params[pnames[ix]]
        if !isnothing(fsumm)  &&  isfinite(parx.unc)
            @gp :- count parx.val .* [1,1] fsumm.fitstat .* [1,1] "w p notit pt 6 ps 3 lc rgb 'red'" :-
            x = range(minimum(x), maximum(x), length=100)
            @gp :- count x (((x .- parx.val) ./ parx.unc).^2) ./ dof(getlh(ppldata)) .+ fsumm.fitstat "w l notit lc rgb 'red'" :-
        end
        @gp :- count "set xtics" xlab=string(pnames[ix]) :-
        @gp :- count ylab=(ix == 1  ?  "Fit stat."  :  "") :-
    end
    @gp
end
