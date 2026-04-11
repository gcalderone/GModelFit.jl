import Gnuplot
import Gnuplot.recipe

Gnuplot.recipe(data::GaussianData{1}) =
    Gnuplot.parseSpecs("set bars 0",
                       coords(domain(data)), values(data), uncerts(data),
                       "with yerr t 'Data' lc rgb 'gray'")

function Gnuplot.recipe(model::AbstractModelEval; keep=Symbol[], skip=Symbol[])
    @assert ndims(domain(model)) == 1
    out = Vector{Gnuplot.AbstractGPSpec}()
    for (ckey, comp) in getcomps(model)
        cname = comp_name(model, ckey)
        (cname in skip)  &&  continue
        if (length(keep) == 0)  ||  (cname in keep)
            append!(out, Gnuplot.parseSpecs(coords(unfolded_domain(getresp(model[()]))), model[cname](), "with lines t '$(cname)'"))
        end
    end
    append!(out, Gnuplot.parseSpecs(coords(domain(model)), model(), "with lines t 'Model' lc rgb 'black' lw 2"))
    return out
end
