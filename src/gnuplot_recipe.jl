import Gnuplot
import Gnuplot.recipe

Gnuplot.recipe(data::Measures{1}) =
    Gnuplot.parseSpecs("set bars 0",
                       coords(domain(data)), values(data), uncerts(data),
                       "with yerr t 'Data' lc rgb 'gray'")

Gnuplot.recipe(meval::ModelEval) = Gnuplot.recipe(GModelFit.ModelSnapshot(model))

function Gnuplot.recipe(model::GModelFit.ModelSnapshot;
                        keep=Symbol[], skip=Symbol[])
    @assert ndims(domain(model)) == 1
    out = Vector{Gnuplot.AbstractGPSpec}()
    for (k,v) in model.buffers
        (k == model.maincomp)  &&  continue
        (k in skip)  &&  continue
        if (length(keep) == 0)  ||  (k in keep)
            #isa(v.comp, GModelFit.FComp)  ||  isa(v.comp, GModelFit.SumReducer)  ||  continue
            append!(out, Gnuplot.parseSpecs(coords(domain(model)), model(k),
                                            "with lines t '$(k)'"))
        end
    end
    append!(out, Gnuplot.parseSpecs(coords(domain(model)), model(),
                                    "with lines t 'Model' lc rgb 'black' lw 2"))
    return out
end
