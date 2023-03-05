import Gnuplot
import Gnuplot.recipe

Gnuplot.recipe(data::Measures{1}) =
    Gnuplot.PlotElement(cmds=["set bars 0"],
                        data=Gnuplot.DatasetBin(coords(domain(data)), values(data), uncerts(data)),
                        plot="with yerr t 'Data' lc rgb 'gray'")

Gnuplot.recipe(model::Model) = Gnuplot.recipe(GModelFit.ModelSnapshot(model))

function Gnuplot.recipe(model::GModelFit.ModelSnapshot)
    @assert ndims(domain(model)) == 1
    out = Vector{Gnuplot.PlotElement}()
    for (k,v) in model.buffers
        (k == model.maincomp)  &&  continue
        #isa(v.comp, GModelFit.FComp)  ||  isa(v.comp, GModelFit.SumReducer)  ||  continue
        push!(out, Gnuplot.PlotElement(
            data=Gnuplot.DatasetBin(coords(domain(model)), model(k)),
            plot="with lines t '$(k)'"))
    end
    push!(out, Gnuplot.PlotElement(
        data=Gnuplot.DatasetBin(coords(domain(model)), model()),
        plot="with lines t 'Model' lc rgb 'black' lw 2"))
    return out
end
