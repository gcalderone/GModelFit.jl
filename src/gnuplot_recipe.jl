import Gnuplot
import Gnuplot.recipe

Gnuplot.recipe(data::Measures{1}) =
    Gnuplot.PlotElement(cmds=["set bars 0"],
                        data=Gnuplot.DatasetBin(coords(domain(data)), values(data), uncerts(data)),
                        plot="with yerr t 'Data' lc rgb 'gray'")

Gnuplot.recipe(model::Model) = Gnuplot.recipe(GFit.ModelSnapshot(model))

function Gnuplot.recipe(model::GFit.ModelSnapshot)
    @assert ndims(domain(model)) == 1
    out = Vector{Gnuplot.PlotElement}()
    for (k,v) in model.buffers
        (k == model.maincomp)  &&  continue
        #isa(v.comp, GFit.FComp)  ||  isa(v.comp, GFit.SumReducer)  ||  continue
        push!(out, Gnuplot.PlotElement(
            data=Gnuplot.DatasetBin(coords(domain(model)), model(k)),
            plot="with lines t '$(k)'"))
    end
    push!(out, Gnuplot.PlotElement(
        data=Gnuplot.DatasetBin(coords(domain(model)), model()),
        plot="with lines t 'Model' lc rgb 'black' lw 2"))
    return out
end
