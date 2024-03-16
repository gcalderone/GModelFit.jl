using Gnuplot, Dates
Gnuplot.quitall()
mkpath("assets")
Gnuplot.options.term = "unknown"
empty!(Gnuplot.options.init)
push!( Gnuplot.options.init, linetypes(:Set1_5, lw=2.5, ps=1.5))
function saveas(file)
    Gnuplot.save(term="pngcairo size 550,350 fontscale 0.8", "assets/$(file).png")
    nothing
end

