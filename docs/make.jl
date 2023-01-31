using Revise, Documenter, GFit, Gnuplot

makedocs(sitename="Gfit.jl",
         authors = "Giorgio Calderone",
         format = Documenter.HTML(prettyurls = false),  # uncomment for local use, comment for deployment
         modules=[GFit],
         pages = [
             "Home" => "index.md",
             "Basic concepts and data types" => "concepts.md",
             "Built-in components" => "builtincomp.md",
             "Custom components" => "customcomp.md",
             "Parameter constraints" => "parameter.md",
             "Multi-dataset fitting" => "multifit.md",
             "Viewers" => "viewers.md",
             "Minimizers" => "minimizers.md",
             "Miscellaneous" => "misc.md",
             "API" => "api.md"
         ])
Gnuplot.quitall()

if !(@isdefined is_compiled)
    is_compiled = true
    error("Re-run with compiled code!")
end
