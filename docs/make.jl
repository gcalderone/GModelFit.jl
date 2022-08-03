using Revise, Documenter, GFit, Gnuplot

makedocs(sitename="Gfit.jl",
         authors = "Giorgio Calderone",
         format = Documenter.HTML(prettyurls = false),  # uncomment for local use, comment for deployment
         modules=[GFit],
         pages = [
             "Home" => "index.md",
             "Basic concepts and data types" => "concepts.md",
             "Built-in components" => "builtincomp.md",
             "Viewers" => "viewers.md",
             "Examples" => "examples.md",
             "API" => "api.md"
         ])
Gnuplot.quitall()

if !(@isdefined is_compiled)
    is_compiled = true
    error("Re-run with compiled code!")
end
