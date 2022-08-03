using Revise, Documenter, GFit

makedocs(sitename="Gfit.jl",
         authors = "Giorgio Calderone",
         format = Documenter.HTML(prettyurls = false),  # uncomment for local use, comment for deployment
         modules=[GFit],
         pages = [
             "Home" => "index.md",
             "Basic concepts and data types" => "concepts.md",
             "Built-in components" => "builtincomp.md",
             "Parameter patching" => "patch.md",
             "Custom minimizers" => "minimizers.md",
             "Multi-dimensional fitting" => "multidim.md",
             "Multi-model fitting" => "multimodel.md",
             "Interactive fitting" => "interactive.md",
             "Viewers" => "viewers.md",
             "Examples" => "examples.md",
             "API" => "api.md"
         ])
