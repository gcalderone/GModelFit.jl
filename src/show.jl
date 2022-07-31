mutable struct ShowSettings
    plain::Bool
    tableformat::TextFormat
    floatformat::String
    showfixed::Bool
    border::Crayon
    header::Crayon
    subheader::Crayon
    fixed::Crayon
    error::Crayon
    highlighted::Crayon
    section::Crayon
    ShowSettings() = new(false, tf_unicode_rounded, "%9.4g", true,
                         crayon"light_blue", crayon"light_blue negative bold",
                         crayon"dark_gray bold", crayon"dark_gray",
                         crayon"light_red blink", crayon"negative", crayon"green bold")
end

const showsettings = ShowSettings()

function printtable(args...; formatters=(), hlines=:none, kw...)
    if showsettings.plain
        c = crayon"default"
        pretty_table(args...; formatters=formatters, alignment=:l, crop=:none, tf=tf_compact, hlines=hlines,
                     border_crayon=c, header_crayon=c, subheader_crayon=c,
                     highlighters=())
    else
        pretty_table(args...; formatters=formatters, alignment=:l, crop=:none, tf=showsettings.tableformat, hlines=hlines,
                     border_crayon=showsettings.border, header_crayon=showsettings.header, subheader_crayon=showsettings.subheader,
                     kw...)
    end
end

function section(io, args...)
    if showsettings.plain
        println(io, args...,)
    else
        println(io, showsettings.section, args..., crayon"default")
    end
end

function left(s::String, maxlen::Int)
    (length(s) <= maxlen)  &&  (return s)
    return s[1:maxlen]
end



function show(io::IO, dom::AbstractDomain)
    section(io, string(typeof(dom)) * " (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    hrule = Vector{Int}()
    push!(hrule, 0, 1, ndims(dom)+1)
    table = Matrix{Union{Int,Float64}}(undef, ndims(dom), 6)
    for i in 1:ndims(dom)
        if isa(dom, CartesianDomain)
            vv = axis(dom, i)
        else
            vv = coords(dom, i)
        end
        steps = 0
        if length(vv) > 1
            steps = vv .- circshift(vv, 1)
            steps = steps[2:end]
        end
        table[i, 1] = i
        table[i, 2] = length(vv)
        table[i, 3:6] = [minimum(vv), maximum(vv), minimum(steps), maximum(steps)]
    end
    printtable(io, table, ["Dim", "Length", "Min val", "Max val", "Min step", "Max step"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 3:6))
end


function show(io::IO, data::AbstractMeasures)
    section(io, typeof(data), ": (length: ", (length(data)), ")")
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    hrule = Vector{Int}()
    push!(hrule, 0, 1)

    names = fieldnames(typeof(data))
    error = Vector{Bool}()
    for i in 1:length(data.labels)
        vv = values(data, i)
        nan = length(findall(isnan.(vv))) + length(findall(isinf.(vv)))
        vv = vv[findall(isfinite.(vv))]
        push!(error, nan > 0)
        table = vcat(table, [data.labels[i] minimum(vv) maximum(vv) mean(vv) median(vv) std(vv) (nan > 0  ?  string(nan)  :  "") ])
    end
    push!(hrule, 0, size(table)[1]+1)
    printtable(io, table, ["", "Min", "Max", "Mean", "Median", "Std. dev.", "Nan/Inf"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 2:6),
               highlighters=(Highlighter((data,i,j) -> error[i], showsettings.error)))
end


function show(io::IO, par::Parameter)
    if par.fixed
        println(io, "Value : ", par.val, "   (FIXED)")
    elseif isnan(par.unc)
        println(io, "Value : ", par.val, "  [", par.low , " : ", par.high, "]")
    else
        if par.val == par.actual
            println(io, "Value : ", par.val, " ± ", par.unc,  "  [", par.low , " : ", par.high, "] ")
        else
            println(io, "Value : ", par.val, " ± ", par.unc,  "  [", par.low , " : ", par.high, "], actual: " , par.actual)
        end
    end
end


function preparetable(comp::AbstractComponent; cname::String="?", cfixed=false)
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    fixed = Vector{Bool}()

    ctype = split(string(typeof(comp)), ".")
    (ctype[1] == "GFit")  &&   (ctype = ctype[2:end])
    ctype = join(ctype, ".")

    for (pname, param) in getparams(comp)
        parname = string(pname)
        parname *= (param.fixed  ?  " (FIXED)"  :  "")
        (!showsettings.showfixed)  &&  param.fixed  &&  continue
        range = strip(@sprintf("%7.2g:%-7.2g", param.low, param.high))
        (range == "-Inf:Inf")  &&  (range = "")
        patch = ""
        isa(param.patch, Symbol)  &&  (patch = string(param.patch))
        isa(param.patch, λFunct)  &&  (patch = param.patch.display)
        isa(param.mpatch,λFunct)  &&  (patch = param.mpatch.display)
        table = vcat(table,
                     [cname * (cfixed  ?  " (FIXED)"  :  "") ctype parname range param.val (patch == ""  ?  ""  :  param.actual) patch])
        push!(fixed, param.fixed)
        if !showsettings.plain
            cname = ""  # delete from following lines within the same component box
            ctype = ""
            cfixed = false
        end
    end
    return (table, fixed)
end


function show(io::IO, comp::AbstractComponent)
    (table, fixed) = preparetable(comp)
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Actual", "Patch"],
               formatters=ft_printf(showsettings.floatformat, 5:7),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]  &&  (j in (3,4,5))), showsettings.fixed)))
end


function show(io::IO, red::λFunct)
    println(io, red.display)
end


function tabledeps(model::Model)
    out0 = tabledeps(model, find_maincomp(model), 0)
    levels = getindex.(out0, 1)
    out = Vector{Tuple}()
    for i in 1:length(out0)
        isamelevel = findfirst(levels[i+1:end] .== levels[i])
        isuper     = findfirst(levels[i+1:end] .<  levels[i])
        final = isnothing(isamelevel)  ||  (!isnothing(isuper)  &&  (isuper < isamelevel))
        push!(out, (final, out0[i]...))
    end

    table = Matrix{Union{String,Int,Float64}}(undef, length(out), 7)
    fixed = Vector{Bool}()
    for i in 1:length(out)
        final = out[i][1]
        level = out[i][2]
        prefix = ""
        (level > 1)  &&  (prefix  = join(fill("  │", level-1)))
        (level > 0)  &&  (prefix *= (final  ?  "  └─ "  :  "  ├─ "))
        table[i, 1] = prefix * out[i][3]
        for j in 2:7
            table[i, j] = out[i][j+2]
        end
        push!(fixed, model.cevals[Symbol(out[i][3])].cfixed)
    end
    return table, fixed
end

function tabledeps(model::Model, cname::Symbol, level::Int)
    out = Vector{Tuple}()
    result = model.cevals[cname].buffer
    v = view(result, findall(isfinite.(result)))
    push!(out, (level, string(cname),
                string(typeof(model[cname])),
                model.cevals[cname].counter,
                minimum(v), maximum(v), mean(v),
                count(isnan.(result)) + count(isinf.(result))))
    deps = dependencies(model[cname])
    for i in 1:length(deps)
        append!(out, tabledeps(model, deps[i], level+1))
    end
    return out
end


function show(io::IO, model::Model)
    println(io)
    section(io, "Components:")
    table, fixed = tabledeps(model)
    printtable(io, table, ["Component", "Type", "Eval. count", "Min", "Max", "Mean", "NaN/Inf"],
               hlines=[0,1, size(table)[1]+1],
               formatters=ft_printf(showsettings.floatformat, 4:6),
               highlighters=Highlighter((data,i,j) -> fixed[i], showsettings.fixed))
    println(io)
    section(io, "Parameters:")
    (length(model.cevals) == 0)  &&  (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    fixed = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, ceval) in model.cevals
        comp = ceval.comp
        (t, f) = preparetable(comp, cname=string(cname), cfixed=ceval.cfixed)
        table = vcat(table, t)
        append!(fixed, f .| ceval.cfixed)
        push!(hrule, length(fixed)+1)
    end
    table = table[:, [1; 3:end]]  # drop the Type column
    printtable(io, table, ["Component", "Param.", "Range", "Value.", "Actual", "Patch"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 4:6),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]   &&  (j in (2,3,4,5))), showsettings.fixed)))
end


function show(io::IO, multi::MultiModel)
    for id in 1:length(multi.models)
        section(io, "\n=====================================================================")
        section(io, "Model $id:")
        show(io, multi.models[id])
    end
    println(io)
end


function show(io::IO, bestfit::Vector{HashHashVector{Parameter}})
    for id in 1:length(bestfit)
        section(io, "\n=====================================================================")
        section(io, "Model $id:")
        show(io, bestfit[id])
    end
    println(io)
end

function show(io::IO, bestfit::HashHashVector{Parameter})
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    fixed = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, hv) in bestfit
        scname = string(cname)
        for (pname, param) in hv
            parname = string(pname)
            parname *= (param.fixed  ?  " (FIXED)"  :  "")
            (!showsettings.showfixed)  &&  param.fixed  &&  continue
            range = strip(@sprintf("%7.2g:%-7.2g", param.low, param.high))
            (range == "-Inf:Inf")  &&  (range = "")
            patch = ""
            isa(param.patch, Symbol)  &&  (patch = string(param.patch))
            isa(param.patch, λFunct)  &&  (patch = param.patch.display)
            isa(param.mpatch,λFunct)  &&  (patch = param.mpatch.display)
            table = vcat(table,
                         [scname parname range param.val param.unc (patch == ""  ?  ""  :  param.actual) patch])
            push!(fixed, param.fixed)
            if !showsettings.plain
                scname = ""  # delete from following lines within the same component box
            end
        end
        push!(hrule, length(fixed)+1)
    end
    printtable(io, table, ["Component", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 4:7),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]   &&  (j in (2,3,4,5))), showsettings.fixed)))
end

function show(io::IO, res::FitResult)
    section(io, "Best fit parameters:")
    show(res.bestfit)

    section(io, "Fit results:")

    println(io, @sprintf("    #Data : %8d              #Free params  : %10d"    , res.ndata, res.nfree))
    println(io, @sprintf("    DOF   : %8d              Red. fit stat.: %10.5g", res.dof, res.fitstat / res.dof))
    print(io,            "    Status: ")
    (crayon, status, message) = as_string(res.status)
    if showsettings.plain
        print(io, @sprintf("%8s", status))
    else
        print(io, crayon, @sprintf("%8s", status), crayon"default")
    end
    println(io, @sprintf("              Elapsed time  : %10.4g s", res.elapsed))

    if message != ""
        println(io, crayon, message, crayon"default")
    end
end
