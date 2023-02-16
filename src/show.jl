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

function printtable(io, table, header, args...; formatters=(), hlines=:none, kw...)
    if showsettings.plain
        c = crayon"default"
        pretty_table(io, table; header=header, formatters=formatters, alignment=:l, crop=:none, tf=tf_compact, hlines=hlines,
                     border_crayon=c, header_crayon=c, subheader_crayon=c,
                     highlighters=())
    else
        pretty_table(io, table; header=header, formatters=formatters, alignment=:l, crop=:none, tf=showsettings.tableformat, hlines=hlines,
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
        vv = data.values[i]
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
    table = Matrix{Union{String,Float64}}(undef, 0, 8)
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
        isa(param.patch, FunctDesc)  &&  (patch = param.patch.display)
        isa(param.mpatch,FunctDesc)  &&  (patch = param.mpatch.display)
        table = vcat(table,
                     [cname * (cfixed  ?  " (FIXED)"  :  "") ctype parname range param.val (isnan(param.unc)  ?  ""  :  param.unc) (patch == ""  ?  ""  :  param.actual) patch])
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
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               formatters=ft_printf(showsettings.floatformat, 5:7),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]  &&  (j in (3,4,5))), showsettings.fixed)))
end


function show(io::IO, red::FunctDesc)
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
        (level > 1)  &&  (prefix  = join(fill(" │", level-1)))
        (level > 0)  &&  (prefix *= (final  ?  " └─ "  :  " ├─ "))
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
    deps = dependencies(model, cname)
    for i in 1:length(deps)
        append!(out, tabledeps(model, deps[i], level+1))
    end
    return out
end


function show(io::IO, model::Model)
    section(io, "Components:")
    table, fixed = tabledeps(model)
    printtable(io, table, ["Component", "Type", "Eval. count", "Min", "Max", "Mean", "NaN/Inf"],
               hlines=[0,1, size(table)[1]+1],
               formatters=ft_printf(showsettings.floatformat, 4:6),
               highlighters=Highlighter((data,i,j) -> fixed[i], showsettings.fixed))
    println(io)
    section(io, "Parameters:")
    (length(model.cevals) == 0)  &&  (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 8)
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
    printtable(io, table, ["Component", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 4:6),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]), showsettings.fixed)))
end


function show(io::IO, multi::MultiModel)
    for id in 1:length(multi.models)
        println(io)
        section(io, join(fill("=", 30)) * "  Model $id  " * join(fill("=", 30)))
        show(io, multi.models[id])
    end
    println(io)
end


function show(io::IO, status::MinimizerStatus)
    print(io, "Status        : ")
    if status.code == MinOK
        ss = (crayon"green", "OK")
    elseif status.code == MinWARN
        ss = (crayon"bold yellow", "WARNING")
    elseif status.code == MinERROR
        ss = (crayon"bold red", "ERROR")
    elseif status.code == MinDRY
        ss = (crayon"bold red", "DRY")
    else
        error("Unsupported minimizer status code: $(status.code)")
    end

    if showsettings.plain
        println(io, @sprintf("%10s", ss[2]))
    else
        println(io, ss[1], @sprintf("%10s", ss[2]), crayon"default")
    end
    if status.message != ""
        println(io, ss[1], ": ", status.message, crayon"default")
    end
end


function show(io::IO, res::FitResult)
    section(io, "Fit results:")
    println(io, @sprintf("    #Data : %8d       Elapsed time  : %10.5g", res.ndata, res.elapsed))
    println(io, @sprintf("    #Free : %8d       Red. fit stat.: %10.5g", res.nfree, res.fitstat))
    print(  io, @sprintf("    DOF   : %8d       ", res.dof))
    show(io, res.status)
end


function show(io::IO, mb::ModelSnapshot)
    println(io, mb.show)
end
