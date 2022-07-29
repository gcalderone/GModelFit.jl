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



show(io::IO, mime::MIME"text/plain", domain::AbstractDomain) = show(io, domain)
function show(io::IO, dom::AbstractDomain)
    section(io, string(typeof(dom)) * " (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    hrule = Vector{Int}()
    push!(hrule, 0, 1, ndims(dom)+1)
    table = Matrix{Union{Int,Float64}}(undef, ndims(dom), 6)
    for i in 1:ndims(dom)
        a = dom[i]
        b = 0
        if length(a) > 1
            b = a .- circshift(a, 1)
            b = b[2:end]
        end
        table[i, 1] = i
        table[i, 2] = length(a)
        table[i, 3:6] = [minimum(a), maximum(a), minimum(b), maximum(b)]
    end
    if isa(dom, Domain)  &&  (ndims(dom) >= 2) # treat Domain{1} as a Cartesian one, despite it is a linear one.
        printtable(io, table[:, 1:4], ["Dim", "Size", "Min val", "Max val"],
                   hlines=hrule, formatters=ft_printf(showsettings.floatformat, 3:4))
    else
        printtable(io, table, ["Dim", "Size", "Min val", "Max val", "Min step", "Max step"],
                   hlines=hrule, formatters=ft_printf(showsettings.floatformat, 3:6))
    end
end


function show(io::IO, data::Measures)
    section(io, typeof(data), ": (length: ", (length(data.val)), ")")
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    hrule = Vector{Int}()
    push!(hrule, 0, 1)

    names = fieldnames(typeof(data))
    error = Vector{Bool}()
    for name in [:val, :unc]
        a = getfield(data, name)
        nan = length(findall(isnan.(a))) + length(findall(isinf.(a)))
        a = a[findall(isfinite.(a))]
        push!(error, nan > 0)
        table = vcat(table, [string(name) minimum(a) maximum(a) mean(a) median(a) std(a) (nan > 0  ?  string(nan)  :  "") ])
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
        if par.val == par.pval
            println(io, "Value : ", par.val, " ± ", par.unc,  "  [", par.low , " : ", par.high, "] ")
        else
            println(io, "Value : ", par.val, " ± ", par.unc,  "  [", par.low , " : ", par.high, "], patched: " , par.pval)
        end
    end
end


function preparetable(comp::AbstractComponent; cname::String="?", cfixed=false)
    table = Matrix{Union{String,Float64}}(undef, 0, 8)
    fixed = Vector{Bool}()
    watch = Vector{Bool}()

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
                     [cname * (cfixed  ?  " (FIXED)"  :  "") ctype parname range param.val param.unc param.pval patch])
        push!(fixed, param.fixed)
        push!(watch, !isnothing(param.patch))
        if !showsettings.plain
            cname = ""  # delete from following lines within the same component box
            ctype = ""
            cfixed = false
        end
    end
    return (table, fixed, watch)
end


function show(io::IO, comp::AbstractComponent)
    (table, fixed, watch) = preparetable(comp)
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               formatters=ft_printf(showsettings.floatformat, 5:7),
               highlighters=(Highlighter((data,i,j) -> (!(watch[i])  &&  (j in (7,8))), showsettings.fixed),
                             Highlighter((data,i,j) -> (fixed[i]  &&  (j in (3,4,5,6))), showsettings.fixed)))
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
               hlines=[0,1, length(model.cevals)+1,  length(model.cevals)+1],
               formatters=ft_printf(showsettings.floatformat, 4:6),
               highlighters=Highlighter((data,i,j) -> fixed[i], showsettings.fixed))
    println(io)
    section(io, "Parameters:")
    (length(model.cevals) == 0)  &&  (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 8)
    fixed = Vector{Bool}()
    watch = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, ceval) in model.cevals
        comp = ceval.comp
        (t, f, w) = preparetable(comp, cname=string(cname), cfixed=ceval.cfixed)
        table = vcat(table, t)
        append!(fixed, f .| ceval.cfixed)
        append!(watch, w)
        push!(hrule, length(watch)+1)
    end
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 5:7),
               highlighters=(Highlighter((data,i,j) -> (!watch[i]  &&  (j in (7,8))), showsettings.fixed),
                             Highlighter((data,i,j) -> (fixed[i]   &&  (j in (3,4,5,6))), showsettings.fixed)))
end


function show(io::IO, multi::MultiModel)
    for id in 1:length(multi.models)
        section(io, "\n=====================================================================")
        section(io, "Model $id:")
        show(io, multi.models[id])
    end
    println(io)
end


function show(io::IO, res::FitResult)
    section(io, "Fit results:")

    println(io, @sprintf("    #Data  : %8d              Fit-stat: %-10.5g", res.ndata, res.fitstat))
    println(io, @sprintf("    #Param : %8d              Red. GOF: %-10.4g", res.nfree, res.gofstat / res.dof))
    println(io, @sprintf("    DOF    : %8d              Prob.   : %-10.4g", res.dof, 10^res.log10testprob))

    print(io,            "    Status : ")
    (crayon, status, message) = as_string(res.status)
    if showsettings.plain
        print(io, @sprintf("%8s", status))
    else
        print(io, crayon, @sprintf("%8s", status), crayon"default")
    end

    println(io, @sprintf("              Elapsed : %-10.4g s", res.elapsed))

    if message != ""
        println(io, crayon, message, crayon"default")
    end
end
