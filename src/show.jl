mutable struct ShowSettings
    plain::Bool
    tableborders::TextTableBorders
    floatformat::String
    border::Crayon
    header::Crayon
    fixed::Crayon
    error::Crayon
    highlighted::Crayon
    section::Crayon
    ShowSettings() = new(false, text_table_borders__unicode_rounded, "%9.4g",
                         crayon"light_blue",
                         crayon"light_blue negative", crayon"dark_gray",
                         crayon"light_red", crayon"negative", crayon"green bold")
end

const showsettings = ShowSettings()

function printtable(io, table, header; formatters=[], hlines=:none, kws...)
    common = pairs((column_labels=header,
                    formatters=formatters,
                    alignment=:l,
                    column_label_alignment=:c,
                    limit_printing=false,
                    compact_printing=false,
                    fit_table_in_display_horizontally=false,
                    fit_table_in_display_vertically=false,
                    maximum_number_of_columns=-1,
                    maximum_number_of_rows=-1))
    if showsettings.plain
        kws = pairs((table_format=TextTableFormat(borders=text_table_borders__compact,
                                                  horizontal_lines_at_data_rows=hlines),
                     common...))
    else
        kws = pairs((table_format=TextTableFormat(borders=showsettings.tableborders,
                                                  horizontal_lines_at_data_rows=hlines),
                     style=TextTableStyle(table_border=showsettings.border,
                                          first_line_column_label=showsettings.header),
                     common..., kws...))
    end
    pretty_table(io, table; kws...)
end

function section(io, args...; newline=true)
    nl = newline  ?  "\n"  :  ""
    if showsettings.plain
        print(io, args..., nl)
    else
        print(io, showsettings.section, args..., nl, crayon"default")
    end
end


function show(io::IO, dom::AbstractDomain)
    section(io, string(typeof(dom)) * " (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    table = Matrix{Union{Int,Float64}}(undef, ndims(dom), 6)
    for i in 1:ndims(dom)
        if isa(dom, CartesianDomain)
            vv = axes(dom, i)
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
               formatters=[fmt__printf(showsettings.floatformat, 3:6)])
end


function show(io::IO, data::Measures)
    section(io, typeof(data), ": (length: ", (length(data)), ")")
    table = Matrix{Union{String,Float64}}(undef, 0, 7)

    names = fieldnames(typeof(data))
    error = Vector{Bool}()
    for vv in [values(data), uncerts(data)]
        nan = length(findall(isnan.(vv))) + length(findall(isinf.(vv)))
        push!(error, nan > 0)
        vv = vv[findall(isfinite.(vv))]
        (length(vv) == 0)  &&  (vv = [NaN])
        table = vcat(table, ["" minimum(vv) maximum(vv) mean(vv) median(vv) std(vv) (nan > 0  ?  string(nan)  :  "") ])
    end
    table[:, 1] .= ["values", "uncerts"]

    if showsettings.plain
        highlighters = nothing
    else
        highlighters = [TextHighlighter((data,i,j) -> error[i], showsettings.error)]
    end
    printtable(io, table, ["", "Min", "Max", "Mean", "Median", "Std. dev.", "Nan/Inf"],
               formatters=[fmt__printf(showsettings.floatformat, 2:6)],
               highlighters=highlighters)
end


function show(io::IO, par::Parameter)
    if par.fixed
        println(io, "Value : ", par.val, "   (fixed)")
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


function show(io::IO, red::FunctDesc)
    println(io, red.display)
end


function prepare_params_table(input::Union{Model, ModelSnapshot, AbstractComponent, ComponentSnapshot})
    table = Matrix{Union{String,Float64}}(undef, 0, 8)
    fixed = Vector{Bool}()
    warns = Vector{Bool}()

    for (key, param) in getparams(input)
        if isa(key, Tuple)
            cname = key[1]
            pname = key[2]
        else
            cname = Symbol("")
            pname = key
        end

        range =  (isfinite(param.low )  ?  strip(@sprintf("%7.2g", param.low ))  :  (param.low  > 0  ?  "+∞"  :  "-∞")) * " : "
        range *= (isfinite(param.high)  ?  strip(@sprintf("%7.2g", param.high))  :  (param.high > 0  ?  "+∞"  :  "-∞"))
        patch = ""
        reparam = ""
        isa(param.patch, Symbol)       &&  (patch = string(param.patch))
        isa(param.patch, FunctDesc)    &&  (patch = param.patch.display)
        isa(param.reparam, FunctDesc)  &&  (reparam = string(param.reparam))
        table = vcat(table,
                     permutedims([string(cname),
                                  string(pname) * (param.fixed  ?  " (fixed) "  :  ""),
                                  range, param.val, param.unc,
                                  ((patch == "")  &&  (reparam == "")  ?  ""  :  param.actual),
                                  patch, reparam]))
        push!(fixed, param.fixed)
        isnan(param.unc)     &&  (table[end, 5] = "")  # avoid displaying NaNs
        isnan(param.actual)  &&  (table[end, 6] = "")  # avoid displaying NaNs
        if !param.fixed  &&  isnan(param.unc)
            push!(warns, true)
        else
            push!(warns, false)
        end
    end
    return (table, fixed, warns)
end


function show(io::IO, input::Union{AbstractComponent, ComponentSnapshot})
    table, fixed, warns = prepare_params_table(input)
    if showsettings.plain
        highlighters = nothing
    else
        highlighters = [TextHighlighter((data,i,j) -> (fixed[i]), showsettings.fixed),
                        TextHighlighter((data,i,j) -> (warns[i]  &&  (j in (3:5))), showsettings.error)]
    end
    printtable(io, table[:, 2:end], ["Param.", "Range", "Value", "Uncert.", "Actual", "Patch", "Reparam"],
               formatters=[fmt__printf(showsettings.floatformat, 3:5)],
               highlighters=highlighters)
end


function show(io::IO, model::Union{Model, ModelSnapshot})
    function tree_prefix(ff::Vector{DependencyNode})
        BRANCH = showsettings.plain  ?  "+ "  :  "├─╴"
        BRCONT = showsettings.plain  ?  "| "  :  "│  "
        BREND  = showsettings.plain  ?  "+ "  :  "└─╴"

        function set_tree_prefix!(prefix, ff::Vector{DependencyNode}, node::DependencyNode)
            cnames = getfield.(ff, :cname)
            irow = findfirst(cnames .== node.cname)
            prefix[irow, node.level] = string(node.cname)
            if length(node.childs) > 0
                i1 = findfirst(cnames .== node.childs[1  ].cname)
                i2 = findfirst(cnames .== node.childs[end].cname)
                prefix[i1:i2, node.level] .= BRCONT

                for ichild in 1:length(node.childs)
                    child = node.childs[ichild]
                    i = findfirst(cnames .== child.cname)
                    if ichild < length(node.childs)
                        prefix[i, node.level] = BRANCH
                    else
                        prefix[i, node.level] = BREND
                    end

                    set_tree_prefix!(prefix, ff, child)
                end
            end
        end

        prefix = fill("", length(ff), maximum(getfield.(ff, :level)))
        set_tree_prefix!(prefix, ff, ff[1])
        prefix[prefix .== ""] .= " "^length(BRANCH)
        return [string(rstrip(join(prefix[i,:]))) for i in 1:length(ff)]
    end
    section(io, "Components:")
    (length(keys(model.comps)) == 0)  &&  (return nothing)

    tree = deptree(model)
    allcomps = flatten(tree)
    prefix = tree_prefix(allcomps)
    table = Matrix{Union{String,Int,Float64}}(undef, length(allcomps), 7)
    fixed = Vector{Bool}()
    for i in 1:length(allcomps)
        cname = allcomps[i].cname
        table[i, 1] = prefix[i]
        table[i, 2] = comptype(model, cname)
        if isa(model, ModelSnapshot)
            table[i, 3] = evalcounter(model, cname)
            result = model(cname)
            v = view(result, findall(isfinite.(result)))
            if length(v) > 0
                table[i, 4:6] .= [minimum(v), maximum(v), mean(v)]
            else
                table[i, 4:6] .= ["", "", ""]
            end
            table[i, 7] = count(isnan.(result)) + count(isinf.(result))
        end
        push!(fixed, isfrozen(model, cname))
    end

    if showsettings.plain
        highlighters = nothing
    else
        highlighters = [TextHighlighter((data,i,j) -> fixed[i], showsettings.fixed)]
    end
    if isa(model, ModelSnapshot)
        printtable(io, table, ["Comp.", "Type", "Eval. count", "Min", "Max", "Mean", "NaN/Inf"],
                   formatters=[fmt__printf(showsettings.floatformat, 4:6)],
                   highlighters=highlighters)
    else
        printtable(io, table[:, 1:2], ["Component", "Type"],
                   highlighters=highlighters)
    end
    println(io)

    # -----------------------
    section(io, "Parameters:")
    table, fixed, warns = prepare_params_table(model)

    hrule = Vector{Int}()
    prevcname = ""
    for i in 1:size(table)[1]
        cname = table[i, 1]
        if prevcname != cname
            (i == 1)  ||  push!(hrule, i-1)
            prevcname = cname
            isfrozen(model, Symbol(cname))  &&  (table[i, 1] *= " (frozen)")
        else
            table[i, 1] = ""
        end
        isfrozen(model, Symbol(cname))  &&  (fixed[i] = true)
    end

    if showsettings.plain
        highlighters = nothing
    else
        highlighters = [TextHighlighter((data,i,j) -> (fixed[i]), showsettings.fixed),
                        TextHighlighter((data,i,j) -> (warns[i]  &&  (j in (4:6))), showsettings.error)]
    end
    printtable(io, table, ["Comp.", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch", "Reparam"],
               hlines=hrule, formatters=[fmt__printf(showsettings.floatformat, 4:6)],
               highlighters=highlighters)
end


function show(io::IO, ms::ModelSet)
    for (mname, model) in ms
        println(io)
        section(io, join(fill("=", 30)) * "  Model $mname  " * join(fill("=", 30)))
        show(io, model)
    end
    println(io)
end


getmessage(status::SolverStatusOK) = crayon"green", "OK"
getmessage(status::SolverStatusWarn) = crayon"bold yellow", "WARN\n" * status.message
getmessage(status::SolverStatusError) = crayon"bold red", "ERROR\n" * status.message

function show(io::IO, status::AbstractSolverStatus)
    print(io, "status: ")
    color, ss = getmessage(status)
    if showsettings.plain
        print(io, @sprintf("%-8s", ss))
    else
        print(io, color, @sprintf("%-8s", ss), crayon"default")
    end
end


function show(io::IO, res::FitSummary)
    section(io, "Fit summary:", newline=false)
    print(io, @sprintf(" #data: %d, #free pars: %d, red. fit stat.: %.5g, ", res.ndata, res.nfree, res.fitstat))
    show(io, res.status)
    println(io)
end
