mutable struct ShowSettings
    plain::Bool
    tableborders::TextTableBorders
    floatformat::String
    showfixed::Bool
    border::Crayon
    header::Crayon
    fixed::Crayon
    error::Crayon
    highlighted::Crayon
    section::Crayon
    ShowSettings() = new(false, text_table_borders__unicode_rounded, "%9.4g", true,
                         crayon"light_blue",
                         crayon"light_blue negative bold", crayon"dark_gray",
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
        vv = vv[findall(isfinite.(vv))]
        push!(error, nan > 0)
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


function preparetable(comp::Union{AbstractComponent, GModelFit.ComponentSnapshot};
                      cname::String="?", ctype="?", cfixed=false)
    table = Matrix{Union{String,Float64}}(undef, 0, 8)
    fixed = Vector{Bool}()
    warns = Vector{Bool}()

    for (pname, param) in getparams(comp)
        (!showsettings.showfixed)  &&  param.fixed  &&  continue
        range = strip(@sprintf("%7.2g:%-7.2g", param.low, param.high))
        # (range == "-Inf:Inf")  &&  (range = "")
        patch = ""
        isa(param.patch, Symbol)     &&  (patch = string(param.patch))
        isa(param.patch, FunctDesc)  &&  (patch = param.patch.display)
        isa(param.mpatch,FunctDesc)  &&  (patch = param.mpatch.display)
        table = vcat(table,
                     permutedims([cname * (cfixed  ?  " (fixed)"  :  ""), ctype,
                                  string(pname), range, param.val,
                                  (param.fixed | cfixed  ?  " (fixed)"  :  param.unc),
                                  (patch == ""  ?  ""  :  param.actual), patch]))
        push!(fixed, param.fixed)
        if !param.fixed  &&  isnan(param.unc)
            push!(warns, true)
            table[end,6] = ""
        else
            push!(warns, false)
        end
        if !showsettings.plain
            cname = ""  # delete from following lines within the same component box
            ctype = ""
            cfixed = false
        end
    end
    return (table, fixed, warns)
end


function show(io::IO, comp::Union{AbstractComponent, GModelFit.ComponentSnapshot})
    ctype = isa(comp, AbstractComponent)  ?  string(typeof(comp))  :  "?"
    (table, fixed, warns) = preparetable(comp, ctype=ctype)
    if showsettings.plain
        highlighters = nothing
    else
        highlighters = [TextHighlighter((data,i,j) -> (fixed[i]  &&  (j in (3,4,5,6))), showsettings.fixed),
                        TextHighlighter((data,i,j) -> (warns[i]  &&  (j in (3,4,5,6))), showsettings.error)]
    end
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               formatters=[fmt__printf(showsettings.floatformat, 5:7)],
               highlighters=highlighters)
end


function show(io::IO, red::FunctDesc)
    println(io, red.display)
end


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

function tabledeps(model::Union{Model, ModelSnapshot})
    tree = deptree(model)
    allcomps = flatten(tree)
    prefix = tree_prefix(allcomps)

    table = Matrix{Union{String,Int,Float64}}(undef, length(allcomps), 8)
    fixed = Vector{Bool}()
    for i in 1:length(allcomps)
        cname = allcomps[i].cname
        table[i, 1] = prefix[i]
        table[i, 2] = comptype(model, cname)
        table[i, 3] = count(getproperty.(values(getparams(model[cname])), :fixed) .== false)
        (table[i, 3] == 0)  &&  (table[i, 3] = "")
        if !isa(model, Model)
            table[i, 4] = evalcounter(model, cname)
            result = model(cname)
            v = view(result, findall(isfinite.(result)))
            if length(v) > 0
                table[i, 5:7] .= [minimum(v), maximum(v), mean(v)]
            else
                table[i, 5:7] .= ["", "", ""]
            end
            table[i, 8] = count(isnan.(result)) + count(isinf.(result))
        end
        push!(fixed, isfreezed(model, cname))
    end
    return table, fixed
end


function show(io::IO, model::Union{Model, ModelSnapshot})
    section(io, "Components:")
    (length(keys(model)) == 0)  &&  (return nothing)
    table, fixed = tabledeps(model)
    if showsettings.plain
        highlighters = nothing
    else
        highlighters = [TextHighlighter((data,i,j) -> fixed[i], showsettings.fixed)]
    end
    if !isa(model, Model)
        printtable(io, table, ["Component", "Type", "#Free", "Eval. count", "Min", "Max", "Mean", "NaN/Inf"],
                   formatters=[fmt__printf(showsettings.floatformat, 5:7)],
                   highlighters=highlighters)
    else
        printtable(io, table[:, 1:3], ["Component", "Type", "#Free"],
                   formatters=[fmt__printf(showsettings.floatformat, 5:7)],
                   highlighters=highlighters)
    end
    println(io)

    section(io, "Parameters:")
    (length(keys(model)) == 0)  &&  (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 8)
    fixed = Vector{Bool}()
    warns = Vector{Bool}()
    hrule = Vector{Int}()
    for cname in keys(model)
        comp = model[cname]
        (t, f, w) = preparetable(comp,
                                 cname=string(cname),
                                 ctype=comptype(model, cname),
                                 cfixed=isfreezed(model, cname))
        table = vcat(table, t)
        append!(fixed, f .| isfreezed(model, cname))
        append!(warns, w)
        (length(fixed) > 0)  &&  push!(hrule, length(fixed))
    end

    if showsettings.plain
        highlighters = nothing
    else
        if !isa(model, Model)
            highlighters = [TextHighlighter((data,i,j) -> (fixed[i]), showsettings.fixed),
                            TextHighlighter((data,i,j) -> (warns[i]  &&  (j in (3,4,5,6))), showsettings.error)]
        else
            highlighters = [TextHighlighter((data,i,j) -> (fixed[i]), showsettings.fixed)]
        end
    end
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               hlines=hrule, formatters=[fmt__printf(showsettings.floatformat, 5:7)],
               highlighters=highlighters)
end


function show(io::IO, multi::Union{Vector{Model}, Vector{ModelSnapshot}})
    for id in 1:length(multi)
        println(io)
        section(io, join(fill("=", 30)) * "  Model $id  " * join(fill("=", 30)))
        show(io, multi[id])
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
