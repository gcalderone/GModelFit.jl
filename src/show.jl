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
                         crayon"light_red", crayon"negative", crayon"green bold")
end

const showsettings = ShowSettings()

function printtable(io, table, header, args...; formatters=(), hlines=:none, kw...)
    if showsettings.plain
        pretty_table(io, table; header=header, formatters=formatters, alignment=:l, crop=:none, tf=tf_compact, hlines=hlines,
                     highlighters=())
    else
        pretty_table(io, table; header=header, formatters=formatters, alignment=:l, crop=:none, tf=showsettings.tableformat, hlines=hlines,
                     border_crayon=showsettings.border, header_crayon=showsettings.header, subheader_crayon=showsettings.subheader,
                     kw...)
    end
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
    if showsettings.plain
        highlighters = nothing
    else
        highlighters = Highlighter((data,i,j) -> error[i], showsettings.error)
    end
    printtable(io, table, ["", "Min", "Max", "Mean", "Median", "Std. dev.", "Nan/Inf"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 2:6),
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


function preparetable(comp::Union{AbstractComponent, GModelFit.PV.PVComp{GModelFit.Parameter}};
                      cname::String="?", ctype="?", cfixed=false)
    table = Matrix{Union{String,Float64}}(undef, 0, 8)
    fixed = Vector{Bool}()
    warns = Vector{Bool}()

    for (pname, param) in getparams(comp)
        (!showsettings.showfixed)  &&  param.fixed  &&  continue
        range = strip(@sprintf("%7.2g:%-7.2g", param.low, param.high))
        # (range == "-Inf:Inf")  &&  (range = "")
        patch = ""
        isa(param.patch, Symbol)  &&  (patch = string(param.patch))
        isa(param.patch, FunctDesc)  &&  (patch = param.patch.display)
        isa(param.mpatch,FunctDesc)  &&  (patch = param.mpatch.display)
        table = vcat(table,
                     permutedims([cname * (cfixed  ?  " (fixed)"  :  ""), ctype,
                                  string(pname), range, param.val,
                                  (param.fixed | cfixed  ?  " (fixed)"  :  param.unc),
                                  (patch == ""  ?  ""  :  param.actual), patch]))
        push!(fixed, param.fixed)
        if !param.fixed  &&  (isnan(param.unc)  ||  (param.unc <= 0.))
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


function show(io::IO, comp::Union{AbstractComponent, GModelFit.PV.PVComp{GModelFit.Parameter}})
    ctype = isa(comp, AbstractComponent)  ?  string(typeof(comp))  :  "?"
    (table, fixed, warns) = preparetable(comp, ctype=ctype)
    if showsettings.plain
        highlighters = nothing
    else
        highlighters = (Highlighter((data,i,j) -> (fixed[i]  &&  (j in (3,4,5,6))), showsettings.fixed),
                        Highlighter((data,i,j) -> (warns[i]  &&  (j in (3,4,5,6))), showsettings.error))
    end
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               formatters=ft_printf(showsettings.floatformat, 5:7),
               highlighters=highlighters)
end


function show(io::IO, red::FunctDesc)
    println(io, red.display)
end


function tabledeps(model::Union{Model, ModelSnapshot})
    function alldeps(model::Union{Model, GModelFit.ModelSnapshot}, cname=nothing, level=0)
        out = Vector{Tuple}()
        if isnothing(cname)
            cname = GModelFit.find_maincomp(model)
            append!(out, alldeps(model, cname, level+1))
        else
            push!(out, (level, cname))
            for d in GModelFit.dependencies(model, cname)
                append!(out, alldeps(model, d, level+1))
            end
        end
        return out
    end

    allcomps = alldeps(model)
    maxdepth = maximum(getindex.(allcomps, 1))
    prefix = fill("", length(allcomps), maxdepth)
    for i in 1:length(allcomps)
        prefix[i, allcomps[i][1]] = string(allcomps[i][2])
    end

    BRANCH = "├─╴"
    BRCONT = "│  "
    BREND  = "└─╴"
    for i in 2:length(allcomps)
        for j in 1:(maxdepth-1)
            if prefix[i,j] == ""                    # empty cell
                if any(prefix[1:(i-1), j] .!= "")   # ...it has a parent
                    if prefix[i,j+1] != ""          # ...it is a branch
                        prefix[i,j] = BRANCH        # Add a branch
                    end
                end
            end
        end
    end
    for i in 2:length(allcomps)
        for j in 1:(maxdepth-1)
            if prefix[i,j] == BRANCH                  # it is a branch
                # Check if this is the last row for this branch
                if all(prefix[(i+1):end,j] .== "")
                    prefix[i,j] = BREND
                elseif (prefix[i+1,j] != "")  &&  (prefix[i+1,j] != BRANCH)
                    prefix[i,j] = BREND
                end
            end
        end
    end
    # Join branches with vertical lines
    for i in 2:length(allcomps)
        for j in 1:(maxdepth-1)
            if (prefix[i,j] == "")  &&  (prefix[i-1,j] in [BRANCH, BRCONT])
                prefix[i,j] = BRCONT
            end
        end
    end
    prefix[prefix .== ""] .= " "^length(BRANCH)

    prefix = [string(rstrip(join(prefix[i,:]))) for i in 1:length(allcomps)]
    if showsettings.plain
        prefix = [replace(p,
                          BRANCH => "+ ",
                          BRCONT => "| ",
                          BREND  => "+ ") for p in prefix]
    end

    table = Matrix{Union{String,Int,Float64}}(undef, length(allcomps), 8)
    fixed = Vector{Bool}()
    for i in 1:length(allcomps)
        cname = allcomps[i][2]
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
        highlighters = Highlighter((data,i,j) -> fixed[i], showsettings.fixed)
    end
    if !isa(model, Model)
        printtable(io, table, ["Component", "Type", "#Free", "Eval. count", "Min", "Max", "Mean", "NaN/Inf"],
                   hlines=[0,1, size(table)[1]+1],
                   formatters=ft_printf(showsettings.floatformat, 5:7),
                   highlighters=highlighters)
    else
        printtable(io, table[:, 1:3], ["Component", "Type", "#Free"],
                   hlines=[0,1, size(table)[1]+1],
                   formatters=ft_printf(showsettings.floatformat, 5:7),
                   highlighters=highlighters)
    end
    println(io)

    section(io, "Parameters:")
    (length(keys(model)) == 0)  &&  (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 8)
    fixed = Vector{Bool}()
    warns = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for cname in keys(model)
        comp = model[cname]
        (t, f, w) = preparetable(comp,
                                 cname=string(cname),
                                 ctype=comptype(model, cname),
                                 cfixed=isfreezed(model, cname))
        table = vcat(table, t)
        append!(fixed, f .| isfreezed(model, cname))
        append!(warns, w)
        push!(hrule, length(fixed)+1)
    end

    if showsettings.plain
        highlighters = nothing
    else
        if !isa(model, Model)
            highlighters = (Highlighter((data,i,j) -> (fixed[i]), showsettings.fixed),
                            Highlighter((data,i,j) -> (warns[i]  &&  (j in (3,4,5,6))), showsettings.error))
        else
            highlighters = Highlighter((data,i,j) -> (fixed[i]), showsettings.fixed)
        end
    end
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Actual", "Patch"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 5:7),
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
    print(io, "Status: ")
    color, ss = getmessage(status)
    if showsettings.plain
        print(io, @sprintf("%-8s", ss))
    else
        print(io, color, @sprintf("%-8s", ss), crayon"default")
    end
end


function show(io::IO, res::FitSummary)
    section(io, "Fit results:", newline=false)
    print(io, @sprintf(" #data: %d, #free pars: %d, red. fit stat.: %10.5g, ", res.ndata, res.nfree, res.fitstat))
    show(io, res.status)
    println(io)
end
