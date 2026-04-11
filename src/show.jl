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

function printtable(io, table::DataFrame; formatters=[], highlighters=TextHighlighter[], hlines=Int[])
    common = pairs((
        formatters=formatters,
        alignment=:l,
        column_label_alignment=:c,
        limit_printing=false,
        compact_printing=false,
        fit_table_in_display_horizontally=false,
        fit_table_in_display_vertically=false,
        maximum_number_of_columns=-1,
        maximum_number_of_rows=-1,
        show_first_column_label_only = true))

    if showsettings.plain
        kws = pairs((table_format=TextTableFormat(borders=text_table_borders__compact,
                                                  horizontal_lines_at_data_rows=hlines),
                     common...))
    else
        kws = pairs((table_format=TextTableFormat(borders=showsettings.tableborders,
                                                  horizontal_lines_at_data_rows=hlines),
                     style=TextTableStyle(table_border=showsettings.border,
                                          first_line_column_label=showsettings.header),
                     common..., highlighters))
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

show(io::IO, dom::AbstractDomain) = print(io, string(typeof(dom)) * " (ndims: ", ndims(dom), ", length: ", length(dom), ")")
function show(io::IO, mime::MIME"text/plain", dom::AbstractDomain)
    section(io, string(typeof(dom)) * " (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    table = DataFrame(:Dim => Int[], :Length => Int[], Symbol("Min val") => Float64[], Symbol("Max val") => Float64[],
                      Symbol("Min step") => Float64[], Symbol("Max step") => Float64[])
    for i in 1:ndims(dom)
        if isa(dom, CartesianDomain)
            vv = gridcoords(dom, i)
        else
            vv = coords(dom, i)
        end
        steps = 0
        if length(vv) > 1
            steps = vv .- circshift(vv, 1)
            steps = steps[2:end]
        end
        push!(table, (i, length(vv), minimum(vv), maximum(vv), minimum(steps), maximum(steps)))
    end
    printtable(io, table, formatters=[fmt__printf(showsettings.floatformat, 3:6)])
end

show(io::IO, data::GaussianData) = print(io, typeof(data), ": (length: ", (length(data)), ")")
function show(io::IO, mime::MIME"text/plain", data::GaussianData)
    section(io, typeof(data), ": (length: ", (length(data)), ")")
    table = DataFrame(Symbol("") => String[], :Min => Float64[], :Max => Float64[], :Mean => Float64[],
                      :Median => Float64[], Symbol("Std. dev.") => Float64[], Symbol("Nan/Inf") => String[])
    names = fieldnames(typeof(data))
    error = Vector{Bool}()
    for vv in [values(data), uncerts(data)]
        nan = length(findall(isnan.(vv))) + length(findall(isinf.(vv)))
        push!(error, nan > 0)
        vv = vv[findall(isfinite.(vv))]
        (length(vv) == 0)  &&  (vv = [NaN])
        push!(table, ("", minimum(vv), maximum(vv), mean(vv), median(vv), std(vv), (nan > 0  ?  string(nan)  :  "")))
    end
    table[:, 1] .= ["values", "uncerts"]

    highlighters = [TextHighlighter((data,i,j) -> error[i], showsettings.error)]
    printtable(io, table,
               formatters=[fmt__printf(showsettings.floatformat, 2:6)],
               highlighters=highlighters)
end

function show(io::IO, red::FunctDesc)
    print(io, red.display)
end

# --------------------------------------------------------------------
function show(io::IO, par::Parameter)
    print(io, "Value: ", par.val)
    if par.fixed
        print(io, "   (fixed)")
    end
    print(io, "  [", par.low , " : ", par.high, "]")
end

function show(io::IO, par::ParameterEval)
    print(io, "Value: ", par.val)
    if par.fixed  ||  par.actually_fixed
        print(io, "   (fixed)")
    else
        print(io, " ± ", par.unc)
    end
    print(io, "  [", par.low , " : ", par.high, "]")
    if par.val != par.mval
        print(io, ", model value: " , par.mval)
    end
    par._dirty  &&  print(io, " (dirty)")
end

function show_comps(io, model::AbstractModel, gkey::Tuple)
    function push_row!(model::Model, comp::CompSlot, key::Tuple, prefix::String)
        push!(table, (prefix * string(key[end]), stype(comp), isfrozen(comp)  ?  "(frozen)"  :  ""))
        push!(meta.fixed, isfrozen(comp))
        push!(meta.warns, false)
    end
    function push_row!(model::AbstractModelEval, comp::AbstractCompSlotEval, key::Tuple, prefix::String)
        v = comp.evalbuf
        push!(table, (prefix * string(key[end]), stype(comp), isfrozen(comp)  ?  "(frozen)"  :  "",
                      comp.counter, minimum(v), maximum(v), mean(v), count(.!isfinite.(v))))
        push!(meta.fixed, isfrozen(comp))
        push!(meta.warns, false)
    end
              
    function _recursive(model::AbstractModel, gkey::GT, cname::Symbol, prefix::String) where GT
        BRANCH = showsettings.plain  ?  "+-"  :  "├╴ "
        BRCONT = showsettings.plain  ?  "| "  :  "│  "
        BREND  = showsettings.plain  ?  "+-"  :  "└╴ "
        BROVER = showsettings.plain  ?  "  "  :  "   "

        cnames = dependencies(getcomps(model)[(gkey..., cname)])
        for i in 1:length(cnames)
            cname = cnames[i]
            key = (gkey..., cname)
            push_row!(model, getcomps(model)[key], key, prefix * (i < length(cnames)  ?  BRANCH  :  BREND))
            _recursive(model, gkey, cname, prefix * (i < length(cnames)  ?  BRCONT  :  BROVER))
        end
    end

    if isa(model, Model)
        table = DataFrame(:Component => String[], :Type => String[], :Frozen => String[])
    else
        table = DataFrame(:Component => String[], :Type => String[], :Frozen => String[],
                          Symbol("Eval. count") => Int[], Symbol("Min") => Float64[], Symbol("Max") => Float64[],
                          Symbol("Mean") => Float64[], Symbol("NaN/Inf") => Int[])
    end
    meta = (fixed=Bool[], warns=Bool[], hrule=Int[])
    if length(getcomps(model)) > 0
        cname = find_sequence(model, gkey)[end]
        key = (gkey..., cname)
        push_row!(model, getcomps(model)[key], key, "")
        _recursive(model, gkey, cname, "")
        printtable(io, table, formatters=[fmt__printf(showsettings.floatformat, 5:7)],
                   highlighters=[TextHighlighter((data,i,j) -> (meta.fixed[i]), showsettings.fixed),
                                 TextHighlighter((data,i,j) -> (meta.warns[i]), showsettings.error)])
    end
end

function show_params(io, model::AbstractModel, gkey::Tuple, select_ckey::Union{Nothing, Tuple}=nothing)
    if isa(model, Model)
        table = DataFrame(:Component => String[], :Parameter => String[], :Range => String[], :Value => Float64[],
                          Symbol("Patch / reparam") => String[])
    else
        table = DataFrame(:Component => String[], :Parameter => String[], :Range => String[], :Value => Float64[],
                          Symbol("Uncert.") => Float64[], Symbol("Model value") => Float64[], Symbol("Patch / reparam") => String[])
    end
    meta = (fixed=Bool[], warns=Bool[], hrule=Int[])
    for (ckey, comp) in getcomps(model)
        if isnothing(select_ckey)
            (group_key(model, ckey) == gkey)  ||  continue
        else
            (ckey == select_ckey)  ||  continue
        end
        scname = string(ckey[end])
        if haspath(getparams(model), ckey)
            for (pname, par) in view(getparams(model), ckey)
                range =  (isfinite(par.low )  ?  strip(@sprintf("%7.2g", par.low ))  :  (par.low  > 0  ?  "+∞"  :  "-∞")) * " : "
                range *= (isfinite(par.high)  ?  strip(@sprintf("%7.2g", par.high))  :  (par.high > 0  ?  "+∞"  :  "-∞"))
                patch = ""
                isa(par.patch, Symbol)       &&  (patch = string(par.patch))
                isa(par.patch, FunctDesc)    &&  (patch = par.patch.display)
                isa(par.reparam, FunctDesc)  &&  (patch = par.reparam.display)

                if isa(par, Parameter)
                    push!(meta.fixed, par.fixed)
                    push!(meta.warns, !isfinite(par.val))
                    push!(table, (scname, string(pname[1]) * (meta.fixed[end]  ?  " (fixed) "  :  ""), range, par.val, patch))
                else
                    push!(meta.fixed, par.fixed  ||  par.actually_fixed)
                    push!(meta.warns, !isfinite(par.val)  ||  par._dirty  ||  !isfinite(par.unc)  ||  !isfinite(par.mval))
                    push!(table, (scname, string(pname[1]) * (meta.fixed[end]  ?  " (fixed) "  :  ""), range, par.val, par.unc, par.mval, patch))
                end
                scname = ""
            end
            (nrow(table) > 0)  &&  push!(meta.hrule, nrow(table))
        end
    end
    if nrow(table) > 0
        printtable(io, table, formatters=[fmt__printf(showsettings.floatformat, 4:6)], hlines=meta.hrule,
                   highlighters=[TextHighlighter((data,i,j) -> (meta.fixed[i]), showsettings.fixed),
                                 TextHighlighter((data,i,j) -> ((j in 4:6)  &&  meta.warns[i]), showsettings.error),
                                 TextHighlighter((data,i,j) -> ((j == 6)  &&  (table[i, Symbol("Patch / reparam")] == "")), showsettings.fixed)])
    end
end

function show(io::IO, model::AbstractModel)
    print(io, typeof(model), ": ($(length(getgroups(model))) groups, $(length(getcomps(model))) components, $(length(getparams(model))) parameters)")
    if typeof(model) <: AbstractLikelihood
        print(io, ", $(nfree(model)) free parameters, $(ndata(model)) data points, red. fit stat.: $(gofstat(model))")
    end
end

function show(io::IO, mime::MIME"text/plain", model::AbstractModel)
    for (gkey, group) in getgroups(model)
        show(io, mime, group)
    end
end

function show(io::IO, mime::MIME"text/plain", lh::AbstractLikelihood)
    for (gkey, group) in getgroups(lh)
        show(io, mime, group)
    end
    println(io)
    section(io, "Fit stat.")
    print(io, @sprintf("Nfree pars: %d, Ndata: %d, fit stat.: %.5g", nfree(lh), ndata(lh), gofstat(lh)))
    n = dof(lh)
    if !isnothing(n)
        print(io, @sprintf(", red. chi-squared: %.5g", gofstat(lh) / n))
    end
end

show(io::IO, group::AbstractGroup) = print(io, typeof(group), ": ($(length(getcomps(group))) components, $(length(getparams(group))) parameters)")
function show(io::IO, mime::MIME"text/plain", group::AbstractGroup)
    gkey = group.key
    section(io, "Components", (gkey != ())  ?  " in group $gkey:"   :  ":")
    show_comps(io, group.rootmodel, group.key)
    println(io)
    section(io, "Parameters", (gkey != ())  ?  " in group $gkey:"   :  ":")
    show_params(io, group.rootmodel, group.key)
    println(io)
    if isa(group, GroupLH)  ||  isa(group, GroupLHSnapshot)
        @info gofstat(group)
        print(io, @sprintf("Goodness of fit%s: %.5g (%s)",
                           ((gkey != ())  ?  " in group $gkey:"   :  ":"),
                           gofstat(group), gofstat_name(group)))
    end
    println(io)
end

show(io::IO, comp::AbstractCompSlot) = print(io, typeof(comp), ": ($(length(getparams(comp))) parameters)")
function show(io::IO, mime::MIME"text/plain", comp::AbstractCompSlot)
    section(io, "Parameters:")
    show_params(io, comp.rootmodel, group_key(comp.rootmodel, comp.key), comp.key)
end

function show(io::IO, mime::MIME"text/plain", t::Tuple{Likelihood, AbstractSolverStatus{T}}) where T
    show(io, mime, t[1])
    println(io)
    show(io, mime, t[2])
end


getmessage(status::SolverStatusOK) = crayon"green", "OK"
getmessage(status::SolverStatusWarn) = crayon"bold yellow", "WARN\n" * status.message
getmessage(status::SolverStatusError) = crayon"bold red", "ERROR\n" * status.message

function show(io::IO, status::AbstractSolverStatus{T}) where {T}
    print(io, "Solver status: ")
    color, ss = getmessage(status)
    if showsettings.plain
        print(io, ss)
    else
        print(io, color, ss, crayon"default")
    end
    print(io, crayon"dark_gray", "   (", string(typeof(status.status)), ")")
    print(io, crayon"reset")
end
