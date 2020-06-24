mutable struct ShowSettings
    tableformat::PrettyTableFormat
    floatformat::String
    border::Crayon
    header::Crayon
    subheader::Crayon
    fixed::Crayon
    error::Crayon
    highlighted::Crayon
    section::Crayon
    fixedpars::Bool
    ShowSettings() = new(unicode_rounded, "%9.4g",
                         crayon"light_blue", crayon"light_blue negative bold",
                         crayon"dark_gray bold", crayon"dark_gray",
                         crayon"light_red blink", crayon"negative", crayon"green bold",
                         true)
end

const showsettings = ShowSettings()

function printtable(args...; kw...)
    pretty_table(args...; tf=showsettings.tableformat, kw...,
                 border_crayon=showsettings.border,
                 header_crayon=showsettings.header,
                 subheader_crayon=showsettings.subheader, crop=:horizontal)
end

section(io, args...) = println(io, showsettings.section, args...)

function left(s::String, maxlen::Int)
    (length(s) <= maxlen)  &&  (return s)
    return s[1:maxlen]
end

function show(io::IO, dom::AbstractCartesianDomain)
    section(io, "Cartesian domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
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
    printtable(io, table, ["Dim", "Size", "Min val", "Max val", "Min step", "Max step"],
               formatters=ft_printf(showsettings.floatformat, 3:6))
end


function show(io::IO, dom::AbstractLinearDomain)
    section(io, "Linear domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
    table = Matrix{Union{Int,Float64}}(undef, ndims(dom), 3)
    for i in 1:ndims(dom)
        table[i, 1] = i
        table[i, 2] = getaxismin(dom, i)
        table[i, 3] = getaxismax(dom, i)
    end
    printtable(io, table, ["Dim", "Min val", "Max val"],
               formatters=ft_printf(showsettings.floatformat, 2:3))
end


# Special case for Domain_1D: treat it as a Cartesian domain, despite it is a Linear one.
function show(io::IO, dom::Domain_1D)
    section(io, "Domain (ndims: ", ndims(dom), ", length: ", length(dom), ")")
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
    printtable(io, table, ["Dim", "Size", "Min val", "Max val", "Min step", "Max step"],
               formatters=ft_printf(showsettings.floatformat, 3:6))
end


function show(io::IO, data::AbstractData)
    section(io, typeof(data), "   length: ", (length(data.val)))
    table = Matrix{Union{String,Float64}}(undef, 0, 7)

    names = fieldnames(typeof(data))
    error = Vector{Bool}()
    for name in names
        a = getfield(data, name)
        nan = length(findall(isnan.(a))) + length(findall(isinf.(a)))
        a = a[findall(isfinite.(a))]
        push!(error, nan > 0)
        table = vcat(table, [string(name) minimum(a) maximum(a) mean(a) median(a) std(a) (nan > 0  ?  string(nan)  :  "") ])
    end
    printtable(io, table, ["", "Min", "Max", "Mean", "Median", "Std. dev.", "Nan/Inf"],
               formatters=ft_printf(showsettings.floatformat, 2:6),
               highlighters=(Highlighter((data,i,j) -> error[i], showsettings.error)))
end


function show(io::IO, par::Parameter)
    if par.free
        println(io, "Value : ", par.val, "  [", par.low , " : ", par.high, "]")
    else
        println(io, "Value : ", par.val, "   (FIXED)")
    end
end


function preparetable(comp::AbstractComponent, cname="")
    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()

    ctype = split(string(typeof(comp)), ".")
    (ctype[1] == "GFit")  &&   (ctype = ctype[2:end])
    ctype = join(ctype, ".")

    for (pname, param) in getparams(comp)
        parname = string(pname[1])
        if pname[2] >= 1
            parname *= "[" * string(pname[2]) * "]"
        end
        parname *= (param.free  ?  ""  :  " (FIXED)")
        (!showsettings.fixedpars)  &&  (!param.free)  &&  continue
        range = strip(@sprintf("%7.2g:%-7.2g", param.low, param.high))
        (range == "-Inf:Inf")  &&  (range = "")
        table = vcat(table, [cname ctype parname param.val range])
        push!(fixed, !param.free)
        push!(error, !(param.low <= param.val <= param.high))
        cname = ""
        ctype = ""
    end
    if length(table) == 0
        table = vcat(table, [cname ctype "" NaN ""])
        push!(fixed, false)
        push!(error, false)
    end
    return (table, fixed, error)
end


show(io::IO, comp::AbstractComponent) =
    show(io, OrderedDict(Symbol("?") => comp), OrderedDict(Symbol("?") => true))


function show(io::IO, dict::OrderedDict{Symbol, T}, cfree::OrderedDict{Symbol, Bool}) where T <: AbstractComponent
    (length(dict) > 0)  ||  (return nothing)
    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, comp) in dict
        (t, f, e) = preparetable(comp, string(cname) .* (cfree[cname]  ?  ""  :  " (FIXED)"))
        table = vcat(table, t)
        append!(fixed, f .| (.!cfree[cname]))
        append!(error, e)
        push!(hrule, length(error)+1)
    end
    printtable(io, table , ["Component" "Type" "Param." "Value" "Range"], alignment=:l,
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, [4]),
               highlighters=(Highlighter((data,i,j) -> fixed[i], showsettings.fixed),
                             Highlighter((data,i,j) -> (error[i] &&  (j in (3,4))), showsettings.error)))
end


show(io::IO, mime::MIME"text/plain", model::Model) = show(io, model)
function show(io::IO, model::Model)
    section(io, "Components:")
    length(model.comps) != 0  || (return nothing)
    show(io, model.comps, model.cfree)

    for i in 1:length(model.preds)
        println(io)
        section(io, "Prediction #$i: " * model.preds[1].label)
        show(io, model.preds[i])
    end
end



function show(io::IO, pred::Prediction)
    (length(pred.cevals) == 0)  &&  (return nothing)
    table = Matrix{Union{String,Int,Float64}}(undef, length(pred.cevals) + 1, 6)
    error = Vector{Bool}()

    i = 1
    for (cname, ceval) in pred.cevals
        result = ceval.eval
        v = view(result, findall(isfinite.(result)))
        (length(v) == 0)  &&  (v = [NaN])
        nan = length(findall(isnan.(result)))
        inf = length(findall(isinf.(result)))
        table[i, 1] = string(cname)
        table[i, 2] = ceval.counter
        table[i, 3:5] = [minimum(v), maximum(v), mean(v)]
        table[i, 6] = (nan+inf > 0  ?  string(nan)  :  "")
        push!(error, (nan+inf > 0))
        i += 1
    end

    result = pred.eval
    v = view(result, findall(isfinite.(result)))
    nan = length(findall(isnan.(result)))
    inf = length(findall(isinf.(result)))
    table[i, 1] = "Reduced"
    table[i, 2] = pred.counter
    table[i, 3:5] = [minimum(v), maximum(v), mean(v)]
    table[i, 6] = (nan+inf > 0  ?  string(nan)  :  "")
    push!(error, (nan+inf > 0))

    printtable(io, table, ["Label", "Counter", "Min", "Max", "Mean", "NaN/Inf"], alignment=:l,
               hlines=[0,1,length(pred.cevals)+1,length(pred.cevals)+2],
               formatters=ft_printf(showsettings.floatformat, 3:5),
               highlighters=(Highlighter((data,i,j) -> (error[i] && j==5), showsettings.error)))
end


show(io::IO, par::BestFitPar) = println(io, par.val, " ± ", par.unc,
                                        (par.val == par.calc  ?  ""  :
                                         " (calculated value: " * string(par.calc) * ")"))


function preparetable(comp::BestFitComp)
    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    watch = Vector{Bool}()

    for (pname, param) in comp
        if isa(param, Vector{BestFitPar})
            for ii in 1:length(param)
                par = param[ii]
                (!showsettings.fixedpars)  &&  (!par.free)  &&  continue
                spname = string(pname) * "[" * string(ii) * "]"
                table = vcat(table, ["" spname par.val par.unc par.calc])
                push!(fixed, !par.free)
                push!(error, !isfinite(par.unc))
                push!(watch, par.val != par.calc)
            end
        else
            par = param
            (!showsettings.fixedpars)  &&  (!par.free)  &&  continue
            spname = string(pname)
            table = vcat(table, ["" spname par.val par.unc par.calc])
            push!(fixed, !par.free)
            push!(error, !isfinite(par.unc))
            push!(watch, par.val != par.calc)
        end
    end
    return (table, fixed, error, watch)
end


function show(io::IO, comp::BestFitComp)
    (table, fixed, error, watch) = preparetable(comp)
    (length(table) == 0)  &&  return
    printtable(io, table , ["Component" "Param." "Value" "Uncert." "Calculated"], alignment=:l,
               hlines=[0,1,size(table)[1]+1], formatters=ft_printf(showsettings.floatformat, [3,4,5]),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]  &&  (j in (2,3,4))), showsettings.fixed),
                             Highlighter((data,i,j) -> (watch[i]  &&  (j==5)), showsettings.highlighted),
                             Highlighter((data,i,j) -> (error[i]  &&  (!fixed[i])  &&  (j==4)), showsettings.error)))
end


function show(io::IO, res::BestFitResult)
    section(io, "Best Fit results:")

    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    watch = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, comp) in res.comps
        if length(comp) > 0
            (t, f, e, w) = preparetable(comp)
            (length(t) > 0)  ||  continue
            t[1,1] = string(cname)
            table = vcat(table, t)
            append!(fixed, f)
            append!(error, e)
            append!(watch, w)
            push!(hrule, length(error)+1)
        end
    end
    printtable(io, table , ["Component" "Param." "Value" "Uncert." "Calculated"], alignment=:l,
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, [3,4,5]),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]  &&  (j in (2,3,4))), showsettings.fixed),
                             Highlighter((data,i,j) -> (watch[i]  &&  (j==5)), showsettings.highlighted),
                             Highlighter((data,i,j) -> (error[i]  &&  (!fixed[i])  &&  (j==4)), showsettings.error)))

    println(io)
    println(io, @sprintf("    #Data  : %10d              Cost : %-10.5g", res.ndata, res.cost))
    println(io, @sprintf("    #Param : %10d              Red. : %-10.4g", res.ndata-res.dof, res.cost / res.dof))
    print(  io, @sprintf("    DOF    : %10d              ", res.dof))
    if res.log10testprob < -3
        println(io, @sprintf("Prob.: 10^%-10.4g", res.log10testprob))
    else
        println(io, @sprintf("Prob.: %10.4g", 10^res.log10testprob))
    end
    printstyled(io, "    Status :  ")
    if res.status == :Optimal
        printstyled(color=:green, io, @sprintf("%-15s", "Optimal"))
    elseif res.status == :NonOptimal
        printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Non Optimal"))
    elseif res.status == :Warn
         printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Warning"))
    elseif res.status == :Error
        printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Error"))
    else
        printstyled(color=printcolorerr(), io, @sprintf("%-15s", "Unknown (" * string(res.status) * "), see fitter output"))
    end
    println(io, @sprintf("        Elapsed: %-10.4g s", res.elapsed))
end
