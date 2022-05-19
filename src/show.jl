mutable struct ShowSettings
    plain::Bool
    tableformat::TextFormat
    floatformat::String
    showfixed::Bool
    showevals::Bool
    border::Crayon
    header::Crayon
    subheader::Crayon
    fixed::Crayon
    error::Crayon
    highlighted::Crayon
    section::Crayon
    ShowSettings() = new(false, tf_unicode_rounded, "%9.4g", true, true,
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
        if isa(dom, Domain{1})
            a = convert(Vector{Float64}, dom)
        else
            a = isa(dom, Domain)  ?  dom[i]  :  axis(dom, i)
        end
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


#=
The following is needed since AbstractData <: AbstractArray:
without it the show(bstractArray) would be invoked.
=#
show(io::IO, mime::MIME"text/plain", data::AbstractData) = show(io, data)

function show(io::IO, data::AbstractData)
    section(io, typeof(data), ": (length: ", (length(data.val)), ")")
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    hrule = Vector{Int}()
    push!(hrule, 0, 1)

    names = fieldnames(typeof(data))
    error = Vector{Bool}()
    for name in names
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
        if par.val == par.patched
            println(io, "Value : ", par.val, " ± ", par.unc,  "  [", par.low , " : ", par.high, "] ")
        else
            println(io, "Value : ", par.val, " ± ", par.unc,  "  [", par.low , " : ", par.high, "], patched: " , par.patched)
        end
    end
end


function preparetable(comp::AbstractComponent; cname::String="?", cfixed=false)
    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    watch = Vector{Bool}()

    ctype = split(string(typeof(comp)), ".")
    (ctype[1] == "GFit")  &&   (ctype = ctype[2:end])
    ctype = join(ctype, ".")

    for (pid, param) in getparams(comp)
        parname = string(pid.name)
        if pid.index >= 1
            parname *= "[" * string(pid.index) * "]"
        end
        parname *= (param.fixed  ?  " (FIXED)"  :  "")
        (!showsettings.showfixed)  &&  param.fixed  &&  continue
        range = strip(@sprintf("%7.2g:%-7.2g", param.low, param.high))
        (range == "-Inf:Inf")  &&  (range = "")
        table = vcat(table,
                     [cname * (cfixed  ?  " (FIXED)"  :  "") ctype parname range param.val param.unc param.patched])
        push!(fixed, param.fixed)
        push!(error, !(param.low <= param.val <= param.high))
        push!(watch, param.val != param.patched)
        if !showsettings.plain
            cname = ""
            ctype = ""
        end
    end
    if length(table) == 0
        table = vcat(table, [cname ctype "" NaN ""])
        push!(fixed, false)
        push!(error, false)
    end
    return (table, fixed, error, watch)
end


function show(io::IO, comp::AbstractComponent)
    (table, fixed, error, watch) = preparetable(comp)
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Patched"],
               formatters=ft_printf(showsettings.floatformat, 5:7),
               highlighters=(Highlighter((data,i,j) -> (watch[i]  &&  (j==7)), showsettings.highlighted),
                             Highlighter((data,i,j) -> (error[i]  && (j in (3,4,5))), showsettings.error),
                             Highlighter((data,i,j) -> fixed[i], showsettings.fixed)))
end


function show(io::IO, red::λFunct)
    println(io, red.display)
end

function show(io::IO, model::Model)
    println(io)
    section(io, "Components:")
    (length(model.cevals) == 0)  &&  (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 7)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    watch = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, ceval) in model.cevals
        comp = ceval.comp
        (t, f, e, w) = preparetable(comp, cname=string(cname), cfixed=(ceval.cfixed >= 1))
        table = vcat(table, t)
        append!(fixed, f .| (ceval.cfixed >= 1))
        append!(error, e)
        append!(watch, w)
        push!(hrule, length(error)+1)
    end
    printtable(io, table, ["Component", "Type", "Param.", "Range", "Value", "Uncert.", "Patched"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, 5:7),
               highlighters=(Highlighter((data,i,j) -> (watch[i]  &&  (j==7)), showsettings.highlighted),
                             Highlighter((data,i,j) -> (error[i] &&  (j in (3,4,5))), showsettings.error),
                             Highlighter((data,i,j) -> fixed[i], showsettings.fixed)))

    i = 1
    error = Vector{Bool}()
    table = Matrix{Union{String,Int,Float64}}(undef,
                                              length(model.cevals) + length(model.revals), 6)
    for (cname, ceval) in model.cevals
        result = ceval.buffer
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

    for (rname, reval) in model.revals
        result = reval.buffer
        v = view(result, findall(isfinite.(result)))
        (length(v) == 0)  &&  (v = [NaN])
        nan = length(findall(isnan.(result)))
        inf = length(findall(isinf.(result)))
        table[i, 1] = string(rname)
        table[i, 2] = reval.counter
        table[i, 3:5] = [minimum(v), maximum(v), mean(v)]
        table[i, 6] = (nan+inf > 0  ?  string(nan)  :  "")
        push!(error, (nan+inf > 0))
        i += 1
    end

    if length(model.patchfuncts) > 0
        section(io, "Patch expressions:")
        for pf in model.patchfuncts
            println(io, string(pf.display))
        end
    end

    println(io)
    section(io, "Reducers:")
    for (rname, reval) in model.revals
        print(rname, ": ")
        show(io, reval.red)
        println(io)
    end

    if showsettings.showevals
        section(io, "Evaluations:")
        printtable(io, table, ["Component", "Eval. count", "Min", "Max", "Mean", "NaN/Inf"],
                   hlines=[0,1, length(model.cevals)+1,  length(model.cevals)+length(model.revals)+1],
                   formatters=ft_printf(showsettings.floatformat, 3:5),
                   highlighters=(Highlighter((data,i,j) -> (error[i] && j==5), showsettings.error)))
    end
end


function show(io::IO, multi::MultiModel)
    for id in 1:length(multi.models)
        section(io, "\n=====================================================================")
        section(io, "Model $id:")
        show(io, multi.models[id])
    end
    println(io)

    if length(multi.patchfuncts) > 0
        section(io, "Multi model patch expressions:")
        for pf in multi.patchfuncts
            println(io, string(pf.display))
        end
    end
end


function show(io::IO, res::FitResult)
    section(io, "Fit results:")

    println(io, @sprintf("    #Data  : %8d              Fit-stat: %-10.5g", res.ndata, res.fitstat))
    println(io, @sprintf("    #Param : %8d              Red. GOF: %-10.4g", res.nfree, res.gofstat / res.dof))
    println(io, @sprintf("    DOF    : %8d              Prob.   : %-10.4g", res.dof, 10^res.log10testprob))

    print(io,            "    Status : ")
    (crayon, status, message) = as_string(res.mzer)
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
