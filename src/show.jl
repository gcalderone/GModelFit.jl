mutable struct ShowSettings
    plain::Bool
    tableformat::TextFormat
    floatformat::String
    border::Crayon
    header::Crayon
    subheader::Crayon
    fixed::Crayon
    error::Crayon
    highlighted::Crayon
    section::Crayon
    showfixed::Bool
    ShowSettings() = new(false, tf_unicode_rounded, "%9.4g",
                         crayon"light_blue", crayon"light_blue negative bold",
                         crayon"dark_gray bold", crayon"dark_gray",
                         crayon"light_red blink", crayon"negative", crayon"green bold",
                         true)
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
        (name == :meta)  &&  continue
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
    else
        println(io, "Value : ", par.val, "  [", par.low , " : ", par.high, "]")
    end
end


function preparetable(comp::AbstractComponent; cname::String="?", cfixed=false)
    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()

    ctype = split(string(typeof(comp)), ".")
    (ctype[1] == "GFit")  &&   (ctype = ctype[2:end])
    ctype = join(ctype, ".")

    for (qpname, param) in getparams(comp)
        parname = string(qpname.name)
        if qpname.index >= 1
            parname *= "[" * string(qpname.index) * "]"
        end
        parname *= (param.fixed  ?  " (FIXED)"  :  "")
        (!showsettings.showfixed)  &&  param.fixed  &&  continue
        range = strip(@sprintf("%7.2g:%-7.2g", param.low, param.high))
        (range == "-Inf:Inf")  &&  (range = "")
        table = vcat(table,
                     [cname * (cfixed  ?  " (FIXED)"  :  "") ctype parname param.val range])
        push!(fixed, param.fixed)
        push!(error, !(param.low <= param.val <= param.high))
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
    return (table, fixed, error)
end


function show(io::IO, comp::AbstractComponent)
    (table, fixed, error) = preparetable(comp)
    printtable(io, table, ["Component", "Type", "Param.", "Value", "Range"],
               formatters=ft_printf(showsettings.floatformat, [4]),
               highlighters=(Highlighter((data,i,j) -> fixed[i], showsettings.fixed),
                             Highlighter((data,i,j) -> (error[i] &&  (j in (3,4))), showsettings.error)))
end


function show(io::IO, pred::PredRef)
    println(io)
    section(io, "Prediction: $(pred ⋄ :id):")
    (length(pred.cevals) == 0)  &&  (return nothing)

    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, ceval) in pred.cevals
        comp = ceval.comp
        (t, f, e) = preparetable(comp, cname=string(cname), cfixed=(ceval.cfixed >= 1))
        table = vcat(table, t)
        append!(fixed, f .| (ceval.cfixed >= 1))
        append!(error, e)
        push!(hrule, length(error)+1)
    end
    printtable(io, table, ["Component", "Type", "Param.", "Value", "Range"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, [4]),
               highlighters=(Highlighter((data,i,j) -> fixed[i], showsettings.fixed),
                             Highlighter((data,i,j) -> (error[i] &&  (j in (3,4))), showsettings.error)))

    i = 1
    error = Vector{Bool}()
    table = Matrix{Union{String,Int,Float64}}(undef,
                                              length(pred.cevals) + length(pred.revals), 6)
    for (cname, ceval) in pred.cevals
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

    for (rname, reval) in pred.revals
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

    printtable(io, table, ["Component", "Eval. count", "Min", "Max", "Mean", "NaN/Inf"],
               hlines=[0,1, length(pred.cevals)+1,  length(pred.cevals)+length(pred.revals)+1],
               formatters=ft_printf(showsettings.floatformat, 3:5),
               highlighters=(Highlighter((data,i,j) -> (error[i] && j==5), showsettings.error)))
end


function show(io::IO, model::Model)
    for id in 1:length(model.preds)
        show(io, model[id])
    end
end


show(io::IO, par::BestFitPar) = println(io, par.val, " ± ", par.unc,
                                        (par.val == par.patched  ?  ""  :
                                         " (patched value: " * string(par.patched) * ")"))


function preparetable(cname::Symbol, comp::BestFitComp)
    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    watch = Vector{Bool}()

    cname = string(cname)
    for (pname, params) in comp
        if isa(params, AbstractArray)
            for ii in 1:length(params)
                par = params[ii]
                (!showsettings.showfixed)  &&  par.fixed  &&  (par.val == par.patched)  &&  continue
                spname = string(pname) * "[" * string(ii) * "]"
                table = vcat(table, [cname spname par.val par.unc par.patched])
                push!(fixed, par.fixed)
                push!(error, !isfinite(par.unc))
                push!(watch, par.val != par.patched)
                showsettings.plain  ||  (cname = "")
            end
        else
            par = params
            (!showsettings.showfixed)  &&  par.fixed  &&  (par.val == par.patched)  &&  continue
            spname = string(pname)
            table = vcat(table, [cname spname par.val par.unc par.patched])
            push!(fixed, par.fixed)
            push!(error, !isfinite(par.unc))
            push!(watch, par.val != par.patched)
        end
        showsettings.plain  ||  (cname = "")
    end
    return (table, fixed, error, watch)
end


function show(io::IO, comp::BestFitComp)
    (table, fixed, error, watch) = preparetable(Symbol("?"), comp)
    (length(table) == 0)  &&  return
    printtable(io, table , ["Component", "Param.", "Value", "Uncert.", "Patched"],
               hlines=[0,1,size(table)[1]+1], formatters=ft_printf(showsettings.floatformat, [3,4,5]),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]  &&  (j in (2,3,4))), showsettings.fixed),
                             Highlighter((data,i,j) -> (watch[i]  &&  (j==5)), showsettings.highlighted),
                             Highlighter((data,i,j) -> (error[i]  &&  (!fixed[i])  &&  (j==4)), showsettings.error)))
end


function show(io::IO, pred::BestFitPredRef)
    println(io)
    section(io, "Prediction: $(pred ⋄ :id):")

    table = Matrix{Union{String,Float64}}(undef, 0, 5)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    watch = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for (cname, comp) in pred
        (t, f, e, w) = preparetable(cname, comp)
        (length(t) > 0)  ||  continue
        table = vcat(table, t)
        append!(fixed, f)
        append!(error, e)
        append!(watch, w)
        push!(hrule, length(error)+1)
    end
    printtable(io, table, ["Component", "Param.", "Value", "Uncert.", "Patched"],
               hlines=hrule, formatters=ft_printf(showsettings.floatformat, [3,4,5]),
               highlighters=(Highlighter((data,i,j) -> (fixed[i]  &&  (j in (2,3,4))), showsettings.fixed),
                             Highlighter((data,i,j) -> (watch[i]  &&  (j==5)), showsettings.highlighted),
                             Highlighter((data,i,j) -> (error[i]  &&  (!fixed[i])  &&  (j==4)), showsettings.error)))
end

function show(io::IO, res::BestFitResult)
    section(io, "Best Fit results:")

    table = Matrix{Union{String,Float64}}(undef, 0, 6)
    fixed = Vector{Bool}()
    error = Vector{Bool}()
    watch = Vector{Bool}()
    hrule = Vector{Int}()
    push!(hrule, 0, 1)
    for id in 1:length(res.preds)
        show(io, res[id])
    end

    println(io)
    println(io, @sprintf("    #Data  : %8d              Cost   : %-10.5g", res.ndata, res.cost))
    println(io, @sprintf("    #Param : %8d              Red.   : %-10.4g", res.ndata-res.dof, res.cost / res.dof))
    print(  io, @sprintf("    DOF    : %8d              ", res.dof))
    println(io, @sprintf("Prob.  : %-10.4g", 10^res.log10testprob))
    print(io, "    Status : ")
    status = @sprintf("%8s", string(res.status))
    if showsettings.plain
        print(io, status)
    else
        if res.status == :OK
            print(io, crayon"green", status, crayon"default")
        elseif res.status == :Warn
            print(io, crayon"bold yellow", status, crayon"default")
        elseif res.status == :Error
            print(io, crayon"bold red", status, crayon"default")
        else
            print(io, crayon"red", status, crayon"default")
        end
    end
    println(io, @sprintf("              Elapsed: %-10.4g s", res.elapsed))
end

#=
function savelog(filename::String, args...; plain=true)
    orig = showsettings.plain
    showsettings.plain = plain
    try
        f = open(filename, "w")
        for arg in args
            show(f, arg)
            println(f)
        end
        close(f)
    catch
    end
    showsettings.plain = orig
end
=#
