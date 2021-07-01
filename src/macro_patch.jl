using MacroTools

macro patch!(_expr)
    expr = prettify(_expr)
    mname = Symbol[]
    check = Expr[]
    fixed = Expr[]
    
    function extract_param_expr(ex::Expr)
        if  (ex.head == :.)         &&
            (length(ex.args) == 2)
            cname = ex.args[1]
            # pname = ex.args[2]  unused...
            
            if  (cname.head == :ref)       &&
                (length(cname.args) == 2)
                if isa(cname.args[1], Symbol)
                    return (cname.args[1], :(haskey($(cname.args[1]), $(cname.args[2]))), :($(ex).fixed = true))
                elseif isa(cname.args[1], Expr)  &&
                       (cname.args[1].head == :ref)       &&
                       (length(cname.args[1].args) == 2)  &&
                       isa(cname.args[1].args[1], Symbol)
                    return (cname.args[1].args[1], :(haskey($(cname.args[1]), $(cname.args[2]))), :($(ex).fixed = true))
                end
            end
        end
        return nothing
    end

    function recursive_search(ex)
        if isa(ex, Expr)
            if  (ex.head == :(=))       &&
                (length(ex.args) == 2)
                dd = extract_param_expr(ex.args[1])
                if !isnothing(dd)
                    push!(mname, dd[1])
                    push!(check, dd[2])
                    push!(fixed, dd[3])
                end
                for i in 2:length(ex.args)
                    recursive_search(ex.args[i])
                end
            else
                dd = extract_param_expr(ex)
                if !isnothing(dd)
                    push!(mname, dd[1])
                    push!(check, dd[2])
                else
                    for i in 1:length(ex.args)
                        recursive_search(ex.args[i])
                    end
                end
            end
        end
    end

    recursive_search(expr)
    @assert length(unique(mname)) == 1 "More than one Model or MultiModel name was found, but only one was expected: $mname"
    mname = mname[1]
    check = unique(check);  # for i in 1:length(check);  println(check[i]);  end
    fixed = unique(fixed);  # for i in 1:length(fixed);  println(fixed[i]);  end
    check = Expr(:&&, check...)
    fixed = Expr(:block, fixed...)
    out = prettify(:(
        if $check
        $fixed
        GFit.patch!($mname, GFit.@exprfunc $mname -> $expr)
        end
    ))
    # println(out)
    return esc(out)
end
