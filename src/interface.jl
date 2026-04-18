function getparams(model::Union{Model, ModelSnapshot})
    out = OrderedDict{NTuple{2, Symbol}, Parameter}()
    for (cname, comp) in model.comps
        for (pname, par) in getparams(comp)
            @assert !isnan(par.low)
            @assert !isnan(par.high)
            @assert isfinite(par.val)
            @assert par.low < par.high
            @assert par.val >= par.low  "Value for [:$(cname), :$(pname)] is smaller than the minimum allowed value"
            @assert par.val <= par.high "Value for [:$(cname), :$(pname)] is larger than the maximum allowed value"
            @assert isnothing(par.patch)  ||  isnothing(par.reparam) "Parameter [:$(cname), :$(pname)] can either be patched to another value, or re-parametrized (reparam), but not both."
            par.actually_fixed = par.fixed  ||  !isnothing(par.patch)  ||  isfrozen(model, cname)
            out[(cname, pname)] = par
        end
    end
    return out
end

function getparams(ms::Union{ModelSet, ModelSetSnapshot})
    out = OrderedDict{NTuple{3, Symbol}, Parameter}()
    for (mname, model) in ms.dict
        for (k, par) in getparams(model)
            out[tuple(mname, k...)] = par
        end
    end
    return out
end

length(  v::Union{AbstractComponent, ComponentSnapshot})                          =  length(getparams(v))
iterate( v::Union{AbstractComponent, ComponentSnapshot}, args...)                 = iterate(getparams(v), args...)
keys(    v::Union{AbstractComponent, ComponentSnapshot})                          =    keys(getparams(v))
haskey(  v::Union{AbstractComponent, ComponentSnapshot}, key::Symbol)             =  haskey(getparams(v), key)
getindex(v::Union{AbstractComponent, ComponentSnapshot}, key::Symbol)             =         getparams(v)[ key]

length(  v::Union{Model            , ModelSnapshot}    )                          =  length(v.comps)
iterate( v::Union{Model            , ModelSnapshot}    , args...)                 = iterate(v.comps, args...)
keys(    v::Union{Model            , ModelSnapshot}    )                          =    keys(v.comps)
haskey(  v::Union{Model            , ModelSnapshot}    , key::Symbol)             =  haskey(v.comps, key)
getindex(v::Union{Model            , ModelSnapshot}    , key::Symbol)             =         v.comps[ key]
getindex(v::Union{Model            , ModelSnapshot}    , key::Vararg{Symbol, 2})  =         v.comps[ key[1]][key[2]]

length(  v::Union{ModelSet         , ModelSetSnapshot} )                          = length( v.dict)
iterate( v::Union{ModelSet         , ModelSetSnapshot} , args...)                 = iterate(v.dict, args...)
keys(    v::Union{ModelSet         , ModelSetSnapshot} )                          =    keys(v.dict)
haskey(  v::Union{ModelSet         , ModelSetSnapshot} , key::Symbol)             =  haskey(v.dict, key)
getindex(v::Union{ModelSet         , ModelSetSnapshot} , key::Symbol)             =         v.dict[ key]
getindex(v::Union{ModelSet         , ModelSetSnapshot} , key::Vararg{Symbol, 2})  =         v.dict[ key[1]][key[2]]
getindex(v::Union{ModelSet         , ModelSetSnapshot} , key::Vararg{Symbol, 3})  =         v.dict[ key[1]][key[2]][key[3]]
