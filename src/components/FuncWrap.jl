struct FuncWrap <: AbstractComponent
    func::λFunct
    list::Vector{Symbol}
    hash::HashVector{Parameter}

    function FuncWrap(f::λFunct, args...)
        list = deepcopy(f.args)
        params = HashVector{Parameter}()
        for i in 1:length(f.optargs)
            @assert f.optargs[i].head == :(=)
            @assert isa(f.optargs[i].args[1], Symbol)
            @assert isa(f.optargs[i].args[2], Number)
            params[f.optargs[i].args[1]] = Parameter(f.optargs[i].args[2])
        end
        return new(f, list, params)
    end
end

dependencies(comp::FuncWrap) = getfield(comp, :list)

getproperty(comp::FuncWrap, key::Symbol) = getproperty(getfield(comp, :hash), key)

function getparams(comp::FuncWrap)
    out = OrderedDict{Symbol, Parameter}()
    for (key, val) in getfield(comp, :hash)
        out[key] = val
    end
    return out
end


function prepare!(comp::FuncWrap, domain::AbstractDomain)
    # Discard as many argumnts as the number of dimensions in the domain
    list = getfield(comp, :list)
    @info typeof(list)
    for i in 1:ndims(domain)
        deleteat!(list, 1)
    end
    fill(NaN, length(domain))
end


function evaluate!(buffer::Vector{Float64}, comp::FuncWrap, domain::AbstractDomain,
                   deps, params...)
    buffer .= getfield(comp, :func)(coords(domain)..., deps..., params...)
end
