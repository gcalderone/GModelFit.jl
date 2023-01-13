struct λComp <: AbstractComponent
    func::λFunct
    list::Vector{Symbol}
    params::OrderedDict{Symbol, Parameter}

    function λComp(f::λFunct, args...)
        list = deepcopy(f.args)
        params = OrderedDict{Symbol, Parameter}()
        for i in 1:length(f.optargs)
            @assert f.optargs[1].head != :... "Splat not allowed in LComp"
            @assert f.optargs[i].head == :(=)
            @assert isa(f.optargs[i].args[1], Symbol)
            @assert isa(f.optargs[i].args[2], Number)
            params[f.optargs[i].args[1]] = Parameter(f.optargs[i].args[2])
        end
        return new(f, list, params)
    end
end

# Allow access to parameters as `comp.parname`
propertynames(comp::λComp) = collect(keys(getfield(comp, :params)))
getproperty(comp::λComp, key::Symbol) = getfield(comp, :params)[key]

dependencies(comp::λComp) = getfield(comp, :list)


function prepare!(comp::λComp, domain::AbstractDomain)
    # Discard as many arguments as the number of dimensions in the domain
    list = getfield(comp, :list)
    for i in 1:ndims(domain)
        deleteat!(list, 1)
    end
    fill(NaN, length(domain))
end


# We need to implement two evaluate! methods, with/without deps argument respectively
function evaluate!(buffer::Vector{Float64}, comp::λComp, domain::AbstractDomain,
                   params...)
    buffer .= getfield(comp, :func)(domain..., params...)
end

function evaluate!(buffer::Vector{Float64}, comp::λComp, domain::AbstractDomain,
                   deps, params...)
    buffer .= getfield(comp, :func)(domain..., deps..., params...)
end
