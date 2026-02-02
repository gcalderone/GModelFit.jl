struct FComp <: AbstractComponent
    func::Function
    deps::Vector{Symbol}
    params::OrderedDict{Symbol, Parameter}

    function FComp(f::FunctDesc)
        deps = deepcopy(f.args)
        @assert length(deps) >= 1 "No argument provided for domain coordinates"
        params = OrderedDict{Symbol, Parameter}()
        for i in 1:length(f.optargs)
            @assert f.optargs[1].head != :... "Splat not allowed"
            @assert f.optargs[i].head == :(=)
            @assert isa(f.optargs[i].args[1], Symbol)
            @assert isa(f.optargs[i].args[2], Number)
            params[f.optargs[i].args[1]] = Parameter(f.optargs[i].args[2])
        end
        return new(f.funct, deps, params)
    end

    function FComp(funct::Function, deps=Symbol[]; kws...)
        @assert length(deps) >= 1 "No argument provided for domain coordinates"
        params = OrderedDict{Symbol, Parameter}()
        for (name, val) in kws
            @assert isa(name, Symbol)
            @assert isa(val , Number)
            params[name] = Parameter(val)
        end
        return new(funct, deps, params)
    end
end

# Allow access to parameters as `comp.parname`
propertynames(comp::FComp) = collect(keys(getfield(comp, :params)))
getproperty(comp::FComp, key::Symbol) = getfield(comp, :params)[key]
dependencies(comp::FComp) = getfield(comp, :deps)


function evaluate!(comp::FComp, domain::Domain, output::Vector, params...)
    f = getfield(comp, :func)
    output .= f(domain.axes..., params...)
end

function evaluate!(comp::FComp, domain::Domain, output::Vector, deps::Vector, params...)
    f = getfield(comp, :func)
    output .= f(domain.axes..., deps..., params...)
end

function evaluate!(comp::FComp, domain::CartesianDomain, output::Array, params...)
    f = getfield(comp, :func)
    for I in CartesianIndices(tuple(length.(domain.axes[2:end])...))
        X = getindex.(domain.axes[2:end], Tuple(I))
        output[:, I] .= f(domain.axes[1], X..., params...)
    end
end

function evaluate!(comp::FComp, domain::CartesianDomain, output::Array, deps::Vector, params...)
    f = getfield(comp, :func)
    for I in CartesianIndices(tuple(length.(domain.axes[2:end])...))
        X = getindex.(domain.axes[2:end], Tuple(I))
        d = [view(deps[i], :, I) for i in 1:length(deps)]
        output[:, I] .= f(domain.axes[1], X..., d..., params...)
    end
end
