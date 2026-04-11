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
            @assert isa(f.optargs[i].args[2], Real)
            params[f.optargs[i].args[1]] = Parameter(f.optargs[i].args[2])
        end
        return new(f.funct, deps, params)
    end

    function FComp(funct::Function, deps=Symbol[]; kws...)
        @assert length(deps) >= 1 "No argument provided for domain coordinates"
        params = OrderedDict{Symbol, Parameter}()
        for (name, val) in kws
            @assert isa(name, Symbol)
            @assert isa(val , Real)
            params[name] = Parameter(val)
        end
        return new(funct, deps, params)
    end
end

dependencies(comp::FComp) = comp.deps
getparams(comp::FComp) = comp.params


# ====================================================================
evaluate!(comp::FComp, domain::Domain, output::Vector, params...) =
    output .= comp.func(domain.axes..., params...)

evaluate!(comp::FComp, domain::Domain, output::Vector, deps::Vector, params...) =
    output .= comp.func(domain.axes..., deps..., params...)

function evaluate!(comp::FComp, domain::CartesianDomain, output::Array, params...)
    for I in CartesianIndices(tuple(length.(domain.axes[2:end])...))
        X = getindex.(domain.axes[2:end], Tuple(I))
        output[:, I] .= comp.func(domain.axes[1], X..., params...)
    end
end

function evaluate!(comp::FComp, domain::CartesianDomain, output::Array, deps::Vector, params...)
    for I in CartesianIndices(tuple(length.(domain.axes[2:end])...))
        X = getindex.(domain.axes[2:end], Tuple(I))
        d = [view(deps[i], :, I) for i in 1:length(deps)]
        output[:, I] .= comp.func(domain.axes[1], X..., d..., params...)
    end
end
