struct FComp <: AbstractComponent
    func::Function
    deps::Vector{Symbol}
    params::OrderedDict{Symbol, Parameter}

    function FComp(f::FunctDesc)
        deps = deepcopy(f.args)
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


function evaluate!(comp::FComp, domain::AbstractDomain, output::Vector,
                   deps, params...)
    if length(deps) > 0
        output .= getfield(comp, :func)(deps..., params...)
    else
        output .= getfield(comp, :func)(params...)
    end
end
