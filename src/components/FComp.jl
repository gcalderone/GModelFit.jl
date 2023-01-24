struct FComp <: AbstractComponent
    func::Function
    deps::Vector{Symbol}
    params::OrderedDict{Symbol, Parameter}

    function FComp(f::λFunct)
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


# We need to implement two evaluate! methods, with/without deps argument respectively
function evaluate!(buffer::Vector{Float64}, comp::FComp, domain::AbstractDomain,
                   params...)
    buffer .= getfield(comp, :func)(params...)
end

function evaluate!(buffer::Vector{Float64}, comp::FComp, domain::AbstractDomain,
                   deps, params...)
    buffer .= getfield(comp, :func)(deps..., params...)
end



# ====================================================================
struct FCompv <: GFit.AbstractComponent
    funct::Function
    deps::Vector{Symbol}
    params::Vector{Parameter}

    FCompv(funct::Function, guess::Vector{T}) where T <: Number =
        new(funct, Symbol[], Parameter.(guess))
    FCompv(funct::Function, deps::Vector{Symbol}, guess::Vector{T}) where T <: Number =
        new(funct, deps, Parameter.(guess))
end

dependencies(comp::FCompv) = getfield(comp, :deps)


# We need to implement two evaluate! methods, with/without deps argument respectively
function evaluate!(buffer::Vector{Float64}, comp::FCompv, x::AbstractDomain,
                   params::Vararg{Float64})
    buffer .= getfield(comp, :funct)([params...])
end

function evaluate!(buffer::Vector{Float64}, comp::FCompv, domain::AbstractDomain,
                   deps, params...)
    buffer .= getfield(comp, :func)(deps..., [params...])
end
