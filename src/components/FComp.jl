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


# We need to implement two evaluate! methods, with/without deps argument respectively
function evaluate!(ceval::CompEval{FComp, <: AbstractDomain},
                   params...)
    ceval.buffer .= getfield(ceval.comp, :func)(params...)
end

function evaluate!(ceval::CompEval{FComp, <: AbstractDomain},
                   deps, params...)
    ceval.buffer .= getfield(ceval.comp, :func)(deps..., params...)
end



# ====================================================================
struct FCompv <: GModelFit.AbstractComponent
    funct::Function
    deps::Vector{Symbol}
    params::OrderedDict{Symbol, Parameter}

    FCompv(funct::Function, guess::Vector{T}) where T <: Number =
        FCompv(funct, Symbol[], guess)

    function FCompv(funct::Function, deps::Vector{Symbol}, guess::Vector{T}) where T <: Number
        params = OrderedDict{Symbol, Parameter}()
        for i in 1:length(guess)
            params[Symbol(:p, i)] = Parameter(guess[i])
        end
        new(funct, deps, params)
    end
end

propertynames(comp::FCompv) = collect(keys(getfield(comp, :params)))
getproperty(comp::FCompv, key::Symbol) = getfield(comp, :params)[key]
dependencies(comp::FCompv) = getfield(comp, :deps)


# We need to implement two evaluate! methods, with/without deps argument respectively
function evaluate!(ceval::CompEval{FCompv, <: AbstractDomain},
                   params::Vararg{Float64})
    ceval.buffer .= getfield(ceval.comp, :funct)([params...])
end

function evaluate!(ceval::CompEval{FCompv, <: AbstractDomain},
                   deps, params...)
    ceval.buffer .= getfield(ceval.comp, :func)(deps..., [params...])
end
