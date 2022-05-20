struct FuncWrap <: AbstractComponent
    func::λFunct
    params::Vector{Parameter}
    hash::HashVector{Parameter}

    function FuncWrap(f::λFunct, args...)
        params = HashVector{Parameter}()
        for i in 1:length(f.optargs)
            @assert f.optargs[i].head == :(=)
            @assert isa(f.optargs[i].args[1], Symbol)
            @assert isa(f.optargs[i].args[2], Number)
            push!(params, f.optargs[i].args[1], Parameter(f.optargs[i].args[2]))
        end
        return new(f, getfield(params, :data), params)
    end
end

getproperty(comp::FuncWrap, key::Symbol) = getproperty(getfield(comp, :hash), key)

function prepare!(comp::FuncWrap, domain::AbstractDomain)
    @assert length(getfield(comp, :func).args) == ndims(domain)
    fill(NaN, length(domain))
end


function evaluate!(buffer, comp::FuncWrap, domain::AbstractDomain,
                   params...)
    buffer .= getfield(comp, :func)(coords(domain)..., params...)
end
