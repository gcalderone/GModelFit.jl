struct FuncWrap <: AbstractComponent
    func::Function
    p::Vector{Parameter}

    function FuncWrap(func::Function, args...)
        params = Vector{Parameter}()
        for i in 1:length(args)
            push!(params, Parameter(args[i]))
        end
        return new(func, params)
    end
end


function evaluate!(buffer, comp::FuncWrap, domain::AbstractDomain,
                   params...)
    buffer .= comp.func(domain..., params...)
end
