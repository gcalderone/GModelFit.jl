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

ceval_data(domain::AbstractDomain, comp::FuncWrap) = (nothing, length(domain))

function evaluate(c::CompEval{T, FuncWrap}, params...) where T <: AbstractLinearDomain
    if ndims(c.domain) == 1
        out = c.comp.func(c.domain[1], params...)
    else
        dom = [c.domain[i] for i in 1:ndims(c.domain)]
        out = c.comp.func(dom..., params...)
    end
    if length(c.eval) == 0
        append!(c.eval, out)
    else
        c.eval .= out
    end
end
