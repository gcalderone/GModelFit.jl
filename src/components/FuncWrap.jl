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

compeval_cdata(comp::FuncWrap, domain::AbstractDomain) = nothing
compeval_array(comp::FuncWrap, domain::AbstractDomain) = fill(NaN, length(domain))

function evaluate(c::CompEval{FuncWrap, T}, params...) where T <: AbstractLinearDomain
    if ndims(c.domain) == 1
        c.buffer .= c.comp.func(c.domain[1], params...)
    else
        dom = [c.domain[i] for i in 1:ndims(c.domain)]
        c.buffer .= c.comp.func(dom..., params...)
    end
    nothing
end
