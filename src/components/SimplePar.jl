struct SimplePar <: AbstractComponent
    par::Parameter
    SimplePar(val::Number) = new(Parameter(val))
end

ceval_data(domain::AbstractDomain, comp::SimplePar) = nothing, 1

function evaluate(c::CompEval{T, SimplePar}, par) where T <: AbstractDomain
    c.eval[1] = par
end
