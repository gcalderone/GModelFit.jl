struct SimplePar <: AbstractComponent
    par::Parameter
    SimplePar(val::Number) = new(Parameter(val))
end

compeval_cdata(comp::SimplePar, domain::AbstractDomain) = nothing
compeval_array(comp::SimplePar, domain::AbstractDomain) = [NaN]

function evaluate(c::CompEval{SimplePar, T}, par) where T <: AbstractLinearDomain
    c.buffer[1] = par
end
