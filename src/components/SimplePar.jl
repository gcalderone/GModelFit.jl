struct SimplePar <: AbstractComponent
    par::Parameter
    SimplePar(val::Number) = new(Parameter(val))
end

compeval_cdata(comp::SimplePar, domain::Domain) = nothing
compeval_array(comp::SimplePar, domain::Domain) = [NaN]

function evaluate(c::CompEval{SimplePar, N}, par) where N
    c.buffer[1] = par
end
