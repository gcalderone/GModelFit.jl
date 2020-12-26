struct SimplePar <: AbstractComponent
    par::Parameter
    SimplePar(val::Number) = new(Parameter(val))
end

compeval_array(comp::SimplePar, domain::AbstractDomain) = [NaN]

function evaluate(buffer, comp::SimplePar, domain::AbstractDomain,
                  par)
    buffer[1] = par
end
