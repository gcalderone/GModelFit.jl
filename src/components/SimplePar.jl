struct SimplePar <: AbstractComponent
    par::Parameter
    SimplePar(val::Number) = new(Parameter(val))
end

prepare!(comp::SimplePar, domain::AbstractDomain) = [NaN]

function evaluate!(buffer::Vector{Float64}, comp::SimplePar, domain::AbstractDomain,
                   par)
    buffer[1] = par
end
