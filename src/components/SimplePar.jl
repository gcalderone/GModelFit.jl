struct SimplePar <: AbstractComponent
    par::Parameter
    SimplePar(val::Number) = new(Parameter(val))
end

prepare!(comp::SimplePar, domain::AbstractDomain) = [NaN]
getproperty( comp::SimplePar, key::Symbol) = getproperty(getfield(comp, :par), key)
setproperty!(comp::SimplePar, key::Symbol, v::Real) = setproperty!(getfield(comp, :par), key, v)

function evaluate!(buffer::Vector{Float64}, comp::SimplePar, domain::AbstractDomain,
                   par)
    buffer[1] = par
end
