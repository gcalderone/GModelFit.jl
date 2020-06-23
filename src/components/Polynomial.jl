mutable struct Polynomial <: AbstractComponent
    coeff::Vector{Parameter}

    function Polynomial(args...)
        coeff = [Parameter(arg) for arg in args]
        return new(coeff)
    end
end

ceval_data(domain::AbstractDomain, comp::Polynomial) = (nothing, length(domain))

function evaluate(c::CompEval{Domain_1D, Polynomial}, coeffs...)
    output .= coeffs[1]
    x = domain[1]
    for deg in 1:length(coeffs)-1
        output .+= x.^deg .* coeffs[deg+1]
    end
    return output
end
