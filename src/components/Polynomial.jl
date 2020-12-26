mutable struct Polynomial <: AbstractComponent
    coeff::Vector{Parameter}

    function Polynomial(args...)
        coeff = [Parameter(arg) for arg in args]
        return new(coeff)
    end
end

function evaluate(buffer, comp::Polynomial, domain::Domain{1},
                  coeffs...)
    buffer .= coeffs[1]
    x = domain[1]
    for deg in 1:length(coeffs)-1
        buffer .+= x.^deg .* coeffs[deg+1]
    end
end
