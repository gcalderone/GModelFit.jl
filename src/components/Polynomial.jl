mutable struct Polynomial <: AbstractComponent
    coeff::Vector{Parameter}

    function Polynomial(args...)
        coeff = [Parameter(arg) for arg in args]
        return new(coeff)
    end
end

function evaluate(c::CompEval{Polynomial, Domain{1}}, coeffs...)
    c.eval .= coeffs[1]
    x = domain[1]
    for deg in 1:length(coeffs)-1
        c.eval .+= x.^deg .* coeffs[deg+1]
    end
end
