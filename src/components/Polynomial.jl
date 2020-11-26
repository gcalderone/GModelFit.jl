mutable struct Polynomial <: AbstractComponent
    coeff::Vector{Parameter}

    function Polynomial(args...)
        coeff = [Parameter(arg) for arg in args]
        return new(coeff)
    end
end

compeval_cdata(comp::Polynomial, domain::AbstractDomain, ) = nothing
compeval_array(comp::Polynomial, domain::AbstractDomain) = fill(NaN, length(domain))

function evaluate(c::CompEval{Polynomial, Domain_1D}, coeffs...)
    c.eval .= coeffs[1]
    x = domain[1]
    for deg in 1:length(coeffs)-1
        c.eval .+= x.^deg .* coeffs[deg+1]
    end
end
