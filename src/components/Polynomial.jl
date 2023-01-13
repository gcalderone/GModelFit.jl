mutable struct Polynomial <: AbstractComponent
    coeff::OrderedDict{Symbol, Parameter}

    function Polynomial(args...)
        coeff = OrderedDict{Symbol, Parameter}()
        for i in 1:length(args)
            coeff[Symbol.(:p, i)] = Parameter(args[i])
        end
        return new(coeff)
    end
end

# Allow access to parameters as `comp.parname`
propertynames(comp::Polynomial) = collect(keys(getfield(comp, :coeff)))
getproperty(comp::Polynomial, key::Symbol) = getfield(comp, :coeff)[key]


function evaluate!(buffer::Vector{Float64}, comp::Polynomial, x::AbstractDomain{1},
                   coeffs...)
    buffer .= coeffs[1]
    for deg in 1:length(coeffs)-1
        buffer .+= coords(x).^deg .* coeffs[deg+1]
    end
end
