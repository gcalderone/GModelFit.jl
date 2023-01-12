mutable struct Polynomial <: AbstractComponent
    coeff::HashVector{Parameter}

    function Polynomial(args...)
        coeff = HashVector{Parameter}()
        for i in 1:length(args)
            coeff[Symbol.(:p, i)] = Parameter(args[i])
        end
        return new(coeff)
    end
end

# Allow access to parameters as `comp.parname`
propertynames(comp::Polynomial) = collect(keys(getfield(getfield(comp, :coeff), :dict)))
getproperty(comp::Polynomial, key::Symbol) = getproperty(getfield(comp, :coeff), key)

function getparams(comp::Polynomial)
    out = OrderedDict{Symbol, Parameter}()
    for (key, val) in getfield(comp, :coeff)
        out[key] = val
    end
    return out
end


function evaluate!(buffer::Vector{Float64}, comp::Polynomial, x::AbstractDomain{1},
                   coeffs...)
    buffer .= coeffs[1]
    for deg in 1:length(coeffs)-1
        buffer .+= coords(x).^deg .* coeffs[deg+1]
    end
end
