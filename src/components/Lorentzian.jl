# ====================================================================
# Component structure
#
mutable struct Lorentzian_1D <: AbstractComponent
    norm::Parameter
    center::Parameter
    fwhm::Parameter

    function Lorentzian_1D(norm::Number, center::Number, fwhm::Number)
        @assert norm  > 0
        @assert fwhm > 0
        
        out = new(Parameter(norm), Parameter(center), Parameter(fwhm))
        out.norm.low = 0
        out.fwhm.low = 0
        return out
    end
end


mutable struct Lorentzian_2D <: AbstractComponent
    norm::Parameter
    centerX::Parameter
    centerY::Parameter
    fwhmX::Parameter
    fwhmY::Parameter

    function Lorentzian_2D(norm::Number, centerX::Number, centerY::Number, fwhmX::Number, fwhmY::Number)
        @assert norm   > 0
        @assert fwhmX > 0
        @assert fwhmY > 0

        out = new(Parameter(norm), Parameter(centerX), Parameter(centerY), Parameter(fwhmX), Parameter(fwhmY))
        out.norm.low = 0
        out.fwhmX.low = 0
        out.fwhmY.low = 0

        return out
    end
end

Lorentzian(norm, center, fwhm) = Lorentzian_1D(norm, center, fwhm)
Lorentzian(norm, centerX, centerY, fwhmX, fwhmY) = Lorentzian_2D(norm, centerX, centerY, fwhmX, fwhmY)


# ====================================================================
# Allocate the output array and the component `cdata`
struct Lorentzian_cdata <: AbstractComponentData; end

ceval_data(domain::Domain_1D, comp::Lorentzian_1D) = (nothing, length(domain))
ceval_data(domain::Domain_2D, comp::Lorentzian_2D) = (nothing, length(domain))


# ====================================================================
# Evaluate component 
function evaluate(c::CompEval{Domain_1D, Lorentzian_1D},
                  norm, center, fwhm)
    x = domain[1]

    @. (output = norm /
        (1. +
         ((x - center) / fwhm)^2.
        ))
    return output
end


function evaluate(c::CompEval{Domain_2D, Lorentzian_2D},
                  norm, centerX, centerY, fwhmX, fwhmY)
    x = domain[1]
    y = domain[2]

    @. (output = norm /
        (1. +
         ((x - centerX) / fwhmX)^2. +
         ((y - centerY) / fwhmY)^2.
         ))

    return output
end
