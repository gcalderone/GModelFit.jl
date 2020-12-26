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
compeval_cdata(comp::Lorentzian_1D, domain::Domain{1}) = nothing
compeval_cdata(comp::Lorentzian_2D, domain::Domain{2}) = nothing
compeval_array(comp::Lorentzian_1D, domain::Domain{1}) = fill(NaN, length(domain))
compeval_array(comp::Lorentzian_2D, domain::Domain{2}) = fill(NaN, length(domain))


# ====================================================================
# Evaluate component 
function evaluate(ce::CompEval{Lorentzian_1D, 1},
                  norm, center, fwhm)
    x = ce.domain[1]

    @. (c.eval = norm /
        (1. +
         ((x - center) / fwhm)^2.
        ))
end


function evaluate(ce::CompEval{Lorentzian_2D, 2},
                  norm, centerX, centerY, fwhmX, fwhmY)
    x = ce.domain[1]
    y = ce.domain[2]

    @. (c.eval = norm /
        (1. +
         ((x - centerX) / fwhmX)^2. +
         ((y - centerY) / fwhmY)^2.
         ))
end
