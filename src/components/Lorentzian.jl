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
# Evaluate component 
function evaluate!(::Lorentzian_1D, domain::AbstractDomain{1}, output,
                   norm, center, fwhm)
    X = coords(domain)
    @. (output = norm /
        (1. +
         ((X - center) / fwhm)^2.
        ))
end


function evaluate!(::Lorentzian_2D, domain::AbstractDomain{2}, output,
                   norm, centerX, centerY, fwhmX, fwhmY)
    x = coords(domain, 1)
    y = coords(domain, 2)

    @. (output = norm /
        (1. +
         ((x - centerX) / fwhmX)^2. +
         ((y - centerY) / fwhmY)^2.
         ))
end
