struct Lorentzian <: AbstractComponent
    params::OrderedDict{Symbol, Parameter}

    function Lorentzian(norm::Real, center::Real, fwhm::Real)
        @assert norm  > 0
        @assert fwhm > 0
        p = OrderedDict{Symbol, Parameter}(
            :norm   => Parameter(norm),
            :center => Parameter(center),
            :fwhm   => Parameter(fwhm))
        p[:norm].low  = 0
        p[:fwhm].low = 0
        return new(p)
    end

    function Lorentzian(norm::Real, centerX::Real, centerY::Real, fwhmX::Real, fwhmY::Real)
        @assert norm  > 0
        @assert fwhmX > 0
        @assert fwhmY > 0
        p = OrderedDict{Symbol, Parameter}(
            :norm    => Parameter(norm),
            :centerX => Parameter(centerX),
            :centerY => Parameter(centerY),
            :fwhmX   => Parameter(fwhmX),
            :fwhmY   => Parameter(fwhmY))
        p[:norm].low  = 0
        p[:fwhmX].low = 0
        p[:fwhmY].low = 0
        return new(p)
    end
end

getparams(comp::Lorentzian) = comp.params

# ====================================================================
function evaluate!(::Lorentzian, domain::Domain{1}, output::Vector,
                   norm, center, fwhm)
    X = coords(domain)
    output .= @. norm / (1. + ((X - center) / fwhm)^2.)
end


function evaluate!(::Lorentzian, domain::Domain{2}, output::Vector,
                   norm, centerX, centerY, fwhmX, fwhmY)
    x = coords(domain, 1)
    y = coords(domain, 2)
    output .= @. norm /
        (1. +
        ((x - centerX) / fwhmX)^2. +
        ((y - centerY) / fwhmY)^2.)
end


function evaluate!(::Lorentzian, domain::CartesianDomain{2}, output::Matrix,
                   norm, centerX, centerY, fwhmX, fwhmY)
    x = axes(domain, 1)
    y = axes(domain, 2)
    for i in 1:length(y)
        output[:, i] .= @. norm /
            (1. +
            ((x    - centerX) / fwhmX)^2. +
            ((y[i] - centerY) / fwhmY)^2.)
    end
end
