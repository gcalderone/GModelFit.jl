struct Gaussian <: AbstractComponent
    params::OrderedDict{Symbol, Parameter}

    function Gaussian(norm::Real, center::Real, sigma::Real)
        @assert norm  > 0
        @assert sigma > 0
        p = OrderedDict{Symbol, Parameter}(
            :norm   => Parameter(norm),
            :center => Parameter(center),
            :sigma  => Parameter(sigma))
        p[:norm].low  = 0
        p[:sigma].low = 0
        return new(p)
    end

    function Gaussian(norm::Real, centerX::Real, centerY::Real, sigma::Real)
        @assert norm   > 0
        @assert sigma  > 0
        p = OrderedDict{Symbol, Parameter}(
            :norm    => Parameter(norm),
            :centerX => Parameter(centerX),
            :centerY => Parameter(centerY),
            :sigma   => Parameter(sigma))
        p[:norm].low  = 0
        p[:sigma].low = 0
        return new(p)
    end

    function Gaussian(norm::Real, centerX::Real, centerY::Real, sigmaX::Real, sigmaY::Real, angle::Real)
        @assert norm   > 0
        @assert sigmaX > 0
        @assert sigmaY > 0
        p = OrderedDict{Symbol, Parameter}(
            :norm    => Parameter(norm),
            :centerX => Parameter(centerX),
            :centerY => Parameter(centerY),
            :sigmaX  => Parameter(sigmaX),
            :sigmaY  => Parameter(sigmaY),
            :angle   => Parameter(angle))
        p[:norm].low  = 0
        p[:sigmaX].low = 0
        p[:sigmaY].low = 0
        return new(p)
    end
end

getparams(comp::Gaussian) = comp.params

# ====================================================================
function evaluate!(::Gaussian, domain::Domain{1}, output::Vector,
                   norm, center, sigma)
    X = coords(domain)
    output .= @. exp( ((X - center) / sigma)^2. / (-2.)) /
        2.5066282746310002 / sigma * norm # sqrt(2pi) = 2.5066282746310002
end


evaluate!(comp::Gaussian, domain::AbstractDomain{2}, output::AbstractArray,
          norm, centerX, centerY, sigma) = evaluate!(comp, domain, output,
                                                     norm, centerX, centerY, sigma, sigma, 0.)

function evaluate!(::Gaussian, domain::Domain{2}, output::Vector,
                   norm, centerX, centerY, sigmaX, sigmaY, angle)
    angle *= -pi / 180.
    a =  (cos(angle) / sigmaX)^2 / 2  +  (sin(angle) / sigmaY)^2 / 2
    b = -sin(2angle) / sigmaX^2  / 2  +  sin(2angle) / sigmaY^2  / 2
    c =  (sin(angle) / sigmaX)^2 / 2  +  (cos(angle) / sigmaY)^2 / 2

    x = coords(domain, 1)
    y = coords(domain, 2)

    output .= @. norm *
        exp(-(
            a * (x - centerX)^2. +
            b * (x - centerX) * (y - centerY) +
            c *                 (y - centerY)^2.
            )) / 6.283185307179586 / sigmaX / sigmaY # 2pi = 6.283185307179586
end


function evaluate!(::Gaussian, domain::CartesianDomain{2}, output::Matrix,
                   norm, centerX, centerY, sigmaX, sigmaY, angle)
    angle *= -pi / 180.
    a =  (cos(angle) / sigmaX)^2 / 2  +  (sin(angle) / sigmaY)^2 / 2
    b = -sin(2angle) / sigmaX^2  / 2  +  sin(2angle) / sigmaY^2  / 2
    c =  (sin(angle) / sigmaX)^2 / 2  +  (cos(angle) / sigmaY)^2 / 2

    x = axes(domain, 1)
    y = axes(domain, 2)

    for i in 1:length(y)
        output[:, i] .= @. norm *
            exp(-(
                a * (x - centerX)^2. +
                b * (x - centerX) * (y[i] - centerY) +
                c *                 (y[i] - centerY)^2.
                )) / 6.283185307179586 / sigmaX / sigmaY # 2pi = 6.283185307179586
    end
end
