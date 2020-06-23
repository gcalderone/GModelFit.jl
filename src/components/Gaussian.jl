# ====================================================================
# Component structure
#
mutable struct Gaussian_1D <: AbstractComponent
    norm::Parameter
    center::Parameter
    sigma::Parameter

    function Gaussian_1D(norm::Number, center::Number, sigma::Number)
        @assert norm  > 0
        @assert sigma > 0
        
        out = new(Parameter(norm), Parameter(center), Parameter(sigma))
        out.norm.low = 0
        out.sigma.low = 0
        return out
    end
end


mutable struct Gaussian_2D <: AbstractComponent
    norm::Parameter
    centerX::Parameter
    centerY::Parameter
    sigmaX::Parameter
    sigmaY::Parameter
    angle::Parameter

    function Gaussian_2D(norm::Number, centerX::Number, centerY::Number, sigmaX::Number, sigmaY::Number, angle::Number)
        @assert norm   > 0
        @assert sigmaX > 0
        @assert sigmaY > 0

        out = new(Parameter(norm), Parameter(centerX), Parameter(centerY), Parameter(sigmaX), Parameter(sigmaY), Parameter(angle))
        out.norm.low = 0
        out.sigmaX.low = 0
        out.sigmaY.low = 0
        return out
    end
end

Gaussian(norm, center, sigma) = Gaussian_1D(norm, center, sigma)
Gaussian(norm, centerX, centerY, sigmaX, sigmaY, angle) = Gaussian_2D(norm, centerX, centerY, sigmaX, sigmaY, angle)
function Gaussian(norm, centerX, centerY, sigma)
    out = Gaussian_2D(norm, centerX, centerY, sigma, sigma, 0.)
    out.sigmaX.free = false
    out.sigmaX.expr = "this_sigmaY"
    out.angle.free = false
    return out
end


# ====================================================================
# Prepare component `cdata`
mutable struct Gaussian_cdata
    ix::Vector{Int}
    iy::Vector{Int}
    Gaussian_cdata() = new(Vector{Int}(), Vector{Int}())
end

ceval_data(domain::Domain_1D, comp::Gaussian_1D) = (Gaussian_cdata(), length(domain))
ceval_data(domain::Domain_2D, comp::Gaussian_2D) = (Gaussian_cdata(), length(domain))


# ====================================================================
# Evaluate component 
function evaluate(c::CompEval{Domain_1D, Gaussian_1D},
                  norm, center, sigma)
    # TODO: optimize using cdata
    x = c.domain[1]
    @. (c.eval = exp( ((x - center) / sigma)^2. / (-2.)) / 
        2.5066282746310002 / sigma * norm) # sqrt(2pi) = 2.5066282746310002
end


function evaluate(c::CompEval{Domain_2D, Gaussian_2D},
                   norm, centerX, centerY, sigmaX, sigmaY, angle)
    angle *= -pi / 180.
    a =  (cos(angle) / sigmaX)^2 / 2  +  (sin(angle) / sigmaY)^2 / 2
    b = -sin(2angle) / sigmaX^2  / 2  +  sin(2angle) / sigmaY^2  / 2
    c =  (sin(angle) / sigmaX)^2 / 2  +  (cos(angle) / sigmaY)^2 / 2

    # TODO: optimize using cdata
    x = c.domain[1]
    y = c.domain[2]
    
    @. (c.eval = norm *
        exp(
            -(
                a * (x - centerX)^2. +
                b * (x - centerX) * (y - centerY) +
                c *                 (y - centerY)^2.
            )
        ) / 6.283185307179586 / sigmaX / sigmaY) # 2pi = 6.283185307179586
end
