# ====================================================================
# Component structure
mutable struct OffsetSlope_1D <: AbstractComponent
    offset::Parameter
    x0::Parameter
    slope::Parameter

    function OffsetSlope_1D(offset::Number, x0::Number, slope::Number)
        out = new(Parameter(offset), Parameter(x0), Parameter(slope))
        out.x0.free = false
        return out
    end
end


mutable struct OffsetSlope_2D <: AbstractComponent
    offset::Parameter
    x0::Parameter
    y0::Parameter
    slopeX::Parameter
    slopeY::Parameter

    function OffsetSlope_2D(offset::Number, x0::Number, y0::Number, slopeX::Number, slopeY::Number)
        out = new(Parameter(offset), Parameter(x0), Parameter(y0), Parameter(slopeX), Parameter(slopeY))
        out.x0.free = false
        out.y0.free = false
        return out
    end
end

OffsetSlope(offset, x0, slope) = OffsetSlope_1D(offset, x0, slope)
OffsetSlope(offset, x0, y0, slopeX, slopeY) = OffsetSlope_2D(offset, x0, y0, slopeX, slopeY)

ceval_data(domain::Domain_1D, comp::OffsetSlope_1D) = (nothing, length(domain))
ceval_data(domain::Domain_2D, comp::OffsetSlope_2D) = (nothing, length(domain))


# ====================================================================
# Evaluate component 
function evaluate(c::CompEval{Domain_1D, OffsetSlope_1D},
                  offset, x0, slope)
    @. (c.eval = slope * (c.domain[1] - x0) + offset)
end


function evaluate(c::CompEval{Domain_2D, OffsetSlope_2D},
                   offset, x0, y0, slopeX, slopeY)
    x = c.domain[1]
    y = c.domain[2]
    @. (c.eval = 
        slopeX * (x - x0) + 
        slopeY * (y - y0) +
        offset)
end
