# ====================================================================
# Component structure
mutable struct OffsetSlope_1D <: AbstractComponent
    offset::Parameter
    x0::Parameter
    slope::Parameter

    function OffsetSlope_1D(offset::Number, x0::Number, slope::Number)
        out = new(Parameter(offset), Parameter(x0), Parameter(slope))
        out.x0.fixed = true
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
        out.x0.fixed = true
        out.y0.fixed = true
        return out
    end
end

OffsetSlope(offset, x0, slope) = OffsetSlope_1D(offset, x0, slope)
OffsetSlope(offset, x0, y0, slopeX, slopeY) = OffsetSlope_2D(offset, x0, y0, slopeX, slopeY)


# ====================================================================
# Evaluate component
function evaluate(buffer, comp::OffsetSlope_1D, domain::Domain{1},
                  offset, x0, slope)
    @. (buffer = slope * (c.domain[1] - x0) + offset)
end


function evaluate(buffer, comp::OffsetSlope_2D, domain::Domain{2},
                  offset, x0, y0, slopeX, slopeY)
    x = domain[1]
    y = domain[2]
    @. (buffer =
        slopeX * (x - x0) +
        slopeY * (y - y0) +
        offset)
end
