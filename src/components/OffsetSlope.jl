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
function evaluate!(::OffsetSlope_1D, domain::Domain{1}, output::Vector,
                   offset, x0, slope)
    X = coords(domain)
    @. (output = slope * (X - x0) + offset)
end


function evaluate!(::OffsetSlope_2D, domain::Domain{2}, output::Vector,
                   offset, x0, y0, slopeX, slopeY)
    x = coords(domain, 1)
    y = coords(domain, 2)
    @. (output =
        slopeX * (x - x0) +
        slopeY * (y - y0) +
        offset)
end


function evaluate!(::OffsetSlope_2D, domain::CartesianDomain{2}, output::Matrix,
                   offset, x0, y0, slopeX, slopeY)
    x = axes(domain, 1)
    y = axes(domain, 2)
    for i in 1:length(y)
        output[:, i] .=
            @. slopeX * (x - x0) +
            slopeY * (y[i] - y0) +
            offset
    end
end
