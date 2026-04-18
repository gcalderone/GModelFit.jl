struct OffsetSlope <: AbstractComponent
    params::OrderedDict{Symbol, Parameter}

    function OffsetSlope(offset::Real, x0::Real, slope::Real)
        p = OrderedDict{Symbol, Parameter}(
            :offset => Parameter(offset),
            :x0     => Parameter(x0),
            :slope  => Parameter(slope))
        p[:x0].fixed = true
        return new(p)
    end

    function OffsetSlope(offset::Real, x0::Real, y0::Real, slopeX::Real, slopeY::Real)
        p = OrderedDict{Symbol, Parameter}(
            :offset => Parameter(offset),
            :x0     => Parameter(x0),
            :y0     => Parameter(y0),
            :slopeX => Parameter(slopeX),
            :slopeY => Parameter(slopeY))
        p[:x0].fixed = true
        p[:y0].fixed = true
        return new(p)
    end

    function OffsetSlope(offset::Real, x0::Real, y0::Real, z0::Real, slopeX::Real, slopeY::Real, slopeZ::Real)
        p = OrderedDict{Symbol, Parameter}(
            :offset => Parameter(offset),
            :x0     => Parameter(x0),
            :y0     => Parameter(y0),
            :z0     => Parameter(z0),
            :slopeX => Parameter(slopeX),
            :slopeY => Parameter(slopeY),
            :slopeZ => Parameter(slopeZ))
        p[:x0].fixed = true
        p[:y0].fixed = true
        p[:z0].fixed = true
        return new(p)
    end
end
getparams(comp::OffsetSlope) = comp.params

# ====================================================================
function evaluate!(::OffsetSlope, domain::Domain{1}, output::Vector,
                   offset, x0, slope)
    X = coords(domain)
    output .= @. slope * (X - x0) + offset
end


function evaluate!(::OffsetSlope, domain::Domain{2}, output::Vector,
                   offset, x0, y0, slopeX, slopeY)
    x = coords(domain, 1)
    y = coords(domain, 2)
    output .= @. slopeX * (x - x0)  +  slopeY * (y - y0)  +  offset
end

function evaluate!(::OffsetSlope, domain::CartesianDomain{2}, output::Matrix,
                   offset, x0, y0, slopeX, slopeY)
    x = axes(domain, 1)
    y = axes(domain, 2)
    for i in 1:length(y)
        output[:, i] .= @. slopeX * (x - x0)  +  slopeY * (y[i] - y0)  +  offset
    end
end

function evaluate!(::OffsetSlope, domain::Domain{3}, output::Vector,
                   offset, x0, y0, z0, slopeX, slopeY, slopeZ)
    x = coords(domain, 1)
    y = coords(domain, 2)
    z = coords(domain, 3)
    output .= @. slopeX * (x - x0)  +  slopeY * (y - y0)  +  slopeZ * (z - z0)  +  offset
end

function evaluate!(::OffsetSlope, domain::CartesianDomain{3}, output::Array,
                   offset, x0, y0, z0, slopeX, slopeY, slopeZ)
    x = axes(domain, 1)
    y = axes(domain, 2)
    z = axes(domain, 3)
    for i in 1:length(z)
        for j in 1:length(y)
            output[:, j, i] .= @. slopeX * (x - x0)  +  slopeY * (y[j] - y0)  +  slopeZ * (z[i] - z0)  +  offset
        end
    end
end
