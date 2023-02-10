mutable struct Polynomial <: AbstractComponent
    params::OrderedDict{Symbol, Parameter}
    function Polynomial(args...)
        params = OrderedDict{Symbol, Parameter}()
        for i in 1:length(args)
            params[Symbol(:p, (i-1))] = Parameter(args[i])
        end
        new(params)
    end
end


function evaluate!(buffer::Vector{Float64}, comp::Polynomial, x::AbstractDomain{1},
                   params...)
    buffer .= params[1]
    for deg in 1:length(params)-1
        buffer .+= coords(x).^deg .* params[deg+1]
    end
end
