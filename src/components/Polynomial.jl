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


function evaluate!(::Polynomial, domain::AbstractDomain{1}, output,
                   params...)
    output .= params[1]
    for deg in 1:length(params)-1
        output .+= coords(domain).^deg .* params[deg+1]
    end
end
