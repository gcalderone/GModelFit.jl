struct Polynomial <: AbstractComponent
    params::OrderedDict{Symbol, Parameter}
    function Polynomial(args...)
        params = OrderedDict{Symbol, Parameter}()
        for i in 1:length(args)
            params[Symbol(:p, (i-1))] = Parameter(args[i])
        end
        new(params)
    end
end
getparams(comp::Polynomial) = comp.params

function evaluate!(::Polynomial, domain::Domain{1}, output::Vector,
                   params...)
    output .= evalpoly.(coords(domain), Ref(params))
end
