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


function evaluate!(ceval::CompEval{Polynomial, <: AbstractDomain{1}},
                   params...)
    ceval.buffer .= params[1]
    for deg in 1:length(params)-1
        ceval.buffer .+= coords(ceval.domain).^deg .* params[deg+1]
    end
end
