struct SumReducer <: AbstractComponent
    list::Vector{Symbol}
    SumReducer() = new(Symbol[])
    SumReducer(args::AbstractSet{Symbol}) = new(collect(args))
    SumReducer(args::Vector{Symbol}) = new(args)
    SumReducer(args::Vararg{Symbol}) = new([args...])
end

dependencies(comp::SumReducer) = comp.list


function evaluate!(ceval::CompEval{SumReducer, <: AbstractDomain})
    ceval.buffer .= 0.
end

function evaluate!(ceval::CompEval{SumReducer, <: AbstractDomain},
                   args)
    ceval.buffer .= 0.
    for i in 1:length(args)
        ceval.buffer .+= args[i]
    end
end
