struct SumReducer <: AbstractComponent
    list::Vector{Symbol}
    SumReducer() = new(Symbol[])
    SumReducer(args::AbstractSet{Symbol}) = new(collect(args))
    SumReducer(args::Vector{Symbol}) = new(args)
    SumReducer(args::Vararg{Symbol}) = new([args...])
end

dependencies(comp::SumReducer) = comp.list


function evaluate!(::SumReducer, ::AbstractDomain, output::Vector, deps)
    output .= 0.
    for i in 1:length(deps)
        output .+= deps[i]
    end
end
