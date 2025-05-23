struct SumReducer <: AbstractComponent
    list::Vector{Symbol}
    SumReducer() = new(Symbol[])
    SumReducer(args::AbstractSet{Symbol}) = new(collect(args))
    SumReducer(args::Vector{Symbol}) = new(args)
    SumReducer(args::Vararg{Symbol}) = new([args...])
end

dependencies(comp::SumReducer) = comp.list


function evaluate!(::SumReducer, ::AbstractDomain, output, deps)
    output .= deps[1]
    for i in 2:length(deps)
        output .+= deps[i]
    end
end
