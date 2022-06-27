struct SumReducer <: AbstractComponent
    list::Vector{Symbol}
    SumReducer() = new(Symbol[])
    SumReducer(args::Vector{Symbol}) = new(args)
    SumReducer(args::Vararg{Symbol}) = new([args...])
end

dependencies(comp::SumReducer) = comp.list


function evaluate!(buffer::Vector{Float64}, comp::SumReducer, domain::AbstractDomain)
    buffer .= 0.
end

function evaluate!(buffer::Vector{Float64}, comp::SumReducer, domain::AbstractDomain,
                   args)
    buffer .= 0.
    for i in 1:length(args)
        buffer .+= args[i]
    end
end
