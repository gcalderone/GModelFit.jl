struct CDomain <: AbstractComponent
    dim::Int
    CDomain(dim::Int=1) = new(dim)
end

compeval_array(comp::CDomain, domain::AbstractDomain) = deepcopy(flatten(domain)[comp.dim])

evaluate(c::CompEval{CDomain, TDomain}) where TDomain = nothing
