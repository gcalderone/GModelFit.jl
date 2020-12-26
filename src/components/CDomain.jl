struct CDomain <: AbstractComponent
    dim::Int
    CDomain(dim::Int=1) = new(dim)
end

compeval_cdata(comp::CDomain, domain::Domain) = nothing
compeval_array(comp::CDomain, domain::Domain) = deepcopy(flatten(domain)[comp.dim])

evaluate(c::CompEval{CDomain, N}) where N = nothing
