struct CDomain <: AbstractComponent
    dim::Int
    CDomain(dim::Int=1) = new(dim)
end

ceval_data(domain::AbstractDomain, comp::CDomain) = true, length(domain)

function evaluate(c::CompEval{T, CDomain}) where T <: AbstractLinearDomain
    if c.cdata
        c.eval .= c.domain[c.comp.dim]
        c.cdata = false
    end
    nothing
end
