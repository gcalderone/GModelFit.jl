struct FComp <: GFit.AbstractComponent
    funct::Function
    params::Vector{Parameter}
    FComp(funct::Function, guess::Vector{T}) where T <: Number =
        new(funct, Parameter.(guess))
end

function prepare!(comp::FComp, domain::AbstractDomain)
    reshape(comp.funct([getfield.(comp.params, :val)...]), :)
end

function evaluate!(buffer::Vector{Float64}, comp::FComp, x::AbstractDomain,
                   params::Vararg{Float64})
    buffer .= getfield(comp, :funct)([params...])
end
