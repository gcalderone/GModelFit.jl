# using DataStructures; import Base: getindex, setindex!
struct PVSet{T}
    icomp::OrderedDict{NTuple{2, Symbol}, Vector{Int}}
    index::OrderedDict{NTuple{3, Symbol}, Int}
    vec::Vector{T}
    allow_overwrite::Bool
    PVSet{T}(; allow_overwrite=false) where T =
        new{T}(OrderedDict{NTuple{2, Symbol}, Vector{Int}}(),
               OrderedDict{NTuple{3, Symbol}, Int}(),
               Vector{T}(), allow_overwrite)
end

getindex(pv::PVSet, key::Vararg{Symbol, 2}) = view(pv.vec, get(pv.icomp, key, Int[]))  # allow for components with no parameters
getindex(pv::PVSet, key::Vararg{Symbol, 3}) = pv.vec[pv.index[key]]

function setindex!(pv::PVSet{T}, v::T, key::Vararg{Symbol, 3}) where T
    i = get(pv.index, key, nothing)
    if isnothing(i)
        haskey(pv.icomp, (key[1], key[2]))  ||  (pv.icomp[(key[1], key[2])] = Vector{Int}())
        push!(pv.icomp[(key[1], key[2])], length(pv.vec) + 1)
        pv.index[key] = length(pv.vec) + 1
        push!(pv.vec, v)
    else
        if pv.allow_overwrite
            pv.vec[i] = v
        else
            error("Entry $key is already present and overwriting is not allowed")
        end
    end
end

# ====================================================================
struct PVModel{T}
    mname::Symbol
    m::PVSet
    PVModel(mname::Symbol, m::PVSet{T}) where T = new{T}(mname, m)
end

getindex( pv::PVModel, mname::Symbol, cname::Symbol, pname::Symbol) = getindex( pv.m   ,    mname, cname, pname)
getindex( pv::PVModel               , cname::Symbol, pname::Symbol) = getindex( pv.m   , pv.mname, cname, pname)
getindex( pv::PVModel               , cname::Symbol               ) = getindex( pv.m   , pv.mname, cname)
setindex!(pv::PVModel, v            , cname::Symbol, pname::Symbol) = setindex!(pv.m, v, pv.mname, cname, pname)
