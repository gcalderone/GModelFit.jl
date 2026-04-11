using TypedJSON
import TypedJSON: lower, reconstruct

TypedJSON.lower(::T) where {T <: Union{Model, ModelEval, Likelihood}} = error("Can't serialize object of type $(T), try taking a snapshot.")


# Can't serialize response object, use RespSnapshot instead
TypedJSON.lower(s::GaussianData) = TypedJSON.JSONDict(GaussianData(RespSnapshot(s.resp), [getfield(s, f) for f in fieldnames(typeof(s))[2:end]]...))

# Can't serialize solver's status, replace it with nothing
TypedJSON.lower(s::CompSnapshot)    = TypedJSON.JSONDict(CompSnapshot(   nothing, (), [getfield(s, f) for f in fieldnames(typeof(s))[3:end]]...))
TypedJSON.lower(s::GroupSnapshot)   = TypedJSON.JSONDict(GroupSnapshot(  nothing, (), [getfield(s, f) for f in fieldnames(typeof(s))[3:end]]...))
TypedJSON.lower(s::GroupLHSnapshot) = TypedJSON.JSONDict(GroupLHSnapshot(nothing, (), [getfield(s, f) for f in fieldnames(typeof(s))[3:end]]...))

function TypedJSON.lower(s::ModelSnapshot)
    out = TypedJSON.JSONDict(s)
    delete!(out.dict, :ifree)  # not needed in serialized JSON
    return out
end

function TypedJSON.lower(s::LikelihoodSnapshot)
    out = TypedJSON.JSONDict(s)
    delete!(out.dict, :ifree)  # not needed in serialized JSON
    return out
end

TypedJSON.reconstruct(::Val{Symbol("GModelFit.Domain")}                   , dict) = Domain(values(dict[:coords])...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.CartesianDomain")}          , dict) = CartesianDomain(                    values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.RespSnapshot")}             , dict) = RespSnapshot(                       values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.GaussianData")}             , dict) = GaussianData(                       values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.ParameterEval")}            , dict) = GModelFit.ParameterEval(            values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.CompSnapshot")}             , dict) = GModelFit.CompSnapshot(             values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.GroupSnapshot")}            , dict) = GModelFit.GroupSnapshot(            values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.GroupLHSnapshot")}          , dict) = GModelFit.GroupLHSnapshot(          values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.ModelSnapshot")}            , dict) = GModelFit.ModelSnapshot(            values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.LikelihoodSnapshot")}       , dict) = GModelFit.LikelihoodSnapshot(       values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.Solvers.SolverStatusOK")}   , dict) = GModelFit.Solvers.SolverStatusOK()
TypedJSON.reconstruct(::Val{Symbol("GModelFit.Solvers.SolverStatusWarn")} , dict) = GModelFit.Solvers.SolverStatusWarn( values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.Solvers.SolverStatusError")}, dict) = GModelFit.Solvers.SolverStatusError(values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.FunctDesc")}                , dict) = GModelFit.FunctDesc((args...) -> nothing, collect(values(dict))[2:end]...)
