using TypedJSON
import TypedJSON: lower, reconstruct

# Can't serialize solver's status, replace it with nothing
TypedJSON.lower(s::FitSummary) =
    TypedJSON.JSONDict(FitSummary(s.start, s.elapsed, s.ndata, s.nfree, s.fitstat, s.status, nothing))

TypedJSON.reconstruct(::Val{Symbol("GModelFit.Domain")}, dict) = Domain(dict[:axis]...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.CartesianDomain")}, dict) = CartesianDomain(dict[:axis]..., roi=dict[:roi])
TypedJSON.reconstruct(::Val{Symbol("GModelFit.Parameter")}, dict) = GModelFit.Parameter(values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.ComponentSnapshot")}, dict) = GModelFit.ComponentSnapshot(values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.ModelSnapshot")}, dict) = GModelFit.ModelSnapshot(values(dict)...)

TypedJSON.reconstruct(::Val{Symbol("GModelFit.FunctDesc")}, dict) =
    GModelFit.FunctDesc((args...) -> nothing,
                        dict[:display], dict[:args], dict[:optargs])

TypedJSON.reconstruct(::Val{Symbol("GModelFit.Solvers.FitSummary")}, dict) = GModelFit.Solvers.FitSummary(values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.Solvers.SolverStatusOK")}, dict) = GModelFit.Solvers.SolverStatusOK()
TypedJSON.reconstruct(::Val{Symbol("GModelFit.Solvers.SolverStatusWarn")}, dict) = GModelFit.Solvers.SolverStatusWarn(values(dict)...)
TypedJSON.reconstruct(::Val{Symbol("GModelFit.Solvers.SolverStatusError")}, dict) = GModelFit.Solvers.SolverStatusError(values(dict)...)

function TypedJSON.reconstruct(::Val{Symbol("GModelFit.Measures")}, dict)
    dom = dict[:domain]
    tmp = dict[:values]
    if isa(dom, CartesianDomain)
        return Measures(dom, reshape(tmp[1], size(dom)), reshape(tmp[2], size(dom)))
    else
        return Measures(dom, tmp[1], tmp[2])
    end
end
