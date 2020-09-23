function todict(param::Parameter)
    out = MDict()
    out[:__meta] = param.meta
    out[:__fixed] = param.fixed
    out[:__value] = param.val
    out[:__low] = param.low
    out[:__high] = param.high
    out[:__error] = !(param.low <= param.val <= param.high)
    return out
end


function todict(comp::AbstractComponent)
    out = MDict()

    ctype = split(string(typeof(comp)), ".")
    (ctype[1] == "GFit")  &&   (ctype = ctype[2:end])
    ctype = join(ctype, ".")
    out[:__type] = ctype

    out[:__params] = MDict()
    for (pname, param) in getparams(comp)
        parname = pname[1]
        if pname[2] >= 1
            parname = Symbol(parname, "[:__ * string(pname[2]) * ]")
        end
        out[:__params][parname] = todict(param)
    end
    return out
end


function todict(ceval::CompEval)
    y = ceval.eval
    i = findall(isfinite.(y))

    out = MDict()
    out[:__counter] = ceval.counter
    out[:__min] = minimum(y[i])
    out[:__max] = maximum(y[i])
    out[:__mean] = mean(y[i])
    out[:__error] = (length(i) == length(y))
    out[:__y] = y
    return out
end


function todict(reval::ReducerEval)
    y = reval.eval
    i = findall(isfinite.(y))

    out = MDict()
    out[:__counter] = reval.counter
    out[:__min] = minimum(y[i])
    out[:__max] = maximum(y[i])
    out[:__mean] = mean(y[i])
    out[:__error] = (length(i) == length(y))
    out[:__y] = y
    return out
end


function todict(pred::Prediction)
    out = MDict()
    out[:__meta] = pred.meta
    out[:__x] = domain(pred)
    out[:__components] = MDict()
    for (cname, ceval) in pred.cevals
        out[:__components][cname] = todict(ceval)
    end
    out[:__reducers] = MDict()
    for (rname, reval) in pred.revals
        out[:__reducers][rname] = todict(reval)
    end
    out[:__main_reducer] = pred.rsel
    return out
end


function todict(pred::Prediction, data::Measures_1D)
    out = MDict()
    out[:__meta] = data.meta
    out[:__y] = data.val
    out[:__err] = data.unc
    out[:__residuals] = (data.val .- geteval(pred)) ./ data.unc
    return out
end


function todict(param::BestFitPar)
    out = MDict()
    out[:__val] = param.val
    out[:__unc] = param.unc
    out[:__fixed] = param.fixed
    out[:__patched] = param.patched
    return out
end

function todict(comp::BestFitComp)
    out = MDict()
    params = getfield(comp, :params)
    for (pname, param) in params
        out[pname] = todict(param)
    end
    return out
end

function todict(res::BestFitResult)
    out = MDict()
    out[:__components] = MDict()
    for (cname, comp) in res.comps
        out[:__components][cname] = todict(comp)
    end
    out[:__ndata] = res.ndata
    out[:__dof] = res.dof
    out[:__cost] = res.cost
    out[:__status] = res.status
    out[:__log10testprob] = res.log10testprob
    out[:__elapsed] = res.elapsed
    return out
end


todict(model::Model, data::T) where T <: AbstractMeasures = todict(model, [data])
todict(model::Model, data::T, bestfit::BestFitResult) where T <: AbstractMeasures = todict(model, [data], bestfit)

function todict(model::Model,
    data::Union{Nothing, Vector{T}}=nothing,
    bestfit::Union{Nothing, BestFitResult}=nothing) where T <: AbstractMeasures

    out = MDict()

    out[:__components] = MDict()
    for (cname, comp) in model.comps
        out[:__components][cname] = todict(comp)
        out[:__components][cname][:__fixed] = model.cfixed[cname]
        out[:__components][cname][:__meta] = metadict(model, cname)
    end

    out[:__predictions] = Vector{MDict}()
    for pred in model.preds
        push!(out[:__predictions], todict(pred))
        for (name, ceval) in pred.cevals
            out[:__predictions][end][:__components][name][:__meta] = metadict(model, name)
        end
        for (name, ceval) in pred.revals
            out[:__predictions][end][:__reducers][name][:__meta] = metadict(model, name)
        end
    end

    if !isnothing(data)
        out[:__data] = Vector{MDict}()
        @assert length(model.preds) == length(data)
        for i in 1:length(data)
            push!(out[:__data], todict(model.preds[i], data[i]))
        end
    end

    if !isnothing(bestfit)
        out[:__bestfit] = todict(bestfit)
    end

    return out
end
