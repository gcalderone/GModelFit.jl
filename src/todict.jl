const todict_opt = Dict(
    :rebin => 1,
    :addcomps => true
)

rebin_data(rebin, v) = rebin_data(rebin, v, ones(eltype(v), length(v)))[1]
function rebin_data(rebin::Int, v, e)
    @assert length(v) == length(e)
    @assert 1 <= rebin <= length(v)
    (rebin == 1)  &&  (return (v, e))
    nin = length(v)
    nout = div(nin, rebin)
    val = zeros(eltype(v), nout)
    unc = zeros(eltype(e), nout)
    w = 1 ./ (e.^2)
    for i in 1:nout-1
        r = (i-1)*rebin+1:i*rebin
        val[i] = sum(w[r] .* v[r]) ./ sum(w[r])
        unc[i] = sqrt(1 / sum(w[r]))
    end
    r = (nout-1)*rebin+1:nin
    val[nout] = sum(w[r] .* v[r]) ./ sum(w[r])
    unc[nout] = sqrt(1 / sum(w[r]))
    return (val, unc)
end


function todict(param::Parameter)
    out = MDict()
    out[:meta] = param.meta
    out[:fixed] = param.fixed
    out[:value] = param.val
    out[:low] = param.low
    out[:high] = param.high
    out[:error] = !(param.low <= param.val <= param.high)
    return out
end


function todict(comp::AbstractComponent)
    out = MDict()

    ctype = split(string(typeof(comp)), ".")
    (ctype[1] == "GFit")  &&   (ctype = ctype[2:end])
    ctype = join(ctype, ".")
    out[:type] = ctype

    out[:params] = MDict()
    for (pname, param) in getparams(comp)
        parname = pname[1]
        if pname[2] >= 1
            parname = Symbol(parname, "[", pname[2], "]")
        end
        out[:params][parname] = todict(param)
    end
    return out
end


function todict(ceval::CompEval)
    y = ceval.buffer
    i = findall(isfinite.(y))

    out = MDict()
    out[:counter] = ceval.counter
    out[:min] = minimum(y[i])
    out[:max] = maximum(y[i])
    out[:mean] = mean(y[i])
    out[:error] = (length(i) == length(y))
    if todict_opt[:addcomps]
        out[:y] = rebin_data(todict_opt[:rebin], y)
    else
        out[:y] = Vector{Float64}()
    end
    return out
end


function todict(reval::ReducerEval)
    y = reval.buffer
    i = findall(isfinite.(y))

    out = MDict()
    out[:counter] = reval.counter
    out[:min] = minimum(y[i])
    out[:max] = maximum(y[i])
    out[:mean] = mean(y[i])
    out[:error] = (length(i) == length(y))
    out[:y] = rebin_data(todict_opt[:rebin], y)
    return out
end


function todict(pred::Prediction)
    out = MDict()
    out[:meta] = pred.meta
    out[:x] = rebin_data(todict_opt[:rebin], domain(pred))
    out[:components] = MDict()
    for (cname, ceval) in pred.cevals
        out[:components][cname] = todict(ceval)
    end
    out[:reducers] = MDict()
    for (rname, reval) in pred.revals
        out[:reducers][rname] = todict(reval)
    end
    out[:main_reducer] = pred.rsel
    return out
end


function todict(pred::Prediction, data::Measures_1D)
    out = MDict()
    out[:meta] = data.meta
    p = rebin_data(todict_opt[:rebin], geteval(pred))
    y, err = rebin_data(todict_opt[:rebin], data.val, data.unc)
    out[:y] = y
    out[:err] = err
    out[:residuals] = (y .- p) ./ err
    return out
end


function todict(param::BestFitPar)
    out = MDict()
    out[:val] = param.val
    out[:unc] = param.unc
    out[:fixed] = param.fixed
    out[:patched] = param.patched
    return out
end

function todict(comp::BestFitComp)
    out = MDict()
    params = getfield(comp, :params)
    for (pname, param) in params
        if isa(param, Vector)
            for i in 1:length(param)
                out[Symbol(pname, "[", i, "]")] = todict(param[i])
            end
        else
            out[pname] = todict(param)
        end
    end
    return out
end

function todict(res::BestFitResult)
    out = MDict()
    out[:components] = MDict()
    for (cname, comp) in res.comps
        out[:components][cname] = todict(comp)
    end
    out[:ndata] = res.ndata
    out[:dof] = res.dof
    out[:cost] = res.cost
    out[:status] = res.status
    out[:log10testprob] = res.log10testprob
    out[:elapsed] = res.elapsed
    return out
end


todict(model::Model, data::T) where T <: AbstractMeasures = todict(model, [data])
todict(model::Model, data::T, bestfit::BestFitResult) where T <: AbstractMeasures = todict(model, [data], bestfit)

function todict(model::Model,
                data::Union{Nothing, Vector{T}}=nothing,
                bestfit::Union{Nothing, BestFitResult}=nothing) where T <: AbstractMeasures

    out = MDict()

    out[:components] = MDict()
    for (cname, comp) in model.comps
        out[:components][cname] = todict(comp)
        out[:components][cname][:fixed] = model.cfixed[cname]
        out[:components][cname][:meta] = metadict(model, cname)
    end

    out[:predictions] = Vector{MDict}()
    for pred in model.preds
        push!(out[:predictions], todict(pred))
        for (name, ceval) in pred.cevals
            out[:predictions][end][:components][name][:meta] = metadict(model, name)
        end
        for (name, ceval) in pred.revals
            out[:predictions][end][:reducers][name][:meta] = metadict(model, name)
        end
    end

    if !isnothing(data)
        out[:data] = Vector{MDict}()
        @assert length(model.preds) == length(data)
        for i in 1:length(data)
            push!(out[:data], todict(model.preds[i], data[i]))
        end
    end

    if !isnothing(bestfit)
        out[:bestfit] = todict(bestfit)
    end

    return out
end
