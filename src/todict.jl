const todict_opt = Dict(
    :rebin => 1,
    :addcomps => true,
    :selcomps => Vector{Symbol}()
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
    for (pid, param) in getparams(comp)
        parname = pid.name
        if pid.index >= 1
            parname = Symbol(parname, "[", pid.index, "]")
        end
        out[:params][parname] = todict(param)
    end
    return out
end


function todict(ceval::CompEval; forceadd=false)
    y = ceval.buffer
    i = findall(isfinite.(y))

    out = MDict()
    out[:counter] = ceval.counter
    out[:min] = minimum(y[i])
    out[:max] = maximum(y[i])
    out[:mean] = mean(y[i])
    out[:error] = (length(i) == length(y))
    if todict_opt[:addcomps] || forceadd
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
    out[:x] = rebin_data(todict_opt[:rebin], pred.domain[1])
    out[:components] = MDict()
    out[:compevals]  = MDict()
    for (cname, ceval) in pred.cevals
        out[:components][cname] = todict(ceval.comp)
        out[:components][cname][:fixed] = (ceval.cfixed >= 1)
        out[:compevals][ cname] = todict(ceval, forceadd = (cname in todict_opt[:selcomps]))
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
    p = rebin_data(todict_opt[:rebin], geteval(pred))
    y, err = rebin_data(todict_opt[:rebin], data.val, data.unc)
    out[:meta] = data.meta
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
    for (pname, params) in comp
        if isa(params, AbstractArray)
            for i in 1:length(params)
                out[Symbol(pname, "[", i, "]")] = todict(params[i])
            end
        else
            out[pname] = todict(params)
        end
    end
    return out
end

function todict(res::BestFitResult)
    out = MDict()
    preds = [MDict(:components => MDict()) for id in 1:length(res.preds)]
    for id in 1:length(res.preds)
        for (cname, comp) in res.preds[id]
            preds[id][:components][cname] = todict(comp)
        end
    end
    out[:predictions] = preds
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

function recursive_copy!(from::MDict, to::MDict)
    for (key, value) in from
        if haskey(to, key)
            @assert isa(value,   MDict)
            @assert isa(to[key], MDict)
            recursive_copy!(value, to[key])
        else
            haskey(to, :meta)  ||  (to[:meta] = MDict())
            to[:meta][key] = value
        end
    end
end

function todict(model::Model,
                data::Union{Nothing, Vector{T}}=nothing,
                bestfit::Union{Nothing, BestFitResult}=nothing) where T <: AbstractMeasures
    out = MDict()

    out[:predictions] = Vector{MDict}()
    for id in 1:length(model.preds)
        push!(out[:predictions], todict(model.preds[id]))
        recursive_copy!(model.meta[id], out[:predictions][end])
    end

    if !isnothing(data)
        out[:data] = Vector{MDict}()
        @assert length(model.preds) == length(data)
        for id in 1:length(data)
            push!(out[:data], todict(model.preds[id], data[id]))
        end
    end

    if !isnothing(bestfit)
        out[:bestfit] = todict(bestfit)
    end

    return out
end
