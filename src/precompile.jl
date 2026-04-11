using PrecompileTools
@compile_workload begin
    io = devnull

    x = [0.1, 1.1, 2.1, 3.1, 4.1]
    dom = Domain(x)
    model = Model(@fd (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))
    data = GModelFit.mock(GaussianData, dom, model, seed=1)
    bestfit, fsumm = fit(data, model)
    gofstat(bestfit)
    TypedJSON.roundtrip(snapshot(bestfit))
    TypedJSON.roundtrip(snapshot(ModelEval(dom, model)))

    x = 0:0.05:6
    dom = Domain(x)
    model = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                  :l2  => GModelFit.Gaussian(1, 3, 0.5),
                  :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                  :main => SumReducer(:l1, :l2, :bkg))
    model[:l2, :norm].reparam = @fd (m, v) -> v + m[:l1, :norm]
    data = GModelFit.mock(GaussianData, dom, model, seed=1)
    bestfit, fsumm = fit(data, model; solver=GModelFit.Solvers.cmpfit())
    gofstat(bestfit)
    TypedJSON.roundtrip(snapshot(bestfit))
    TypedJSON.roundtrip(snapshot(ModelEval(dom, model)))
    export_params(bestfit)
    export_folded_eval(bestfit)
    export_unfolded_eval(bestfit)
    show(io, dom)
    show(io, bestfit)
    show(io, fsumm)

    x = 0:0.05:6
    model = Model{Tuple{Symbol}}()
    model[:a, :l1]   = GModelFit.Gaussian(1, 2, 0.2)
    model[:a, :l2]   = GModelFit.Gaussian(1, 3, 0.5)
    model[:a, :bkg]  = GModelFit.OffsetSlope(0.5, 1, 0.1)
    model[:a, :main] = SumReducer(:l1, :l2, :bkg)

    model[:b, :l1]   = GModelFit.Gaussian(0.8, 2.1, 0.1)
    model[:b, :l2]   = GModelFit.Gaussian(1.2, 2.5, 0.4)
    model[:b, :bkg]  = GModelFit.OffsetSlope(0.5, 1, 0.1)
    model[:b, :main] = SumReducer(:l1, :l2, :bkg)
    freeze!(model[:a, :bkg])
    freeze!(model[:b, :bkg])
    thaw!(  model[:a, :bkg])
    thaw!(  model[:b, :bkg])

    model[:b, :bkg, :offset].patch = @fd m -> m[:a, :bkg, :offset]
    model[:b, :bkg, :slope].patch  = @fd m -> m[:a, :bkg, :slope]
    model[:a, :l2 , :center].patch = @fd m -> m[:b, :l2 , :center]

    dom = Dict((:a,) => Domain(x), (:b,) => Domain(x))
    data = GModelFit.mock(GaussianData, dom, model, seed=1)
    bestfit, fsumm = fit(data, model)
    gofstat(bestfit)
    TypedJSON.roundtrip(snapshot(bestfit))
    TypedJSON.roundtrip(snapshot(ModelEval(dom, model)))
    export_params(bestfit[(:a,)])
    export_folded_eval(bestfit[(:a,)])
    export_unfolded_eval(bestfit[(:a,)])
    show(io, bestfit)
    show(io, fsumm)
end
