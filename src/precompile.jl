using PrecompileTools
@compile_workload begin
    let
        io = IOBuffer()

        x = [0.1, 1.1, 2.1, 3.1, 4.1]
        domain = Domain(x)
        model = Model(@fd (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))
        data = GModelFit.mock(Measures, model, domain, seed=1)
        bestfit, fsumm = fit(model, data)

        x = 0:0.05:6
        model = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                      :l2  => GModelFit.Gaussian(1, 3, 0.5),
                      :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                      :main => SumReducer(:l1, :l2, :bkg));
        model[:l2, :norm].cast = @fd (m, v) -> v + m[:l1, :norm]
        data = GModelFit.mock(Measures, model, Domain(x), seed=1)
        bestfit, fsumm = fit(model, data, GModelFit.Solvers.cmpfit())
        show(io, bestfit)
        show(io, stat)

        x = 0:0.05:6
        model1 = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                       :l2  => GModelFit.Gaussian(1, 3, 0.5),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg));

        model2 = Model(:l1  => GModelFit.Gaussian(0.8, 2.1, 0.1),
                       :l2  => GModelFit.Gaussian(1.2, 2.5, 0.4),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg));

        ms = ModelSet(:a => model1, :b => model2)
        freeze!(model1, :bkg);
        freeze!(model2, :bkg);
        thaw!(model1, :bkg);
        thaw!(model2, :bkg);

        model2[:bkg, :offset].patch = @fd m -> m[:a, :bkg, :offset]
        model2[:bkg, :slope].patch  = @fd m -> m[:a, :bkg, :slope]
        model1[:l2 , :center].patch = @fd m -> m[:b, :l2 , :center]

        data = GModelFit.mock(Measures, ms, [Domain(x), Domain(x)], seed=1)
        bestfit, fsumm = fit(ms, data)
        show(io, bestfit)
        show(io, stat)
    end
end
