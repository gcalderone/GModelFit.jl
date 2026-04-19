using PrecompileTools
@compile_workload begin
    let
        io = IOBuffer()

        x = [0.1, 1.1, 2.1, 3.1, 4.1]
        domain = Domain(x)
        model = Model(@fd (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))
        data = GModelFit.mock(Measures, model, domain, seed=1)
        fitstat(model, data)
        bestfit, fsumm = fit(model, data)

        x = 0:0.05:6
        model = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                      :l2  => GModelFit.Gaussian(1, 3, 0.5),
                      :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                      :main => SumReducer(:l1, :l2, :bkg));
        model[:l2, :norm].reparam = @fd (m, v) -> v + m[:l1, :norm]
        data = GModelFit.mock(Measures, model, Domain(x), seed=1)
        bestfit, fsumm = fit(model, data, GModelFit.Solvers.cmpfit())
        show(io, bestfit)
        show(io, stat)

        x = 0:0.05:6
        ms = ModelSet()
        ms[:a] = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                       :l2  => GModelFit.Gaussian(1, 3, 0.5),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg));

        ms[:b] = Model(:l1  => GModelFit.Gaussian(0.8, 2.1, 0.1),
                       :l2  => GModelFit.Gaussian(1.2, 2.5, 0.4),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg));
        freeze!(ms[:a], :bkg);
        freeze!(ms[:b], :bkg);
        thaw!(ms[:a], :bkg);
        thaw!(ms[:b], :bkg);

        ms[:b][:bkg, :offset].patch = @fd m -> m[:a, :bkg, :offset]
        ms[:b][:bkg, :slope].patch  = @fd m -> m[:a, :bkg, :slope]
        ms[:a][:l2 , :center].patch = @fd m -> m[:b, :l2 , :center]

        data = GModelFit.mock(Measures, ms, Dict(:a => Domain(x), :b => Domain(x)), seed=1)
        fitstat(ms, data)
        bestfit, fsumm = fit(ms, data)
        show(io, bestfit)
        show(io, stat)
    end
end
