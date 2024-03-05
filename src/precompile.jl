#=

Note: I can't compile ModelSnapshot since it involves a show()
call, which in turn involve PrettyTables, which is not currently
possible to precompile (an error is raised).

=#

using PrecompileTools
@compile_workload begin
    let
        x = [0.1, 1.1, 2.1, 3.1, 4.1]
        domain = Domain(x)
        model = Model(@λ (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))
        data = GModelFit.mock(Measures, model, domain, seed=1)
        status = fit(model, data)

        
        x = 0:0.05:6
        model = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                      :l2  => GModelFit.Gaussian(1, 3, 0.5),
                      :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                      :main => SumReducer(:l1, :l2, :bkg));
        model[:l2].norm.patch = :l1
        model[:l2].norm.patch = @λ (m, v) -> v + m[:l1].norm
        data = GModelFit.mock(Measures, model, Domain(x), seed=1)
        status = fit(model, data, minimizer=cmpfit())


        x = 0:0.05:6
        model1 = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                       :l2  => GModelFit.Gaussian(1, 3, 0.5),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg));

        model2 = Model(:l1  => GModelFit.Gaussian(0.8, 2.1, 0.1),
                       :l2  => GModelFit.Gaussian(1.2, 2.5, 0.4),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg));

        model = [model1, model2]
        freeze!(model[1], :bkg);
        freeze!(model[2], :bkg);
        thaw!(model[1], :bkg);
        thaw!(model[2], :bkg);

        model[2][:bkg].offset.mpatch = @λ m -> m[1][:bkg].offset
        model[2][:bkg].slope.mpatch  = @λ m -> m[1][:bkg].slope
        model[1][:l2].center.mpatch = @λ m -> m[2][:l2].center

        data = GModelFit.mock(Measures, model, [Domain(x), Domain(x)], seed=1)
        status = fit(model, data)
    end
end
