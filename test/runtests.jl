using Random, Test, TypedJSON, DataStructures, Gnuplot
using GModelFit
import GModelFit: PVModel, PVSet

@testset "GModelFit.jl Comprehensive Test Suite" begin

    # ====================================================================
    @testset "1. Data Structures & PVSet" begin
        m = PVSet{Int64}() # Default allow_overwrite is false
        a = PVModel(:M1, m)
        b = PVModel(:M2, m)

        m[:M1, :ca, :pa1] = 1
        m[:M1, :ca, :pa2] = 2
        m[:M2, :cb, :pb1] = 3
        
        # Test overwrite protection
        @test_throws Exception (m[:M1, :ca, :pa1] = 1)

        @test a[:ca, :pa1] == 1
        @test a[:ca, :pa2] == 2
        @test b[:cb, :pb1] == 3
        @test a[:M1, :ca, :pa1] == 1
        @test a[:M1, :ca, :pa2] == 2
        @test a[:M2, :cb, :pb1] == 3
        @test b[:M1, :ca, :pa1] == 1
        @test b[:M1, :ca, :pa2] == 2
        @test b[:M2, :cb, :pb1] == 3

        @test all(a[:ca] .== [1, 2])
        @test all(b[:cb] .== [3])
        @test length(b[:nonexisting]) == 0
        @test all(m[:M1, :ca] .== [1, 2])
        @test all(m[:M2, :cb] .== [3])
        @test length(m[:M1, :nonexisting]) == 0
        
        # Explicit allow_overwrite test
        m_allow = PVSet{Int64}(allow_overwrite=true)
        m_allow[:M1, :ca, :pa1] = 1
        m_allow[:M1, :ca, :pa1] = 2
        @test m_allow[:M1, :ca, :pa1] == 2
    end

    # ====================================================================
    @testset "2. Core Logic & Components" begin
        @testset "Patch and Reparam Mutually Exclusive" begin
            comp = GModelFit.Polynomial(1.0, 2.0)
            comp.params[:p0].patch = :p1
            comp.params[:p0].reparam = @fd (p, v) -> v * 2
            m = Model(comp)
            @test_throws AssertionError GModelFit.getparams(m)
        end

        @testset "Parameter Bounds Enforcement" begin
            comp = GModelFit.Polynomial(1.0)
            comp.params[:p0].low = 5.0 # Set lower bound higher than current val
            m = Model(comp)
            @test_throws AssertionError GModelFit.getparams(m)
            
            comp.params[:p0].low = -Inf
            comp.params[:p0].high = 0.0 # Set upper bound lower than current val
            @test_throws AssertionError GModelFit.getparams(m)
            
            # Physical bounds on standard components
            g_comp = GModelFit.Gaussian(1.0, 0.0, 1.0)
            g_comp[:norm].val = -1.0 # Force negative norm
            @test_throws AssertionError GModelFit.getparams(Model(g_comp))
        end

        @testset "Circular Dependency Detection" begin
            c1 = GModelFit.FComp((x, deps...) -> x, [:comp2])
            c2 = GModelFit.FComp((x, deps...) -> x, [:comp1])
            m = Model(:comp1 => c1, :comp2 => c2)
            @test_throws AssertionError GModelFit.deptree(m)
        end
        
        @testset "Component Specifics" begin
            # 3D z0 Parameter Assignment mapping
            offset, x0, y0, z0, slopeX, slopeY, slopeZ = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0
            comp = GModelFit.OffsetSlope(offset, x0, y0, z0, slopeX, slopeY, slopeZ)
            @test comp[:z0].val == 4.0
            
            # FComp Splatting Rejection
            @test_throws AssertionError GModelFit.FComp(@fd (x, args...) -> x)
        end
    end

    # ====================================================================
    @testset "3. Mocking, IO, & Display" begin
        @testset "Mock Dataset Function Purity" begin
            m = Model(GModelFit.Gaussian(1.0, 0.0, 1.0))
            dom = Domain(collect(1.0:10.0))
            baseline_eval = m(dom)
            eval_copy = copy(baseline_eval)
            mock_data = GModelFit.mock(Measures, m, dom)
            # Ensure internal cache wasn't modified in-place
            @test m(dom) == eval_copy 
        end

        @testset "Base.show Handles NaN-Only Measures Gracefully" begin
            dom = Domain(10)
            meas = Measures(dom, fill(NaN, 10), fill(NaN, 10))
            io = IOBuffer()
            @test (show(io, meas); true)
        end

        @testset "Gnuplot Recipe Accepts Serialized Snapshots" begin
            comp = GModelFit.Polynomial(1.0, 2.0)
            m = Model(:poly => comp)
            dom = Domain(collect(1.0:10.0))
            m(dom) # Prime evaluation
            mseval = GModelFit.ModelSetEval{Float64}(ModelSet(:_ => m), [dom])
            GModelFit.update_eval!(mseval)
            snapshot = GModelFit.ModelSnapshot(mseval.vec[1], GModelFit.getparams(m))
            recipe_out = Gnuplot.recipe(snapshot)
            @test recipe_out isa Vector{Gnuplot.AbstractGPSpec}
        end
    end

    # ====================================================================
    @testset "4. 1D Fitting & Multi-component Models" begin
        @testset "Simple 1D Polynomial" begin
            x = 0:0.1:5
            model = Model(:main => @fd (x, a2=1, a1=1, a0=5) -> @. (a2 * x^2  +  a1 * x  + a0))
            data = GModelFit.mock(Measures, model, Domain(x), seed=1)
            bestfit, fsumm = fit(model, data)
            
            @test isapprox(bestfit[:main, :a2].val, 1, atol=3 * bestfit[:main, :a2].unc)
            @test isapprox(bestfit[:main, :a1].val, 1, atol=3 * bestfit[:main, :a1].unc)
            @test isapprox(bestfit[:main, :a0].val, 5, atol=3 * bestfit[:main, :a0].unc)
            @test fsumm.ndata == 51
            @test fsumm.nfree == 3
            @test isapprox(fsumm.fitstat, 1.7, atol=0.2)
        end

        @testset "Complex 1D Function" begin
            f = @fd (x, p1=1, p2=1.e-3, p3=1e-6, p4=4, p5=5) -> @. (p1 + p2 * x + p3 * x^2 + p4 * sin(p5 * x))  * cos(x)
            x = 1.:50:10000
            model = Model(f)
            data = GModelFit.mock(Measures, model, Domain(x), seed=1)
            bestfit, fsumm = fit(model, data)
            
            @test abs(bestfit[:main, :p1].val - 1)    / bestfit[:main, :p1].unc < 3
            @test abs(bestfit[:main, :p2].val - 1e-3) / bestfit[:main, :p2].unc < 3
            @test abs(bestfit[:main, :p3].val - 1e-6) / bestfit[:main, :p3].unc < 3
            @test abs(bestfit[:main, :p4].val - 4   ) / bestfit[:main, :p4].unc < 3
            @test abs(bestfit[:main, :p5].val - 5   ) / bestfit[:main, :p5].unc < 3
            @test fsumm.ndata == 200
            @test fsumm.nfree == 5
            @test isapprox(fsumm.fitstat, 1, atol=0.2)
        end

        @testset "Component Combination & Dependency Injection" begin
            f1 = @fd (x, p1=1, p2=1e-3, p3=1e-6) -> @.  p1  +  p2 * x  +  p3 * x^2
            f2 = @fd (x, p4=4, p5=5) -> @. p4 * sin(p5 * x)
            f3 = @fd (x) -> cos.(x)

            x = 1.:50:10000
            
            # Setup 1 (All at once)
            model_a = Model(:f1 => f1, :f2 => f2, :f3 => f3,
                            :main => @fd (x, f1, f2, f3) -> (f1 .+ f2) .* f3)
            data_a = GModelFit.mock(Measures, model_a, Domain(x), seed=1)
            bestfit_a, fsumm_a = fit(model_a, data_a)
            
            # Setup 2 (Incremental construction)
            model_b = Model(:f1 => f1)
            model_b[:f2] = f2
            model_b[:f3] = f3
            model_b[:main] = @fd (x, f1, f2, f3) -> (f1 .+ f2) .* f3
            data_b = GModelFit.mock(Measures, model_b, Domain(x), seed=1)
            bestfit_b, fsumm_b = fit(model_b, data_b)
            
            for (bestfit, fsumm) in [(bestfit_a, fsumm_a), (bestfit_b, fsumm_b)]
                @test abs(bestfit[:f1, :p1].val - 1)    / bestfit[:f1, :p1].unc < 3
                @test abs(bestfit[:f1, :p2].val - 1e-3) / bestfit[:f1, :p2].unc < 3
                @test abs(bestfit[:f1, :p3].val - 1e-6) / bestfit[:f1, :p3].unc < 3
                @test abs(bestfit[:f2, :p4].val - 4   ) / bestfit[:f2, :p4].unc < 3
                @test abs(bestfit[:f2, :p5].val - 5   ) / bestfit[:f2, :p5].unc < 3
                @test fsumm.ndata == 200
                @test fsumm.nfree == 5
                @test isapprox(fsumm.fitstat, 1, atol=0.2)
            end
        end
    end

    # ====================================================================
    @testset "5. 2D & 3D Fitting (Cartesian Domains)" begin
        @testset "1D Distance Baseline" begin
            dom = Domain(-5:0.1:5)
            model = Model(:main => @fd (x, cx=0.) -> @. x-cx)
            data = GModelFit.mock(Measures, model, dom)
            bestfit, fsumm = fit(model, data)
            @test fsumm.fitstat < 1.5
        end

        @testset "2D Domain Evaluation" begin
            dom_base = CartesianDomain(-5:0.1:5, -5:0.1:5)
            for dom in [dom_base, Domain(dom_base)]
                local model = Model(:main => @fd (x, y, cx=0., cy=0.) -> @. sqrt((x-cx)^2 + (y-cy)^2))
                local data = GModelFit.mock(Measures, model, dom)
                local bestfit, fsumm = fit(model, data)
                @test fsumm.fitstat < 1.5

                local model_bkg = Model(:bkg => GModelFit.OffsetSlope(0, 0, 0, 1, 1),
                                        :main => @fd (x, y, bkg, cx=0., cy=0.) -> @. bkg + sqrt((x-cx)^2 + (y-cy)^2))
                local data_bkg = GModelFit.mock(Measures, model_bkg, dom)
                local bestfit_bkg, fsumm_bkg = fit(model_bkg, data_bkg)
                @test fsumm_bkg.fitstat < 1.5
            end
        end

        @testset "3D Domain Evaluation" begin
            dom_base = CartesianDomain(-5:0.1:5, -5:0.1:5, -5:0.1:5)
            for dom in [dom_base, Domain(dom_base)]
                local model = Model(:main => @fd (x, y, z, cx=0., cy=0., cz=0.) -> @. sqrt((x-cx)^2 + (y-cy)^2 + (z-cz)^2))
                local data = GModelFit.mock(Measures, model, dom)
                local bestfit, fsumm = fit(model, data)
                @test fsumm.fitstat < 1.5

                local model_bkg = Model(:bkg => GModelFit.OffsetSlope(0, 0, 0, 0, 1, 1, 1),
                                        :main => @fd (x, y, z, bkg, cx=0., cy=0., cz=0.) -> @. bkg + sqrt((x-cx)^2 + (y-cy)^2 + (z-cz)^2))
                local data_bkg = GModelFit.mock(Measures, model_bkg, dom)
                local bestfit_bkg, fsumm_bkg = fit(model_bkg, data_bkg)
                @test fsumm_bkg.fitstat < 1.5
            end
        end
    end

    # ====================================================================
    @testset "6. Parameter Patches & Reparams" begin
        x = 0:0.05:6
        model = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                      :l2  => GModelFit.Gaussian(1, 3, 0.5),
                      :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                      :main => SumReducer(:l1, :l2, :bkg))
        
        @testset "Baseline Multi-Component" begin
            data = GModelFit.mock(Measures, model, Domain(x), seed=1)
            bestfit, fsumm = fit(model, data)
            @test abs(bestfit[:l1, :norm].val    - 1)   / bestfit[:l1, :norm].unc < 3
            @test abs(bestfit[:l1, :center].val  - 2)   / bestfit[:l1, :center].unc < 3
            @test abs(bestfit[:l1, :sigma].val   - 0.2) / bestfit[:l1, :sigma].unc < 3
            @test abs(bestfit[:l2, :norm].val    - 1)   / bestfit[:l2, :norm].unc < 3
            @test abs(bestfit[:l2, :center].val  - 3)   / bestfit[:l2, :center].unc < 3
            @test abs(bestfit[:l2, :sigma].val   - 0.5) / bestfit[:l2, :sigma].unc < 3
            @test abs(bestfit[:bkg, :offset].val - 0.5) / bestfit[:bkg, :offset].unc < 3
            @test abs(bestfit[:bkg, :slope].val  - 0.1) / bestfit[:bkg, :slope].unc < 3
            @test isnan(bestfit[:bkg, :x0].unc)
            @test fsumm.ndata == 121
            @test fsumm.nfree == 8
            @test isapprox(fsumm.fitstat, 1.3, atol=0.2)
        end

        @testset "Parameter Patching" begin
            model[:l2, :norm].patch = :l1 # Tie l2 norm to l1
            data = GModelFit.mock(Measures, model, Domain(x), seed=1)
            bestfit, fsumm = fit(model, data)
            @test abs(bestfit[:l1, :norm].val    - 1)   / bestfit[:l1, :norm].unc < 3
            @test abs(bestfit[:l1, :center].val  - 2)   / bestfit[:l1, :center].unc < 3
            @test abs(bestfit[:l1, :sigma].val   - 0.2) / bestfit[:l1, :sigma].unc < 3
            @test isnan(bestfit[:l2, :norm].unc) # Patched parameter should have NaN uncertainty
            @test abs(bestfit[:l2, :center].val  - 3)   / bestfit[:l2, :center].unc < 3
            @test abs(bestfit[:l2, :sigma].val   - 0.5) / bestfit[:l2, :sigma].unc < 3
            @test abs(bestfit[:bkg, :offset].val - 0.5) / bestfit[:bkg, :offset].unc < 3
            @test abs(bestfit[:bkg, :slope].val  - 0.1) / bestfit[:bkg, :slope].unc < 3
            @test isnan(bestfit[:bkg, :x0].unc)
            @test fsumm.ndata == 121
            @test fsumm.nfree == 7
            @test isapprox(fsumm.fitstat, 1.3, atol=0.2)
        end

        @testset "Parameter Reparametrizing (Lambda Functions)" begin
            model[:l2, :norm].patch = nothing
            model[:l2, :norm].reparam = @fd (m, v) -> v + m[:l1, :norm]
            data = GModelFit.mock(Measures, model, Domain(x), seed=1)
            bestfit, fsumm = fit(model, data)
            @test abs(bestfit[:l1, :norm].val    - 1)   / bestfit[:l1, :norm].unc < 3
            @test abs(bestfit[:l1, :center].val  - 2)   / bestfit[:l1, :center].unc < 3
            @test abs(bestfit[:l1, :sigma].val   - 0.2) / bestfit[:l1, :sigma].unc < 3
            @test abs(bestfit[:l2, :norm].actual - 2)   < 0.2
            @test abs(bestfit[:l2, :center].val  - 3)   / bestfit[:l2, :center].unc < 3
            @test abs(bestfit[:l2, :sigma].val   - 0.5) / bestfit[:l2, :sigma].unc < 3
            @test abs(bestfit[:bkg, :offset].val - 0.5) / bestfit[:bkg, :offset].unc < 3
            @test abs(bestfit[:bkg, :slope].val  - 0.1) / bestfit[:bkg, :slope].unc < 3
            @test isnan(bestfit[:bkg, :x0].unc)
            @test fsumm.ndata == 121
            @test fsumm.nfree == 8
            @test isapprox(fsumm.fitstat, 1.3, atol=0.2)
        end
    end

    # ====================================================================
    @testset "7. Multi-Model Sets (ModelSet)" begin
        x = 0:0.05:6
        model1 = Model(:l1  => GModelFit.Gaussian(1, 2, 0.2),
                       :l2  => GModelFit.Gaussian(1, 3, 0.5),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg))

        model2 = Model(:l1  => GModelFit.Gaussian(0.8, 2.1, 0.1),
                       :l2  => GModelFit.Gaussian(1.2, 2.5, 0.4),
                       :bkg => GModelFit.OffsetSlope(0.5, 1, 0.1),
                       :main => SumReducer(:l1, :l2, :bkg))

        ms = ModelSet(:a => model1, :b => model2)
        
        @testset "Freezing Components Across Models" begin
            freeze!(model1, :bkg)
            freeze!(model2, :bkg)
            data = GModelFit.mock(Measures, ms, [Domain(x), Domain(x)], seed=1)
            bestfit, fsumm = fit(ms, data, GModelFit.Solvers.cmpfit())

            @test abs(bestfit[:a, :l1, :norm].val    - 1)   / bestfit[:a, :l1, :norm].unc < 3
            @test abs(bestfit[:a, :l1, :center].val  - 2)   / bestfit[:a, :l1, :center].unc < 3
            @test abs(bestfit[:a, :l1, :sigma].val   - 0.2) / bestfit[:a, :l1, :sigma].unc < 3
            @test abs(bestfit[:a, :l2, :norm].val    - 1)   / bestfit[:a, :l2, :norm].unc < 3
            @test abs(bestfit[:a, :l2, :center].val  - 3)   / bestfit[:a, :l2, :center].unc < 3
            @test abs(bestfit[:a, :l2, :sigma].val   - 0.5) / bestfit[:a, :l2, :sigma].unc < 3
            @test isnan(bestfit[:a, :bkg, :offset].unc)
            @test isnan(bestfit[:a, :bkg, :x0].unc)
            @test isnan(bestfit[:a, :bkg, :slope].unc)
            
            @test abs(bestfit[:b, :l1, :norm].val    - 0.8) / bestfit[:b, :l1, :norm].unc < 3
            @test abs(bestfit[:b, :l1, :center].val  - 2.1) / bestfit[:b, :l1, :center].unc < 3
            @test abs(bestfit[:b, :l1, :sigma].val   - 0.1) / bestfit[:b, :l1, :sigma].unc < 3
            @test abs(bestfit[:b, :l2, :norm].val    - 1.2) / bestfit[:b, :l2, :norm].unc < 3
            @test abs(bestfit[:b, :l2, :center].val  - 2.5) / bestfit[:b, :l2, :center].unc < 3
            @test abs(bestfit[:b, :l2, :sigma].val   - 0.4) / bestfit[:b, :l2, :sigma].unc < 3
            @test isnan(bestfit[:b, :bkg, :offset].unc)
            @test isnan(bestfit[:b, :bkg, :x0].unc)
            @test isnan(bestfit[:b, :bkg, :slope].unc)
            
            @test fsumm.ndata == 242
            @test fsumm.nfree == 12
            @test isapprox(fsumm.fitstat, 1., atol=0.2)
        end

        @testset "Patching Across Models" begin
            thaw!(model1, :bkg)
            thaw!(model2, :bkg)
            ms[:b, :bkg, :offset].patch = @fd m -> m[:a, :bkg, :offset]
            ms[:b, :bkg, :slope].patch  = @fd m -> m[:a, :bkg, :slope]
            ms[:a, :l2 , :center].patch = @fd m -> m[:b, :l2, :center]
            
            data = GModelFit.mock(Measures, ms, [Domain(x), Domain(x)], seed=1)
            bestfit, fsumm = fit(ms, data)

            @test abs(bestfit[:a, :l1, :norm].val    - 1)   / bestfit[:a, :l1, :norm].unc < 3
            @test abs(bestfit[:a, :l1, :center].val  - 2)   / bestfit[:a, :l1, :center].unc < 3
            @test abs(bestfit[:a, :l1, :sigma].val   - 0.2) / bestfit[:a, :l1, :sigma].unc < 3
            @test abs(bestfit[:a, :l2, :norm].val    - 1)   / bestfit[:a, :l2, :norm].unc < 3
            @test isnan(bestfit[:a, :l2, :center].unc)
            @test abs(bestfit[:a, :l2, :sigma].val   - 0.5) / bestfit[:a, :l2, :sigma].unc < 3
            @test abs(bestfit[:a, :bkg, :offset].val - 0.5) / bestfit[:a, :bkg, :offset].unc < 3
            @test isnan(bestfit[:a, :bkg, :x0].unc)
            @test abs(bestfit[:a, :bkg, :slope].val  - 0.1) / bestfit[:a, :bkg, :slope].unc < 3
            
            @test abs(bestfit[:b, :l1, :norm].val    - 0.8) / bestfit[:b, :l1, :norm].unc < 3
            @test abs(bestfit[:b, :l1, :center].val  - 2.1) / bestfit[:b, :l1, :center].unc < 3
            @test abs(bestfit[:b, :l1, :sigma].val   - 0.1) / bestfit[:b, :l1, :sigma].unc < 3
            @test abs(bestfit[:b, :l2, :norm].val    - 1.2) / bestfit[:b, :l2, :norm].unc < 3
            @test abs(bestfit[:b, :l2, :sigma].val   - 0.4) / bestfit[:b, :l2, :sigma].unc < 3
            @test abs(bestfit[:b, :bkg, :offset].actual - 0.5) < 0.2
            @test isnan(bestfit[:b, :bkg, :offset].unc)
            @test isnan(bestfit[:b, :bkg, :x0].unc)
            @test abs(bestfit[:b, :bkg, :slope].actual - 0.1) < 0.2
            @test isnan(bestfit[:b, :bkg, :slope].unc)
            
            @test fsumm.ndata == 242
            @test fsumm.nfree == 13
            @test isapprox(fsumm.fitstat, 1., atol=0.2)
        end
    end

    # ====================================================================
    @testset "8. Serialization & Deserialization" begin
        # Generate some quick data to serialize from the previous tests context
        x = 0:0.05:6
        model = Model(:l1 => GModelFit.Gaussian(1, 2, 0.2))
        data = GModelFit.mock(Measures, model, Domain(x), seed=1)
        bestfit, fsumm = fit(model, data)

        json_file = "test.json"
        TypedJSON.serialize(json_file, (bestfit, fsumm, data))
        bestfit2, fsumm2, data2 = TypedJSON.deserialize(json_file)

        for f in fieldnames(typeof(fsumm))
            (f == :solver_retval)  &&  continue
            @test getfield(fsumm, f) == getfield(fsumm2, f)
        end

        for (key, par) in GModelFit.getparams(bestfit)
            @test isapprox(bestfit[key...].val, bestfit2[key...].val)
        end

        # Cleanup
        rm(json_file, force=true)
    end
end
