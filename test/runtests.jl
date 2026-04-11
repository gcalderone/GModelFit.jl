using Test
using GModelFit
using GModelFit.Solvers
using Random
using DataFrames
using Statistics

@testset "GModelFit.jl Full Test Suite" begin

    @testset "1. Domains, GaussianData & Response API" begin
        # --- AbstractDomain, Domain, CartesianDomain ---
        x1 = [0.1, 1.1, 2.1, 3.1]
        dom1 = Domain(x1)
        @test dom1 isa AbstractDomain
        @test ndims(dom1) == 1
        @test length(dom1) == 4
        @test coords(dom1, 1) == x1

        x2, y2 = 1:5, 1:10
        dom2 = CartesianDomain(x2, y2)
        @test dom2 isa AbstractDomain
        @test ndims(dom2) == 2
        @test size(dom2) == (5, 10)
        @test gridcoords(dom2, 1) == collect(x2)
        @test coords(Domain(dom2), 1) == repeat(collect(float.(x2)), 10) # Testing the fallback cast

        # --- GaussianData, uncerts, domain, response ---
        vals = rand(4)
        uncs = fill(0.1, 4)
        meas = GaussianData(dom1, vals, uncs)

        @test values(meas) == vals
        @test uncerts(meas) == uncs
        @test length(meas) == 4
        @test domain(meas) === dom1
        @test getresp(meas) isa GModelFit.AbstractResponse
    end

    @testset "2. Components, @fd, & Dependencies" begin
        # --- @fd, FunctDesc ---
        f = @fd (x, a=1, b=2) -> a .* x .+ b
        @test f(10.0, 2.0, 5.0) == 25.0

        # --- SumReducer, stype, dependencies ---
        m_dep = Model{Tuple{}}()
        m_dep[:compA] = GModelFit.FComp(@fd (x, a=1) -> a)
        m_dep[:compB] = GModelFit.FComp(@fd (x, b=1) -> b)
        m_dep[:main]  = SumReducer(:compA, :compB)

        compA = getcomps(m_dep)[(:compA,)]
        main_comp = getcomps(m_dep)[(:main,)]

        # stype and getcomp
        @test getcomp(compA) isa GModelFit.FComp
        @test typeof(stype(compA)) == String
        @test stype(main_comp) == "SumReducer"

        # dependencies
        @test :compA in dependencies(main_comp)
        @test :compB in dependencies(main_comp)
        @test isempty(dependencies(compA))
    end

    @testset "3. Model Structure & Accessors" begin
        # --- Model construction ---
        m = Model(@fd (x, slope=2, int=1) -> slope .* x .+ int)
        m[:bkg] = GModelFit.FComp(@fd (x, bg=0.1) -> bg)

        # --- getgroups, getcomps, getparams, getfreepars ---
        @test length(getgroups(m)) == 1
        @test length(getcomps(m)) == 2
        @test length(getparams(m)) == 3

        meval = ModelEval(Domain(1:5), m)
        @test length(getfreepars(meval)) == 3 # slope, int, bg are all free by default

        # --- flattenkeys ---
        d = Dict((:group1, :comp1) => 42)
        flat = flattenkeys(d, "_")
        @test flat[:group1_comp1] == 42
    end

    @testset "4. Model Manipulation (Freeze/Thaw)" begin
        m = Model(@fd (x, slope=2, int=1) -> slope .* x .+ int)
        comp = getcomps(m)[(:main,)]

        # --- isfrozen, freeze!, thaw! ---
        @test !isfrozen(comp)
        freeze!(comp)
        @test isfrozen(comp)

        meval = ModelEval(Domain(1:5), m)
        @test length(getfreepars(meval)) == 0 # Everything frozen at component level

        thaw!(comp)
        @test !isfrozen(comp)
        meval2 = ModelEval(Domain(1:5), m)
        @test length(getfreepars(meval2)) == 2 # Thawed, parameters are free again
    end

    @testset "5. Likelihood, Stats & Fitting Routines" begin
        # Generate Ground Truth & Mock Data
        x = 1.0:1.0:10.0
        dom = Domain(x)
        truth = Model(@fd (x, slope=2.0, int=5.0) -> slope .* x .+ int)
        data = GModelFit.mock(GaussianData, dom, truth, seed=42, properr=0.01)

        # --- Likelihood ---
        guess = Model(@fd (x, slope=1.0, int=1.0) -> slope .* x .+ int)
        lh = Likelihood(data, guess)

        # --- loglh, gofstat, dof ---
        initial_loglh = loglh(lh)
        initial_gof = gofstat(lh)

        @test initial_loglh < 0 # log-likelihood is negative
        @test initial_gof > 0   # Chi-squared is positive
        @test dof(lh) == 10 - 2 # 10 data points, 2 free parameters

        # --- fit vs fit! ---
        # fit returns a new likelihood, fit! modifies in place
        bestfit_lh, fsumm = fit(data, guess)
        @test gofstat(bestfit_lh) < initial_gof # Fit should improve statistic

        # Test fit! directly on the likelihood
        fsumm_inplace = fit!(lh)
        @test isapprox(gofstat(lh), gofstat(bestfit_lh), rtol=1e-3)
    end

    @testset "6. Solvers Integration" begin
        x = 1.0:1.0:10.0
        dom = Domain(x)
        truth = Model(@fd (x, slope=2.0) -> slope .* x)
        data = GModelFit.mock(GaussianData, dom, truth, seed=1)

        function get_fresh_lh()
            guess = Model(@fd (x, slope=1.0) -> slope .* x)
            return Likelihood(data, guess, autodiff=true)
        end

        # --- lsqfit (Default) ---
        lh_lsq = get_fresh_lh()
        stat_lsq = maximize!(lh_lsq, Solvers.lsqfit())
        @test stat_lsq isa SolverStatusOK
        @test isapprox(getparams(lh_lsq)[(:main, :slope)].val, 2.0, atol=0.1)

        # --- cmpfit ---
        lh_cmp = get_fresh_lh()
        stat_cmp = maximize!(lh_cmp, Solvers.cmpfit())
        # Depending on tolerances, it might return OK or Warn(Status 5/2), both valid fit states
        @test stat_cmp isa SolverStatusOK || stat_cmp isa SolverStatusWarn
        @test isapprox(getparams(lh_cmp)[(:main, :slope)].val, 2.0, atol=0.1)

        # --- curvefit ---
        lh_curv = get_fresh_lh()
        stat_curv = maximize!(lh_curv, Solvers.curvefit())
        @test stat_curv isa SolverStatusOK
        @test isapprox(getparams(lh_curv)[(:main, :slope)].val, 2.0, atol=0.1)

        # --- dry solver ---
        lh_dry = get_fresh_lh()
        stat_dry = maximize!(lh_dry, Solvers.dry())
        @test stat_dry isa SolverStatusWarn # Dry always warns
        @test getparams(lh_dry)[(:main, :slope)].val == 1.0 # Values should not change
    end

    @testset "7. Snapshots & DataFrame Exports" begin
        # --- snapshot ---
        m = ModelEval(Domain(1:5), Model(@fd (x, a=1.0) -> a .* x))
        snap = snapshot(m)
        @test snap isa GModelFit.ModelSnapshot
        @test getcomps(snap)[(:main,)].evalbuf == [1.0, 2.0, 3.0, 4.0, 5.0]

        # Verify detachment (modifying original shouldn't modify snapshot)
        getparams(m)[(:main, :a)].val = 99.0
        @test getparams(snap)[(:main, :a)].val == 1.0

        # --- DataFrame Exports ---
        truth = Model(@fd (x, slope=2.0) -> slope .* x)
        data = GModelFit.mock(GaussianData, Domain(1.0:5.0), truth, seed=1)
        lh = Likelihood(data, truth)

        df_params = export_params(lh)
        @test df_params isa DataFrame
        @test nrow(df_params) == 1
        @test df_params.pname[1] == "(:main, :slope)"

        df_data = GModelFit.export_data(data)
        @test df_data isa DataFrame
        @test :data_vals in propertynames(df_data)

        df_folded = export_folded_eval(lh)
        @test df_folded isa DataFrame
        @test :model in propertynames(df_folded)
        @test :resid in propertynames(df_folded)

        df_unfolded = export_unfolded_eval(lh)
        @test df_unfolded isa DataFrame
        @test :main in propertynames(df_unfolded)
    end
end

@testset "PoissonCounts and Mock Data" begin
    @testset "Constructor and Type Handling" begin
        # Test that the constructor correctly handles Arrays of Integers
        # and safely converts them to Float64 internally.
        counts_int = [10, 0, 5, 2]
        dom = Domain(1:4)

        data = PoissonCounts(dom, counts_int)

        @test data isa PoissonCounts
        @test values(data) isa Vector{Float64}
        @test values(data) == [10.0, 0.0, 5.0, 2.0]
    end

    @testset "Mock Data Generation" begin
        dom = Domain(1:500)
        # Flat model predicting exactly 20.0 counts everywhere
        truth_model = Model(:main => @fd((x, bg=20.0) -> bg))
        meval = ModelEval(dom, truth_model)

        # Generate mock data
        data = GModelFit.mock(PoissonCounts, meval; seed=123)

        @test data isa PoissonCounts
        @test length(values(data)) == 500

        # Check that generated values are mathematically integers (even if typed as Float64)
        @test all(isinteger.(values(data)))

        # For a Poisson distribution with lambda = 20, the mean and variance should both be ~20.
        # We test with a generous tolerance to avoid random test failures (flakiness).
        v = values(data)
        @test isapprox(Statistics.mean(v), 20.0, atol=1.0)
        @test isapprox(var(v), 20.0, atol=2.5)
    end

    @testset "Cash Statistic (C-stat) and Log-Likelihood" begin
        # We will test against analytically known C-stat values.
        # Formula: C_i = 2 * (M_i - D_i + D_i * ln(D_i / M_i)) (if D_i > 0)
        #          C_i = 2 * M_i                               (if D_i == 0)

        # Model evaluates to M = [1.0, 2.0, 0.5]
        model = Model(:main => @fd((x) -> x))
        dom = Domain([1.0, 2.0, 0.5])

        # Data observed is D = [0, 2, 1]
        data = PoissonCounts(dom, [0, 2, 1])
        lh = Likelihood(data, model)

        # Analytical calculation:
        # Bin 1 (D=0, M=1.0): C_1 = 2 * (1.0) = 2.0
        # Bin 2 (D=2, M=2.0): C_2 = 2 * (2.0 - 2.0 + 2.0 * ln(2.0/2.0)) = 0.0
        # Bin 3 (D=1, M=0.5): C_3 = 2 * (0.5 - 1.0 + 1.0 * ln(1.0/0.5)) = 2 * (-0.5 + 0.693147) ≈ 0.386294
        # Total C-stat ≈ 2.386294

        expected_cstat = 2.0 + 0.0 + 2.0 * (-0.5 + log(1.0 / 0.5))

        @test gofstat(lh) ≈ expected_cstat

        # Log-likelihood should exactly map to -0.5 * gofstat
        @test loglh(lh) ≈ -0.5 * expected_cstat
    end

    @testset "Anscombe Residuals" begin
        # Re-using the analytical setup from above
        model = Model(:main => @fd((x) -> x))
        dom = Domain([1.0, 2.0, 0.5])
        data = PoissonCounts(dom, [0, 2, 1])
        lh = Likelihood(data, model)

        # Manually trigger the residual evaluation on the group
        group = lh[()]
        GModelFit.evaluate_resid!(group)
        res = group.residuals

        # Residual formula: sign(D - M) * sqrt(C_i)

        # Bin 1 (D=0, M=1.0): C_1 = 2.0. sign(0 - 1.0) = -1.0. Res = -sqrt(2.0)
        @test res[1] ≈ -sqrt(2.0)

        # Bin 2 (D=2, M=2.0): C_2 = 0.0. sign(2 - 2.0) = 0.0. Res = 0.0
        @test res[2] ≈ 0.0

        # Bin 3 (D=1, M=0.5): C_3 ≈ 0.386. sign(1 - 0.5) = +1.0. Res = +sqrt(C_3)
        c_3 = 2.0 * (0.5 - 1.0 + 1.0 * log(1.0 / 0.5))
        @test res[3] ≈ sqrt(c_3)
    end
end
