using Random, Test, GFit, .PV

# Test Model
mm = PVModel{Float64}()
@assert typeof(mm) == PVModel{Float64}

mm[:comp] # empty component
mm[:comp3][:par] = 3.1
mm[:comp1][:par] = 1.1
mm[:comp1][:alt] = 1.2
mm[:comp3].par = 99

@assert PV.internal_vector(mm) == [99, 1.1, 1.2]

@assert mm[:comp1].par == 1.1
@assert mm[:comp1].alt == 1.2
@assert mm[:comp3].par == 99

@assert keys(mm) == [:comp, :comp3, :comp1]

@assert propertynames(mm[:comp])  == Symbol[]
@assert propertynames(mm[:comp1]) == [:par, :alt]
@assert propertynames(mm[:comp3]) == [:par]

@assert items(mm) == [99, 1.1, 1.2]
@assert items(mm[:comp]) == Float64[]
@assert items(mm[:comp1]) == [1.1, 1.2]
@assert items(mm[:comp3]) == [99]

for (cname, comp) in mm
    println(cname)
    for (pname, par) in comp
        println("  ", pname, ": ", par)
        @assert getproperty(mm[cname], pname) == par
    end
end


# Test multimodel
MM = PVMulti{Float64}()
push!(MM, mm)
push!(MM, deepcopy(mm))
@assert typeof(MM[1]) == PVModel{Float64}
@assert items(MM) == [99., 1.1, 1.2, 99., 1.1, 1.2]

MM[1][:comp1].par = 12345
@assert items(MM) == [99., 12345., 1.2, 99., 1.1, 1.2]

set_items!(MM, [1,2,3,4,5,6])
@assert items(MM) == [1,2,3,4,5,6]

for i in 1:length(MM)
    for (cname, comp) in MM[i]
        println(cname)
        for (pname, par) in comp
            println("  ", pname, ": ", par)
            @assert getproperty(MM[i][cname], pname) == par
        end
    end
end


# ====================================================================
x    = [0.1, 1.1, 2.1, 3.1, 4.1]
meas = [6.29, 7.27, 10.41, 18.67, 25.3]
unc  = [1.1, 1.1, 1.1, 1.2, 1.2]

domain = Domain(x)
data = Measures(domain, meas, unc)
model = Model(domain, @λ (x, a2=1, a1=1, a0=5) -> (a2 .* x.^2  .+  a1 .* x  .+  a0))
res = fit!(model, data)
# @gp data model


# ====================================================================
x = 0:0.1:5
model = Model(Domain(x), :a2 => 1, :a1 => 1, :a0 => 5)
model[:parabola] = @λ (x, a2, a1, a0) -> @. (a2 * x^2  +  a1 * x  + a0)
res = fit!(model, GFit.mock(Measures, model, seed=1))
# @gp x y "w l t 'True model'" x values(data) uncerts(data) "w yerr t 'Data'" x model() "w l t 'Best fit'"
# viewer(model, data, res)


# ====================================================================
x = 0:0.1:5
model = Model(Domain(x), @λ (x, a2=1, a1=1, a0=5) -> @. (a2 * x^2  +  a1 * x  + a0))
res = fit!(model, GFit.mock(Measures, model, seed=1))


# ====================================================================
f = @λ (x, p1=1, p2=1.e-3, p3=1e-6, p4=4, p5=5) ->
@. (p1 + p2 * x + p3 * x^2 + p4 * sin(p5 * x))  *  cos(x)
x = 1.:50:10000
model = Model(Domain(x), f)
res = fit!(model, GFit.mock(Measures, model, seed=1))


# ====================================================================
f1 = @λ (x, p1=1, p2=1e-3, p3=1e-6) -> @.  p1  +  p2 * x  +  p3 * x^2
f2 = @λ (x, p4=4, p5=5) -> @. p4 * sin(p5 * x)
f3 = @λ (x) -> cos.(x)

x = 1.:50:10000
model = Model(Domain(x),
              :f1 => f1, 
              :f2 => f2,
              :f3 => f3,
              :main => @λ (x, f1, f2, f3) -> (f1 .+ f2) .* f3)
res = fit!(model, GFit.mock(Measures, model, seed=1))

# Same results with
model = Model(Domain(x), :f1 => f1)
model[:f2] = f2
model[:f3] = f3
model[:main] = @λ (x, f1, f2, f3) -> (f1 .+ f2) .* f3
res = fit!(model, GFit.mock(Measures, model, seed=1))
#=
GFitViewer.save_json(model, data, res, filename="test1.json")
GFitViewer.save_json(model, data     , filename="test1a.json")
GFitViewer.save_json(model           , filename="test1b.json")
=#


# ====================================================================
x = 0:0.05:6
model = Model(Domain(x),
              :l1  => GFit.Gaussian(1, 2, 0.2),
              :l2  => GFit.Gaussian(1, 3, 0.5),
              :bkg => GFit.OffsetSlope(0.5, 1, 0.1),
              :main => SumReducer(:l1, :l2, :bkg));
res = fit!(model, GFit.mock(Measures, model, seed=1))


# Tie two parameters
model[:l2].norm.patch = :l1
res = fit!(model, GFit.mock(Measures, model, seed=1))

# Patch one parameter to another via a λ function
model[:l2].norm.patch = @λ (v, m) -> v + m[:l1].norm
res = fit!(model, GFit.mock(Measures, model, seed=1))
#=
GFitViewer.save_json(model, data, res, filename="test2.json")
GFitViewer.save_json(model, data     , filename="test2a.json")
GFitViewer.save_json(model           , filename="test2b.json")
=#


# ====================================================================
x = 0:0.05:6
model1 = Model(Domain(x),
               :l1  => GFit.Gaussian(1, 2, 0.2),
               :l2  => GFit.Gaussian(1, 3, 0.5),
               :bkg => GFit.OffsetSlope(0.5, 1, 0.1),
               :main => SumReducer(:l1, :l2, :bkg));

model2 = Model(Domain(x),
               :l1  => GFit.Gaussian(0.8, 2.1, 0.1),
               :l2  => GFit.Gaussian(1.2, 2.5, 0.4),
               :bkg => GFit.OffsetSlope(0.5, 1, 0.1),
               :main => SumReducer(:l1, :l2, :bkg));

model = MultiModel(model1, model2)
freeze!(model[1], :bkg);
freeze!(model[2], :bkg);
data = GFit.mock(Measures, model, seed=1)
res = fit!(model, data, GFit.cmpfit())
GFit.print_param_covariance(res, sort=true, select=["[2][l1].norm"])



thaw!(model[1], :bkg);
thaw!(model[2], :bkg);

model[2][:bkg].offset.mpatch = @λ m -> m[1][:bkg].offset
model[2][:bkg].slope.mpatch  = @λ m -> m[1][:bkg].slope


model[1][:l2].center.mpatch = @λ m -> m[2][:l2].center

@time res = fit!(model, data)

#=
@gp x y1 "w l t 'True model'" x values(data1) uncerts(data1) "w yerr t 'Data'" x model1() "w l t 'Best fit'"
@gp x y2 "w l t 'True model'" x values(data2) uncerts(data2) "w yerr t 'Data'" x model2() "w l t 'Best fit'"
viewer(model, [data1, data2], res)

GFitViewer.save_json(model, [data1, data2], res, filename="test3.json")
GFitViewer.save_json(model, [data1, data2]     , filename="test3a.json")
GFitViewer.save_json(model                     , filename="test3b.json")
=#
