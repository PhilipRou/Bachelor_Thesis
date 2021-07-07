using Plots
using LsqFit

Z = [1.678, 1.29, 1.086, 1.033, 1.018, 1.011]
uncer = [0.034, 0.015, 0.003, 0.001, 0.001, 0.001]
N = [12, 16, 24, 36, 48, 60]
a = [1/12, 1/16, 1/24, 1/36, 1/48, 1/60]


# function f(x, Qmax)
#     k = 2*pi/(x^2)
#     return 2 * (sin(k*Qmax) - k*Qmax*cos(k*Qmax)) / (k*Qmax - sin(k*Qmax)*cos(k*Qmax))
# end
#
# yvals = []
#
# for i = 1:6
#     push!(yvals, f(N[i], 5))
# end
#
# image2 = plot(N, Z, seriestype = :scatter, size = (750,600), yerror = uncer)
# image2 = plot!(N, yvals, xticks = 12:12:60)


# Fit erstellen:
model(x, p) = p[1] .+ p[2] .*x .+ p[3] .* x.^2 .+ p[4] .* x.^3
p0 = [1.0, 1.0, 1.0, 1.0]
fit = curve_fit(model, a, Z, 1 ./uncer.^2, p0)
χ_sq = round(sum(fit.resid.^2), digits = 3)
degof = dof(fit)
χ_sq_dof = round(sum((fit.resid.^2))/dof(fit), digits = 3)      # χ²/dof

fit_x = []
fit_vals = []
P = fit.param
for i = 0:9
    x = i/100
    push!(fit_x, x)
    push!(fit_vals, model(x,P))
end

A = round(P[1], digits=3)
B = round(P[2], digits=2)
C = round(P[3], digits=1)
D = round(P[4], digits=1)
A_err = round(stderror(fit)[1], digits = 3)
B_err = round(stderror(fit)[2], digits = 2)
C_err = round(stderror(fit)[3], digits = 1)
D_err = round(stderror(fit)[4], digits = 1)


image = plot(a, Z, yerror = uncer, label = "Z-values", seriestype = :scatter)
image = plot!(fit_x, fit_vals, label = "fit: f(x) = A + Bx + Cx³ + Dx⁵")
image = plot!(
title = "χ²/#dof = $χ_sq / $degof = $χ_sq_dof
A = $A ± $A_err,   B = $B ± $B_err,
C = $C ± $C_err,   D = $D ± $D_err",
size=(750,600),
# xticks = int_min:int_max,
# yticks = int_min:int_max,
xlabel = "lattice spacing a ",
# ylabel = "Continuous Q",
legend = :topleft,
labelfontsize = 16,
legendfontsize = 12,
titlefontsize = 16
)


# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\topup_ZoverA_polyfit_odd")
