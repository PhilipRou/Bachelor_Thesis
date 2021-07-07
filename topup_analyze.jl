using Plots
using Statistics
using DelimitedFiles
using OffsetArrays
using LsqFit

# Wir arbeiten ohnehin mit quadratischen Gittern:
N = 24
N_t = N
N_x = N

# Daten der Simulation einlesen:
# Hinter dem letzten \\ steht der Name des in topup_sim.jl erstellten
# Datei, davor ist der Dateipfad der Datei angegeben (welcher gleich dem
# von topup_sim.jl sein muss). Hierbei müssen \\ statt \ verwendet werden.
file = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\Julia-stuff\\compare_topup_$N_t.x.$N_x.txt")
nb_meas = size(file, 1)

# Arrays mit allen diskreten und kont. Ladungen
int_charges = round.(Int, file[:,1])
cont_charges = file[:,2]


##### Nun zum Sortieren der kont. Ladungen nach ihren jeweiligen diskr. Ladungen:

# Es wird ein Array sorted erstellt, welches für jede diskr. Ladung einen
# Array beherbergt, in dem die zugehörigen kont. Ladungen enthalten sind.
# Mit OffsetArray werden diese Arrays einen Index haben, der gleich der
# jeweils zugehörigen diskr. Ladung ist
# (Beispiel: sorted[-3] ist der Array, der alle kontinuierlichen Ladungen
# beinhaltet, deren diskr. Ladung -3 ist)

int_min = minimum(int_charges)
int_max = maximum(int_charges)
int_total = abs(int_min) + abs(int_max) +1

sorted = []
for k = 1:int_total
    push!(sorted, [])
end
sorted = OffsetArray(sorted, int_min:int_max)

for k = 1:nb_meas
    for q = int_min:int_max
        if file[k,1] == q
            push!(sorted[q], file[k,2])
        end
    end
end



##### Bestimmung der Renormierungskonstante über Fits

# Im Folgenden

means = []
uncer = []

for q = int_min:int_max
    push!(means, mean(sorted[q]))
    push!(uncer, std(sorted[q])/sqrt(length(sorted[q])))
end


# Fit erstellen:
model(x, p) = p[1]*x
p0 = [1.0]
fit = curve_fit(model, int_min:int_max, means, 1 ./uncer.^2, p0)
χ_sq = round(sum(fit.resid.^2), digits = 3)
degof = dof(fit)
χ_sq_dof = round(sum((fit.resid.^2))/dof(fit), digits = 3)      # χ²/dof
α = fit.param[1]
α_err = stderror(fit)[1]
Z = 1/α
Z_err = Z^2*α_err

α = round(α, digits = 3)
α_err = round(α_err, digits = 3)
Z = round(Z, digits = 3)
Z_err = round(Z_err, digits = 3)

image = plot(int_min:1:int_max, means, yerror = uncer, seriestype = :scatter)
image = plot!([int_min,int_max], [model(int_min, α), model(int_max, α)])
image = plot!(
title = "$N_t×$N_x lattice,   χ²/#dof = $χ_sq / $degof = $χ_sq_dof
slope = $α ± $α_err   ⇒   Z = $Z ± $Z_err",
size=(750,600),
xticks = int_min:int_max,
yticks = int_min:int_max,
xlabel = "Discrete Q",
ylabel = "Continuous Q",
legend = :false,
labelfontsize = 16,
legendfontsize = 12,
titlefontsize = 16
)

display(image)


# Der folgende Befehl speichert die aktuelle Grafik im Dateipfad
# C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis
# unter dem Namen topup_analyze_N_t.x.N_x. wobei N_t und N_x durch die
# jeweiligen Werte ersetzt werden (also z.B. "24.x.24.")

savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\topup_analyze_$N_t.x.$N_x.")
