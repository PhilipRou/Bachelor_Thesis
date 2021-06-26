using Optim
using DelimitedFiles


# Der Hauptteil des Programms befindet sich in einer for-Schleife, sodass
# bei einmaligem Ausführen die Renormierungskonstanten aller Gitter
# von N = 12 bis N = 60 bestimmt werden. Die Ergebnisse werden einfach
# im REPL ausgegeben (es sind ja nur 5 Werte)
for i = 1:5
# Wir arbeiten ohnehin mit quadratischen Gittern:
N = 12*i
N_t = N
N_x = N

# Daten der Simulation einlesen:
# Hinter dem letzten \\ steht der Name des in topup_sim.jl erstellten
# Datei, davor ist der Dateipfad der Datei angegeben (welcher gleich dem
# von topup_sim.jl sein muss). Hierbei müssen \\ statt \ verwendet werden.
file = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\Julia-stuff\\compare_topup_$N_t.x.$N_x.txt")
cont_charges = file[:,2]


# Die Minimumsstelle dieser Funktion ist der Wert der Renormierungskonstante
function chi_sq(Z::Vector)
    sum((Z[1]*cont_charges - round.(Z[1]*cont_charges)).^2)
end

# Bestimmung der Miminumsstelle (minimizer)
lower = [1.0]
upper = [2.0]
initial = [1.5]
inner_optimizer = GradientDescent()
result = optimize(chi_sq, lower, upper, initial, Fminbox(inner_optimizer))
Z_LSQ = Optim.minimizer(result)

println("","N = $N:    Z = ", round(Z_LSQ[1], digits=3))



end
