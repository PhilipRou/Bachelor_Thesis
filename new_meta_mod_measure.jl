using Plots
using Statistics
using DelimitedFiles
#=
Die Phasen von Eichfeld-Konfigurationen Uμ( ⃗n) werden in einem Array
config[μ][nₜ][nₛ] gespeichert, wobei ⃗n = (nₜ, nₛ). Im weiteren soll
μ = 1 für die zeitliche und μ = 2 für die räumliche Dimension stehen.
(es wird auch μ ≡ 1 für temp. und ν ≡ 2 für räuml. synonym verwendet)
Phasen der Eichtransformationen Ω( ⃗n) werden als trans[nₜ, nₛ] gespeichert.

In For-Schleifen steht der Index i allermeistens für einen zeitlichen und
der Index j für einen räumlichen Index.
=#
N_t = 32
N_x = N_t
β = (N_t*N_x)/80
r = 2*pi*0.04
sweeps = 500000
sweeps_by_ten = Int64(0.1*sweeps)
sweeps_count = 0
cut = 62500
nb_accepted = 0
nb_metro = 1
nb_snapshots = Int64(sweeps/5)

δq = 0.01
w = 0.0001
Q_max = 15             # ⭕⭕⭕ muss von 15 auf 20 erhöht werden für N ≥ 36 ⭕⭕⭕
Q_thrs = 5
# counting = 80474137
boundary_count = 0

bb = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\Julia-stuff\\new_meta_bb_$N_t.txt")
bias_potential = bb[size(bb,1),:] .- bb[size(bb,1),1]

# rest_potential = zeros(size(bb,2))
# for i = 1:length(rest_potential)
#     Q_x = abs(-Q_max + (i-1)*δq)
#     if Q_x > Q_thrs
#         rest_potential[i] = (w/(δq^2))*(Q_x-Q_thrs)^2
#     end
# end



# Funktion um eine n-Instanton-Konfig zu bauen
function instanton(n)
new_config = [[],[]]
for i = 1:N_t
    push!(new_config[1], zeros(N_x))
    push!(new_config[2], n*i*2*pi/(N_t*N_x).*ones(N_x))
end
for j = 1:N_x
    new_config[1][N_t][j] = -n*j*2*pi/N_x
end
new_config
end

function P12(nt, nx, config)
NT = nt%N_t +1
NX = nx%N_x +1
return config[1][nt][nx] + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
end

function action(config)
a = 0.0
for i = 1:N_t
    for j = 1:N_x
        plaq = exp(im*P12(i,j,config))
        a += β*real(1-plaq)
    end
end
return a
end

# Funktion, um den adjungierten Staple in Richtung mu bei ⃗n = (nₜ,nₓ) einer
# gegebenen Konfig zu berechnen
function staple_dagger(mu, nt, nx, config)
a = 0
b = 0

NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
NX = nx%N_x +1  # nₛ+1 mit periodische Randbed.
TN = (nt + N_t -2)%N_t +1   # nₜ-1 mit periodische Randbed.
XN = (nx + N_x -2)%N_x +1   # nₛ-1 mit periodische Randbed.

if mu == 1
    a = -config[2][nt][nx] - config[1][nt][NX] + config[2][NT][nx]
    b = config[2][nt][XN] - config[1][nt][XN] - config[2][NT][XN]
elseif mu == 2
    a = -config[1][nt][nx] - config[2][NT][nx] + config[1][nt][NX]
    b = config[1][TN][nx] - config[2][TN][nx] - config[1][TN][NX]
end

exp(im*a) + exp(im*b)
end

# Funktion, um den Wirkungsunterschied einer gegebenen Konfig zu messen,
# deren Link U_old = config[mu][nt][nx] zu U_new geändert wurde
function Delta_S(U_old, U_new, mu, nt, nx, config)
stap = staple_dagger(mu, nt, nx, config)
-β * real((exp(im*U_new)- exp(im*U_old)) * stap)
end

function cont_charge(config)
Q = 0.0
for i = 1:N_t
    for j = 1:N_x
        Q += (1/(2*pi)) * sin(P12(i,j,config))
    end
end
return Q
end

function int_charge(config)
Q = 0.0
for i = 1:N_t
    for j = 1:N_x
        plaq = exp(im*P12(i,j,config))
        Q += (1/(2*pi)) * imag(log(plaq))
    end
end
return Q
end

# Funktion um den Unterschied in der kont. Ladung einer gegebenen Konfig
# zu bestimmen, nachdem U[mu][nt][nx] geändert wurde (U_old -> U_new)
function Delta_Q(U_old, U_new, mu, nt, nx, config)
NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
NX = nx%N_x +1  # nₛ+1 mit periodische Randbed.
TN = (nt + N_t -2)%N_t +1   # nₜ-1 mit periodische Randbed.
XN = (nx + N_x -2)%N_x +1   # nₛ-1 mit periodische Randbed.

if mu == 1
    a = U_new + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
    b = config[1][nt][XN] + config[2][NT][XN] - U_new - config[2][nt][XN]
    c = U_old + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
    d = config[1][nt][XN] + config[2][NT][XN] - U_old - config[2][nt][XN]
elseif mu == 2
    a = U_new - config[1][TN][NX] - config[2][TN][nx] + config[1][TN][nx]
    b = config[2][NT][nx] - config[1][nt][NX] - U_new + config[1][nt][nx]
    c = U_old - config[1][TN][NX] - config[2][TN][nx] + config[1][TN][nx]
    d = config[2][NT][nx] - config[1][nt][NX] - U_old + config[1][nt][nx]
end

#⭕  hier muss man die Realteile nicht mitschleppen
#(imag(exp(im*a) + exp(im*b) - exp(im*c) - exp(im*d))) / (2*pi)
(sin(a) + sin(b) - sin(c) - sin(d)) / (2*pi)
end

function grid_ind(q)
grid_index = (q+Q_max)/δq + 0.5000000001
return round(Int, grid_index)
end

function return_potential(q)
grid_index = grid_ind(q)
grid_q = -Q_max + (grid_index-1)*δq

(1-(q-grid_q)/δq)*bias_potential[grid_index] + ((q-grid_q)/δq)*bias_potential[grid_index+1]
end



##########   Thermalisierung   ###########



config_1 = []
for μ = 1:2
push!(config_1, [])
for i = 1:N_t
    rand_phase = 2*pi*(rand(N_x) .- 0.5)
    push!(config_1[μ], rand_phase)
end
end

raw_configs = [config_1]
insta = instanton(1)

for k = 1:cut
test = deepcopy(raw_configs[k])

for eo = 0:1
    for i = 1:N_t
        for j = 1:N_x>>1
            J = 2*j -eo
            old_link = deepcopy(test[1][i][J])
            proposal = deepcopy(test[1][i][J])
            proposal += r*2*(rand()-0.5)
            ΔS = Delta_S(old_link, proposal, 1, i, J, test)

            if rand() < exp(-ΔS)
                test[1][i][J] = proposal
                nb_accepted += 1
            end
        end
    end
    for j = 1:N_x
        for i = 1:N_t>>1
            I = 2*i -eo
            old_link = deepcopy(test[2][I][j])
            proposal = deepcopy(test[2][I][j])
            proposal += r*2*(rand()-0.5)
            ΔS = Delta_S(old_link, proposal, 2, I, j, test)

            if rand() < exp(-ΔS)
                test[2][I][j] = proposal
                nb_accepted += 1
            end
        end
    end
end
push!(raw_configs, test)
end

# Topologische Ladung auf 0 setzen
charge_therm = cont_charge(raw_configs[cut+1])
raw_configs[cut+1] .+= -instanton(charge_therm)

# Zweite Thermalisierung (nun im Null-Sektor)
for k = cut+1:2*cut
test = deepcopy(raw_configs[k])

for eo = 0:1
    for i = 1:N_t
        for j = 1:N_x>>1
            J = 2*j -eo
            old_link = deepcopy(test[1][i][J])
            proposal = deepcopy(test[1][i][J])
            proposal += r*2*(rand()-0.5)
            ΔS = Delta_S(old_link, proposal, 1, i, J, test)


            if rand() < exp(-ΔS)
                test[1][i][J] = proposal
                nb_accepted += 1
            end
        end
    end
    for j = 1:N_x
        for i = 1:N_t>>1
            I = 2*i -eo
            old_link = deepcopy(test[2][I][j])
            proposal = deepcopy(test[2][I][j])
            proposal += r*2*(rand()-0.5)
            ΔS = Delta_S(old_link, proposal, 2, I, j, test)
            # ΔS = Delta_S(old_link, proposal, 2, I, j, test) + bias_potential[ind_prop] - bias_potential[old_ind]


            if rand() < exp(-ΔS)
                test[2][I][j] = proposal
                nb_accepted += 1
            end
        end
    end
end
push!(raw_configs, test)
end



acceptance_therm = round(100*nb_accepted/(2*cut*2*N_t*N_x), digits = 2)
println("","Acceptance during thermalization: $acceptance_therm%")



#           ##########   Simulation   ###########



test = deepcopy(raw_configs[length(raw_configs)])
Q = cont_charge(test)
int_charges = []# ["discrete charges for N = $N_t"]
cont_charges = []#["continuous charges"]
# weights = ["weights exp(V(Q_cont))"]
actions = []#["actions"]
biases = []#["bias potential V at respective cont. charge"]
# avg_weights = []
plaquettes =[]# ["plaquettes"]
for k = 1:sweeps

Q = cont_charge(test)

for eo = 0:1
    for i = 1:N_t
        for j = 1:N_x>>1
            J = 2*j -eo
            old_link = deepcopy(test[1][i][J])
            proposal = deepcopy(test[1][i][J])
            proposal += r*2*(rand()-0.5)

            Q_prop = Q + Delta_Q(old_link, proposal, 1, i, J, test)
            old_ind = grid_ind(Q)
            new_ind = grid_ind(Q_prop)
            ΔS = Delta_S(old_link, proposal, 1, i, J, test) + bias_potential[new_ind] - bias_potential[old_ind]


            # if rand() < exp(rest_potential[old_ind]-rest_potential[new_ind])
            if rand() < exp(-ΔS)
                test[1][i][J] = proposal
                nb_accepted += 1
                Q = Q_prop
            end
        end
    end
end
#⭕  äusserste Schleife über Richtungen
for eo = 0:1
    for j = 1:N_x
        for i = 1:N_t>>1
            I = 2*i -eo
            old_link = deepcopy(test[2][I][j])
            proposal = deepcopy(test[2][I][j])
            proposal += r*2*(rand()-0.5)

            Q_prop = Q + Delta_Q(old_link, proposal, 2, I, j, test)
            old_ind = grid_ind(Q)
            new_ind = grid_ind(Q_prop)
            ΔS = Delta_S(old_link, proposal, 2, I, j, test) + bias_potential[new_ind] - bias_potential[old_ind]


            # if rand() < exp(rest_potential[old_ind]-rest_potential[new_ind])
            if rand() < exp(-ΔS)
                test[2][I][j] = proposal
                nb_accepted += 1
                Q = Q_prop
            end
        end
    end
end
push!(int_charges, int_charge(test))                    #⭕
push!(cont_charges, cont_charge(test))                  #⭕
push!(actions, action(test))                            #⭕
# push!(plaquettes, P12(test))                            #⭕
push!(biases, bias_potential[grid_ind(Q)])              #⭕


# top. Update:
if k%nb_metro == 0
    proposal = deepcopy(test)
    P = rand()
    proposal .+= sign(P-0.5) .* insta

    Q_prop = cont_charge(proposal)
    old_ind = grid_ind(Q)
    ind_prop = grid_ind(Q_prop)
    ΔS = action(proposal) - action(test) + bias_potential[ind_prop] - bias_potential[old_ind]

    if rand() < exp(-ΔS)
        test = proposal
        Q = Q_prop
        # KEIN UPDATE DES BIAS POTENTIAL MEHR
    end
end

# for the long runs, so you know how much time has passed:
if (k)%(sweeps_by_ten) == 0
    sweeps_count += 10
    println("We're already $sweeps_count% deep in the simulation!")
end
end



#               ##########   Berechnen   ###########



weights = exp.(biases)

# Funktion um gewichtetes Mittel eines Arrays x zu bestimmen
function get_weighted_mean(x)
return sum(x .* weights )/sum(weights)
end

function bootstrap(K, observable)
bootstraps = []                     # Array für die K Mittelwerte
for i = 1:K
    boots = []
    for j = 1:length(observable)    # length(observable) = N aus Kap. 4.1.3
        push!(boots, observable[rand(1:length(observable))])
    end
    push!(bootstraps, mean(boots))
end
O_tilda = mean(bootstraps)
σ = sqrt(mean((bootstraps .- O_tilda).^2))
return [O_tilda, σ]
end

function bootstrap_weighted(K, observable)
bootstraps = []                     # Array für die K Mittelwerte
for i = 1:K
    boots = []
    boot_weights = []
    for j = 1:length(observable)    # length(observable) = N aus Kap. 4.1.3
        rand_num = rand(1:length(observable))
        push!(boot_weights, weights[rand_num])
        push!(boots, observable[rand_num]*weights[rand_num])
    end
    push!(bootstraps, sum(boots)/sum(boot_weights))
end
O_tilda = mean(bootstraps)
σ = sqrt(mean((bootstraps .- O_tilda).^2))
return [O_tilda, σ]
end



acceptance = round(100*nb_accepted/((sweeps+2*cut)*2*N_t*N_x), digits = 2)
int_charges = round.(Int, int_charges)
q_sq_mean = get_weighted_mean(int_charges.^2)
q_sq_boot = bootstrap_weighted(50, int_charges.^2)
println("Acceptance total: $acceptance%")
println("   ", "q_sq_mean = ", round(q_sq_mean, digits=3))
println("   ","q_sq_boot = ", round(q_sq_boot[1], digits=3), " ± ", round(q_sq_boot[2], digits=3))
# println("Boundary count: $boundary_count")
println("\n")



#               ##########   Plotten   ###########



histo_int = histogram(
int_charges, bins = -30.25:0.5:31, c = :salmon, background_color = :grey,
bar_edges = :false)

histo_int = plot!(
size = (750,600), title = "Discrete top. charges on frozen bias potential
$N_t×$N_x, β = $β, $sweeps sweeps",
xticks = -26:2:26,
xlabel = "Discrete Q",
legend = :false,
labelfontsize = 16,
# legendfontsize = 12,
titlefontsize = 16
)
display(histo_int)


histo_cont = histogram(
cont_charges, bins = -20:0.1:20, c = :salmon, background_color = :grey,
bar_edges = :false)
histo_cont = plot!(
size = (750,600), title = "Continuous top. charges on frozen bias potential
$N_t×$N_x, β = $β, $sweeps sweeps",
xticks = -20:2:20,
xlabel = "Continuous Q",
legend = :false,
labelfontsize = 16,
# legendfontsize = 12,
titlefontsize = 16
)

display(histo_cont)

#=
histo_weights = histogram(
weights, #= bins = -8:0.1:8, =# c = :salmon, background_color = :grey,
bar_edges = :false)

histo_int = plot!(
size = (750,600), title = "Weights: e^(-V)
$N_t×$N_x, β = $β, $sweeps sweeps",
# xticks = -8:1:8,
xlabel = "Weights",
legend = :false,
labelfontsize = 16,
# legendfontsize = 12,
titlefontsize = 16
)

display(histo_weights)


histo_weights_inv = histogram(
weights.^-1, #= bins = -8:0.1:8, =# c = :salmon, background_color = :grey,
bar_edges = :false)

histo_int = plot!(
size = (750,600), title = "Inverted Weights: e^(V)
$N_t×$N_x, β = $β, $sweeps sweeps",
# xticks = -8:1:8,
xlabel = "Inverted Weights",
legend = :false,
labelfontsize = 16,
# legendfontsize = 12,
titlefontsize = 16
)

display(histo_weights_inv)
=#

wab = deepcopy(weights) ./ sum(weights)

histo_int_weighted = histogram(
int_charges, bins = -30.25:0.5:31, c = :salmon, background_color = :grey,
bar_edges = :false, weights = Float64.(wab))
histo_int_weighted = plot!(
size = (750,600), title = "Discrete top. charges on frozen bias potential, WEIGHTED
$N_t×$N_x, β = $β, $sweeps sweeps",
xticks = -26:2:26,
xlabel = "Discrete Q",
legend = :false,
labelfontsize = 16,
# legendfontsize = 12,
titlefontsize = 16
)
display(histo_int_weighted)



plot(
cont_charges, weights,
seriestype = :scatter,
legend = :false,
xticks = -8:2:8,
size = (750,600)
)

plot(1:length(cont_charges), cont_charges, size=(750,600))
plot(1:length(cont_charges), int_charges, size=(750,600))

savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\new_meta_mod_measure_$N_t.cont.$N_t.")
