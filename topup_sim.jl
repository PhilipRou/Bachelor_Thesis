# Die Phasen von Eichfeld-Konfigurationen Uμ( ⃗n) werden in einem Array
# config[μ][nₜ][nₛ] gespeichert, wobei ⃗n = (nₜ, nₛ). Im weiteren soll
# μ = 1 für die zeitliche und μ = 2 für die räumliche Dimension stehen.
# (es wird auch μ ≡ 1 für temp. und ν ≡ 2 für räuml. synonym verwendet)
# Phasen der Eichtransformationen Ω( ⃗n) werden als trans[nₜ, nₛ] gespeichert.

using Plots
using Statistics
using DelimitedFiles


N_t = 12        # Anzahl temporaler Gitterpunkte
N_x = 12    	# Anzahl räumlicher Gitterpunkte
β = (N_x*N_t)/80      # 1/g², wobei g die Eich-Kopplung ist
r = 2*pi*0.1      # Schrittweite im Metropolis (nur in pos. Richtung, s.u.)
nb_sim = 4000   # Anzahl der Update-Durchläufe. Pro Durchlauf werden:
nb_metro = 3    # nb_metro Metropolis-Updates und
nb_top = 1      # nb_top topologische-Updates durchgeführt.
Nb_updates = nb_sim*(nb_metro+nb_top)   # Anzahl insg. durchgeführter Updates
cut = Int64(0.6*Nb_updates) # Thermalization
N_skip = 20         # Nur jede N_skip-te Konfig wird gemessen
link_accepted = 0   # Misst die Akzeptanz-Rate der Link-Updates im Metropolis
R = r/(2*pi)


# Die erste Feldkonfig. wird per Zufall generiert (hot start), mit Werten ∈ [0,2π]
config_1 = []
for mu = 1:2
    push!(config_1, [])
    for i = 1:N_t
        a = (2*pi).*rand(N_x)
        push!(config_1[mu], a)
    end
end

# Funktion um einen adjungierten Staple in mu-Richtung an der Stelle
# ⃗n = (nₜ,nₛ) zu bauen
function staple_dagger(mu, nt, nx, config)
    a = 0
    b = 0

    NT = nt%N_t +1
    NX = nx%N_x +1
    TN = (nt + N_t -2)%N_t +1
    XN = (nx + N_x -2)%N_x +1

    if mu == 1
        a = -config[2][nt][nx] - config[1][nt][NX] + config[2][NT][nx]
        b = config[2][nt][XN] - config[1][nt][XN] - config[2][NT][XN]
    elseif mu == 2
        a = -config[1][nt][nx] - config[2][NT][nx] + config[1][nt][NX]
        b = config[1][TN][nx] - config[2][TN][nx] - config[1][TN][NX]
    end

    exp(im*a) + exp(im*b)
end


# Metropolis einer Konfiguration "config"
# Das Ergebnis ist ein Array, dessen erster Eintrag die geupdatete Konfig
# und dessen zweiter Eintrag die Anzahl der akzeptierten Link-Updates ist
function metro_update(config)
    nb_updates = 0
    test = deepcopy(config)
    for eo = 0:1
        for i = 1:N_t
            for j = 1:Int64(0.5*N_x)
                J = 2*j -eo
                a = test[1][i][J]
                b = staple_dagger(1,i,J,test)
                proposal = deepcopy(a)
                q = 2*(rand(1)[1] -0.5)
                proposal += r*q

                Delta = -β * real((exp(im*proposal)- exp(im*a)) * b)
                P = rand(1)[1]
                if P < exp(-Delta)
                    test[1][i][J] = proposal
                    nb_updates += 1
                end
            end
        end
        for j = 1:N_x
            for i = 1:Int64(0.5*N_t)
                I = 2*i -eo
                a = test[2][I][j]
                b = staple_dagger(2,I,j,test)
                proposal = deepcopy(a)
                q = 2*(rand(1)[1] -0.5)
                proposal += r*q

                Delta = -β * real((exp(im*proposal)- exp(im*a)) * b)
                P = rand(1)[1]
                if P < exp(-Delta)
                    test[2][I][j] = proposal
                    nb_updates += 1
                end
            end
        end
    end
    result = [test, nb_updates]
    result
end

# Funktion, um die PHASE einer Plaquette Pμν( ⃗n) zu berechnen
function P12(nt, nx, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NX = nx%N_x +1  # nₛ+1 mit periodische Randbed.
    config[1][nt][nx] + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
end

# Funktion um Wirkung einer gegebenen Konfig zu messen
function action(config)
    S = 0
        for i = 1:N_t
            for j = 1:N_x
                plaq = exp(im*P12(i,j,config))
                S += β*real(1-plaq)
            end
        end
    S
end

# Funktion um die topologische Ladung einer Konfig zu messen
function int_charge(config)
    Q = 0
    for i = 1:N_t
        for j = 1:N_x
            plaq = exp(im*P12(i,j,config))
            Q += (imag(log(plaq))) / (2*pi)
        end
    end
    Q
end

# Funktion um die "kontinuierliche" topologische Ladung einer Konfig zu messen
function cont_charge(config)
    Q = 0
    for i = 1:N_t
        for j = 1:N_x
            plaq = exp(im*P12(i,j,config))
            Q += (imag(plaq)) / (2*pi)
        end
    end
    Q
end

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


# Funktion für ein Update der topologischen Ladung
function top_update(config)
    if rand() < 0.5
        insta_1 = instanton(1)
    else
        insta_1 = -instanton(1)
    end
    new_config = config.+insta_1

    Delta = action(new_config) - action(config)
    if rand() < exp(-Delta)
        return new_config
    else
        return config
    end
end


#                  ##### Die eigentliche Simulation #####
# Zunächst werden die "rohen" Konfigs durch Metropolis- und top. Updates erstellt
raw_configs = [config_1]
for i = 1:nb_sim
    for j = 1:nb_metro
        J = (nb_metro+nb_top)*(i-1) + j
        a = metro_update(raw_configs[J])
        link_accepted += a[2]
        push!(raw_configs, a[1])
    end
    for j = 1:nb_top
        J = (nb_metro+nb_top)*(i-1) + nb_metro + j
        push!(raw_configs, top_update(raw_configs[J]))
    end
end

acceptance = round(link_accepted/(nb_sim*nb_metro*N_t*N_x*2), digits = 3)

# Gemessen wird nur jede N_skip-te Konfig, nach der Thermalization
configs = []
for i = cut:size(raw_configs, 1)
    if i%N_skip == 1
        push!(configs, raw_configs[i])
    end
end
nb_configs = size(configs, 1)


# Plots
actions = []
int_charges = []
cont_charges = []
x = []
for i = 1:nb_configs
    push!(actions, action(configs[i]))
    push!(int_charges, int_charge(configs[i]))
    push!(cont_charges, cont_charge(configs[i]))
    push!(x, i)
end


# Die ganzzahligen werden auf ganze Zahlen gerundet (begrenzte
# Maschinenpräzision führt zu Abweichungen der Ordnung e-10 und kleiner)
int_charges = round.(Int, int_charges)
# Ganzzahlige Ladungen werden für das Histogram so verschoben, dass die
# Säulen auf den Zahlen liegen (s.u.)
#int_charges = int_charges.-0.25


# Das 2D-Histogram zum Vergleich der diskreten und kont. Ladungen
histo2d = histogram2d(
int_charges, cont_charges, bins = 100, c = :reds, background_color = :grey)

histo2d = plot!(
size = (750,600),
title = "Different top. charges: $nb_metro metrop., $nb_top top. update
$N_t×$N_x,   β = $β,   ($Nb_updates - $cut) sweeps,   N_skip = $N_skip",
xticks = -6:1:6,
yticks = -6:1:6,
xlabel = "Discrete Q",
ylabel = "Continuous Q",
labelfontsize = 16,
titlefontsize = 16
)

display(histo2d)


# Folgender Befehl erstellt eine txt-Datei in dem Ordner, in dem sich
# dieses Programm befindet, mit dem Namen compare_topup_N_t.x.N_x.txt
# (wobei N_t und N_x im Namen durch die Zahlenwerte der aktuellen Simulation
# ersetzt werden), die die Daten der kontinuierlichen und diskreten top.
# Ladungen in zwei Spalten beinhalten

writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\Julia-stuff\\compare_topup_$N_t.x.$N_x.txt", [int_charges cont_charges])


# Folgender Befehl speichert die aktuelle Grafik unter dem Dateipfad
# C:\Users\proue\OneDrive\Desktop\Physik Uni\Thesis unter dem Namen
# compare_topup_24x24 ab

savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\topup_sim_$N_t.x.$N_x.")


# histo1 = histogram(actions, bins = 35,
# label = "Wilson Gauge Actions",
# title = "Nₓ×Nₜ = $N_x × $N_t,   β = $β,   nb_sim⋅(nb_metro+nb_top) = $nb_sim⋅($nb_metro+$nb_top),
# cut = $cut,   N_skip = $N_skip,   r/2π = $R,   acceptance = $acceptance " )
# histo2 = histogram(int_charges, bins = 130,
# label = "Topological int. Charges")
# image_histo = plot(histo1, histo2, layout=(2,1), size=(800,600))
# display(image_histo)
#
# plot1 = plot(x, int_charges)

#savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\
#Lattice Gauge Theory Varnhorst\\topup_32x32")
