using Plots
using Statistics
using OffsetArrays
#=
Die Phasen von Eichfeld-Konfigurationen Uμ( ⃗n) werden in einem Array
config[μ][nₜ][nₛ] gespeichert, wobei ⃗n = (nₜ, nₛ). Im weiteren soll
μ = 1 für die zeitliche und μ = 2 für die räumliche Dimension stehen.
(es wird auch μ ≡ 1 für temp. und ν ≡ 2 für räuml. synonym verwendet)
Phasen der Eichtransformationen Ω( ⃗n) werden als trans[nₜ, nₛ] gespeichert.
=#
N_t = 24       # Anzahl temporaler Gitterpunkte
N_x = 24    	# Anzahl räumlicher Gitterpunkte
β = (N_x*N_t)/80      # 1/a²g², wobei g die Eich-Kopplung ist
sweeps = 200000          # Updates im Metropolis
r = 2*pi*0.03      # Schrittweite im Metropolis (nur in pos. Richtung, s.u.)
cut = 1200 #Int64(0.6*sweeps)  # Thermalization
N_skip = 20         # Nur jede N_skip-te Konfig wird gemessen
δq = 0.005
w = 0.005
Q_thrs = 13.0


# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\meta_8")



# Die erste Feldkonfig. wird per Zufall generiert (hot start), mit Werten ∈ [0,2π]
config_1 = []
for mu = 1:2
    push!(config_1, [])
    for i = 1:N_t
        a = (2*pi).*rand(N_x)
        push!(config_1[mu], a)
    end
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

# Funktion um den Unterschied in der kont. Ladung einer gegebenen Konfig
# zu bestimmen, nachdem U[mu][nt][nx] geändert wurde (U_old -> U_new)
function DeltaQ(U_old, U_new, mu, nt, nx, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NX = nx%N_x +1  # nₛ+1 mit periodische Randbed.
    TN = (nt + N_t -2)%N_t +1   # nₜ-1 mit periodische Randbed.
    XN = (nx + N_x -2)%N_x +1   # nₛ-1 mit periodische Randbed.

    if mu == 1
        a = U_old + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
        b = config[1][nt][XN] + config[2][NT][XN] - U_old - config[2][nt][XN]
        c = U_new + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
        d = config[1][nt][XN] + config[2][NT][XN] - U_new - config[2][nt][XN]
    elseif mu == 2
        a = U_old - config[1][TN][NX] - config[2][TN][nx] + config[1][TN][nx]
        b = config[2][NT][nx] - config[1][nt][NX] - U_old + config[1][nt][nx]
        c = U_new - config[1][TN][NX] - config[2][TN][nx] + config[1][TN][nx]
        d = config[2][NT][nx] - config[1][nt][NX] - U_new + config[1][nt][nx]
    end

    (imag(exp(im*a) + exp(im*b) - exp(im*c) - exp(im*d))) / (2*pi)
end


# Funktion, um die PHASE einer Plaquette Pμν( ⃗n) zu berechnen
function P12(nt, nx, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NX = nx%N_x +1  # nₛ+1 mit periodische Randbed.
    config[1][nt][nx] + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
end

# Funktion um die ganzzahlige topologische Ladung einer Konfig zu messen
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

# Das Bias-Potential wird auf einem Grid gespeichert
gridq = []
len_gridq = Int64(2*Q_thrs/δq)

for k = 0:len_gridq+2
    a = -Q_thrs + (k-1) *δq
    push!(gridq, a)
end

gridq = OffsetArray(gridq, 0:len_gridq+2)


Vbias = zeros(len_gridq+3)
Vbias = OffsetArray(Vbias, 0:len_gridq+2)


# Die + 0.5000000001 anstelle von +0.5 ist wichtig, da in Julia
# round(Int, 4.5) = 4 ist, aber round(Int, 4.5000000001) = 5
function gridint(q)
    a = (q-gridq[1])/δq + 0.5000000001
    round(Int, a)
end


# Funktion um das Bias-Potential zu aktualisieren
function update_bias(q)
    i = gridint(q)
    Vbias[i] += w*(1 - (q-gridq[i])/δq )

    if i != len_gridq +2
        Vbias[i+1] += w*(q-gridq[i])/δq
    end
end




#           ***** Hauptteil der Simulation *****



# Thermalisierung:
# Es werden "cut-viele" Konfigs mithilfe des Metropolis-Algorithmus konstruiert,
# um das System zu thermalisieren. Diese raw_raw_configs finden keine
# weitere Verwendung

link_accepted = 0
raw_raw_configs = [config_1]
for k = 1:cut
    test = deepcopy(raw_raw_configs[k])

    for eo = 0:1  # even-odd
        for i = 1:N_t
            for j = 1:Int64(0.5*N_x)
                J = 2*j -eo
                a = test[1][i][J]
                proposal = deepcopy(a)
                step = 2*(rand() -0.5)
                proposal += r*step

                b = staple_dagger(1,i,J,test)
                Δ = -β * real((exp(im*proposal)- exp(im*a)) * b)
                if rand() < exp(-Δ)
                    test[1][i][J] = proposal
                    link_accepted += 1
                end
            end
        end
        for j = 1:N_x
            for i = 1:Int64(0.5*N_t)
                I = 2*i -eo
                a = test[2][I][j]
                proposal = deepcopy(a)
                step = 2*(rand() -0.5)
                proposal += r*step

                b = staple_dagger(2,I,j,test)
                Δ = -β * real((exp(im*proposal)- exp(im*a)) * b)
                if rand() < exp(-Δ)
                    test[2][I][j] = proposal
                    link_accepted += 1
                end
            end
        end
    end

    push!(raw_raw_configs, test)
end
acceptance = round(100*link_accepted/(2*cut*N_t*N_x), digits = 3)
println("", "Acceptance during thermalization: $acceptance%")

# Metadynamics:
# Die eigentliche Simulation, in der Metadynamics verwendet werden. Dies
# geschieht ähnlich wie in der Thermalisierung, welche hier wieder verwendet
# und um neue Zeilen ergänzt wurde, die mit einem ❌ markiert werden, wobei
# modifizierte Zeilen mit einem ⭕ gekennzeichnet werden

raw_configs = [raw_raw_configs[cut]]
first_charge = cont_charge(raw_configs[1])
update_bias(first_charge)
for k = 1:sweeps
    test = deepcopy(raw_configs[k])
    Q = cont_charge(test)   #❌ new line

    for eo = 0:1  # even-odd
        for i = 1:N_t
            for j = 1:Int64(0.5*N_x)
                J = 2*j -eo
                a = test[1][i][J]
                proposal = deepcopy(a)
                step = 2*(rand() -0.5)
                proposal += r*step

                Q_prop = Q + DeltaQ(a, proposal, 1, i, J, test) #❌ new line
                ind = gridint(Q)    #❌ new line
                ind_prop = gridint(Q_prop)  #❌ new line
                b = staple_dagger(1,i,J,test)
                Δ = -β * real((exp(im*proposal)- exp(im*a)) * b) + Vbias[ind_prop] -Vbias[ind] #⭕ modified line

                if rand() < exp(-Δ)
                    test[1][i][J] = proposal
                    Q = Q_prop  #❌ new line
                    link_accepted += 1
                end
                update_bias(Q) #❌ new line
            end
        end
        for j = 1:N_x
            for i = 1:Int64(0.5*N_t)
                I = 2*i -eo
                a = test[2][I][j]
                proposal = deepcopy(a)
                step = 2*(rand() -0.5)
                proposal += r*step

                Q_prop = Q + DeltaQ(a, proposal, 2, I, j, test) #❌ new line
                ind = gridint(Q)    #❌ new line
                ind_prop = gridint(Q_prop)  #❌ new line
                b = staple_dagger(2,I,j,test)
                Δ = -β * real((exp(im*proposal)- exp(im*a)) * b) + Vbias[ind_prop] -Vbias[ind] #⭕ modified line

                if rand() < exp(-Δ)
                    test[2][I][j] = proposal
                    Q = Q_prop  #❌ new line
                    link_accepted += 1
                end
                update_bias(Q)  #❌ new line
            end
        end
    end

    push!(raw_configs, test)
end

acceptance = round(100*link_accepted/(2*(cut+sweeps)*N_t*N_x), digits = 3)
println("", "Acceptance total: $acceptance%")
println("\n ")



nb_raw = length(raw_configs)
cont_charges = []
for k = 1:nb_raw
    push!(cont_charges, cont_charge(raw_configs[k]))
end


# image1 = plot(1:nb_raw, cont_charges, legend = :false)


xvals = []
yvals = []

for k = 0:len_gridq+2
    push!(xvals, gridq[k])
    push!(yvals, Vbias[k])
end

image2 = plot(xvals, yvals)
image2 = plot!(
title = "Metadynamics:   $N_t×$N_x lattice,   $sweeps sweeps, \n w = $w,   δq = $δq",
size=(750,600),
xticks = -Q_thrs:Q_thrs,
xlabel = "q grid",
ylabel = "bias potential",
legend = :false,
labelfontsize = 16,
legendfontsize = 12,
titlefontsize = 16
)
display(image2)
