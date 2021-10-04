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
N_t = 24
N_x = 24
β = (N_t*N_x)/80
r = 2*pi*0.05
sweeps = 80000
sweeps_by_ten = Int64(0.1*sweeps)
sweeps_count = 0
cut = 10000
nb_accepted = 0
nb_metro = 1
nb_snapshots = Int64(sweeps/5)

δq = 0.01
w = 0.0001
Q_max = 15
Q_thrs = 5
boundary_count = 0

bb = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\Julia-stuff\\new_meta_bb.txt")
bias_potential = bb[size(bb,1),:] .- bb[size(bb,1),1]



function P12(nt, nx, config)
    NT = nt%N_t +1
    NX = nx%N_x +1
    return config[1][nt][nx] + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
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

function cont_charge(config)
    Q = 0.0
    for i = 1:N_t
        for j = 1:N_x
            Q += (1/(2*pi)) * sin(P12(i,j,config))
        end
    end
    return Q
end

function grid_ind(q)
    grid_index = (q+Q_max)/δq + 0.5000000001
    return round(Int, grid_index)
end



rand_charges = []
rand_weights = []
for k = 1:sweeps
    config = []
    for μ = 1:2
        push!(config, [])
        for i = 1:N_t
            rand_phase = 2*pi*(rand(N_x) .- 0.5)
            push!(config[μ], rand_phase)
        end
    end
    if abs(cont_charge(config)) < 10
        push!(rand_charges, int_charge(config))
        push!(rand_weights, exp(-bias_potential[grid_ind(Q)]))
    end
end

χ_mean = sum((rand_charges.^2) .* rand_weights)/sum(rand_weights)
