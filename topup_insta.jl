# Es werden 5 Plots erstellt: die kont. und diskreten top. Ladungen von
# n-Instanton-Konfigurationen werden in Abhängigkeit von n für quadratische
# Gitter der Kantenlänge N = 12, 24, ..., 60 aufgezeichnet. Da ein
# sinusförmiger Verlauf der kont. Ladung entsteht, dessen Frequenz
# bekannt ist, wird außerdem eine Datei erstellt, die die Amplituden
# beinhaltet, sodass dieser Verlauf eindeutig bestimmt ist.
#
# Die in der Thesis verwendeten Abbildungen wurden mit einer leicht abgewandelten
# Version dieses Programms erstellt, welche verlangt, dass jeder Plot manuell
# angepasst werden muss. Um dies zu vermeiden, wurde diese Version vorgelegt
# (hier müssen nur die entsprechenden Dateipfade eingesetzt werden, s.u.)

using Plots
using DelimitedFiles

text1 = []
amplitudes = []
text2 = []
anomals = []

# Da direkt alle 5 Plots erstellt werden, ist der Hauptteil des Programms in
# einer großen for-Schleife
for i = 1:5

N = 12*i
N_t = N
N_x = N

# Funktion, um die PHASE einer Plaquette Pμν( ⃗n) zu berechnen
function P12(nt, nx, config)
    NT = nt%N_t +1  # nₜ+1 mit periodische Randbed.
    NX = nx%N_x +1  # nₛ+1 mit periodische Randbed.
    config[1][nt][nx] + config[2][NT][nx] - config[1][nt][NX] - config[2][nt][nx]
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


##### Generieren der Plots und des Textfiles mit Amplituden


insta_int_charges = []
insta_cont_charges = []
c = 0.5*N^2+1
for i = -c:c
    push!(insta_int_charges,int_charge(instanton(i)))
    push!(insta_cont_charges,cont_charge(instanton(i)))
end

IMAGE = plot(-c:c,insta_int_charges, label="discr. top. charges")
IMAGE = plot!(-c:c, insta_cont_charges, label="cont. top. charges")
IMAGE = plot!(
title = "Discrete VS Continuous Top. Charges of n-Instanton-Configs:
$N_t×$N_x-Lattice",
size=(750,600),
xticks = -c+1:4*N:c-1,
yticks = -c+1:4*N:c-1,
xlabel = "n of the n-Instanton-Configuration",
legend = :bottomright,
labelfontsize = 16,
legendfontsize = 12,
titlefontsize = 16
)

display(IMAGE)

# Der folgende Befehl speichert die erstellten Grafiken in dem Dateipfad
# C:\Users\proue\OneDrive\Desktop\Physik Uni\Thesis\Julia-stuff
# unter dem Namen compare_insta_N_t.x.N_x. wobei für N_t und N_x die jeweiligen
# Werte eingesetzt werden:

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\compare_insta_$N_t.x.$N_x.")

push!(text1, "Amplitude für N=$N:")
push!(amplitudes, cont_charge(instanton(0.25*N^2)))
push!(text2, "diskr. Lad. bei n = 0.5*N^2:")
aus = abs(round(Int, insta_int_charges[2]))
push!(anomals, string("± ", aus))


end

# Der folgende Befehl erstellt die Datei compare_insta.txt in dem Dateipfad
# C:\Users\proue\OneDrive\Desktop\Physik Uni\Thesis\Julia-stuff

# writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Thesis\\Julia-stuff\\compare_insta.txt", [text1 amplitudes text2 anomals])
