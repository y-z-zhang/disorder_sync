# In[]
using DelimitedFiles
using Statistics
cd("/home/yuanzhao/random_hetero/experiments/")
pwd()
current_homo = readdlm("oc091319_6.txt");
current_homo = current_homo[:,1:16]
shift = Array{Float64}(undef,16)
for i in 1:16
    shift[i] = (maximum(current_homo[:,i])+minimum(current_homo[:,i]))/2
end
current_homo = current_homo' .- shift

# In[]
using Hilbert
z = hilbert(current_homo)
amp = abs.(z)
phase = angle.(z);

# In[]
using Plots; pyplot()
t = range(0,stop=100,length=20001)
plot(t,amp[1,1:20001],xlims=(50,100))#,label="Amplitude of analytic signal")
plot!(t,real(z),label="Shifted experimental data")
plot!(t,imag(z),xlims=(68,90),ylims=(-.2,.3),
    xlabel="Time(s)",ylabel="Current (mA)",
    xtickfont=font(12),ytickfont = font(12),guidefont=font(15),legendfont=font(10),
    label="Hilbert transform of shifted experimental data",legend=:topright)

# In[]
plot(t,angle.(z),xlims=(70,90),label="Phase of analytic signal")

# In[]
using LinearAlgebra
using Plots
using Statistics
#?eigvals
M=diagm(0=>[.01,-.1,-.1])
ϵ=1e-3
e=zeros(100)
for i in 1:100
    ΔM=ϵ*randn(3,3)
    ΔM=ΔM
    e[i]=real(eigvals(M+ΔM)[3])
end
mean(e)
