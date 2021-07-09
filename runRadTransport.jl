include("settings.jl")
include("DLRASolver.jl")

using PyPlot
using DelimitedFiles

#close("all")

s = Settings(802,"LineSource");
recordFNorm = true;
############################
# rank 10
cflHigh = 1.0
s.r = 10;
dx = s.dx;
s.dt = cflHigh*dx;
s.stabilization = 0;
solver = DLRSolver(s);
@time tEnd, uUIr10,timeUIr10,FNormUIr10 = SolveUnconventional(solver,recordFNorm);

s.stabilization = 2; # use stable treatment for projector splitting
solver = DLRSolver(s);
@time tEnd, uPSstabler10,timePSstabler10,FNormPSstabler10 = SolveProjectorSplitting(solver,recordFNorm);

s.stabilization = 0; # discretize first, DLR second
solver = DLRSolver(s);
@time tEnd, uPSr10,timePSr10,FNormPSr10 = SolveProjectorSplitting(solver,recordFNorm);

############################
# rank 15
cflHigh = 1.0
s.r = 15;
dx = s.dx;
s.dt = cflHigh*dx;
s.stabilization = 0;
solver = DLRSolver(s);
@time tEnd, uUI,timeUI,FNormUI = SolveUnconventional(solver,recordFNorm);

s.stabilization = 2; # use stable treatment for projector splitting
solver = DLRSolver(s);
@time tEnd, uPSstable,timePSstable,FNormPSstable = SolveProjectorSplitting(solver,recordFNorm);

s.stabilization = 0; # discretize first, DLR second
solver = DLRSolver(s);
@time tEnd, uPS,timePS,FNormPS = SolveProjectorSplitting(solver,recordFNorm);

s.tEnd = tEnd;

# store reference solution
v = readdlm("PlaneSourceRaw", ',')
uEx = zeros(length(v));
for i = 1:length(v)
    if v[i] == ""
        uEx[i] = 0.0;
    else
        uEx[i] = Float64(v[i])
    end
end
x = collect(range(-1.5,1.5,length=(2*length(v)-1)));
uEx = [uEx[end:-1:2];uEx]

#####################################
##           plot r = 10           ##
#####################################

# plot stable versions with CFL=1
fig = figure("Figure1a",figsize=(15, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](x,uEx, "k-", linewidth=2, label="exact", alpha=0.8)
ax[:plot](s.xMid,uPSr10[:,1], "g--", linewidth=2, label="projector-splitting (discretize first)", alpha=0.6)
ax[:plot](s.xMid,uPSstabler10[:,1], "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](s.xMid,uUIr10[:,1], "b-.", linewidth=2, label="unconventional", alpha=0.8)
#ax[:plot](s.xMid,u[:,1], "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="upper left", fontsize=20)
ylabel(L"$\phi$", fontsize=25)
xlabel(L"$x$", fontsize=25)
ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=25) 
tight_layout()
show()

## investigate time evolution of Frobenius norm

fig = figure("Figure2a",figsize=(10, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](timePSr10,FNormPSr10, "g--", linewidth=2, label="projector-splitting (discretize first)", alpha=0.8)
ax[:plot](timePSstabler10,FNormPSstabler10, "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](timeUIr10,FNormUIr10, "b-.", linewidth=2, label="unconventional", alpha=0.8)
#ax[:plot](s.xMid,u[:,1], "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="upper left", fontsize=20)
ylabel(L"$\Vert\psi\Vert_F$", fontsize=20)
xlabel(L"$t$", fontsize=20)
ax.set_xlim([0.0,s.tEnd])
ax.set_ylim([min(minimum(FNormUIr10),minimum(FNormPSstabler10)),max(maximum(FNormUIr10),maximum(FNormPSstabler10),maximum(FNormPSr10))])
ax.tick_params("both",labelsize=20) 
tight_layout()
show()

#####################################
##           plot r = 15           ##
#####################################

# plot stable versions with CFL=1
fig = figure("Figure1b",figsize=(15, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](x,uEx, "k-", linewidth=2, label="exact", alpha=0.8)
ax[:plot](s.xMid,uPSstable[:,1], "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](s.xMid,uUI[:,1], "b-.", linewidth=2, label="unconventional", alpha=0.8)
#ax[:plot](s.xMid,uPS[:,1], "g--", linewidth=2, label="projector-splitting", alpha=0.6)
#ax[:plot](s.xMid,u[:,1], "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="upper left", fontsize=20)
ylabel(L"$\phi$", fontsize=25)
xlabel(L"$x$", fontsize=25)
ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=25) 
tight_layout()
show()

## investigate time evolution of Frobenius norm

fig = figure("Figure2b",figsize=(10, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](timePS,FNormPS, "g--", linewidth=2, label="projector-splitting (discretize first)", alpha=0.8)
ax[:plot](timePSstable,FNormPSstable, "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](timeUI,FNormUI, "b-.", linewidth=2, label="unconventional", alpha=0.8)
#ax[:plot](s.xMid,u[:,1], "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="upper right", fontsize=20)
ylabel(L"$\Vert\psi\Vert_F$", fontsize=20)
xlabel(L"$t$", fontsize=20)
ax.set_xlim([0.0,s.tEnd])
#ax.set_ylim([min(minimum(FNormUI),minimum(FNormPSstable)),max(maximum(FNormUI),maximum(FNormPSstable),minimum(FNormPS)+5.0)])
ax.set_ylim([min(minimum(FNormUIr10),minimum(FNormPSstabler10)),max(maximum(FNormUIr10),maximum(FNormPSstabler10),maximum(FNormPSr10))])
ax.tick_params("both",labelsize=20) 
tight_layout()
show()
