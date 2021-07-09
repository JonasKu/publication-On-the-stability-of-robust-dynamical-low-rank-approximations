include("settings.jl")
include("DLRASolver.jl")

using DelimitedFiles
using PyPlot

#close("all")
Nx = 2002;
s = Settings(Nx,"UQ");

recordFNorm = true;
############################
# rank 5
s.r = 5;
dx = s.dx;
s.cfl = 1.0;
s.dt = s.cfl*dx;
s.stabilization = 0;
solver = DLRSolver(s);
@time tEnd, uUIr5,timeUIr5,FNormUIr5 = SolveUnconventional(solver,recordFNorm);

s.stabilization = 2; # use stable treatment for projector splitting
solver = DLRSolver(s);
@time tEnd, uPSstabler5,timePSstabler5,FNormPSstabler5 = SolveProjectorSplitting(solver,recordFNorm);

s.stabilization = 0; # discretize first, DLR second
solver = DLRSolver(s);
@time tEnd, uPSr5,timePSr5,FNormPSr5 = SolveProjectorSplitting(solver,recordFNorm);

############################
# rank 10
s.r = 10;
dx = s.dx;
s.cfl = 1.0;
s.dt = s.cfl*dx;
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
x = s.xMid;
Nx = length(x);
uEx = zeros(length(x));
VarEx = zeros(length(x));
N = s.nPN;
q = Quadrature(10*N,"Gauss");
b = Basis(q,s);
aXi = ((1.0 - s.sigmaS) .+s.sigmaS*q.xi).^3;
for j = 1:s.NCells
    uEx[j] = IntegralVec(q, IC(s,x[j].-tEnd*aXi).*b.PhiQuad[:,1]*0.5,-1.0,1.0);
    VarEx[j] = IntegralVec(q, (IC(s,x[j].-tEnd*aXi).-uEx[j]).^2 .*b.PhiQuad[:,1]*0.5,-1.0,1.0);
end
VarPS = zeros(Nx);
VarUI = zeros(Nx);
VarPSstable = zeros(Nx);
VarPSr5 = zeros(Nx);
VarUIr5 = zeros(Nx);
VarPSstabler5 = zeros(Nx);

for j = 1:Nx
    VarPS[j] = uPS[j,2:end]'uPS[j,2:end];
    VarUI[j] = uUI[j,2:end]'uUI[j,2:end];
    VarPSstable[j] = uPSstable[j,2:end]'uPSstable[j,2:end];
    VarPSr5[j] = uPSr5[j,2:end]'uPSr5[j,2:end];
    VarUIr5[j] = uUIr5[j,2:end]'uUIr5[j,2:end];
    VarPSstabler5[j] = uPSstabler5[j,2:end]'uPSstabler5[j,2:end];
end


#####################################
##           plot r = 5           ##
#####################################

# plot stable versions with CFL=1
fig = figure("Figure4a",figsize=(15, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](x,uEx, "k-", linewidth=2, label="exact", alpha=0.8)
ax[:plot](s.xMid,uPSr5[:,1], "g--", linewidth=2, label="projector-splitting (discretize first)", alpha=0.6)
ax[:plot](s.xMid,uPSstabler5[:,1], "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](s.xMid,uUIr5[:,1], "b-.", linewidth=2, label="unconventional", alpha=0.8)
#ax[:plot](s.xMid,u[:,1], "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="upper right", fontsize=20)
ylabel(L"E[$u$]", fontsize=25)
xlabel(L"$x$", fontsize=25)
ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=25) 
tight_layout()
show()

fig = figure("Figure5a",figsize=(15, 10), dpi=100)
ax = gca()
ax[:plot](x,sqrt.(VarEx), "k-", linewidth=2, label="exact", alpha=0.8)
ax[:plot](s.xMid,sqrt.(VarPSr5), "g--", linewidth=2, label="projector-splitting (discretize first)", alpha=0.8)
ax[:plot](s.xMid,sqrt.(VarPSstabler5), "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](s.xMid,sqrt.(VarUIr5), "b-.", linewidth=2, label="unconventional", alpha=0.8)
#x[:plot](s.xMid,sqrt.(Var), "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="upper right", fontsize=20)
ylabel(L"\sigma_u", fontsize=25)
xlabel(L"$x$", fontsize=25)
#ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=25) 
tight_layout()
show()

## investigate time evolution of Frobenius norm

fig = figure("Figure6a",figsize=(10, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](timePSr5,FNormPSr5, "g--", linewidth=2, label="projector-splitting (discretize first)", alpha=0.8)
ax[:plot](timePSstabler5,FNormPSstabler5, "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](timeUIr5,FNormUIr5, "b-.", linewidth=2, label="unconventional", alpha=0.8)
#ax[:plot](s.xMid,u[:,1], "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="lower left", fontsize=20)
ylabel(L"$\Vert u\Vert_F$", fontsize=20)
xlabel(L"$t$", fontsize=20)
ax.set_xlim([0.0,s.tEnd])
ax.set_ylim([min(minimum(FNormUIr5),minimum(FNormPSstabler5)),max(maximum(FNormUIr5),maximum(FNormPSstabler5),minimum(FNormPSr5)+0.5)])
ax.tick_params("both",labelsize=20) 
tight_layout()
show()

#####################################
##           plot r = 10           ##
#####################################

# plot stable versions with CFL=1
fig = figure("Figure4b",figsize=(15, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](x,uEx, "k-", linewidth=2, label="exact", alpha=0.8)
ax[:plot](s.xMid,uPSstable[:,1], "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](s.xMid,uUI[:,1], "b-.", linewidth=2, label="unconventional", alpha=0.8)
ax[:legend](loc="upper right", fontsize=20)
ylabel(L"E[$u$]", fontsize=25)
xlabel(L"$x$", fontsize=25)
ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=25) 
tight_layout()
show()

fig = figure("Figure5b",figsize=(15, 10), dpi=100)
ax = gca()
ax[:plot](x,sqrt.(VarEx), "k-", linewidth=2, label="exact", alpha=0.8)
ax[:plot](s.xMid,sqrt.(VarPSstable), "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](s.xMid,sqrt.(VarUI), "b-.", linewidth=2, label="unconventional", alpha=0.8)
ax[:legend](loc="upper right", fontsize=20)
ylabel(L"\sigma_u", fontsize=25)
xlabel(L"$x$", fontsize=25)
#ax.set_ylim([0.0,1.2*maximum(uUI[:,1])])
ax.set_xlim([s.a,s.b])
ax.tick_params("both",labelsize=25) 
tight_layout()
show()
## investigate time evolution of Frobenius norm

fig = figure("Figure6b",figsize=(10, 10), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
ax = gca()
ax[:plot](timePS,FNormPS, "g--", linewidth=2, label="projector-splitting (discretize first)", alpha=0.8)
ax[:plot](timePSstable,FNormPSstable, "r:", linewidth=2, label="projector-splitting (DLRA first)", alpha=0.8)
ax[:plot](timeUI,FNormUI, "b-.", linewidth=2, label="unconventional", alpha=0.8)
#ax[:plot](s.xMid,u[:,1], "g:", linewidth=2, label=L"$P_N$", alpha=0.6)
ax[:legend](loc="upper right", fontsize=20)
ylabel(L"$\Vert u\Vert_F$", fontsize=20)
xlabel(L"$t$", fontsize=20)
ax.set_xlim([0.0,s.tEnd])
#ax.set_ylim([min(minimum(FNormUI),minimum(FNormPSstable)),max(maximum(FNormUI),maximum(FNormPSstable),minimum(FNormPS)+5.0)])
ax.set_ylim([min(minimum(FNormUIr5),minimum(FNormPSstabler5)),max(maximum(FNormUIr5),maximum(FNormPSstabler5),minimum(FNormPSr5)+0.5)])
ax.tick_params("both",labelsize=20) 
tight_layout()
show()
