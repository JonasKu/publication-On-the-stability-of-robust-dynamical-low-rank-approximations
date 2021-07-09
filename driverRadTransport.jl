include("settings.jl")
include("DLRASolver.jl")

using NPZ;
using PyPlot;
using DelimitedFiles;

#close("all");

local u, solver, tEnd;

Nx = 801;
s = Settings(Nx+1);
N = s.nPN;
dx = s.dx;

if Nx != s.NCells
    println("ERROR: NCells missmatch")
end

# store exact solution
if s.ICType == "LS"
    v = readdlm("PlaneSourceRawNx500", ',')
    uEx = zeros(length(v));
    for i = 1:length(v)
        if v[i] == ""
            uEx[i] = 0.0;
        else
            println(v[i])
            uEx[i] = Float64(v[i])
        end
    end
    x = collect(range(-1.5,1.5,length=(2*length(v)-1)));
    uEx = [uEx[end:-1:2];uEx]
elseif s.ICType == "shock"
    x = s.xMid;
    Nx = length(x);
    uEx = zeros(length(x));
    VarEx = zeros(length(x));
    N = s.nPN;
    q = Quadrature(2*N,"Gauss");
    b = Basis(q,s);
    aXi = (1.0 - s.sigmaS) .+s.sigmaS*q.xi;
    for j = 1:s.NCells
        uEx[j] = IntegralVec(q, IC(s,x[j].-s.tEnd*aXi).*b.PhiQuad[:,1]*0.5,-1.0,1.0);
        VarEx[j] = IntegralVec(q, (IC(s,x[j].-s.tEnd*aXi).-uEx[j]).^2 .*b.PhiQuad[:,1]*0.5,-1.0,1.0);
    end
end

#r = [2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13; 14; 15; 16]
r = [10; 15];
cfl = [0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8; 0.9; 1.0;]
errorPS = zeros(length(r),length(cfl));
errorPS1 = zeros(length(r),length(cfl));
errorPS2 = zeros(length(r),length(cfl));
errorUI = zeros(length(r),length(cfl));
errorPSVar = zeros(length(r),length(cfl));
errorPS1Var = zeros(length(r),length(cfl));
errorPS2Var = zeros(length(r),length(cfl));
errorUIVar = zeros(length(r),length(cfl));
Var = zeros(Nx);

s.stabilization = 0;

for k = 1:length(r)
    for i = 1:length(cfl)
        s.r = r[k];
        s.cfl = cfl[i];
        s.dt = cfl[i]*dx;
        solver = DLRSolver(s);
        @time tEnd, u = SolveProjectorSplitting(solver);
        errorPS[k,i] = s.dx*sqrt(sum((u[:,1].-uEx).^2));
        println("rank ",r[k], ", cfl ",cfl[i]," : ",errorPS[k,i])
        if s.problem == "UQ"
            for j = 1:Nx
                Var[j] = u[j,2:end]'u[j,2:end];
            end
            errorPSVar[k,i] = s.dx*sqrt(sum((sqrt.(Var).-sqrt.(VarEx)).^2));
            println("Var: rank ",r[k], ", cfl ",cfl[i]," : ",errorPSVar[k,i])
        end
    end
end
println("-> DLR Projector Splitting Standard DONE.")

s.stabilization = 1;

for k = 1:length(r)
    for i = 1:length(cfl)
        s.r = r[k];
        s.cfl = cfl[i];
        s.dt = cfl[i]*dx;
        solver = DLRSolver(s);
        @time tEnd, u = SolveProjectorSplitting(solver);
        errorPS1[k,i] = s.dx*sqrt(sum((u[:,1].-uEx).^2));
        println("rank ",r[k], ", cfl ",cfl[i]," : ",errorPS1[k,i])
        if s.problem == "UQ"
            for j = 1:Nx
                Var[j] = u[j,2:end]'u[j,2:end];
            end
            errorPS1Var[k,i] = s.dx*sqrt(sum((sqrt.(Var).-sqrt.(VarEx)).^2));
            println("Var: rank ",r[k], ", cfl ",cfl[i]," : ",errorPS1Var[k,i])
        end
    end
end
println("-> DLR Projector Splitting 1 DONE.")

s.stabilization = 2;

for k = 1:length(r)
    for i = 1:length(cfl)
        s.r = r[k];
        s.cfl = cfl[i];
        s.dt = cfl[i]*dx;
        solver = DLRSolver(s);
        @time tEnd, u = SolveProjectorSplitting(solver);
        errorPS2[k,i] = s.dx*sqrt(sum((u[:,1].-uEx).^2));
        println("rank ",r[k], ", cfl ",cfl[i]," : ",errorPS2[k,i])
        if s.problem == "UQ"
            for j = 1:Nx
                Var[j] = u[j,2:end]'u[j,2:end];
            end
            errorPS2Var[k,i] = s.dx*sqrt(sum((sqrt.(Var).-sqrt.(VarEx)).^2));
            println("Var: rank ",r[k], ", cfl ",cfl[i]," : ",errorPS2Var[k,i])
        end
    end
end
println("-> DLR Projector Splitting 2 DONE.")

for k = 1:length(r)
    for i = 1:length(cfl)
        s.r = r[k];
        s.cfl = cfl[i];
        s.dt = cfl[i]*dx;
        solver = DLRSolver(s);
        @time tEnd, u = SolveUnconventional(solver);
        errorUI[k,i] = s.dx*sqrt(sum((u[:,1].-uEx).^2));
        println("rank ",r[k], ", cfl ",cfl[i]," : ",errorUI[k,i])
        if s.problem == "UQ"
            for j = 1:Nx
                Var[j] = u[j,2:end]'u[j,2:end];
            end
            errorUIVar[k,i] = s.dx*sqrt(sum((sqrt.(Var).-sqrt.(VarEx)).^2));
            println("Var: rank ",r[k], ", cfl ",cfl[i]," : ",errorUIVar[k,i])
        end
    end
end
println("-> DLR Unconventional Integrator DONE.")

maxVal = 0.01;
minVal = min(minimum(errorUI),minimum(errorPS2))*0.995;

for k = 1:length(r)
    if r[k] == 10
        fig = figure("Figure3a",figsize=(10, 7), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    elseif r[k] == 15 
        fig = figure("Figure3b",figsize=(10, 7), dpi=100)#, facecolor='w', edgecolor='k') # dpi Aufloesung
    else
        continue;
    end
    ax = gca()
    if r[k] == 10
        ax.plot(cfl,errorPS[k,:], "g-->", linewidth=2, label="projector-splitting (discretize first)", alpha=1.0)
    end
    ax.plot(cfl,errorPS2[k,:], "r:+", linewidth=2, label="projector-splitting (DLRA first)", alpha=1.0)
    ax.plot(cfl,errorUI[k,:], "b-.o", linewidth=2, label="unconventional", alpha=1.0)
    
    ylabel(L"error $\phi$", fontsize=20)
    xlabel("CFL", fontsize=20)
    ax.set_xlim([cfl[1],cfl[end]])
    #maxVal = max(maximum(errorUI[k,:]),maximum(errorPS2[k,:]))*1.005;
    #minVal = min(minimum(errorUI[k,:]),minimum(errorPS2[k,:]))*0.995;
    if r[k] == 10
        ax.set_ylim([minVal,maxVal])
        ax.legend(loc="upper left", fontsize=20)
    else
        ax.set_ylim([minVal,maxVal])
        ax.legend(loc="upper right", fontsize=20)
    end
    ax.tick_params("both",labelsize=20) 
    fig.canvas.draw() # Update the figure
    tight_layout()
    PyPlot.savefig("studyCFLRank$(r[k]).png")
end

for k = 1:length(r)
    npzwrite("data/ErrorsCFLStudyExpUnconventionalDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)r$(r[k]).npy", errorUI[k,:])
    npzwrite("data/ErrorsCFLStudyExpProjectorSplittingDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)Stabilization0r$(r[k]).npy", errorPS[k,:])
    npzwrite("data/ErrorsCFLStudyExpProjectorSplittingDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)Stabilization1r$(r[k]).npy", errorPS1[k,:])
    npzwrite("data/ErrorsCFLStudyExpProjectorSplittingDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)Stabilization2r$(r[k]).npy", errorPS2[k,:])

    if s.problem == "UQ"
        npzwrite("data/ErrorsVarCFLStudyExpUnconventionalDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)r$(r[k]).npy", errorUIVar[k,:])
        npzwrite("data/ErrorsVarCFLStudyExpProjectorSplittingDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)Stabilization0r$(r[k]).npy", errorPSVar[k,:])
        npzwrite("data/ErrorsVarCFLStudyExpProjectorSplittingDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)Stabilization1r$(r[k]).npy", errorPS1Var[k,:])
        npzwrite("data/ErrorsVarCFLStudyExpProjectorSplittingDLRNx$(s.Nx)N$(N)tEnd$(s.tEnd)Stabilization2r$(r[k]).npy", errorPS2Var[k,:])
    end
end
