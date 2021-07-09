__precompile__
include("quadrature.jl")
include("basis.jl")

using ProgressMeter
using LinearAlgebra

struct DLRSolver
    # spatial grid of cell interfaces
    x::Array{Float64,1};

    # Solver settings
    settings::Settings;

    # preallocate memory for performance
    outRhs::Array{Float64,2};
    
    # squared L2 norms of Legendre coeffs
    gamma::Array{Float64,1};
    # flux matrix PN system
    A::Array{Float64,2};

    # physical parameters
    sigmaT::Float64;
    sigmaS::Float64;

    # constructor
    function DLRSolver(settings)
        x = settings.x;

        outRhs = zeros(settings.NCells,settings.nPN);

        # setup gamma vector
        gamma = zeros(settings.nPN);
        for i = 1:settings.nPN
            n = i-1;
            gamma[i] = 2/(2*n+1);
        end
        
        # setup flux matrix
        A = zeros(settings.nPN,settings.nPN)
        if settings.problem == "UQ" # UQ for advection equation
            N = settings.nPN;
            q = Quadrature(2*N,"Gauss");
            b = Basis(q,settings);
            aXi = ((1.0 - settings.sigmaS) .+settings.sigmaS*q.xi).^3;
            for i = 1:N
                for j = 1:N
                    A[i,j] = IntegralVec(q, aXi.*b.PhiQuad[:,i].*b.PhiQuad[:,j]*0.5,-1.0,1.0)
                end
            end
        else # radiative transfer
            for i = 1:(settings.nPN-1)
                n = i-1;
                A[i,i+1] = (n+1)/(2*n+1)*sqrt(gamma[i+1])/sqrt(gamma[i]);
            end

            for i = 2:settings.nPN
                n = i-1;
                A[i,i-1] = n/(2*n+1)*sqrt(gamma[i-1])/sqrt(gamma[i]);
            end
        end

        # update dt with correct maximal speed lmax
        lmax = maximum(abs.(eigvals(A)));
        settings.dt = settings.dx*settings.cfl/lmax;

        new(x,settings,outRhs,gamma,A,settings.sigmaT,settings.sigmaS);
    end
end

# source term depends only on x now
function Q(obj::DLRSolver,t::Float64,x::Float64)
    y = zeros(obj.settings.nPN)
    return y;
    if obj.settings.problem == "ManufacturedSolution"
        y[1] = 2*exp(t)*cos.(x)/sqrt(obj.gamma[1]);
        y[2] = -2/3 * exp(t)*sin.(x)/sqrt(obj.gamma[2]);
    elseif obj.settings.problem == "ManufacturedSolutionLinear"
        y[1] = 2.0*cos.(x)/sqrt(obj.gamma[1]);
        y[2] = -2/3 * (t+1.0)*sin.(x)/sqrt(obj.gamma[2]);
    elseif obj.settings.problem == "ManufacturedSolutionSteady"
        y[2] = -2/3 * sin.(x)/sqrt(obj.gamma[2]);
    end
    return y;
end

function SetupIC(obj::DLRSolver)
    u = zeros(obj.settings.NCells,obj.settings.nPN); # Nx interfaces, means we have Nx - 1 spatial cells
    u[:,1] = 2.0/sqrt(obj.gamma[1])*IC(obj.settings,obj.settings.xMid);
    return u;
end

function Rhs(obj::DLRSolver,u::Array{Float64,2},t::Float64=0.0)   
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    for j = 2:(obj.settings.NCells-1) # leave out ghost cells
        obj.outRhs[j,:] = (u[j+1,:]-2*u[j,:]+u[j-1,:])/dt/2 - obj.A*(u[j+1,:]-u[j-1,:])/dx/2;

        # add source term and scattering
        obj.outRhs[j,:] += Q(obj,t,obj.settings.xMid[j]) .- obj.sigmaT*u[j,:];
        obj.outRhs[j,1] += obj.sigmaS*u[j,1]
    end

    return obj.outRhs;
end

function RhsStreaming(obj::DLRSolver,u::Array{Float64,2},t::Float64=0.0)   
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    for j = 2:(obj.settings.NCells-1) # leave out ghost cells
        obj.outRhs[j,:] = (u[j+1,:]-2*u[j,:]+u[j-1,:])/dt/2 - obj.A*(u[j+1,:]-u[j-1,:])/dx/2;
    end

    return obj.outRhs;
end

function RhsStreamingLW(obj::DLRSolver,u::Array{Float64,2},t::Float64=0.0)   
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    for j = 2:(obj.settings.NCells-1) # leave out ghost cells
        obj.outRhs[j,:] =  - obj.A*(u[j+1,:]-u[j-1,:])/dx/2;
    end

    return obj.outRhs;
end

function RhsScattering(obj::DLRSolver,u::Array{Float64,2},t::Float64=0.0)   
    dt = obj.settings.dt;
    dx = obj.settings.dx;
    obj.outRhs .= zeros(size(obj.outRhs));
    for j = 2:(obj.settings.NCells-1) # leave out ghost cells
        # add source term and scattering
        obj.outRhs[j,:] += Q(obj,t,obj.settings.xMid[j]) .- obj.sigmaT*u[j,:];
        obj.outRhs[j,1] += obj.sigmaS*u[j,1]
    end

    return obj.outRhs;
end

function Solve(obj::DLRSolver)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;

    # Set up initial condition
    u = SetupIC(obj);
    
    # time loop
    while t < tEnd
        println("t = ",t);

        # Update time by dt
        yU = RhsStreaming(obj,u);
        u = u .+ dt*yU;
        if obj.settings.problem != "UQ"
            yU = RhsScattering(obj,u);
            u = u .+ dt*yU;
        end

        #println("Norm ",norm(u,2))
        
        t = t+dt;
    end

    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*u;

end

function SolveProjectorSplitting(obj::DLRSolver,recordFNorm::Bool=false)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;
    dx = obj.x[2]-obj.x[1];
    Nx = obj.settings.NCells;
    r = obj.settings.r; # DLR rank
    N = obj.settings.nPN; # here, N is the number of quadrature points

    # Set up initial condition
    u = SetupIC(obj);

    # Low-rank approx of init data:
    X,S,W = svd(u); 
    
    # rank-r truncation:
    X = X[:,1:r]; 
    W = W[:,1:r];
    S = Diagonal(S);
    S = S[1:r, 1:r]; 

    K = zeros(Nx,r);
    KNew = X*S;
    L = zeros(N,r);
    
    Nt = Integer(round(tEnd/dt));

    prog = Progress(Nt,1)

    # stabilization matrix
    L1 = zeros(Nx,Nx);
    for j = 1:Nx
        if obj.settings.stabilization != 1
            if j > 1
                L1[j,j-1] = 0.5;
            end
            if j < Nx
                L1[j,j+1] = 0.5;
            end
        else
            L1[j,j] = 1;
        end
    end

    yK = zeros(Nx,r);
    X1 = zeros(r,r);
    X2 = zeros(r,r);

    timeVec = zeros(Nt);
    FNormVec = zeros(Nt);
    
    # time loop
    #@showprogress 0.1 "Progress " 
    #@gif 
    for n = 1:Nt

        ############## Streaming ##############

        ###### K-step ######
        K .= X*S;

        WAW = W'*obj.A*W;
        for j = 2:(obj.settings.NCells-1) # leave out ghost cells
            yK[j,:] = (K[j+1,:]-2*K[j,:]+K[j-1,:])/dt/2- (WAW*(K[j+1,:]-K[j-1,:])/dx/2);
        end        

        K .= K .+ dt*yK;

        X,S = qr(K); # optimize by choosing XFull, SFull
        X = X[:, 1:obj.settings.r]; 
        S = S[1:obj.settings.r, 1:obj.settings.r];

        ###### S-step ######

        for k = 1:r
            for l = 1:r
                X1[k,l] = 0.0;
                X2[k,l] = 0.0;
                for j = 2:(obj.settings.NCells-1)
                    X1[k,l] = X1[k,l] + X[j,k].*(X[j+1,l]-2*X[j,l]+X[j-1,l])/dt/2;
                    X2[k,l] = X2[k,l] + X[j,k].*(X[j+1,l]-X[j-1,l])/dx/2;
                end
            end
        end

        if obj.settings.stabilization == 0
            S .= S .- dt*X1*S .+ dt*X2*S*WAW;
        elseif obj.settings.stabilization == 1
            S .= S .+ dt*X2*S*WAW;
        elseif obj.settings.stabilization == 2
            S .= X'*L1*X*S .+ dt*X2*S*WAW;
        end

        ###### L-step ######
        L .= W*S';
        
        if obj.settings.stabilization == 0 || obj.settings.stabilization == 2
            L .= L .+ dt*L*X1' .- dt*(obj.A'*L*X2');#  .+ dt*(X'*RhsStreamingLW(obj,X*L'))';
        elseif obj.settings.stabilization == 1
            L .= L .- dt*obj.A'*L*X2';
        end

        ############## Scattering ##############
        if obj.settings.problem != "UQ"
            ###### L-step ######
            L .= L .+ dt*(X'*RhsScattering(obj,X*L'))';
        end
                
        W,S = qr(L);
        W = W[:, 1:obj.settings.r];
        S = S[1:obj.settings.r, 1:obj.settings.r];

        S .= S';
        
        next!(prog) # update progress bar
        timeVec[n] = t;
        if recordFNorm
            FNormVec[n] = norm(0.5*sqrt(obj.gamma[1])*S,2);
        end
        t = t+dt;
    end

    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*X*S*W',timeVec,FNormVec;
end

function SolveUnconventional(obj::DLRSolver,recordFNorm::Bool=false)
    t = 0.0;
    dt = obj.settings.dt;
    tEnd = obj.settings.tEnd;
    dx = obj.x[2]-obj.x[1];
    Nx = obj.settings.NCells;
    r = obj.settings.r; # DLR rank
    N = obj.settings.nPN; # here, N is the number of quadrature points


    # Set up initial condition
    u = SetupIC(obj);

    # Low-rank approx of init data:
    X,S,W = svd(u); 
    
    # rank-r truncation:
    X = X[:,1:r]; 
    W = W[:,1:r];
    S = Diagonal(S);
    S = S[1:r, 1:r]; 

    K = zeros(Nx,r);
    KNew = X*S;
    L = zeros(N,r);
    
    Nt = Integer(round(tEnd/dt));

    yK = zeros(Nx,r);
    X1 = zeros(r,r);
    X2 = zeros(r,r);

    timeVec = zeros(Nt);
    FNormVec = zeros(Nt);

    prog = Progress(Nt,1)
    for n = 1:Nt

        ############## Streaming ##############

        ###### K-step ######
        K .= X*S;

        WAW = W'*obj.A*W;
        for j = 2:(obj.settings.NCells-1) # leave out ghost cells
            yK[j,:] = (K[j+1,:]-2*K[j,:]+K[j-1,:])/dt/2- (WAW*(K[j+1,:]-K[j-1,:])/dx/2);
        end        

        K .= K .+ dt*yK;

        XNew,STmp = qr(K); # optimize bei choosing XFull, SFull
        XNew = XNew[:, 1:obj.settings.r]; 

        MUp = XNew' * X;

        ###### L-step ######
        L = W*S';

        for k = 1:r
            for l = 1:r
                X1[k,l] = 0.0;
                X2[k,l] = 0.0;
                for j = 2:(obj.settings.NCells-1)
                    X1[k,l] = X1[k,l] + X[j,k].*(X[j+1,l]-2*X[j,l]+X[j-1,l])/dt/2;
                    X2[k,l] = X2[k,l] + X[j,k].*(X[j+1,l]-X[j-1,l])/dx/2;
                end
            end
        end

        L .= L .+ dt*L*X1' .- dt*(obj.A'*L*X2');
                
        WNew,STmp = qr(L);
        WNew = WNew[:, 1:obj.settings.r]; 

        NUp = WNew' * W;
        W .= WNew;
        X .= XNew;

        ################## S-step ##################
        S .= MUp*S*(NUp')

        WAW = W'*obj.A*W;
        for k = 1:r
            for l = 1:r
                X1[k,l] = 0.0;
                X2[k,l] = 0.0;
                for j = 2:(obj.settings.NCells-1)
                    X1[k,l] = X1[k,l] + X[j,k].*(X[j+1,l]-2*X[j,l]+X[j-1,l])/dt/2;
                    X2[k,l] = X2[k,l] + X[j,k].*(X[j+1,l]-X[j-1,l])/dx/2;
                end
            end
        end

        S .= S .+ dt*X1*S .- dt*X2*S*WAW;

        ############## Scattering ##############
        if obj.settings.problem != "UQ"
        ###### K-step ######
        K .= X*S;

        K .= K .+ dt*RhsScattering(obj,K*W')*W;

        XNew,STmp = qr(K); # optimize bei choosing XFull, SFull
        XNew = XNew[:, 1:obj.settings.r]; 

        MUp = XNew' * X;

        ###### L-step ######
        L = W*S';

        L .= L .+ dt*(X'*RhsScattering(obj,X*L'))';
                
        WNew,STmp = qr(L);
        WNew = WNew[:, 1:obj.settings.r]; 

        NUp = WNew' * W;
        W .= WNew;
        X .= XNew;

        ################## S-step ##################
        S .= MUp*S*(NUp')

        S .= S .+ dt.*X'*RhsScattering(obj,X*S*W')*W;

        end

        timeVec[n] = t;
        if recordFNorm
            FNormVec[n] = norm(0.5*sqrt(obj.gamma[1])*S,2);
        end
        
        next!(prog) # update progress bar

        t = t+dt;
    end

    # return end time and solution
    return t, 0.5*sqrt(obj.gamma[1])*X*S*W',timeVec,FNormVec;

end
