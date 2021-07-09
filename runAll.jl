# runs simulations for radiation transport with CFL = 1
include("runRadTransport.jl")

# runs simulations for UQ with CFL = 1
include("runUQ.jl")

# runs CFL study for radiation transport
include("driverRadTransport.jl")

# runs CFL study for UQ
include("driverUQ.jl")