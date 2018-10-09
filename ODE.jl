# ------------------------------------------------------------------------------ LIBRARIES
using DifferentialEquations
using Plots


# ------------------------------------------------------------------------------ MODEL DEFINITION
function system(du, u, p, t)

# PARAMETERS
    alpha     =p[1]       # growth of labor productivity
    beta      =p[2]       # growth of population
    delta_m   =p[3]       # depreciation rate
    delta_e   =p[4]       #
    f_m       =p[5]       # allowing surplus
    f_e       =p[6]       #
    gamma_m   =p[7]       # energy efficiency
    gamma_e   =p[8]       # (EROI)
    sigma_m   =p[9]       # throughput to capital
    sigma_e   =p[10]      #
    tau_m     =p[11]      # adjustment speed of production
    tau_e     =p[12]      #
    eta_c     =p[13]      # energy efficiency in consumption
    r_b       =p[14]      # interest rate
    phi       =p[15]      # a_e / a_m
    gamma     =p[16]      # money illusion (inflation)
    kappa_0    =p[17]      # output share
    kappa_1    =p[18]      #
    kappa_2    =p[19]      #
    kappa_3    =p[20]      #
    uetoile   =p[21]      #


# FUNCTIONS
    #
    pi_m = sigma_m-r_b*u[4]-sigma_m*(u[3]+u[1]*gamma_m)

    #
    f_tau_m = kappa_0+kappa_1*atan(kappa_2*pi_m+kappa_3)

    #
    util = (gamma_e/gamma_e-1)*(eta_c+1/gamma_m)*((u[2]*sigma_m)/sigma_e)

    #
    if(util <= uetoile)
        f_tau_e = 0.0
    else
        #f_tau_e = -3 +(4*util)
        f_tau_e = max((util-uetoile)/(1-uetoile))
    end

    #
    f_u = (gamma_e/(gamma_e-1))*(eta_c*(1-f_tau_m))*(1-f_tau_e+1/gamma_m)*((u[2]*sigma_m)/sigma_e)

    #
    q_e = ((eta_c*(1-f_tau_m)*(1-f_tau_e+gamma_m)*u[2]*sigma_m)/((1-gamma_e)*sigma_e*f_u))

    #
    f_upsilon = gamma/(1-f_u)-gamma

    #
    i = ((u[1]/(u[1]+eta_c))*f_upsilon+(eta_c/(u[1]+eta_c))*tau_m*(f_m*(u[3]+gamma_m*u[1]))-1)

    #
    f_phi = phi/((1-u[5])^2)-phi


# DYNAMIC
    # rho
    du[1]= u[1]* (f_upsilon-tau_m*(f_m*(u[3])-1))
    # kappa
    du[2]= u[2]* (sigma_m*f_tau_m-u[2]*sigma_e*f_tau_e-delta_m+delta_e)
    # omega
    du[3]= u[3]* (f_phi+gamma*i-alpha-tau_m*(f_m*(u[3]+gamma_m*u[1]))-1)
    # debt
    du[4]= -u[4]* (eta_c*(f_m*(u[3]+u[1]*gamma_m)-1)+delta_m*f_tau_m-delta_m)+
                    (delta_m*f_tau_m-pi_m)
    # lambda
    du[5]= u[5]* ((((sigma_e*sigma_m)/(sigma_m*u[2]+sigma_e*f_u*phi))*
                  (((1/sigma_e)*u[2]*(sigma_m*f_tau_m-delta_m))+
                   ((1/sigma_m)*f_u*phi*tau_e*(f_e*q_e-1)))) -alpha -beta)

 end


# ------------------------------------------------------------------------------ INITIALISATION
# Time
ending = 3.0
tspan=(0.0, ending)

# Initial state
rho0= 1.0               # relative price of energy
kappa0= 0.9             # Km / Ke
omega0= 0.78            # wage share
debt0= 0.1              # debt
lambda0= 0.75           # employment rate

# Parameters
alpha     = 0.1
beta      = 0.1
delta_m   = 0.0521
delta_e   = 0.0521
f_m       = 1.3
f_e       = 1.3
gamma_m   = 5.0
gamma_e   = 50.0
sigma_m   = 0.33
sigma_e   = 0.33
tau_m     = 4.0
tau_e     = 4.0
eta_c     = 1.0
r_b       = 0.03
phi       = 1.0
gamma     = 0.7
kappa_0   = -0.0065
kappa_1   = exp(-5)
kappa_2   = 20.0
kappa_3   = 11.991
uetoile   = 0.75

#
init=[rho0; kappa0; omega0; debt0; lambda0]
parameters = [alpha; beta; delta_m; delta_e; f_m; f_e; gamma_m; gamma_e;
              sigma_m; sigma_e; tau_m; tau_e; eta_c; r_b; phi; gamma;
              kappa_0; kappa_1; kappa_2; kappa_3; uetoile]

#
problem = ODEProblem(system, init, tspan, parameters)


# ------------------------------------------------------------------------------ SOLVE
#result = solve(problem, alg_hints = [:stiff])
#result = solve(problem, Rosenbrock32())
result = solve(problem)


# ------------------------------------------------------------------------------ PLOT
plot(result, ylims=(-0.5,1.5))
