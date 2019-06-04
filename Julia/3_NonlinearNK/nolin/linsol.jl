include("gensys.jl")

function model_solution(para)

    # Theta = [tau kappa psi1 psi2 rA piA gammaQ rho_R rho_g rho_z sigma_R sigma_g sigma_z]
    tau     = para[1]
    kappa   = para[2]
    psi1    = para[3]
    psi2    = para[4]
    rA      = para[5]
    piA     = para[6]
    gammaQ  = para[7]
    rho_R   = para[8]
    rho_g   = para[9]
    rho_z   = para[10]
    sigma_R = para[11]
    sigma_g = para[12]
    sigma_z = para[13]
    bet = 1/(1+rA/400)

    # Equation indices

    eq_1   = 1  #** (2.1) on \hat{y}(t) **/
    eq_2   = 2  #** (2.1) on \hat{pi}(t) **/
    eq_3   = 3  #** (2.1) on \hat{R}(t) **/
    eq_4   = 4  #** \hat{y}(t-1) **/
    eq_5   = 5  #** \hat{g} process **/
    eq_6   = 6  #** \hat{z} process **/
    eq_7   = 7  #** \hat{y} expectational error **/
    eq_8   = 8  #** \hat{pi} expectational error **/

    # Variable indices

    y_t    = 1
    pi_t   = 2
    R_t    = 3
    y1_t   = 4
    g_t    = 5
    z_t    = 6
    Ey_t1  = 7
    Epi_t1 = 8

    # Expectation error indices (eta)

    ey_sh  = 1
    epi_sh = 2

    # Shock indices (eps)

    z_sh = 1
    g_sh = 2
    R_sh = 3

    # SUMMARY

    neq  = 8
    neta = 2
    neps = 3

    # /** initialize matrices **/

    GAM0 = zeros(neq,neq)
    GAM1 = zeros(neq,neq)
       C = zeros(neq,1)
     PSI = zeros(neq,neps)
     PPI = zeros(neq,neta)


    # =========================================================================
    #                 EQUILIBRIUM CONDITIONS: CANONICAL SYSTEM
    # =========================================================================

    # =========================================================================
    #          1.
    # =========================================================================

    GAM0[eq_1,y_t] =  1
    GAM0[eq_1,R_t] =  1/tau
    GAM0[eq_1,g_t] = -(1-rho_g)
    GAM0[eq_1,z_t] = -rho_z/tau
    GAM0[eq_1,Ey_t1] = -1
    GAM0[eq_1,Epi_t1] = -1/tau

    # =========================================================================
    #          2.
    # =========================================================================

    GAM0[eq_2,y_t] = -kappa;
    GAM0[eq_2,pi_t] = 1;
    GAM0[eq_2,g_t] =  kappa;
    GAM0[eq_2,Epi_t1] = -bet;

    # =========================================================================
    #          3.
    # =========================================================================

    GAM0[eq_3,y_t] = -(1-rho_R)*psi2
    GAM0[eq_3,pi_t] = -(1-rho_R)*psi1
    GAM0[eq_3,R_t] = 1
    GAM0[eq_3,g_t] = (1-rho_R)*psi2
    GAM1[eq_3,R_t] = rho_R
    PSI[eq_3,R_sh] = 1

    # =========================================================================
    #          4.
    # =========================================================================

    GAM0[eq_4,y1_t] = 1
    GAM1[eq_4,y_t] = 1

    # =========================================================================
    #          5.
    # =========================================================================

    GAM0[eq_5,g_t] = 1
    GAM1[eq_5,g_t] = rho_g
    PSI[eq_5,g_sh] = 1

    # =========================================================================
    #          6.
    # =========================================================================

    GAM0[eq_6,z_t] = 1
    GAM1[eq_6,z_t] = rho_z
    PSI[eq_6,z_sh] = 1

    # =========================================================================
    #          7.
    # =========================================================================

    GAM0[eq_7,y_t] = 1
    GAM1[eq_7,Ey_t1] = 1
    PPI[eq_7,ey_sh] = 1

    # =========================================================================
    #          8.
    # =========================================================================

    GAM0[eq_8,pi_t] = 1
    GAM1[eq_8,Epi_t1] = 1
    PPI[eq_8,epi_sh] = 1


    # =========================================================================
    #            QZ(generalized Schur) decomposition by GENSYS
    # =========================================================================

    T1, TC, T0, fmat, fwt, ywt, gev, eu, loose = gensys(GAM0,GAM1,C,PSI,PPI)

    return T1, T0, eu

end

function sysmat(T1,T0,para)

    # Theta = [tau kappa psi1 psi2 rA piA gammaQ rho_R rho_g rho_z sigma_R sigma_g sigma_z]
    tau     = para[1]
    kappa   = para[2]
    psi1    = para[3]
    psi2    = para[4]
    rA      = para[5]
    piA     = para[6]
    gammaQ  = para[7]
    rho_R   = para[8]
    rho_g   = para[9]
    rho_z   = para[10]
    sigma_R = para[11]
    sigma_g = para[12]
    sigma_z = para[13]

    eq_y = 1
    eq_pi = 2
    eq_ffr = 3

    # /** number of observation variables **/

    ny = 3

    # /** model variable indices **/

    y_t    = 1
    pi_t   = 2
    R_t    = 3
    y1_t   = 4
    g_t    = 5
    z_t    = 6
    Ey_t1  = 7
    Epi_t1 = 8

    # /** shock indices **/

    z_sh = 1
    g_sh = 2
    R_sh = 3

    # =========================================================================
    #                           TRANSITION EQUATION
    #
    #            s(t) = Phi*s(t-1) + R*e(t)
    #            e(t) ~ iid N(0,Se)
    #
    # =========================================================================

    nep = size(T0,2)

    Phi = T1

    R   = T0

    Se  = zeros(nep,nep)

    Se[z_sh,z_sh] = (sigma_z)^2
    Se[g_sh,g_sh] = (sigma_g)^2
    Se[R_sh,R_sh] = (sigma_R)^2

    # =========================================================================
    #                           MEASUREMENT EQUATION
    #
    #            y(t) = a + b*s(t) + u(t)
    #            u(t) ~ N(0,HH)
    #
    # =========================================================================

    A           = zeros(ny,1)
    A[eq_y,1]   = gammaQ
    A[eq_pi,1]  = piA
    A[eq_ffr,1] = piA+rA+4*gammaQ

    nstate = size(Phi,2)

    B = zeros(ny,nstate)

    B[eq_y,y_t]   = 1
    B[eq_y,y1_t]  = -1
    B[eq_y, z_t]  = 1
    B[eq_pi,pi_t] = 4
    B[eq_ffr,R_t] = 4

    H = zeros(ny,ny)
    # with measurement errors (from dsge1_me.yaml)
    H[eq_y,y_t] = (0.20*0.579923)^2
    H[eq_pi,pi_t] = (0.20*1.470832)^2
    H[eq_ffr,R_t] = (0.20*2.237937)^2

    return A,B,H,R,Se,Phi

end

function g(Phi,sqrtS2,old_states,shocks)

    new_states = Phi*old_states + sqrtS2*shocks

    return new_states

end

function pdfy(para,yt,new_states,old_states)

    rA      = para[5]
    piA     = para[6]
    gammaQ  = para[7]

    eq_y = 1
    eq_pi = 2
    eq_ffr = 3

    # /** model variable indices **/

    y_t    = 1
    pi_t   = 2
    R_t    = 3
    y1_t   = 4
    g_t    = 5
    z_t    = 6
    Ey_t1  = 7
    Epi_t1 = 8

    me = zeros(3)
    # with measurement errors (from dsge1_me.yaml)
    me[eq_y] = 0.20*0.579923
    me[eq_pi] = 0.20*1.470832
    me[eq_ffr] = 0.20*2.237937

	temp = 0.0
    pdf1 = 1.0

	# log gdp growth
    temp = new_states[y_t] - old_states[y_t] + gammaQ + new_states[z_t]
    temp =(1.0/me[eq_y])*(yt[eq_y]-temp)
    pdf1 = pdf1*(1.0/me[eq_y])*exp(-0.5*temp*temp)

	# log inflation (annual)
    temp = 4*new_states[pi_t] + piA
    temp = (1.0/me[eq_pi])*(yt[eq_pi]-temp)
    pdf1 = pdf1*(1.0/me[eq_pi])*exp(-0.5*temp*temp)

	# log nominal rate (annual)
    temp = 4*new_states[R_t] + piA+rA+4*gammaQ
    temp = (1.0/me[eq_ffr])*(yt[eq_ffr]-temp)
    pdf1 = pdf1*(1.0/me[eq_ffr])*exp(-0.5*temp*temp)

    pdf1 = pdf1*1.0/(sqrt(2.0*pi)^3)

    return pdf1

end

function KF(A,B,H,R,Se,Phi,y)

    # Initialize the State Vector at the Stationary Distribution
    T,l    = size(y)
    n,n    = size(Phi)
    s        = zeros(T+1,n)
    P        = zeros(T+1,n,n)
    s[1,:]   = zeros(n,1)'

    a = inv(Matrix(1.0I,n*n,n*n) - kron(Phi,Phi))*reshape(R*Se*R',n*n,1)
    P[1,:,:] = reshape(a,n,n)

    # Kalman Filter Recursion
    sprime             = zeros(n,1)
    Pprime             = zeros(n,n)
    errorprediction    = ones(T,l)
    Varerrorprediction = ones(T,l,l)
    liki               = ones(T,1)
    measurepredi       = ones(T,l)

    for i=1:T

        # Updating Step

        sprime = Phi*s[i,:]
        Pprime = Phi*P[i,:,:]*Phi' + R*Se*R'

        # Prediction Step

        yprediction = A + B*sprime

        v = y[i,:] - yprediction

        F = B*Pprime*B' + H

        kgain    = Pprime*B'*inv(F)
        s[i+1,:] = (sprime + kgain*v)'
        P[i+1,:,:] = Pprime - kgain*B*Pprime
        errorprediction[i,:] = v'
        Varerrorprediction[i,:,:] = F
        temp = 0.5*v'*inv(F)*v
        liki[i] = -0.5*l*log(2*pi) - 0.5*log(det(F)) - temp[1]
        measurepredi[i,:] = y[i,:]-v

    end


    statepredi = s
    varstatepredi = P

    return liki,measurepredi,statepredi,varstatepredi

end
