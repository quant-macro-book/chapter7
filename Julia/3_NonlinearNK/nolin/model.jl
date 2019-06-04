mutable struct model

    np::Int64
    ns::Int64
    nexog::Int64
    ngh::Int64
    nvars::Int64
    nv::Int64
    nz::Int64
    pctrn::Float64
    mg::Float64
    mz::Float64
    mr::Float64
    damp::Float64
    tol::Float64
    zlbflag::Bool

#     yy::Array{Float64,2}
#     T::Int64
#     me::Array{Float64,1}

    xgrid::Array{Float64,2}
    bbt::Array{Float64,2}
    bbtinv::Array{Float64,2}
    slopecon::Array{Float64,2}
    coeffcn::Array{Float64,1}
    coeffpn::Array{Float64,1}
    coeffcb::Array{Float64,1}
    coeffpb::Array{Float64,1}

    model() = new()

end

function new_model()

    m = model()
    # number of variables
    m.np = 2
    m.ns = 4
    m.nexog = 3
    m.ngh = 3
    m.nvars = m.ns+m.nexog
    m.nv = 1+m.np*m.ns
    m.nz = m.ngh^m.nexog

    # m.pctrn = 0.05*100
    m.pctrn = 0.05
    m.mg = 2.0
    m.mz = 2.0
    m.mr = 2.0

    m.damp = 0.7
    m.tol = 1e-5
    m.zlbflag = 1

#     m.yy = readdlm("./us.txt")
#     m.T = size(m.yy)[1]
#     m.me = zeros(3)
#     m.me[1] = std(m.yy[:,1]/100)*0.2
#     m.me[2] = std(m.yy[:,2]/100)*0.2
#     m.me[3] = std(m.yy[:,3]/100)*0.2

    return m

end

function qnwnorm1(n::Int64,mu::Float64,sigma::Float64)
# Originally from Fackler and Miranda's Compecon toolbox for MATLAB
# Copyright (c) 1997-2010, Paul L. Fackler & Mario J. Miranda

# Based on an algorithm in W.H. Press, S.A. Teukolsky, W.T. Vetterling
# and B.P. Flannery, "Numerical Recipes in FORTRAN", 2nd ed.  Cambridge
# University Press, 1992.
    maxit::Int64 = 100
    pim4::Float64 = 1/pi^0.25
    m::Int64 = floor((n+1)/2)
    x = zeros(Float64,n)
    w = zeros(Float64,n)
    z::Float64 = 0.0
    for i::Int64=1:m
        # Reasonable starting values
        if (i==1)
            z = sqrt(2*n+1)-1.85575*((2*n+1)^(-1/6))
        elseif (i==2)
            z = z-1.14*(n^0.426)/z
        elseif (i==3)
            z = 1.86*z+0.86*x[1]
        elseif (i==4)
            z = 1.91*z+0.91*x[2]
        else
            z = 2*z+x[i-2]
        end
        # root finding iterations
        its::Int64 = 0
        pp::Float64 = 0
        # type for p1,p2,p3,j,z1?
        # println([its maxit])
        while (its<maxit)
            its = its+1
            p1 = pim4
            p2 = 0
            for j=1:n
                p3 = p2
                p2 = p1
                p1 = z*sqrt(2/j)*p2-sqrt((j-1)/j)*p3
            end
            pp = sqrt(2*n)*p2
            # println([its pp])
            z1 = z
            z  = z1-p1/pp
            if (abs(z-z1)<1e-14)
                break
            end
        end
        # if its>=maxit
        #    error('failure to converge in qnwnorm1')
        # end
        x[n+1-i] = z
        x[i] = -z
        w[i] = 2/(pp*pp)
        w[n+1-i] = w[i]
    end
    # println(size(x))
    w = w/sqrt(pi)
    x = x*sqrt(2)

    x = mu*ones(Float64,n) + sigma*x

    return x, w

end

function poly2s(xx,np,ns)

    nv = 1+np*ns
    f = zeros(1,nv)
    f[1] = 1.0

    for is = 1:ns

        x = xx[is]
        f[2*is] = x
        f[2*is+1] = 2*x^2-1

    end

    return f
end

function makegrid(np,ns)

    nv = 1+np*ns
    xgrid = zeros(nv,ns)
    x = [0 -1 1]

    iv = 2

    for is=1:ns

        for ip=1:np

            xgrid[iv,is] = x[ip+1]
            iv = iv+1

        end

    end

    return xgrid
end

function steadystate(para)

    # Theta = [tau kappa psi1 psi2 rA piA gammaQ rho_R rho_g rho_z sigma_R sigma_g sigma_z]
    invnu   = 6.0
    gyss    = 0.2
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
    sigma_R = para[11]/100
    sigma_g = para[12]/100
    sigma_z = para[13]/100

    gbar = 1.0/(1.0-gyss)
    bet  = 1/(1+rA/400)
    piss = 1+piA/400
    ass  = 1+gammaQ/100
    rnss = ass*piss/bet
    phi  = tau*(invnu-1)/piss^2/kappa # to be checked

    # variables rn0, c0, pi0, y0, rnot0, fc0, fp0, gnow, znow, rnow
    css  = (1.0-1.0/invnu)^(1.0/tau)
    yss  = gbar*css
    fcss = bet*css^(-tau)/ass/piss
    fpss = bet*phi*css^(-tau)*yss*(piss-piss)*piss

    ss = [rnss css piss yss rnss fcss fpss 0.0 0.0 0.0]

    return ss
end

function pf(fc0,fp0,piss,invnu,tau,phi)

    try
        global c0 = fc0^(-1/tau)
    catch
        global c0 = 1.0 # css
        conv = false
    end
    # solve quadratic equation for pi
    a0 = .5*phi*piss^2*invnu + (1-invnu) + invnu*c0^tau + fp0*c0^tau
    a1 = -phi*piss*(1-invnu)
    a2 = -phi*(1-.5*invnu)
    try
        global pi0 = (a1-sqrt(a1^2-4*a0*a2))/2/a2
    catch
        global pi0 = piss
        conv = false
    end

    return c0,pi0

end

function decr(endogvarm,shocks,para,slopecon,coeffcn,coeffpn,coeffcb,coeffpb,np,ns,zlbflag)

    # Theta = [tau kappa psi1 psi2 rA piA gammaQ rho_R rho_g rho_z sigma_R sigma_g sigma_z]
    invnu   = 6.0
    gyss    = 0.2
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
    sigma_R = para[11]/100
    sigma_g = para[12]/100
    sigma_z = para[13]/100

    gbar = 1/(1-gyss)
    bet  = 1/(1+rA/400)
    piss = 1+piA/400
    ass  = 1+gammaQ/100
    rnss = ass*piss/bet
    phi  = tau*(invnu-1)/piss^2/kappa # to be checked

    rn0 = endogvarm[1]
    gnow = endogvarm[8]
    znow = endogvarm[9]
    # gp  = rho_g*gnow + sigma_g*shocks[1]
    # zp  = rho_z*znow + sigma_z*shocks[2]
    # rp  = sigma_R*shocks[3]
    gp  = rho_g*gnow + shocks[1]
    zp  = rho_z*znow + shocks[2]
    rp  = shocks[3]
    ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow)

    xrn = slopecon[1,1]*rn0 + slopecon[1,2]
    xgp = slopecon[2,1]*gp + slopecon[2,2]
    xzp = slopecon[3,1]*zp + slopecon[3,2]
    xrp = slopecon[4,1]*rp + slopecon[4,2]

    # first assume the ZLB is not binding, and use coeffcn and coeffpn
    temp = poly2s([xrn xgp xzp xrp],np,ns)*coeffcn
    fc1 = temp[1]
    temp = poly2s([xrn xgp xzp xrp],np,ns)*coeffpn
    fp1 = temp[1]
    # next period's c and pi (obtained by next period's fc and fp)
    c1,pi1 = pf(fc1,fp1,piss,invnu,tau,phi)
    y1 = c1/(1/gbar/exp(gp) - phi/2*(pi1-piss)^2)
    try
        global rn1 = rn0^rho_R*( rnss*(pi1/piss)^psi1*(y1/ystar)^psi2 )^(1-rho_R)*exp(rp)
    catch
        global rn1 = rnss
        conv = false
    end

    endogvarp = [rn1 c1 pi1 y1 rn1 fc1 fp1 gp zp rp]

    # then check if the ZLB is violated by rn1, and use coeffcb and coeffpb instead
    if (rn1<1.0)

        temp = poly2s([xrn xgp xzp xrp],np,ns)*coeffcb
        fc1 = temp[1]
        temp = poly2s([xrn xgp xzp xrp],np,ns)*coeffpb
        fp1 = temp[1]

        # next period's c and pi (obtained by next period's fc and fp)
        c1,pi1 = pf(fc1,fp1,piss,invnu,tau,phi)
        y1 = c1/(1/gbar/exp(gp) - phi/2*(pi1-piss)^2)
        try
            global rn1 = rn0^rho_R*( rnss*(pi1/piss)^psi1*(y1/ystar)^psi2 )^(1-rho_R)*exp(rp)
        catch
            global rn1 = rnss
            conv = false
        end

        if (zlbflag)
            endogvarp = [rn1 c1 pi1 y1 1.0 fc1 fp1 gp zp rp]
        else
            endogvarp = [rn1 c1 pi1 y1 rn1 fc1 fp1 gp zp rp]
        end

    end

    return endogvarp

end

function solve(m::model,para)

    invnu   = 6.0
    gyss    = 0.2
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
    sigma_R = para[11]/100
    sigma_g = para[12]/100
    sigma_z = para[13]/100

    gbar = 1/(1-gyss)
    bet  = 1/(1+rA/400)
    piss = 1+piA/400
    ass  = 1+gammaQ/100
    rnss = ass*piss/bet
    phi  = tau*(invnu-1)/piss^2/kappa # to be checked

    # steady state values
    ss = steadystate(para)
    rnss = ss[1]
    css  = ss[2]
    piss = ss[3]
    yss  = ss[4]
    rnss = ss[5]
    fcss = ss[6]
    fpss = ss[7]

    # set up grid points
    m.xgrid = makegrid(m.np,m.ns)
    m.bbt = zeros(m.nv,m.nv)
    for iv = 1:m.nv

        m.bbt[iv,:] = poly2s(m.xgrid[iv,:],m.np,m.ns)

    end
    m.bbtinv = inv(m.bbt)

    # set bounds
    rnmin = (1-m.pctrn)*rnss
    rnmax = (1+m.pctrn)*rnss
    gmin = -m.mg*sigma_g/sqrt(1-rho_g^2)
    gmax = m.mg*sigma_g/sqrt(1-rho_g^2)
    zmin = -m.mz*sigma_z/sqrt(1-rho_z^2)
    zmax = m.mz*sigma_z/sqrt(1-rho_z^2)
    rmin = -m.mr*sigma_R
    rmax = m.mr*sigma_R

    m.slopecon = zeros(m.ns,2)
    m.slopecon[1,1] = 2/(rnmax-rnmin)
    m.slopecon[1,2] = -(rnmax+rnmin)/(rnmax-rnmin)
    m.slopecon[2,1] = 2/(gmax-gmin)
    m.slopecon[2,2] = -(gmax+gmin)/(gmax-gmin)
    m.slopecon[3,1] = 2/(zmax-zmin)
    m.slopecon[3,2] = -(zmax+zmin)/(zmax-zmin)
    m.slopecon[4,1] = 2/(rmax-rmin)
    m.slopecon[4,2] = -(rmax+rmin)/(rmax-rmin)

    # gh nodes and weights
    x,w = qnwnorm1(m.ngh,0.0,1.0)
    ghweights = zeros(m.nz)
    ghnodes = zeros(m.nz,m.nexog)
    for ighr = 1:m.ngh

        for ighz = 1:m.ngh

            for ighg = 1:m.ngh

                index = m.ngh^2*(ighr-1)+m.ngh*(ighz-1)+ighg
                ghweights[index] = w[ighg]*w[ighz]*w[ighr]
                ghnodes[index,:] = [x[ighg] x[ighz] x[ighr]]

            end

        end

    end


    # initial values
    # NOTE: We use an index-function approach with a pair of policy functions.
    # One assumes the ZLB always binds and the other assumes the ZLB never
    # binds. The next period's policy function is given by a weighted average
    # of the policy function in the ZLB regime and the policy function in the
    # non-ZLB regime with an indicator function. The value of the indicator
    # function is one when the notional rate is greater than the ZLB,
    # otherwise zero.
    cvec0n  = css*ones(m.nv)
    pivec0n = piss*ones(m.nv)
    rnvec0n = rnss*ones(m.nv)
    yvec0n  = yss*ones(m.nv)
    cvec0b  = css*ones(m.nv)
    pivec0b = piss*ones(m.nv)
    rnvec0b = rnss*ones(m.nv)
    yvec0b  = yss*ones(m.nv)
    cvec1n  = zeros(m.nv)
    pivec1n = zeros(m.nv)
    rnvec1n = zeros(m.nv)
    yvec1n  = zeros(m.nv)
    cvec1b  = zeros(m.nv)
    pivec1b = zeros(m.nv)
    rnvec1b = zeros(m.nv)
    yvec1b  = zeros(m.nv)

    fcvec0n = fcss*ones(m.nv)
    fpvec0n = fpss*ones(m.nv)
#     fcvec0b = bet*css^(-tau)/ass/piss*ones(m.nv) # assumes rn0 = 1.0
    fcvec0b = fcss*ones(m.nv) # assumes rn0 = 1.0
    fpvec0b = fpss*ones(m.nv)
    fcvec1n = zeros(m.nv)
    fpvec1n = zeros(m.nv)
    fcvec1b = zeros(m.nv)
    fpvec1b = zeros(m.nv)

    diff = 1e+4
    iter = 0
    global conv = true

    while ((diff>m.tol) && (iter<1000))

        # fitting polynomials
        m.coeffcn = m.bbtinv*fcvec0n
        m.coeffpn = m.bbtinv*fpvec0n
        m.coeffcb = m.bbtinv*fcvec0b
        m.coeffpb = m.bbtinv*fpvec0b

        for iv = 1:m.nv

            # variables rnot0, c0, pi0, y0, rn0, fc0, fp0, gnow, znow, rnow
            rnpast = (rnmax-rnmin)/2*m.xgrid[iv,1] + (rnmax+rnmin)/2
            gnow = (gmax-gmin)/2*m.xgrid[iv,2] + (gmax+gmin)/2
            znow = (zmax-zmin)/2*m.xgrid[iv,3] + (zmax+zmin)/2
            rnow = (rmax-rmin)/2*m.xgrid[iv,4] + (rmax+rmin)/2
            # natural output
            # BUG: it should be inside the expectation loop too
            ystar = (1-1/invnu)^(1/tau)*gbar*exp(gnow)

            # current period's c and pi (obtained by current period's fc and fp)
            # in the non-ZLB regime
            fc0n = fcvec0n[iv]
            fp0n = fpvec0n[iv]
            c0n,pi0n = pf(fc0n,fp0n,piss,invnu,tau,phi)
            y0n = c0n/(1/gbar/exp(gnow) - phi/2*(pi0n-piss)^2)
            try
                global rn0n = rnpast^rho_R*( rnss*(pi0n/piss)^psi1*(y0n/ystar)^psi2 )^(1-rho_R)*exp(rnow)
            catch
                global rn0n = rnss
                conv = false
            end
            # endogvarn = [rn0n c0n pi0n y0n rn0n fc0n fp0n gnow znow rnow]

            # in the ZLB regime
            fc0b = fcvec0b[iv]
            fp0b = fpvec0b[iv]
            c0b,pi0b = pf(fc0b,fp0b,piss,invnu,tau,phi)
            y0b = c0b/(1/gbar/exp(gnow) - phi/2*(pi0b-piss)^2)
            try
                global rn0b = rnpast^rho_R*( rnss*(pi0b/piss)^psi1*(y0b/ystar)^psi2 )^(1-rho_R)*exp(rnow)
            catch
                global rn0b = rnss
                conv = false
            end
            # endogvarb = [rn0b c0b pi0b y0b rn0b fc0b fp0b gnow znow rnow]

            # update the expectation terms fc and fp with interpolation
            shocks = zeros(3)
            fc0n = 0.0
            fp0n = 0.0
            fc0b = 0.0
            fp0b = 0.0
            for iz = 1:m.nz

                shocks[1] = sigma_g*ghnodes[iz,1]
                shocks[2] = sigma_z*ghnodes[iz,2]
                shocks[3] = sigma_R*ghnodes[iz,3]

                # in the non-ZLB regime
                endogvarn = [rn0n c0n pi0n y0n rn0n fc0n fp0n gnow znow rnow]
                endogvarp = decr(endogvarn,shocks,para,m.slopecon,m.coeffcn,m.coeffpn,m.coeffcb,m.coeffpb,m.np,m.ns,m.zlbflag)
                c1n = endogvarp[2]
                pi1n = endogvarp[3]
                y1n = endogvarp[4]
                zp = endogvarp[9]
                fcxn = c1n^(-tau)*rn0n/(ass*exp(zp))/pi1n # NOTE: max operator is needed?
                fpxn = c1n^(-tau)*y1n*(pi1n-piss)*pi1n/y0n

                # in the ZLB regime
                endogvarb = [rn0b c0b pi0b y0b 1.0 fc0b fp0b gnow znow rnow]
                endogvarp = decr(endogvarb,shocks,para,m.slopecon,m.coeffcn,m.coeffpn,m.coeffcb,m.coeffpb,m.np,m.ns,m.zlbflag)
                c1b = endogvarp[2]
                pi1b = endogvarp[3]
                y1b = endogvarp[4]
                zp = endogvarp[9]
#                 fcxb = c1b^(-tau)*1.0/(ass*exp(zp))/pi1b
                fcxb = c1b^(-tau)*rn0b/(ass*exp(zp))/pi1b
                fpxb = c1b^(-tau)*y1b*(pi1b-piss)*pi1b/y0b

                fc0n = fc0n + ghweights[iz]*bet*fcxn
                fp0n = fp0n + ghweights[iz]*bet*phi*fpxn
                fc0b = fc0b + ghweights[iz]*bet*fcxb
                fp0b = fp0b + ghweights[iz]*bet*phi*fpxb

            end # for iz

            cvec1n[iv]  = c0n
            pivec1n[iv] = pi0n
            yvec1n[iv]  = y0n
            rnvec1n[iv] = rn0n
            fcvec1n[iv] = fc0n
            fpvec1n[iv] = fp0n

            cvec1b[iv]  = c0b
            pivec1b[iv] = pi0b
            yvec1b[iv]  = y0b
            rnvec1b[iv] = rn0b
            fcvec1b[iv] = fc0b
            fpvec1b[iv] = fp0b

        end # for iv

        # calculate the norm between the old and new policy functions
        diffcn = maximum(abs.(cvec1n-cvec0n))
        diffpn = maximum(abs.(pivec1n-pivec0n))
        diffrn = maximum(abs.(rnvec1n-rnvec0n))
        diffyn = maximum(abs.(yvec1n-yvec0n))
        diffn  = maximum([diffcn diffpn diffrn diffyn])

        diffcb = maximum(abs.(cvec1b-cvec0b))
        diffpb = maximum(abs.(pivec1b-pivec0b))
        diffrb = maximum(abs.(rnvec1b-rnvec0b))
        diffyb = maximum(abs.(yvec1b-yvec0b))
        diffb  = maximum([diffcb diffpb diffrb diffyb])

        diff = maximum([diffn diffb])

        # update the policy functions
        cvec0n  = m.damp*cvec0n  + (1.0-m.damp)*cvec1n
        pivec0n = m.damp*pivec0n + (1.0-m.damp)*pivec1n
        rnvec0n = m.damp*rnvec0n + (1.0-m.damp)*rnvec1n
        yvec0n  = m.damp*yvec0n  + (1.0-m.damp)*yvec1n
        fcvec0n = m.damp*fcvec0n + (1.0-m.damp)*fcvec1n
        fpvec0n = m.damp*fpvec0n + (1.0-m.damp)*fpvec1n

        cvec0b  = m.damp*cvec0b  + (1.0-m.damp)*cvec1b
        pivec0b = m.damp*pivec0b + (1.0-m.damp)*pivec1b
        rnvec0b = m.damp*rnvec0b + (1.0-m.damp)*rnvec1b
        yvec0b  = m.damp*yvec0b  + (1.0-m.damp)*yvec1b
        fcvec0b = m.damp*fcvec0b + (1.0-m.damp)*fcvec1b
        fpvec0b = m.damp*fpvec0b + (1.0-m.damp)*fpvec1b

        # counter for iterations
        iter = iter + 1

        println([iter diffn diffb])

        if (any(isnan,[fcvec0n fpvec0n fcvec0b fpvec0b])==true)
            conv = false
        end

        if (typeof([fcvec0n fpvec0n fcvec0b fpvec0b])==Array{Complex{Float64},2})
            conv = false
        end

        if (conv==false) break end

    end # while

    if (iter>=1000) conv = false end

#     println(ghnodes)

    return conv
end

# function g(m,para,old_states,shocks)
#
#     new_states = decr(old_states,shocks,para,m.slopecon,m.coeffcn,m.coeffpn,m.coeffcb,m.coeffpb,m.np,m.ns,m.zlbflag)
#
#     return new_states
#
# end
#
# function pdfy(para,yt,new_states,old_states)
#
#     rA      = para[5]
#     piA     = para[6]
#     gammaQ  = para[7]
#
#     eq_y = 1
#     eq_pi = 2
#     eq_ffr = 3
#
#     # /** model variable indices **/
#
#     y_t    = 4
#     pi_t   = 3
#     R_t    = 5
#     g_t    = 8
#     z_t    = 9
#
#     me = zeros(3)
#     # with measurement errors (from dsge1_me.yaml)
#     me[eq_y] = 0.20*0.579923
#     me[eq_pi] = 0.20*1.470832
#     me[eq_ffr] = 0.20*2.237937
#
# 	temp = 0.0
#     pdf1 = 1.0
#
# 	# log gdp growth
#     temp = new_states[y_t] - old_states[y_t] + gammaQ + new_states[z_t]
#     temp =(1.0/me[eq_y])*(yt[eq_y]-temp)
#     pdf1 = pdf1*(1.0/me[eq_y])*exp(-0.5*temp*temp)
#
# 	# log inflation (annual)
#     temp = 4*new_states[pi_t] + piA
#     temp = (1.0/me[eq_pi])*(yt[eq_pi]-temp)
#     pdf1 = pdf1*(1.0/me[eq_pi])*exp(-0.5*temp*temp)
#
# 	# log nominal rate (annual)
#     temp = 4*new_states[R_t] + piA+rA+4*gammaQ
#     temp = (1.0/me[eq_ffr])*(yt[eq_ffr]-temp)
#     pdf1 = pdf1*(1.0/me[eq_ffr])*exp(-0.5*temp*temp)
#
#     pdf1 = pdf1*1.0/(sqrt(2.0*pi)^3)
#
#     return pdf1
#
# end
