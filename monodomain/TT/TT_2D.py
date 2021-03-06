#Solving monodomain equation in Firedrake
#Justin Crum
#4/16/21

#To do:
#1.  Implement iIon as in Ten Tusscher
#3.  Need to set up initial conditions including what the activiation in
#    a corner will look like.

from firedrake import *
import argparse
import numpy as np
from irksome import GaussLegendre, Dt, TimeStepper, RadauIIA , LobattoIIIC


# parser = argparse.ArgumentParser(
#         description="Allows for input of order and mesh refinement.")
# parser.add_argument("-O", "--Order", 
#         type=int, help="Input the order of the polynomials.")
# args = parser.parse_args()

mesh = RectangleMesh(100, 100, 70, 70, quadrilateral=True)
#polyOrder = args.Order
polyOrder = 1
outfile = File("TTmonodomain2D.pvd")

#Set up the function space and test/trial functions.
V = FunctionSpace(mesh, "Q", polyOrder)
Z = V * V * V * V * V * V * V * V * V * V * V * V * V * V * V * V * V * V * V

x, y = SpatialCoordinate(mesh)

InitialPotential = conditional(Or(x < 3.5, And(And(31<= x, x < 39), And(0 <= y, y < 35))), Constant(2.0), Constant(-86.709))

uu = Function(Z)
vu, vnaic, vcaic, vkic, vcasr, vcass, vrbar, vxr1, vxr2, vxs, vm, vh, vj, vd, vf, vf2, vfcass, vs, vr = TestFunctions(Z)

uu.sub(0).interpolate(InitialPotential) #voltage
uu.sub(1).interpolate(Constant(8.604))  #NaIC
uu.sub(2).interpolate(Constant(0.000126)) #CaICT
uu.sub(3).interpolate(Constant(136.89)) #KIC
uu.sub(4).interpolate(Constant(3.64))  #CaSRT
uu.sub(5).interpolate(Constant(0.00036)) #CaSST
uu.sub(6).interpolate(Constant(0.9073)) #rBar
uu.sub(7).interpolate(Constant(0.00621)) #xR1
uu.sub(8).interpolate(Constant(0.4712)) #xR2
uu.sub(9).interpolate(Constant(0.0095)) #xS
uu.sub(10).interpolate(Constant(0.00172)) #m
uu.sub(11).interpolate(Constant(0.7444)) #h
uu.sub(12).interpolate(Constant(0.7045)) #j
uu.sub(13).interpolate(Constant(3.375e-5)) #d
uu.sub(14).interpolate(Constant(0.7888)) #f
uu.sub(15).interpolate(Constant(0.9755)) #f2
uu.sub(16).interpolate(Constant(0.9953)) #fCass
uu.sub(17).interpolate(Constant(0.999998)) #s
uu.sub(18).interpolate(Constant(2.42e-8)) #r

(uCurr, NaIC, CaICT, KIC, CaSRT, CaSST, rBar, xR1, xR2, xS, m, h, j, d, f, f2, fCass, s, r) = split(uu)

butcher_tableau = LobattoIIIC(2) #RadauIIA(2) #RadauIIA(1)  #LobattoIIIC(2)
ns = butcher_tableau.num_stages

x, y = SpatialCoordinate(mesh)
dt = Constant(0.01)
t = Constant(0.0)

#Model parameters, taken from the literature.
chi = Constant(1.)    
capacitance = Constant(1.)
sigma11 = Constant(1.0)
sigma22 = Constant(0.25)
sigma12 = Constant(0.5)
sigma21 = Constant(0.5)
sigma =  as_matrix([[sigma11, sigma12], [sigma21, sigma22]])

#Constants I will need.
R = Constant(8.3143) #gas constant, J K**-1 mol**-1
T = Constant(310) #K, temperature
F = Constant(96.4867) #C/mmol
zCa = Constant(2.)
zK = Constant(1.)
zNa = Constant(1.)

gamma = Constant(0.35) #voltage dependance parameter for iNaCa
CaO = Constant(2.) #mM
NaO = Constant(140.) #mM
kO = Constant(5.4) #mM

alpha = Constant(2.5)

vC = Constant(16.404) #mu m^3
vSR = Constant(1.094) #mu m^3
vSS = Constant(0.05468) #mu m^3

gNa = Constant(14.838) #nS/pF 
gTo = Constant(0.294) #nS/pF, we use the **epicardial** model, same as the Niederer paper.
gKr = Constant(0.096) #units are nS/pF.
gK1 = Constant(5.405) #nS/pF.
gpK = Constant(0.0146) #nS/pF
gbNa = Constant(0.00029) #nS/pF
gbCa = Constant(0.000592) #nS/pF
gCaL = Constant(3.98 * 10**(-5)) #cm / (mS microF)
gKs = Constant(0.392) #nS/pF
gpCa = Constant(0.025) #nS/pF

pNaK = Constant(1.362) #pA/pF
pKNa = Constant(0.03)

kmK = Constant(1.0) #mM
kmNa = Constant(40.) #mM
kpCa = Constant(0.0005) #mM
kmNai = Constant(87.5) #mM
kmCa = Constant(1.38) #mM
kSat = Constant(0.1)

vLeak = Constant(0.00036) #mM/ms
vMaxup = Constant(0.006375) #mM/ms
kUp = Constant(0.00025) #mM
vRel = Constant(40.8) #mM/ms
vXfer = Constant(0.0038) #mM/ms

k1prime = Constant(0.15)
k2prime = Constant(0.045)
k3 = Constant(0.60)
k4 = Constant(0.000015)

maxSR = Constant(2.5)
minSR = Constant(1.)

EC = Constant(1.5)
BUFC = Constant(0.2)
kBUFC = Constant(0.001)
BUFSR = Constant(10.)
kBUFSR = Constant(0.3)
BUFSS = Constant(0.4)
kBUFSS = Constant(0.00025)



CaIC = CaICT- BUFC
CaSS = CaSST - BUFSS
CaSR = CaSRT - BUFSR

eCa = R * T / (zCa * F) * ln(CaO / CaIC)
eNa = R * T / (zNa * F) * ln(NaO/NaIC)
eK = R * T / (zK * F) * ln(kO/KIC)
eKs = R * T / F * ln((kO + pKNa * NaO) / (KIC + pKNa * NaIC))

mInf = 1 / (1 + exp((-56.86 - uCurr)/9.03))**2
hInf = 1 / (1 + exp( (uCurr + 71.55)/7.43))**2
jInf = 1 / (1 + exp( (uCurr + 71.55)/7.43))**2

alphaM = 1 / (1 + exp((-60-uCurr)/5))
betaM = 0.1 / (1 + exp((uCurr+35)/5)) + 0.1 / (1 + exp((uCurr-50)/200))
tauM = alphaM * betaM

alphaHLT =  0.057 * exp(-(uCurr + 80)/6.8)
alphaHGT = 0
betaHGT =  0.77 / (0.13 * (1 + exp(-(uCurr + 10.66)/11.1)))
betaHLT = 2.7 * exp(0.079 * uCurr) + 3.1e5 * exp(0.3485 * uCurr)
alphaH = conditional(uCurr >= -40, alphaHGT, alphaHLT)
betaH = conditional(uCurr >= -40, betaHGT, betaHLT)
tauH = 1 / (alphaH + betaH)

alphaJGT = 0
alphaJLT = (uCurr + 37.78) * (-2.5428e4*exp(0.2444*uCurr) - 6.948e-6*exp(-0.04391*uCurr)) / (1 + exp(0.311 * (uCurr + 79.23)))
betaJGT = (0.6 * exp(0.057 * uCurr)) / (1 + exp(-0.1*(uCurr+32)))
betaJLT = (0.02424 * exp(-0.01052 * uCurr)) / (1 + exp(-0.1378 * (uCurr + 40.14)))
alphaJ = conditional(uCurr >= - 40, alphaJGT, alphaJLT)
betaJ = conditional(uCurr >= - 40, betaJGT, betaJLT)
tauJ = 1 / (alphaJ + betaJ)

#iNa = gNa * mInf**3 * hInf * jInf * (uCurr - eNa)
iNa = gNa * m**3 * h * j * (uCurr - eNa)

#The transient outward current, iTo.
rInf = 1 / (1 + exp((20 - uCurr)/6))
sInf = 1 / (1 + exp((uCurr + 20)/5))
tauR = 9.5 * exp(-(uCurr + 40)**2/1800) + 0.8
tauS = 85 * exp(-(uCurr+45)**2 / 320) + 5 / (1 + exp((uCurr - 20) / 5)) + 3

iTo = gTo * rInf * sInf * (uCurr - eK)

#Rapid delayed rectifier current
xR1Inf = 1 / (1 + exp((-26-uCurr)/7))
xR2Inf = 1 / (1 + exp((uCurr+88)/24))

alphaxR1 = 450 / (1 + exp((-45-uCurr)/10))
betaxR1 = 6 / (1 + exp((uCurr+30)/11.5))
tauxR1 = alphaxR1 * betaxR1

alphaxR2 = 3 / (1 + exp((-60-uCurr)/20))
betaxR2 = 1.12 / (1 + exp((uCurr - 60)/20))
tauxR2 = alphaxR2 * betaxR2

iKr = gKr * sqrt(kO / 5.4) * xR1Inf * xR2Inf * (uCurr - eK)

#The inward rectifier current.
alphaK1 = 0.1 / (1 + exp(0.06 * (uCurr - eK - 200)))
betaK1 = (3 * exp(0.0002 * (uCurr - eK +100)) + exp(0.1 * (uCurr - eK -10))) / (1 + exp(-0.5 * (uCurr - eK)))
xK1Inf = alphaK1 / (alphaK1 + betaK1)

iK1 = gK1 * sqrt(kO / 5.4) * xK1Inf * (uCurr - eK)

#The Na+/Ca2+ pump current.
iNaK = pNaK * kO * NaIC / ((kO + kmK) * (NaIC + kmNa) * (1 + 0.1245 * exp(-0.1 * uCurr * F / (R * T)) + 0.0353 * exp(-uCurr * F /(R*T))))

#The Na+/Ca2+ exchanger current.
iNaCaTop = exp(gamma * uCurr * F / (R * T)) * NaIC**3 * CaO  - exp((gamma - 1) * uCurr * F /(R * T)) * NaO**3 * CaIC * alpha
iNaCaBottom = (kmNai**3 + NaO**3) * (kmCa + CaO) * (1 + kSat * exp((gamma - 1) * uCurr * F / (R * T)))
iNaCa = iNaCaTop / iNaCaBottom

#The current IpCa.
ipCa = gpCa * (CaIC / (kpCa + CaIC))

#The current iPk.
iPk = gpK * ((uCurr - eK) / (1 + exp((25-uCurr)/5.98)))

#Background Currents.
ibNa = gbNa * (uCurr - eNa)
ibCa = gbCa * (uCurr - eCa)

#Now for the newer currents.
#L-Type Ca2+ Current.
dInf = 1 / (1 + exp((-8-uCurr)/7.5)) 
alphaD = 1.4 / (1 + exp((-35 - uCurr)/13)) + 0.25
betaD = 1.4 / (1 + exp((uCurr + 5)/5))
gammaD = 1 / (1 + exp((50 - uCurr)/20))
tauD = alphaD * betaD + gammaD

fInf = 1 / (1 + exp((uCurr + 20)/ 7))
alphaF = 1102.5 * exp(-((uCurr+27)/15)**2)
betaF = 200 / (1 + exp((13 - uCurr)/10))
gammaF = 180 / (1 + exp((uCurr+30)/10)) + 20
tauF = alphaF + betaF + gammaF

f2Inf = 0.67 / (1 + exp((uCurr + 35)/ 7)) + 0.33
alphaF2 = 600 * exp(-(uCurr+25)**2/170)
betaF2 = 31 / (1 + exp((25-uCurr)/10))
gammaF2 = 16 / (1 + exp((uCurr+30)/10))
tauF2 = alphaF2 + betaF2 + gammaF2

fCassInf = 0.6 / (1 + (CaSS/0.05)**2) + 0.4
taufCass = 80 / (1 + (CaSS / 0.05)**2) + 2

iCaL = 4 * gCaL * d * f * f2 * fCass * (uCurr-15)*F**2/(R * T) * (0.25 * CaSS * exp(2 * (uCurr - 15) * F / (R * T)) - CaO) / (exp(2 * (uCurr - 15) * F / (R * T)) - 1)

#Slow delayed rectifier current.
xSInf = 1 / (1 + exp((-5-uCurr)/14))
alphaxS = 1400 / sqrt(1 + exp((5-uCurr)/6))
betaxS = 1 / (1 + exp((uCurr - 35)/15))
tauxS = alphaxS * betaxS + 80

iKs = gKs * xSInf**2 * (uCurr - eKs)


#Calcium Dynamics

kCasr = maxSR - (maxSR - minSR) / (1 + (EC / CaSR)**2)
k1 = k1prime / kCasr
k2 = k2prime / kCasr
O = (k1 * CaSS**2 * rBar) / (k3 + k1 * CaSS**2)
iLeak = vLeak * (CaSR - CaIC)
iUp = vMaxup / (1 + kUp**2/CaIC**2)
iRel = vRel * O * (CaSR - CaSS)
iXfer = vXfer * (CaSS - CaIC)

CaIBUFC = CaIC * BUFC / (CaIC + kBUFC)
CaSRBUFSR = CaSR * BUFSR / (CaSR + kBUFSR)
CaSSBUFSS = CaSS * BUFSS / (CaSS + kBUFSS)

O = (k1 * CaSS**2 * rBar) / (k3 + k1 * CaSS**2)


###These may need to change for benchmark testing.
iStim = Constant(0.0)
iAx = Constant(0.0)
iIon = iNa + iK1 + iTo + iKr + iKs + iCaL + iNaCa + iNaK + ipCa + iPk + ibCa + ibNa

Fu = inner(chi * capacitance * Dt(uCurr), vu)*dx + inner(grad(uCurr), sigma * grad(vu))*dx + inner(chi * iIon, vu)*dx
FNaIC = inner(Dt(NaIC), vnaic)*dx - inner((NaIC + ibNa + 3 * iNaK + 3 * iNaCa) / (vC * F), vnaic)*dx
FCaITotal = inner(Dt(CaICT), vcaic)*dx + inner((ibCa + ipCa - 2 * iNaCa) / (2 * vC * F), vcaic)*dx - inner(vSR / vC * (iLeak - iUp) + iXfer, vcaic)*dx
FKIC = inner(Dt(KIC), vkic)*dx + inner((iK1 + iTo + iKr + iKs - 2 * iNaK + iPk + iStim - iAx) / (vC * F), vkic)*dx
FCaSRT = inner(Dt(CaSRT), vcasr)*dx - inner(iUp - iLeak - iRel, vcasr)*dx
FCaSST = inner(Dt(CaSST), vcass)*dx + inner((iCaL / (2 * vSS * F) + (vSR / vSS) * iRel - (vC / vSS) * iXfer), vcass)*dx
FrBar = inner(Dt(rBar), vrbar)*dx - inner(-k2 * CaSS * rBar + k4 * (1 - rBar), vrbar)*dx
FxR1 = inner(Dt(xR1), vxr1)*dx - inner(alphaxR1 * (1 - xR1) + betaxR1 * xR1, vxr1)*dx
FxR2 = inner(Dt(xR2), vxr2)*dx - inner(alphaxR2 * (1 - xR2) + betaxR2 * xR2, vxr2)*dx
FxS = inner(Dt(xS), vxs)*dx - inner(alphaxS * (1 - xS) + betaxS * xS, vxs)*dx
Fm = inner(Dt(m), vm)*dx - inner(alphaM * (1 - m) + betaM * m, vm)*dx
Fh = inner(Dt(h), vh)*dx - inner(alphaH * (1 - h) + betaH * h, vh)*dx
Fj = inner(Dt(j), vj)*dx - inner(alphaJ * (1 - j) + betaJ * j, vj)*dx
Fd = inner(Dt(d), vd)*dx - inner(alphaD * (1 - d) + betaD * d, vd)*dx
Ff = inner(Dt(f), vf)*dx - inner(alphaF * (1 - f) + betaF * f, vf)*dx
Ff2 = inner(Dt(f2), vf2)*dx - inner(alphaF2 * (1 - f2) + betaF2 * f2, vf2)*dx
FfCass = inner(Dt(fCass), vfcass)*dx - inner((fCassInf - fCass) / taufCass, vfcass)*dx
Fr = inner(Dt(r), vr)*dx - inner((rInf - r) / tauR, vr)*dx
Fs = inner(Dt(s), vs)*dx - inner((sInf - s) / tauS, vs)*dx

F = Fu + FNaIC + FCaITotal + FKIC + FCaSRT + FCaSST + FrBar + FxR1 + FxR2 + FxS + Fm + Fh + Fj + Fd + Ff + Ff2 + FfCass + Fr + Fs

params = {"mat_type": "aij",
           "ksp_type": "cg", #gmres or cg preonly
           "pc_type": "lu", #gamg, lu
           #"snes_lag_jacobian": -2,
            #"snes_lag_preconditioner": -2,
           "snes_atol" : 1e-10,
           #"snes_monitor" : None
           }

stepper = TimeStepper(F, butcher_tableau, t, dt, uu,
                      solver_parameters=params)

uFinal, NaICF, CaICF, KICF, CaSRF, CaSSF, rBarF, xR1F, xR2F, xSF, mF, hF, jF, dF, fF, f2F, fCassF, sF, rF = uu.split()
outfile.write(uFinal, time=0)

for j in range(1):
    stepper.advance()
    t.assign(float(t) + float(dt))
    if (j % 5  == 0):
        print("Time step", j)
        outfile.write(uFinal, time=j * float(dt))


#     #Update parameter values for the cell model for the next iteration.
#     CaIC = CaICT- BUFC
#     CaSS = CaSST - BUFSS
#     CaSR = CaSRT - BUFSR


#     #First create iNa.
#     eCa = R * T / (zCa * F) * log(CaO / CaIC)
#     eNa = R * T / (zNa * F) * log(NaO/NaIC)
#     eK = R * T / (zK * F) * log(kO/KIC)
#     eKs = R * T / F * log((kO + pKNa * NaO) / (KIC + pKNa * NaIC))

#     mInf = 1 / (1 + exp((-56.86 - uCurr)/9.03))**2
#     hInf = 1 / (1 + exp( (uCurr + 71.55)/7.43))**2
#     jInf = 1 / (1 + exp( (uCurr + 71.55)/7.43))**2

#     alphaM = 1 / (1 + exp((-60-V)/5))
#     betaM = 0.1 / (1 + exp((V+35)/5)) + 0.1 / (1 + exp((V-50)/200))
#     tauM = alphaM * betaM

#     alphaHLT =  0.057 * exp(-(uCurr + 80)/6.8)
#     alphaHGT = 0
#     betaHGT =  0.77 / (0.13 * (1 + exp(-(uCurr + 10.66)/11.1)))
#     betaHLT = 2.7 *  exp(0.079 * uCurr) + 3.1e5 * exp(0.3485 * uCurr)
#     alphaH = Function(V).assign(conditional(uCurr >= -40, alphaHGT, alphaHLT))
#     betaH = Function(V).assign(conditional(uCurr >= -40, betaHGT, betaHLT))
#     tauH = 1 / (alphaH + betaH)

#     alphaJGT = 0
#     alphaJLT = (uCurr + 37.78) * (-2.5428e4*exp(0.2444*uCurr) - 6.948e-6*exp(-0.04391*uCurr)) / (1 + exp(0.311 * (uCurr + 79.23)))
#     betaJGT = (0.6 * exp(0.057 * uCurr)) / (1 + exp(-0.1*(uCurr+32)))
#     betaJLT = (0.02424 * exp(-0.01052 * uCurr)) / (1 + exp(-0.1378 * (uCurr + 40.14)))
#     alphaJ = Function(V).assign(conditional(uCurr >= - 40, alphaJGT, alphaJLT))
#     betaJ = Function(V).assign(conditional(uCurr >= - 40, betaJGT, betaJLT))
#     tauJ = 1 / (alphaJ + betaJ)

#     iNa = gNa * mInf**3 * hInf * jInf * (uCurr - eNa)
#     #iNa = gNa * m**3 * h * j * (uCurr - eNa)

#     #The transient outward current, iTo.
#     rInf = 1 / (1 + exp((20 - uCurr)/6))
#     sInf = 1 / (1 + exp((uCurr + 20)/5))
#     tauR = 9.5 * exp(-(uCurr + 40)**2/1800) + 0.8
#     tauS = 85 * exp(-(uCurr+45)**2 / 320) + 5 / (1 + exp((uCurr - 20) / 5)) + 3

#     iTo = gTo * rInf * sInf * (uCurr - eK)

#     #Rapid delayed rectifier current
#     xR1Inf = 1 / (1 + exp((-26-uCurr)/7))
#     xR2Inf = 1 / (1 + exp((uCurr+88)/24))

#     alphaxR1 = 450 / (1 + exp((-45-uCurr)/10))
#     betaxR1 = 6 / (1 + exp((uCurr+30)/11.5))
#     tauxR1 = alphaxR1 * betaxR1

#     alphaxR2 = 3 / (1 + exp((-60-uCurr)/20))
#     betaxR2 = 1.12 / (1 + exp((uCurr - 60)/20))
#     tauxR2 = alphaxR2 * betaxR2

#     iKr = gKr * sqrt(kO / 5.4) * xR1Inf * xR2Inf * (uCurr - eK)

#     #The inward rectifier current.
#     alphaK1 = 0.1 / (1 + exp(0.06 * (uCurr - eK - 200)))
#     betaK1 = (3 * exp(0.0002 * (uCurr - eK +100)) + exp(0.1 * (uCurr - eK -10))) / (1 + exp(-0.5 * (uCurr - eK)))
#     xK1Inf = alphaK1 / (alphaK1 + betaK1)

#     iK1 = gK1 * sqrt(kO / 5.4) * xK1Inf * (uCurr - eK)

#     #The Na+/Ca2+ pump current.
#     iNaK = pNaK * kO * NaIC / ((kO + kmK) * (NaIC + kmNa) * (1 + 0.1245 * exp(-0.1 * uCurr * F / (R * T)) + 0.0353 * exp(-uCurr * F /(R*T))))

#     #The Na+/Ca2+ exchanger current.
#     iNaCaTop = exp(gamma * uCurr * F / (R * T)) * NaIC**3 * CaO  - exp((gamma - 1) * uCurr * F /(R * T)) * NaO**3 * CaIC * alpha
#     iNaCaBottom = (kmNai**3 + NaO**3) * (kmCa + CaO) * (1 + kSat * exp((gamma - 1) * uCurr * F / (R * T)))
#     iNaCa = iNaCaTop / iNaCaBottom

#     #The current IpCa.
#     ipCa = gpCa * (CaIC / (kpCa + CaIC))

#     #The current iPk.
#     iPk = gpK * ((uCurr - eK) / (1 + exp((25-V)/5.98)))

#     #Background Currents.
#     ibNa = gbNa * (uCurr - eNa)
#     ibCa = gbCa * (uCurr - eCa)

#     #Now for the newer currents.
#     #L-Type Ca2+ Current.
#     dInf = 1 / (1 + exp((-8-uCurr)/7.5)) 
#     alphaD = 1.4 / (1 + exp((-35 - uCurr)/13)) + 0.25
#     betaD = 1.4 / (1 + exp((uCurr + 5)/5))
#     gammaD = 1 / (1 + exp((50 - uCurr)/20))
#     tauD = alphaD * betaD + gammaD

#     fInf = 1 / (1 + exp((uCurr + 20)/ 7))
#     alphaF = 1102.5 * exp(-((uCurr+27)/15)**2)
#     betaF = 200 / (1 + exp((13 - uCurr)/10))
#     gammaF = 180 / (1 + exp((uCurr+30)/10)) + 20
#     tauF = alphaF + betaF + gammaF

#     f2Inf = 0.67 / (1 + exp((V + 35)/ 7)) + 0.33
#     alphaF2 = 600 * exp(-(uCurr+25)**2/170)
#     betaF2 = 31 / (1 + exp((25-uCurr)/10))
#     gammaF2 = 16 / (1 + exp((uCurr+30)/10))
#     tauF2 = alphaF2 + betaF2 + gammaF2

#     fCassInf = 0.6 / (1 + (CaSS/0.05)**2) + 0.4
#     taufCass = 80 / (1 + (CaSS / 0.05)**2) + 2

#     iCaL = 4 * gCaL * dInf * fInf * f2Inf * fCassInf * (V-15)*F**2/(R * T) * (0.25 * CaSS * exp(2 * (uCurr - 15) * F / (R * T)) - CaO) / (exp(2 * (uCurr - 15) * F / (R * T)) - 1)
#     #iCaL = 4 * g * d * f * f2 * fCass * (V-15)*F**2/(R * T) * (0.25 * CaSS * exp(2 * (uCurr - 15) * F / (R * T)) - CaO) / (exp(2 * (uCurr - 15) * F / (R * T)) - 1)

#     #Slow delayed rectifier current.
#     xSInf = 1 / (1 + exp((-5-uCurr)/14))
#     alphaxS = 1400 / sqrt(1 + exp((5-uCurr)/6))
#     betaxS = 1 / (1 + exp((uCurr - 35)/15))
#     tauxS = alphaxS * betaxS + 80

#     iKs = gKs * xSInf**2 * (uCurr - eKs)


#     #Calcium Dynamics
#     iLeak = vLeak * (CaSR - CaIC)
#     iUp = vMaxup / (1 + kUp**2/CaIC**2)
#     iRel = vRel * O * (CaSR - CaSS)
#     iXfer = vXfer * (CaSS - CaIC)

#     kCasr = maxSR - (maxSR - minSR) / (1 + (EC / CaSR)**2)
#     k1 = k1prime / kCasr
#     k2 = k2prime / kCasr

#     CaIBUFC = CaIC * BUFC / (CaIC + kBUFC)
#     CaSRBUFSR = CaSR * BUFSR / (CaSR + kBUFSR)
#     CaSSBUFSS = CaSS * BUFSS / (CaSS + kBUFSS)

#     O = (k1 * CaSS**2 * rBar) / (k3 + k1 * CaSS**2)




