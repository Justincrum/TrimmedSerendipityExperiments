#Solving monodomain equation in Firedrake
#Justin Crum
#4/16/21

#To do:
#1.  Implement iIon as in Ten Tusscher
#3.  Need to set up initial conditions including what the activiation in
#    a corner will look like.

from firedrake import *
import argparse

parser = argparse.ArgumentParser(
        description="Allows for input of order and mesh refinement.")
parser.add_argument("-O", "--Order", 
        type=int, help="Input the order of the polynomials.")
args = parser.parse_args()

#Create an extruded mesh with geometry 20x7x3, with cube elements of size 0.1x0.1x0.1
msh = RectangleMesh(200, 70, 20, 7, quadrilateral=True)
mesh = ExtrudedMesh(msh, layers=30, layer_height=0.1)
polyOrder = args.Order

#Set up the function space and test/trial functions.
h1Space = FunctionSpace(mesh, "Q", polyOrder)
uCurr = Function(h1Space, name="potential")

u = TrialFunction(h1Space)
v = TestFunction(h1Space)

uCurr.assign(-85.23)

#Set the initial conditions.
x, y, z = SpatialCoordinate(mesh)

#Set up an output file.  I think I need to set up initial conditions before this
#if I want to use to IC's other than 0.
outfile = File("monodomain.pvd")
outfile.write(uCurr)

#Model parameters, taken from the literature.
#A = Constant(0.04)
#vRest = Constant(-85.23)
#vPeak = Constant(40.)
#vTh = Constant(-65.)
chi = Constant(140.)    
#iIon = A**2 * (uCurr - vRest) * (uCurr - vTh) * (uCurr - vPeak)
capacitance = Constant(0.01)
sigma =  as_matrix([[0.133, 0.0, 0.0], [0.0, 0.0176, 0.0], [0.0, 0.0, 0.0176]])



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

#Initial conditions!
NaIC = Function(h1Space, name="NaIC").assign(8.604)
CaIC = Function(h1Space, name="CaIC").assign(0.000126)
KIC = Function(h1Space, name="KIC").assign(136.89)
CaSR = Function(h1Space).assign(3.64)
CaSS = Function(h1Space).assign(0.00036)
rBar = Function(h1Space).assign(0.9073)
xR1 = Function(h1Space).assign(0.00621)
xR2 = Function(h1Space).assign(0.4712)
xS = Function(h1Space).assign(0.0095)
m = Function(h1Space).assign(0.00172)
h = Function(h1Space).assign(0.7444)
j = Function(h1Space).assign(0.7045)
d = Function(h1Space).assign(3.375e-5)
f = Function(h1Space).assign(0.7888)
f2 = Function(h1Space).assign(0.9755)
fCass = Function(h1Space).assign(0.9953)
s = Function(h1Space).assign(0.999998)
r = Function(h1Space).assign(2.42e-8)

iLeak = vLeak * (CaSR - CaIC) 
iUp = vMaxup / (1 + kUp**2/CaIC**2)
iRel = vRel * O * (CaSR - CaSS)
iXfer = vXfer * (CaSS - CaIC)

eCa = R * T / (z * F) * np.log(CaO / CaIC)
eNa = R * T / (z * F) * np.log(NaO/NaIC)
eK = R * T / (z * F) * np.log(kO/KIC)
eKs = R * T / F * np.log((kO + pKNa * NaO) / (KIC + pKNa * NaIC))

mInf = 1 / (1 + e**((-56.86 - uCurr)/9.03))**2
hInf = 1 / (1 + e**( (uCurr + 71.55)/7.43))**2
jInf = 1 / (1 + e**( (uCurr + 71.55)/7.43))**2

alphaM = 1 / (1 + e**((-60-uCurr)/5))
betaM = 0.1 / (1 + e**((uCurr+35)/5)) + 0.1 / (1 + e**((uCurr-50)/200))
tauM = alphaM * betaM

alphaHLT =  0.057 * e**(-(uCurr + 80)/6.8)
alphaHGT = 0
betaHGT =  0.77 / (0.13 * (1 + e**(-(uCurr + 10.66)/11.1)))
betaHLT = 2.7 * e ** (0.079 * uCurr) + 3.1e5 * e**(0.3485 * uCurr)
alphaH = Function(h1Space).assign(conditional(uCurr >= -40, alphaHGT, alphaHLT))
betaH = Function(h1Space).assign(conditional(uCurr >= -40, betaHGT, betaHLT))
tauH = 1 / (alphaH + betaH)

alphaJGT = 0
alphaJLT = (uCurr + 37.78) * (-2.5428e4*e**(0.2444*uCurr) - 6.948e-6*e**(-0.04391*uCurr)) / (1 + e**(0.311 * (uCurr + 79.23)))
betaJGT = (0.6 * e**(0.057 * uCurr)) / (1 + e**(-0.1*(uCurr+32)))
betaJLT = (0.02424 * e**(-0.01052 * uCurr)) / (1 + e**(-0.1378 * (uCurr + 40.14)))
alphaJ = Function(h1Space).assign(conditional(uCurr >= - 40, alphaJGT, alphaJLT))
betaJ = Function(h1Space).assign(conditional(uCurr >= - 40, betaJGT, betaJLT))
tauJ = 1 / (alphaJ + betaJ)

#iNa = gNa * mInf**3 * hInf * jInf * (uCurr - eNa)
iNa = gNa * m**3 * h * j * (uCurr - eNa)

#The transient outward current, iTo.
rInf = 1 / (1 + e**((20 - uCurr)/6))
sInf = 1 / (1 + e**((uCurr + 20)/5))
tauR = 9.5 * e**(-(uCurr + 40)**2/1800) + 0.8
tauS = 85 * e**(-(uCurr+45)**2 / 320) + 5 / (1 + e**((uCurr - 20) / 5)) + 3

iTo = gTo * rInf * sInf * (uCurr - eK)

#Rapid delayed rectifier current
xR1Inf = 1 / (1 + e**((-26-uCurr)/7))
xR2Inf = 1 / (1 + e**((uCurr+88)/24))

alphaxR1 = 450 / (1 + e**((-45-uCurr)/10))
betaxR1 = 6 / (1 + e**((uCurr+30)/11.5))
tauxR1 = alphaxR1 * betaxR1

alphaxR2 = 3 / (1 + e**((-60-uCurr)/20))
betaxR2 = 1.12 / (1 + e**((uCurr - 60)/20))
tauxR2 = alphaxR2 * betaxR2

iKr = gKr * np.sqrt(kO / 5.4) * xR1Inf * xR2Inf * (uCurr - eK)

#The inward rectifier current.
alphaK1 = 0.1 / (1 + e**(0.06 * (uCurr - eK - 200)))
betaK1 = (3 * e**(0.0002 * (uCurr - eK +100)) + e**(0.1 * (uCurr - eK -10))) / (1 + e**(-0.5 * (uCurr - eK)))
xK1Inf = alphaK1 / (alphaK1 + betaK1)

iK1 = gK1 * np.sqrt(kO / 5.4) * xK1Inf * (uCurr - eK)

#The Na+/Ca2+ pump current.
iNaK = pNaK * kO * NaIC / ((kO + kmK) * (NaIC + kmNa) * (1 + 0.1245 * e**(-0.1 * uCurr * F / (R * T)) + 0.0353 * e**(-uCurr * F /(R*T))))

#The Na+/Ca2+ exchanger current.
iNaCaTop = e**(gamma * uCurr * F / (R * T)) * NaIC**3 * CaO  - e**((gamma - 1) * uCurr * F /(R * T)) * NaO**3 * CaIC * alpha
iNaCaBottom = (kmNai**3 + NaO**3) * (kmCa + CaO) * (1 + kSat * e**((gamma - 1) * uCurr * F / (R * T)))
iNaCa = iNaCaTop / iNaCaBottom

#The current IpCa.
ipCa = gpCa * (CaIC / (kpCa + CaIC))

#The current ipK.
ipK = gpK * ((uCurr - eK) / (1 + e**((25-V)/5.98)))

#Background Currents.
ibNa = gbNa * (uCurr - eNa)
ibCa = gbCa * (uCurr - eCa)

#Now for the newer currents.
#L-Type Ca2+ Current.
dInf = 1 / (1 + e**((-8-uCurr)/7.5)) 
alphaD = 1.4 / (1 e**((-35 - uCurr)/13)) + 0.25
betaD = 1.4 / (1 e**((uCurr + 5)/5))
gammaD = 1 / (1 e**((50 - uCurr)/20))
tauD = alphaD * betaD + gammaD

fInf = 1 / (1 + e**((uCurr + 20)/ 7))
alphaF = 1102.5 * e**(-((uCurr+27)/15)**2)
betaF = 200 / (1 + e**((13 - uCurr)/10))
gammaF = 180 / (1 + e**((uCurr+30)/10)) + 20
tauF = alphaF + betaF + gammaF

f2Inf = 0.67 / (1 + e**((V + 35)/ 7)) + 0.33
alphaF2 = 600 * e**(-(uCurr+25)**2/170)
betaF2 = 31 / (1 + e**((25-uCurr)/10))
gammaF2 = 16 / (1 + e**((uCurr+30)/10))
tauF2 = alphaF2 + betaF2 + gammaF2

fCassInf = 0.6 / (1 + (CaSS/0.05)**2) + 0.4
taufCass = 80 / (1 + (CaSS / 0.05)**2) + 2

iCaL = 4 * g * d * f * f2 * fCass * (V-15)*F**2/(R * T) * (0.25 * CaSS * e**(2 * (uCurr - 15) * F / (R * T)) - CaO) / (e**(2 * (uCurr - 15) * F / (R * T)) - 1)

#Slow delayed rectifier current.
xSInf = 1 / (1 + e**((-5-uCurr)/14))
alphaxS = 1400 / np.sqrt(1 + e**((5-uCurr)/6))
betaxS = 1 / (1 + e**((uCurr - 35)/15))
tauxS = alphaxS * betaxS + 80

iKs = gKs * xSInf**2 * (uCurr - eKs)


#Calcium Dynamics
iLeak = vLeak * (CaSR - CaIC)
iUp = vMaxup / (1 + kUp**2/CaIC**2)
iRel = vRel * O * (CaSR - CaSS)
iXfer = vXfer * (CaSS - CaIC)

kCasr = maxSR - (maxSR - minSR) / (1 + (EC / CaSR)**2)
k1 = k1prime / kCasr
k2 = k2prime / kCasr

CaIBUFC = CaIC * BUFC / (CaIC + kBUFC)
CaSRBUFSR = CaSR * BUFSR / (CaSR + kBUFSR)
CaSSBUFSS = CaSS * BUFSS / (CaSS + kBUFSS)

O = (k1 * CaSS**2 * rBar) / (k3 + k1 * CaSS**2)

#Time step for dU/dT -> dUt
tFinal = 2.0
dt = 0.05
t = 0.0
step = 0

while t <= tFinal:
        step += 1

        #Compute the ODEs to get new values for all the concentrations.
        #These will be used in determining the new values of all of the variables below.
        #After that we can then compute the new values of the currents, and from there
        #we can compute the function I_ion(u).
        NaIC = dt * (-1) * (NaIC + ibNa + 3 * iNaK + 3 * iNaCa) / (vC * F) + NaIC 
        KIC = dt * (-1) * (iK1 + iTo + iKr + iKs - 2 * iNaK + iPk + iStim - iAx) / (vC * F) + KIC
        rBar = dt * (-k2 * CaSS * rBar + k4 * (1 - rBar)) + rBar
        CalciumITotal = dt * ((-1) * (ibCa + ipCa - 2 * iNaCa) / (2 * vC * F) + vSR / vC * (iLeak - iUp) + iXfer) + CalciumITotal
        CalciumSR = dt * (iUp - iLeak - iRel) + CalciumSR
        CalciumSS = - dt * (iCaL / (2 * vSS * F) + (vSR / vSS) * iRel - (vC / vSS) * iXfer) + CalciumSS
        d = dt * (alphaD * (1 - d) + betaD * (d)) + d
        f = dt * (alphaF * (1 - f) + betaF * (f)) + f
        f2 = dt * (alphaF2 * (1 - f) + betaF2 * (f2)) + f2
        fCass = dt * (fCassInf - fCass) / taufCass + fCass 
        xS = dt * (alphaxS * (1 - xS) + betaxS * xS) + xS
        m = dt * (alphaM * (1 - m) + betaM * m) + m
        h = dt * (alphaH * (1 - h) + betaH * h) + H
        j = dt * (alphaJ * (1 - j) + betaJ * j) + j
        r = dt * (rInf - r) / tauR + r
        s = dt * (sInf - s) / tauS + s
        xR1 = dt * (alphaxR1 * (1 - xR1) + betxR1 * xR1) + xR1
        xR2 = dt * (alphaxR2 * (1 - xR2) + betaxR2 * xR2) + xR2

        #Update calcium concentrations for next iteration.
        CaIC = CalciumITotal - BUFC
        CaSS = CalciumSS - BUFSS
        CaSR = CalciumSR - BUFSR


        #First create iNa.  TO DO:  FIGURE OUT HOW TO WRITE H AND J FUNCTIONS.
        eCa = R * T / (z * F) * np.log(CaO / CaIC)
        eNa = R * T / (z * F) * np.log(NaO/NaIC)
        eK = R * T / (z * F) * np.log(kO/KIC)
        eKs = R * T / F * np.log((kO + pKNa * NaO) / (KIC + pKNa * NaIC))

        mInf = 1 / (1 + e**((-56.86 - uCurr)/9.03))**2
        hInf = 1 / (1 + e**( (uCurr + 71.55)/7.43))**2
        jInf = 1 / (1 + e**( (uCurr + 71.55)/7.43))**2

        alphaM = 1 / (1 + e**((-60-V)/5))
        betaM = 0.1 / (1 + e**((V+35)/5)) + 0.1 / (1 + e**((V-50)/200))
        tauM = alphaM * betaM

        alphaHLT =  0.057 * e**(-(uCurr + 80)/6.8)
        alphaHGT = 0
        betaHGT =  0.77 / (0.13 * (1 + e**(-(uCurr + 10.66)/11.1)))
        betaHLT = 2.7 * e ** (0.079 * uCurr) + 3.1e5 * e**(0.3485 * uCurr)
        alphaH = Function(h1Space).assign(conditional(uCurr >= -40, alphaHGT, alphaHLT))
        betaH = Function(h1Space).assign(conditional(uCurr >= -40, betaHGT, betaHLT))
        tauH = 1 / (alphaH + betaH)

        alphaJGT = 0
        alphaJLT = (uCurr + 37.78) * (-2.5428e4*e**(0.2444*uCurr) - 6.948e-6*e**(-0.04391*uCurr)) / (1 + e**(0.311 * (uCurr + 79.23)))
        betaJGT = (0.6 * e**(0.057 * uCurr)) / (1 + e**(-0.1*(uCurr+32)))
        betaJLT = (0.02424 * e**(-0.01052 * uCurr)) / (1 + e**(-0.1378 * (uCurr + 40.14)))
        alphaJ = Function(h1Space).assign(conditional(uCurr >= - 40, alphaJGT, alphaJLT))
        betaJ = Function(h1Space).assign(conditional(uCurr >= - 40, betaJGT, betaJLT))
        tauJ = 1 / (alphaJ + betaJ)

        iNa = gNa * mInf**3 * hInf * jInf * (uCurr - eNa)
        #iNa = gNa * m**3 * h * j * (uCurr - eNa)

        #The transient outward current, iTo.
        rInf = 1 / (1 + e**((20 - uCurr)/6))
        sInf = 1 / (1 + e**((uCurr + 20)/5))
        tauR = 9.5 * e**(-(uCurr + 40)**2/1800) + 0.8
        tauS = 85 * e**(-(uCurr+45)**2 / 320) + 5 / (1 + e**((uCurr - 20) / 5)) + 3

        iTo = gTo * rInf * sInf * (uCurr - eK)
        #iTo = gTo * r * s * (uCurr - eK)

        #Rapid delayed rectifier current
        xR1Inf = 1 / (1 + e**((-26-uCurr)/7))
        xR2Inf = 1 / (1 + e**((uCurr+88)/24))

        alphaxR1 = 450 / (1 + e**((-45-uCurr)/10))
        betaxR1 = 6 / (1 + e**((uCurr+30)/11.5))
        tauxR1 = alphaxR1 * betaxR1

        alphaxR2 = 3 / (1 + e**((-60-uCurr)/20))
        betaxR2 = 1.12 / (1 + e**((uCurr - 60)/20))
        tauxR2 = alphaxR2 * betaxR2

        iKr = gKr * np.sqrt(kO / 5.4) * xR1Inf * xR2Inf * (uCurr - eK)

        #The inward rectifier current.
        alphaK1 = 0.1 / (1 + e**(0.06 * (uCurr - eK - 200)))
        betaK1 = (3 * e**(0.0002 * (uCurr - eK +100)) + e**(0.1 * (uCurr - eK -10))) / (1 + e**(-0.5 * (uCurr - eK)))
        xK1Inf = alphaK1 / (alphaK1 + betaK1)

        iK1 = gK1 * np.sqrt(kO / 5.4) * xK1Inf * (uCurr - eK)

        #The Na+/Ca2+ pump current.
        iNaK = pNaK * kO * NaIC / ((kO + kmK) * (NaIC + kmNa) * (1 + 0.1245 * e**(-0.1 * uCurr * F / (R * T)) + 0.0353 * e**(-uCurr * F /(R*T))))

        #The Na+/Ca2+ exchanger current.
        iNaCaTop = e**(gamma * uCurr * F / (R * T)) * NaIC**3 * CaO  - e**((gamma - 1) * uCurr * F /(R * T)) * NaO**3 * CaIC * alpha
        iNaCaBottom = (kmNai**3 + NaO**3) * (kmCa + CaO) * (1 + kSat * e**((gamma - 1) * uCurr * F / (R * T)))
        iNaCa = iNaCaTop / iNaCaBottom

        #The current IpCa.
        ipCa = gpCa * (CaIC / (kpCa + CaIC))

        #The current ipK.
        ipK = gpK * ((uCurr - eK) / (1 + e**((25-V)/5.98)))

        #Background Currents.
        ibNa = gbNa * (uCurr - eNa)
        ibCa = gbCa * (uCurr - eCa)

        #Now for the newer currents.
        #L-Type Ca2+ Current.
        dInf = 1 / (1 + e**((-8-uCurr)/7.5)) 
        alphaD = 1.4 / (1 e**((-35 - uCurr)/13)) + 0.25
        betaD = 1.4 / (1 e**((uCurr + 5)/5))
        gammaD = 1 / (1 e**((50 - uCurr)/20))
        tauD = alphaD * betaD + gammaD

        fInf = 1 / (1 + e**((uCurr + 20)/ 7))
        alphaF = 1102.5 * e**(-((uCurr+27)/15)**2)
        betaF = 200 / (1 + e**((13 - uCurr)/10))
        gammaF = 180 / (1 + e**((uCurr+30)/10)) + 20
        tauF = alphaF + betaF + gammaF

        f2Inf = 0.67 / (1 + e**((V + 35)/ 7)) + 0.33
        alphaF2 = 600 * e**(-(uCurr+25)**2/170)
        betaF2 = 31 / (1 + e**((25-uCurr)/10))
        gammaF2 = 16 / (1 + e**((uCurr+30)/10))
        tauF2 = alphaF2 + betaF2 + gammaF2

        fCassInf = 0.6 / (1 + (CaSS/0.05)**2) + 0.4
        taufCass = 80 / (1 + (CaSS / 0.05)**2) + 2

        iCaL = 4 * gCaL * dInf * fInf * f2Inf * fCassInf * (V-15)*F**2/(R * T) * (0.25 * CaSS * e**(2 * (uCurr - 15) * F / (R * T)) - CaO) / (e**(2 * (uCurr - 15) * F / (R * T)) - 1)
        #iCaL = 4 * g * d * f * f2 * fCass * (V-15)*F**2/(R * T) * (0.25 * CaSS * e**(2 * (uCurr - 15) * F / (R * T)) - CaO) / (e**(2 * (uCurr - 15) * F / (R * T)) - 1)

        #Slow delayed rectifier current.
        xSInf = 1 / (1 + e**((-5-uCurr)/14))
        alphaxS = 1400 / np.sqrt(1 + e**((5-uCurr)/6))
        betaxS = 1 / (1 + e**((uCurr - 35)/15))
        tauxS = alphaxS * betaxS + 80

        iKs = gKs * xSInf**2 * (uCurr - eKs)


        #Calcium Dynamics
        iLeak = vLeak * (CaSR - CaIC)
        iUp = vMaxup / (1 + kUp**2/CaIC**2)
        iRel = vRel * O * (CaSR - CaSS)
        iXfer = vXfer * (CaSS - CaIC)

        kCasr = maxSR - (maxSR - minSR) / (1 + (EC / CaSR)**2)
        k1 = k1prime / kCasr
        k2 = k2prime / kCasr

        CaIBUFC = CaIC * BUFC / (CaIC + kBUFC)
        CaSRBUFSR = CaSR * BUFSR / (CaSR + kBUFSR)
        CaSSBUFSS = CaSS * BUFSS / (CaSS + kBUFSS)

        O = (k1 * CaSS**2 * rBar) / (k3 + k1 * CaSS**2)

        #Problem setup.
        a = (v * u)*dx + (dt / (chi * capacitance) * dot(grad(v), sigma * grad(u)))*dx
        L = (v * uCurr)*dx - (dt / capacitance * iIon * v)*dx
        t += dt

        #No need to set up explicit boundary conditions.  We have already
        #applied the Neumann BC of \grad u \cdot n = 0 by getting rid of 
        #the boundary integral including a \grad u \cdot n factor.

        #Set up solver parameters and solve.
        params = {"snes_type": "newtonls",
                "snes_linesearch_type": "basic",
                "snes_monitor": None,
                "snes_converged_reason": None,
                "mat_type": "aij",
                "snes_max_it": 10,
                "snes_lag_jacobian": -2,
                "snes_lag_preconditioner": -2,
                "ksp_type": "preonly",
                "ksp_converged_reason": None,
                "ksp_monitor_true_residual": None,
                "pc_type": "lu",
                "snes_rtol": 1e-12,
                "snes_atol": 1e-20,
                "pc_factor_mat_solver_type": "mumps",
                "mat_mumps_icntl_14": "200",
                "mat_mumps_icntl_11": "2"}

        solve(a == L, uCurr, solver_parameters=params)
        #t += dt
        if step % 5 == 0:
                outfile.write(uCurr, time=t)



