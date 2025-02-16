from sympy import symbols, Matrix, simplify, I, sqrt, im , conjugate, expand
from sympy.physics.quantum import Dagger
from enum import Enum

import ModelGenerator as ModelGenerator
import argparse


# SM paramters

Cg = symbols('C_g',real=True)
Cgs = symbols('C_gs',real=True)
sigma0 = Matrix([[1,0],[0,1]])
sigma1 = Matrix([[0,1],[1,0]])
sigma2 = Matrix([[0,-I],[I,0]])
sigma3 = Matrix([[1,0],[0,-1]]) 
m_electron = symbols('C_MassElectron',real=True)
m_mu = symbols('C_MassMu',real=True)
m_tau = symbols('C_MassTau',real=True)
m_up = symbols('C_MassUp',real=True)
m_charm = symbols('C_MassCharm',real=True)
m_top = symbols('C_MassTop',real=True)
m_down = symbols('C_MassDown',real=True)
m_strange = symbols('C_MassStrange',real=True)
m_bottom = symbols('C_MassBottom',real=True)



# CKM Matrix
Vud,Vus,Vub,Vcd,Vcs,Vcb,Vtd,Vts,Vtb = symbols('Vud Vus Vub Vcd Vcs Vcb Vtd Vts Vtb')
VCKM = Matrix([[Vud,Vus,Vub],[Vcd,Vcs,Vcb],[Vtd,Vts,Vtb]])

#parameters
musq = symbols('musq',real=True)
lam = symbols('lam',real=True)
a1 = symbols('a1',real=True)
a2 = symbols('a2',real=True)
b2 = symbols('b2',real=True)
b3 = symbols('b3',real=True)
b4 = symbols('b4',real=True)

params=[musq,lam,a1,a2,b2,b3,b4]
# params=[musq,lam,b2,a1,a2,b3,b4] # change the order so the tadpole equations are used to replace musq and b2

#CT params
dmusq = symbols('dmusq',real=True)
dlam = symbols('dlam',real=True)
da1 = symbols('da1',real=True)
da2 = symbols('da2',real=True)
db2 = symbols('db2',real=True)
db3 = symbols('db3',real=True)
db4 = symbols('db4',real=True)

dparams=[dmusq,dlam,da1,da2,db2,db3,db4]
# dparams=[dmusq,dlam,db2,da1,da2,db3,db4] # change the order to correspond with params above

#VEVs
vH = symbols('vH', real=True)
vS = symbols('vS', real=True)
wH = symbols('wH', real=True)
wS = symbols('wS', real=True)

#Higgsfields
rho1,eta1,zeta1,psi1 = symbols('rho1 eta1 zeta1 psi1', real=True)
zetaS = symbols('zetaS', real=True)
# Higgsfields=[rho1,eta1,zeta1,psi1,zetaS]
Higgsfields=[rho1,eta1,psi1,zeta1,zetaS] # order: both zetas are here at the end! To match with my Maple implementation.
CTTadpoles = symbols('dT1:{}'.format(len(Higgsfields)+1),real=True)

#doublet
# phi = Matrix([[Higgsfields[0]+I*Higgsfields[1]], [Higgsfields[2]+I*Higgsfields[3]]]) * 1/sqrt(2)
phi = Matrix([[Higgsfields[0]+I*Higgsfields[1]], [Higgsfields[3]+I*Higgsfields[2]]]) * 1/sqrt(2) # order, see above

#singlet
phiS = Higgsfields[4]

#replacements
# higgsvevAtZeroTemp = [0, 0, vH, 0, vS]
# higgsVEVAtFiniteTemp = [0, 0, wH, 0, wS]
higgsvevAtZeroTemp = [0, 0, 0, vH, vS] # order, see above
higgsVEVAtFiniteTemp = [0, 0, 0, wH, wS] # order, see above
zeroTempVEV = [(Higgsfields[i],higgsvevAtZeroTemp[i]) for i in range(len(Higgsfields)) ]
finiteTempVEV = [(Higgsfields[i],higgsVEVAtFiniteTemp[i]) for i in range(len(Higgsfields)) ]
fieldsZero = [(x,0) for x in Higgsfields]

phiSq = simplify((Dagger(phi)*phi)[0])

VHiggsNC = -musq*phiSq + lam*phiSq**2 + a1/2*phiSq*phiS + a2/2*phiSq*phiS**2 + b2/2*phiS**2 + b3/3*phiS**3 + b4/4*phiS**4
VHiggs = simplify(VHiggsNC)

if im(VHiggs) != 0:
    raise Exception("Higgs potential has an imaginary part with " + str(im(VHiggs)))

# Generate the model
RxSM = ModelGenerator.ModelGenerator(params,dparams,CTTadpoles,Higgsfields,VHiggs,zeroTempVEV, finiteTempVEV)




# Set Gauge fields
W1, W2, W3, B0 = symbols('W1 W2 W3 B0',real=True)


Dmu = -I*Cg/2 * (sigma1*W1 + sigma2 * W2 + sigma3*W3) -I*Cgs/2 * sigma0 * B0
VGauge = simplify(Dagger(Dmu*phi)*(Dmu*phi))[0,0]

# Generate Lepton Potentials
NuL = symbols('veL vmuL vtauL',real=True)
ER = symbols('eR muR tauR', real=True)
EL = symbols('eL muL tauL', real=True)
LepBase = NuL + ER + EL

ye = sqrt(2)*m_electron/vH
ymu = sqrt(2)*m_mu/vH
ytau = sqrt(2)*m_tau/vH

PiLep = Matrix([[ye,0,0],[0,ymu,0],[0,0,ytau]])

VFLep = 0
for i in range(len(NuL)):
    for j in range(len(ER)):
        VFLep += (NuL[i] * PiLep[i,j] * ER[j])*phi[0]

for i in range(len(EL)):
    for j in range(len(ER)):
        VFLep += (EL[i] * PiLep[i,j] * ER[j])*phi[1]

VFLep = simplify(VFLep)


# Generate Quark Potentials
UL = symbols('uL cL tL', real=True)
DL = symbols('dL sL bL', real=True)
UR = symbols('uR cR tR', real=True)
DR = symbols('dR sR bR', real=True)
QuarkBase = UR + DR + UL + DL

yb = sqrt(2)*m_bottom/vH
yc = sqrt(2)*m_charm/vH
yd = sqrt(2)*m_down/vH
ys = sqrt(2)*m_strange/vH
yt = sqrt(2)*m_top/vH
yu = sqrt(2)*m_up/vH

DownCoupling = Matrix([[yd,0,0],[0,ys,0],[0,0,yb]])
UpCoupling = Matrix([[yu,0,0],[0,yc,0],[0,0,yt]])

ULVector = Matrix([[x] for x in UL])
DLVector = Matrix([[x] for x in DL])
URVector = Matrix([[x] for x in UR])
DRVector = Matrix([[x] for x in DR])

VQuark = ULVector.transpose() * VCKM * DownCoupling * DRVector * phi[0]
VQuark+= DLVector.transpose() * DownCoupling * DRVector * phi[1] 

VQuark += ULVector.transpose() * UpCoupling * URVector * phi[1].conjugate()
VQuark+= -DLVector.transpose() * Dagger(VCKM)*UpCoupling*URVector*phi[0].conjugate()
VQuark = simplify(VQuark[0,0])

# Get the tesnors
RxSM.setGauge([W1,W2,W3,B0],VGauge)
RxSM.setLepton(LepBase, VFLep)
RxSM.setQuark(QuarkBase, VQuark)


def setAdditionalCTEquations():
    # additional equations to define a unique CT solution point 
    additionaEquations = []
    # additionaEquations.append(dmusq)
    # additionaEquations.append(db2)
    # additionaEquations.append(db3) # choose b3 as MSbar parameter: CT of b3 is zero
    # additionaEquations.append(db4) # choose b4 as MSbar parameter: CT of b4 is zero
    additionaEquations.append(da2)
    additionaEquations.append(db3)
    # additionaEquations.append(db4)
    additionaEquations.append(CTTadpoles[-1]) # set singlet tadpole CT to zero
    return additionaEquations



parser = argparse.ArgumentParser()
parser.add_argument('-s','--show',choices=['ct','tensor','treeSimpl','CTSimpl'],required=True,help='The part of the model to be printed')

if __name__ == "__main__":
    args = parser.parse_args()
    method = args.show

    printCT = method == 'ct'
    printTensors = method == 'tensor'

    if printCT:
        print("//Begin CT Calculation")
        RxSM.printCTForCPP(setAdditionalCTEquations())
        print("//End CT Calculation")

        print("//Begin CT Order for set_CT_Pot_Par")
        RxSM.printCTOrder()
        print("//End CT Order")

    if printTensors:
        RxSM.printModelToCPP()

    if method == 'treeSimpl':
        RxSM.printTreeSimplified()

    if method == 'CTSimpl':
        RxSM.printVCTSimplified()
