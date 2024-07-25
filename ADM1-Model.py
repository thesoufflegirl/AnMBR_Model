from gekko import GEKKO
import pandas as pd

m = GEKKO(remote=False) # create GEKKO model


######################################################
#Initial Concentrations (in mg/L), These can be altered based on concentration of system. X are the non dissolved fractions, S are dissolved fractions

CODo = 1000

Ssuo = 0 #Dissolved Sugar
Saao = 0 #amino acids
Sfao = 0 #long chain fatty acids
Svao = 0 #valerate
Sbuo = 0.013 #butyrate
Sproo = 0.233 #propionate
Saco = 0.157 #acetate
Sh2o = 0 #hydrogen gas
Sc4o = 0 #methane gas
SICo = 0 #inorganic carbon
SINo = 0.4 #inorganic nitrogen
SIo = 0 #Soluble inerts
Xco = 0 #composites - these are the dead biomass basically - they break down into lipids, proteins, and carbohydrates during disintegration
Xcho = CODo*0.848 #carbohydrates
Xproo = CODo*0.229 #Proteins 
Xlio = CODo*0.07 #lipids
#Bacteria in influent
Xsuo = 0 #sugar degraders
Xaao = 0 #amino acid degradrs
Xfao = 0 #fatty acid degraders
Xc4o = 0 #methane degraders
Xpr_o = 0 #propionate degraders
Xaco = 0 #acetate degraders
Xh2o = 0 #hydrogen degraders
XIo = 0 #Particulate(insoluble) Inerts
Spho = 40 #Phosphorus

VFAon = 0 

#########################################################
#Constants

kdis = 0.5 #Dissociation constant for composites
kch = 10 #dissociation constants for carbohydrates
kpr = 10 #dissociation constant for proteins
kli = 10 #dissociation constant for lipids
kdec = 0.3 #decay constant for bacteria (assumes all degraders have same decay constant)


Q = 1000  #flow
V = 500 #volume (flow and volume can be altered in ratios to affect HRT)
Vliq = 0.8*V #volume of reactor that's liquid
Vgas = 0.2*V #volume of reactor that's gas
T = 298 #temp K
R = 0.08315 #Gas law constant
patm = 1 #pressure in atmospheres


######################################################### ANMBR physics

#Yield From Catabolism (f)
#these are constants in  a matrix e.g. f_Su,aa would be row 1 column 2

d = {'Ssu':[0,0,0,0,0.13,0.27,0.41,0.19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Saa':[0,0,0,0.23,0.26,0.05,0.04,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Sfa':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.95,0,0,0,0,0,0,0,0],
                'Sva':[0,0.23,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Sbu':[0.13,0.26,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Spro':[0.27,0.05,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Sac':[0.41,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Sh2':[0.19,0.06,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Sc4':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'SIC':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'SIN':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'SI':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xc':[0,0,0,0,0,0,0,0,0,0,0,0.1,0,0.2,0.3,0,0,0,0,0,0,0,0,0.1],
                'Xch':[0,0,0,0,0,0,0,0,0,0,0,0,0.2,0,0,0,0,0,0,0,0,0,0,0],
                'Xpr':[0,0,0,0,0,0,0,0,0,0,0,0,0.3,0,0,0,0,0,0,0,0,0,0,0],
                'Xli':[0,0,0.95,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xsu':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xaa':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xfa':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xc4':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xpro':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xac':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'Xh2':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
                'XI':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] }


f = pd.DataFrame(data=d,index = ['Ssu','Saa','Sfa', 'Sva','Sbu','Spro','Sac','Sh2','Sc4','SIC','SIN','SI','Xc','Xch','Xpr','Xli','Xsu','Xaa','Xfa','Xc4','Xpro','Xac','Xh2','XI'])

#Microbial Biomass Yield
Y = {'su':0.1,'aa':0.08,'fa':0.06,'c4':0.06,'pro':0.04,'ac':0.05,'h2':0.06,'bu':0.08}

#y.get to get htings 
#Carbon content 
C = {'Ssu':0.0313,
                'Saa':0.03,
                'Sfa':0.0217,
                'Sva':0.024,
                'Sbu':0.025,
                'Spro':0.0268,
                'Sac':0.0313,
                'Sh2':0,
                'Sc4':0.0156,
                'SIC':1,
                'SIN':0,
                'SI':0.03,
                'Xc':0.02786,
                'Xch':0.0313,
                'Xpr':0.03,
                'Xli':0.022,
                'Xsu':0.0313,
                'Xaa':0.03,
                'Xfa':0.0217,
                'Xc4':0.0156,
                'Xpro':0.0268,
                'Xac':0.0313,
                'XI':0.03,
                'h2':0,
                'c4':0.0156} 

#add in arrhenius constants to monod stuff          
#Nitrogen Content
N = {'xc':0.0027,'I':0.0043,'bac':0.0057,'aa':0.007}
#Monod Constants
km = {'su':30,'aa':50,'fa':6,'c4':20,'pro':13,'ac':8,'h2':35} #max uptake rate for bacteria
#Half saturation constnats for bacteria
KS = {'su':0.5,'aa':0.3,'fa':0.4,'c4':0.2,'pro':0.1,'ac':0.15,'h2':7*10**-6,'IN':0.0001}

#
#Biochemical Rate Coefficients 
#These shouldn't need to be changed - you can change the values yields up above

d = {           'p1':[0,0,0,0,0.13,0.27,0.41,0.19,0,0,0,f.at['SI','Xc'],-1,f.at['Xch','Xc'],f.at['Xpr','Xc'],f.at['Xli','Xc'],0,0,0,0,0,0,0,f.at['XI','Xc']],
                'p2':[1,0,0,0,0,0.,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],
                'p3':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],
                'p4':[1-f.at['Sfa','Xli'],0,1-f.at['Sfa','Xli'],0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],
                'p5':[-1,0,0,0,(1-Y.get('su'))*f.at['Sbu','Ssu'],(1-Y.get('su'))*f.at['Spro','Ssu'],(1-Y.get('su'))*f.at['Sac','Ssu'],(1-Y.get('su'))*f.at['Sh2','Ssu'],
                        0,-C.get('Ssu')*Y.get('su') ,-Y.get('su')*N.get('bac'),0,0,0,0,0,Y.get('su'),0,0,0,0,0,0,0],
                'p6':[0,-1,0,(1-Y.get('aa'))*f.at['Sva','Sfa'],(1-Y.get('aa'))*f.at['Sbu','Sfa'],(1-Y.get('aa'))*f.at['Spro','Sfa'],(1-Y.get('aa'))*f.at['Sac','Sfa'],(1-Y.get('aa'))*f.at['Sh2','Sfa'],0,-C.get('Saa')*Y.get('aa'),-Y.get('su')*N.get('bac'),0,0,0,0,0,0,Y.get('aa'),0,0,0,0,0,0],
                'p7':[0,0,-1,0,0,0,(1-Y.get('fa'))*0.7,(1-Y.get('fa'))*0.3,0,0,-Y.get('fa')*N.get('bac'),0,0,0,0,0,0,Y.get('aa'),0,0,0,0,0,0],
                'p8':[0,0,0,-1,0,(1-Y.get('c4'))*0.54,(1-Y.get('c4'))*0.31,(1-Y.get('c4'))*0.15,0,0,-Y.get('c4')*N.get('bac'),0,0,0,0,0,0,0,0,Y.get('c4'),0,0,0,0],
                'p9':[0,0,0,0,-1,0,(1-Y.get('c4'))*0.8,(1-Y.get('c4'))*0.2,0,0,-Y.get('c4')*N.get('bac'),0,0,0,0,0,0,0,0,Y.get('c4'),0,0,0,0],
                'p10':[0,0,0,0,0,-1,(1-Y.get('pro'))*0.57,(1-Y.get('pro'))*0.43,0,-C.get('Spro')*Y.get('pro'),-Y.get('pro')*N.get('bac'),0,0,0,0,0,0,0,0,0,Y.get('pro'),0,0,0],
                'p11':[0,0,0,0,0,0,-1,0,(1-Y.get('ac')),-C.get('Sac')*Y.get('ac'),-Y.get('ac')*N.get('bac'),0,0,0,0,0,0,0,0,0,0,Y.get('ac'),0,0],
                'p12':[0,0,0,0,0,0,0,-1,1-Y.get('h2'),0,-Y.get('h2')*N.get('bac'),0,0,0,0,0,0,0,0,0,0,0,Y.get('h2'),0],
                'p13':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0,0,0,0,0],
                'p14':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1,0,0,0,0,0,0],
                'p15':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0, 0, 0,-1,0,0,0,0,0],
                'p16':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0, 0, 0, 0,-1,0,0,0,0],
                'p17':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0, 0, 0,0,-1,0,0,0],
                'p18':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0, 0, 0,0,0,-1,0,0],
                'p19':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0, 0, 0, 0,0,0,0,-1,0],}

v = pd.DataFrame(data=d,index = ['Ssu','Saa','Sfa', 'Sva','Sbu','Spro','Sac','Sh2','Sc4','SIC','SIN','SI','Xc','Xch','Xpr','Xli','Xsu','Xaa','Xfa','Xc4','Xpro','Xac','Xh2','XI'])
print(v)

a = 0.8



#Create variables for concentrations - If model is not converging - These can be adjusted to have better initial guesses
#Better initial guesses will make model converge in fewer iterations
Ssu = m.Var(Ssuo*a,lb = 0, ub = 10000) # create GEKKO variable
Saa = m.Var(Saao*a,lb = 0, ub = 10000) 
Sfa = m.Var(Sfao*a,lb = 0, ub = 10000)
Sva = m.Var(Svao*a,lb = 0, ub = 10000)
Sbu = m.Var(Sbuo*a,lb = 0, ub = 10000)
Spro = m.Var(Sproo*a,lb = 0, ub = 10000)
Sac = m.Var(Saco,lb = 0, ub = 10000)
Sh2 = m.Var(Sh2o,lb = 0, ub = 10000)
Sc4 = m.Var(Sc4o,lb = 0, ub = 10000)
SIC = m.Var(SICo,lb = 0, ub = 10000) 
SIN = m.Var(SINo,lb = 0, ub = 10000)
SI = m.Var(SIo,lb = 0, ub = 10000) 
Xc = m.Var(0,lb = 0, ub = 10000)
Xch = m.Var(0,lb = 0, ub = 10000)
Xpro = m.Var(0,lb = 0, ub = 10000)
Xli = m.Var(0,lb = 0, ub = 10000)
Xsu = m.Var(0,lb = 0, ub = 10000)
Xaa = m.Var(Xaao,lb = 0, ub = 10000)
Xfa = m.Var(Xfao,lb = 0, ub = 10000)
Xc4 = m.Var(Xc4o,lb = 0, ub = 10000)
Xpr = m.Var(Xpr_o,lb = 0, ub = 10000)
Xac = m.Var(Xaco,lb = 0, ub = 10000)
Xh2 = m.Var(Xh2o,lb = 0, ub = 10000)
XI = m.Var(0,lb = 0, ub = 10000)
Sh = m.Var(10**(-8),0,10)
Shgas = m.Var(0,lb = 0, ub = 10000)
Sc4gas = m.Var(0,lb = 0, ub = 10000)
Sco2gas = m.Var(0,lb = 0, ub = 10000)
ph2 = m.Var(0,lb = 0, ub = 10000)
pc4 = m.Var(0,lb = 0, ub = 10000)
pco2 = m.Var(0,lb = 0, ub = 10000)
pgas = m.Var(0,lb = 0, ub = 10000)
qgas = m.Var(0,lb = 0, ub = 10000)
temp = m.Var(0,lb = 0)
P_anmbr = m.Var(0)

conc = {'Ssu':Ssu,
                'Saa':Saa,
                'Sfa':Sfa,
                'Sva':Sva,
                'Sbu':Sbu,
                'Spro':Spro,
                'Sac':Sac,
                'Sh2':Sh2,
                'Sc4':Sc4,
                'SIC':SIC,
                'SIN':SIN,
                'SI':SI,
                'Xc':Xc,
                'Xch':Xch,
                'Xpr':Xpr,
                'Xli':Xli,
                'Xsu':Xsu,
                'Xaa':Xaa,
                'Xfa':Xfa,
                'Xc4':Xc4,
                'Xpro':Xpro,
                'Xac':Xac,
                'Xh2':Xh2,
                'XI':XI} 
initConc = {'Ssu':Ssu,
                'Saa':Saao,
                'Sfa':Sfao,
                'Sva':Svao,
                'Sbu':Sbuo,
                'Spro':Sproo,
                'Sac':Saco,
                'Sh2':Sh2o,
                'Sc4':Sc4o,
                'SIC':SICo,
                'SIN':SINo,
                'SI':SIo,
                'Xc':Xco,
                'Xch':Xcho,
                'Xpr':Xpro,
                'Xli':Xlio,
                'Xsu':Xsuo,
                'Xaa':Xaao,
                'Xfa':Xfao,
                'Xc4':Xc4o,
                'Xpro':Xproo,
                'Xac':Xaco,
                'Xh2':Xh2o,
                'XI':XIo} 
    


#Inhibition Constants
#For the first pH function, pHUL and pHLL are upper and lower limits where the group of organisms is 50%
ph_UL_aa = 5.5
ph_LL_aa = 4.0
ph_UL_ac = 7.0
ph_LL_ac = 6.0
ph_UL_h2 = 6.0
ph_LL_h2 = 5.0
n_aa = 3/(ph_UL_aa-ph_LL_aa) 
n_ac = 3/(ph_UL_ac-ph_LL_ac) 
n_h2 = 3/(ph_UL_h2-ph_LL_h2)
K_ph_aa = 10 - (ph_UL_aa+ph_LL_aa)/2
K_ph_ac = 10 - (ph_UL_ac+ph_LL_ac)/2
K_ph_h2 = 10 - (ph_UL_h2+ph_LL_h2)/2
K_I_h2_fa = 5*10**-6
K_I_h2_pro = 3.5*10**-6
K_I_nh3 = 0.0018
K_I_h2_c4 = 1*10**-5

#Henry's Law Constants
KhH2 = 0.00072
Khco2 = 0.025
Khc4 = 0.0011
kla = 200
kp = 5*10^4

I1  = {'ph_aa':K_ph_aa**n_aa/(K_ph_aa**n_aa+Sh**n_ac),
       'IN_lim':SIN/(SIN+KS.get('IN')),
       'h2_fa':K_I_h2_fa /(K_I_h2_fa +Sh2),
        'h2_pro':K_I_h2_pro/(K_I_h2_pro+Sh2),
        'ph_ac': K_ph_ac**n_ac/(K_ph_ac**n_ac+Sh**n_ac),
        'h2_c4': K_I_h2_c4/(K_I_h2_c4+Sh2),
        'nh3':K_I_nh3/(K_I_nh3+SIN),
        'ph_h2':K_ph_h2**(n_h2)/(K_ph_h2**(n_h2)+Sh**n_h2)
}

I = {'su':I1.get('ph_aa')*I1.get('IN_lim'),
     'aa':I1.get('ph_aa')*I1.get('IN_lim'),
     'ph_aa':K_ph_aa**n_aa/(K_ph_aa**n_aa+Sh**n_aa),
     'IN_lim':SIN/(SIN+KS.get('IN')),
     'fa':I1.get('ph_aa')*I1.get('IN_lim')*I1.get('h2_fa'),
     'h2_fa':K_I_h2_fa /(K_I_h2_fa +Sh2),
     'c4':I1.get('ph_aa')*I1.get('IN_lim')*I1.get('h2_c4'),
     'h2_c4': K_I_h2_c4/(K_I_h2_c4+Sh2),
     'pro':I1.get('ph_aa')*I1.get('IN_lim')*I1.get('h2_pro'),
     'h2_pro':K_I_h2_pro/(K_I_h2_pro+Sh2),
     'ac':I1.get('ph_ac')*I1.get('IN_lim')*I1.get('nh3'),
     'ph_ac': K_ph_ac**n_ac/(K_ph_ac**n_ac+Sh**n_ac),
     'h2':I1.get('ph_h2')*I1.get('IN_lim'),
     'ph_h2':K_ph_h2**(n_h2)/(K_ph_h2**(n_h2)+Sh**n_h2),
     'nh3':K_I_nh3/(K_I_nh3+SIN)}


#Speed of reactions
rho = {'p1':kdis*Xc,
    'p2':kch*Xch,
    'p3':kpr*Xpr,
    'p4':kli*Xli,
    'p5':km.get('su')*Ssu/(KS.get('su')+Ssu)*Xsu*(1/(1+Ssu/I.get('su'))),
    'p6':km.get('aa')*Saa/(KS.get('aa')+Saa)*Xaa*(1/(1+Saa/I.get('aa'))),
    'p7':km.get('fa')*Sfa/(KS.get('fa')+Sfa)*Xfa*(1/(1+Sfa/I.get('fa'))),
    'p8':km.get('c4')*Sva/(KS.get('c4')+Sva)*Xc4*Sva/(Sva+Sbu)*I.get('c4'),
    'p9':km.get('c4')*Sbu/(KS.get('c4')+Sbu)*Xc4*Sbu/(Sbu+Sva)*I.get('c4'),
    'p10':km.get('pro')*Spro/(KS.get('pro')+Spro)*Xpro*I.get('pro'),
    'p11':km.get('ac')*Sac/(KS.get('ac')+Sac)*Xac*I.get('ac'),
    'p12':km.get('h2')*Sh2/(KS.get('h2')+Sh2)*Xh2*I.get('h2'),
    'p13':kdec*Xsu,
    'p14':kdec*Xaa,
    'p15':kdec*Xfa,
    'p16':kdec*Xc4,
    'p17':kdec*Xpro,
    'p18':kdec*Xac,
    'p19':kdec*Xh2o}



m.Equation(Ssu.dt() == Ssuo+rho.get('p2') + (1-f.at['Sfa','Xli'])*rho.get('p4')-rho.get('p5'))
m.Equation(Saa.dt() == Saao+ rho.get('p3')-rho.get('p6')) #2
m.Equation(Sfa.dt() == Sfao+ (1-f.at['Sfa','Xli'])*rho.get('p4')-rho.get('p7')) #3
m.Equation(Sva.dt() == Svao+ (1-Y.get('aa')) * f.at['Sva','Sfa'] *rho.get('p6') -rho.get('p8')) #4
m.Equation(Sh2.dt() == Sh2o+(1-Y.get('su'))*f.at['Sh2','Ssu']*rho.get('p5') + (1-Y.get('aa'))*f.at['Sh2','Saa']*rho.get('p6')+ (1-Y.get('fa'))*0.3*rho.get('p7')+(1-Y.get('c4'))*0.15*rho.get('p8')+(1-Y.get('c4'))*0.2*rho.get('p9')+(1-Y.get('pro'))*0.43*rho.get('p10')-rho.get('p12'))
m.Equation(Sc4.dt() == Sc4o+(1-Y.get('ac'))*rho.get('p11')+(1-Y.get('h2'))*rho.get('p12'))
m.Equation(SIC.dt() == SICo+((C.get('SIN')*Y.get('su')*N.get('bac')-C.get('Xsu')*Y.get('su')*rho.get('p5') -(-C.get('Saa')+C.get('Sva'))*(1-Y.get('aa'))*f.at['Sva','Saa']+C.get('Sbu')*(1-Y.get('bu'))*f.at['Sbu','Saa']+C.get('Spro')*(1-Y.get('pro'))*f.at['Spro','Saa']+C.get('Sac')*(1-Y.get('ac'))*f.at['Sac','Saa']+C.get('h2')*(1-Y.get('h2'))*f.at['Sh2','Saa']+N.get('aa')*(-Y.get('aa'))*N.get('bac')*C.get('SIN') +Y.get('aa')*C.get('Saa'))*rho.get('p6') -(-C.get('Spro')+(1-Y.get('pro'))*0.57*C.get('Sac') +(1-Y.get('pro'))*0.43*C.get('h2')-Y.get('pro')*N.get('bac')*C.get('SIN') +Y.get('pro')*C.get('Spro'))*rho.get('p10') -(-C.get('Sac')+(1-Y.get('ac')*C.get('Sac')-Y.get('ac')*N.get('bac') )+Y.get('ac')*C.get('Sac'))-(-C.get('h2')+C.get('c4')*(1-Y.get('h2'))-Y.get('h2')*N.get('bac')*C.get('SIN')+Y.get('h2')*C.get('Sh2')))*rho.get('p10'))
m.Equation(SIN.dt() == SINo -Y.get('su')*N.get('bac')*rho.get('p5')+N.get('aa')*(-Y.get('aa'))*N.get('bac')*rho.get('p6')-Y.get('fa')*N.get('bac')*rho.get('p7')-Y.get('c4')*N.get('bac')*rho.get('p8')-Y.get('c4')*N.get('bac')*rho.get('p9')-Y.get('pro')*N.get('bac')*rho.get('p10')-Y.get('ac')*N.get('bac')*rho.get('p11')-Y.get('h2')*N.get('bac')*rho.get('p12'))
m.Equation(SI.dt() == SIo+f.at['SI','Xc']*rho.get('p1'))
m.Equation(Xc.dt() == Xco-rho.get('p1')+rho.get('p13')+rho.get('p14')+rho.get('p15')+rho.get('p16')+rho.get('p17')+rho.get('p18')+rho.get('p19'))
m.Equation(Xpr.dt() == Xpro+f.at['Xpr','Xc']*rho.get('p1')-rho.get('p3'))
m.Equation(Xsu.dt() == Xsuo+Y.get('su')*rho.get('p5')-rho.get('p13'))
m.Equation(Xaa.dt() == Xaao+ Y.get('aa')*rho.get('p6')-rho.get('p14'))
m.Equation(Xfa.dt() == Xfao+Y.get('fa')*rho.get('p7')-rho.get('p15'))
m.Equation(Xc4.dt() == Xc4o+Y.get('c4')*rho.get('p8')+Y.get('c4')*rho.get('p9')-rho.get('p16'))
m.Equation(Xpro.dt() == Xproo+Y.get('pro')*rho.get('p9')-rho.get('p17'))
m.Equation(Xac.dt() == Xaco+Y.get('ac')*rho.get('p10')-rho.get('p18'))
m.Equation(Xh2.dt() == Xh2o+Y.get('h2')*rho.get('p11')-rho.get('p19'))
m.Equation(XI.dt() ==XIo+ f.at['XI','Xc']*rho.get('p1'))


#Gas transfer Equations
m.Equation(Shgas.dt()  == Vliq/Vgas*kla*(Sh2-16*KhH2*ph2))
m.Equation(Sc4gas.dt() == Vliq/Vgas*kla*(Sc4-64*Khc4*pc4))
m.Equation(Sco2gas.dt() == Vliq/Vgas*kla*(Sc4-Khco2*pco2))
m.Equation(pc4 == Sc4gas*R*T/64)
m.Equation(pco2 == Sco2gas*R*T)
m.Equation(ph2 == Shgas*R*T/16)
m.Equation(pgas == pc4+pco2+ph2)
m.Equation(qgas == kp*(pgas-patm)*pgas/patm)


#########################################################################End ANMBR Physics
###########################################Start Degasser Physics
#######Solve Henry's law problem first to find pout
#Assume everything is at steady state

P_degas = m.Var(0,0) # create GEKKO variable
Sc4gas_degas = m.Var(0,0) #effluent concentration coming out
pin_degas  = m.Var(0,0)
pout_degas = m.Var(0,0)
Wc4gas_degas = m.Var(0,0)
T_degas = T #Temperature
#k = kp from above heat capacity ratio, continue calling kp in this section

#Influent concentration = effluent concentration from above, partial pressure may change from temp, so recalculate
m.Equation(pin_degas == Sc4gas*R*T_degas/64 + Sco2gas*R*T_degas*(Shgas*R*T_degas/16))
m.Equation(Sc4gas_degas == kp*(pout_degas-patm)*pout_degas/patm)
m.Equation(Wc4gas_degas == Sc4gas_degas*Q) #molar flow of gas in (assume all gas in from influent)
m.Equation(P_degas == (kp*Wc4gas_degas*R*T_degas)/(kp-1)*((pout_degas/pin_degas)*(kp-1)/kp-1))
#coming out of this equation we have power needed, concentration of effluent methane, molar flow of methane gas
#one degree of freedom, minimized later through LCA


m.Minimize(Ssu+Saa+Sfa+Sva+Sbu+Spro+Sac+SIN)

m.time = [0,10]      # time points if doing non-steady state (these are in hours, adding more time points will make model run slower)
m.options.IMODE = 2 # simulation mode 

m.options.NODES = 2
m.options.CSV_WRITE=2  
m.options.MAX_ITER = 100  #Increase if not converging 
     
                        #   with internal nodes
stop = False

k = 1
while (stop == True or k <10):

    try:
     m.options.SOLVER= 2           
     m.solve(disp= False)     # solve
     m.Equation(Sbu.dt() == Sbuo+(1-Y.get('su'))*f.at['Sbu','Ssu']*rho.get('p5')+(1-Y.get('aa'))*f.at['Sbu','Saa']*rho.get('p6')-rho.get('p9')) #5
     m.Equation(Xch.dt() == Xcho+f.at['Xch','Xc']*rho.get('p1')-rho.get('p2'))
     m.Equation(Xli.dt() == Xlio+f.at['Xli','Xc']*rho.get('p1')-rho.get('p4'))
     m.Equation(Spro.dt() == Sproo+(1-Y.get('su'))*f.at['Spro','Ssu']*rho.get('p5')+(1-Y.get('aa'))*f.at['Spro','Saa']*rho.get('p6')+(1-Y.get('c4'))*0.54*rho.get('p8'))
     m.Equation(Sac.dt() == Saco+(1-Y.get('su'))*f.at['Sac','Ssu']*rho.get('p5') + (1-Y.get('aa'))*f.at['Sac','Saa']*rho.get('p6')+ (1-Y.get('fa'))*0.7*rho.get('p7')+(1-Y.get('c4'))*0.31*rho.get('p8')+(1-Y.get('c4'))*0.8*rho.get('p9')+(1-Y.get('pro'))*0.57*rho.get('p10')-rho.get('p11'))
     m.Equation(Ssu.dt() == Ssuo + rho.get('p2') + (1-f.at['Sfa','Xli'])*rho.get('p4')-rho.get('p5')) #1
     m.solve(disp= False )

     #Gas transfer Equations
     m.Equation(Shgas.dt()  == Vliq/Vgas*kla*(Sh2-16*KhH2*ph2))
     m.Equation(Sc4gas.dt() == Vliq/Vgas*kla*(Sc4-64*Khc4*pc4))
     m.Equation(Sco2gas.dt() == Vliq/Vgas*kla*(Sc4-Khco2*pco2))
     m.Equation(pc4 == Sc4gas*R*T/64)
     m.Equation(pco2 == Sco2gas*R*T)
     m.Equation(ph2 == Shgas*R*T/16)
     m.Equation(pgas == pc4+pco2+ph2)
     m.Equation(qgas == kp*(pgas-patm)*pgas/patm)
    

     stop = True
    except:
        a = a+1
        k = 1 +k
        m.open_folder()
        Ssu = m.Var(Ssuo*a,lb = 0, ub = 10000) # create GEKKO variable
        Saa = m.Var(Saao*a,lb = 0, ub = 10000) 
        Sfa = m.Var(Sfao*a,lb = 0, ub = 10000)
        Sva = m.Var(Svao*a,lb = 0, ub = 10000)
        Sbu = m.Var(Sbuo*a,lb = 0, ub = 10000)
        Spro = m.Var(Sproo*a,lb = 0, ub = 10000)
        Sac = m.Var(Saco,lb = 0, ub = 10000)
        Sh2 = m.Var(Sh2o,lb = 0, ub = 10000)
        Sc4 = m.Var(Sc4o,lb = 0, ub = 10000)
        SIC = m.Var(SICo,lb = 0, ub = 10000) 
        SIN = m.Var(SINo,lb = 0, ub = 10000)
        SI = m.Var(SIo,lb = 0, ub = 10000) 
        Xc = m.Var(0,lb = 0, ub = 10000)
        Xch = m.Var(0,lb = 0, ub = 10000)
        Xpro = m.Var(0,lb = 0, ub = 10000)
        Xli = m.Var(0,lb = 0, ub = 10000)
        Xsu = m.Var(0,lb = 0, ub = 10000)
        Xaa = m.Var(Xaao,lb = 0, ub = 10000)
        Xfa = m.Var(Xfao,lb = 0, ub = 10000)
        Xc4 = m.Var(Xc4o,lb = 0, ub = 10000)
        Xpr = m.Var(Xpr_o,lb = 0, ub = 10000)
        Xac = m.Var(Xaco,lb = 0, ub = 10000)
        Xh2 = m.Var(Xh2o,lb = 0, ub = 10000)
        XI = m.Var(0,lb = 0, ub = 10000)
        Sh = m.Var(10**(-8),0,10)
        Shgas = m.Var(0,lb = 0, ub = 10000)
        Sc4gas = m.Var(0,lb = 0, ub = 10000)
        Sco2gas = m.Var(0,lb = 0, ub = 10000)
        ph2 = m.Var(0,lb = 0, ub = 10000)
        pc4 = m.Var(0,lb = 0, ub = 10000)
        pco2 = m.Var(0,lb = 0, ub = 10000)
        pgas = m.Var(0,lb = 0, ub = 10000)
        qgas = m.Var(0,lb = 0, ub = 10000)

print(Ssu.value)
print(Saa.value)
print(Sfa.value)
print(Sva.value)
print(Sbu.value)
print(Spro.value)
print(Sc4.value)
print(Sac.value)
print(SIN.value)
print(Xch.value)


#After the model is solved you can print the final value of any variables in this area. Model is set up to minimize total soluble COD (What would be in the effluent after membrane separation)
#To print a variable use the format "print(varName.value)"

