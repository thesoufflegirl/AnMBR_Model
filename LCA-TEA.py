from gekko import GEKKO
import pandas as pd

m = GEKKO(remote=False) # create GEKKO model

scenario = 1
VFAon = 0

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
m.Equation(Sbu.dt() == Sbuo +(1-Y.get('su')) * f.at['Sbu','Ssu'] *rho.get('p5')+(1-Y.get('aa'))* f.at['Sbu','Saa'] *rho.get('p6') -rho.get('p9'))
m.Equation(Spro.dt() == Sproo +(1-Y.get('su')) * f.at['Spro','Ssu'] *rho.get('p5')+(1-Y.get('aa'))* f.at['Spro','Saa'] *rho.get('p6')+(1-Y.get('c4')*0.54)*rho.get('p8') -rho.get('p10'))
m.Equation(Sac.dt() == Saco +(1-Y.get('su')) * f.at['Sac','Ssu'] *rho.get('p5')+(1-Y.get('aa'))* f.at['Sac','Saa'] *rho.get('p6')+(1-Y.get('c4')*0.31)*rho.get('p8')+(1-Y.get('c4')*0.8)*rho.get('p9')+(1-Y.get('pro')*0.57*rho.get('p10')) -rho.get('p11'))
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
beta = m.Var(0,lb=0,ub = 1)
T_degas = T #Temperature
#k = kp from above heat capacity ratio, continue calling kp in this section

#Influent concentration = effluent concentration from above, partial pressure may change from temp, so recalculate
m.Equation(pin_degas == Sc4gas*R*T_degas/64 + Sco2gas*R*T_degas*(Shgas*R*T_degas/16))
Q = (1-beta)*Q
bfreq = 454609 #change if touching influent conc or flow
cfreq = 1/(2*14*Q) 
citric_acid = 2500*Q/1892705*1000/(1000*1000)
bleach= 2500*Q/1892705*5000/(1000*1000)


#VFAs
####################################Check external distillation file for energy needs vs energy produced and cost of equipment
#If feasible turn VFAon = 1
#Please note - THIS IS FOR MINIMUM ENERGY FOR SEPARATION, it uses an infinite plate approximation, the capital
#cost to build this is infinity. This is here to ensure there's not a situation where this becomes 
#environmentally viable. Also check if cost of operation is less than cost of VFAs produced, in which case capital expenditures
#on column will need to be included, and additional work to do the tradeoff between capex/energy needs to be evalauted 
#if that's the case there's examples of writing those equations in Gross's from CMU's group
#Their group's journal article is here: http://egon.cheme.cmu.edu/Papers/Caballero_Optimization_of_Distillation_Processes.pdf
#This proces has integer variables in it that I couldn't just if/and my way around, so it can slow down the model/require more iterations

col_nums = 1
minVFARemoval = 0.5 #this only is applicable as the max allowable in the bototms
B = m.Var(Q,lb = 0) 
D = m.Var(0,lb = 0)
xb = m.Var(1,lb=0,ub = 1)
xd = m.Var(0,lb=0,ub=1)
Rmin = m.Var(0,lb=0)
prev = m.Var(0, 0, 1, True)
prev2 = m.Var(0, 0, 1, True)
prev3 = m.Var(0, 0, 1, True)
kvfa = m.Var(0,lb = 0,ub = 1)
qcond = 0
qreb = 0
if VFAon == 1:
    l = Sva/109+Sac/60.05+Spro/74+Sbu/88.1
    #VFAs in tops water in bottoms
    K = 1.8e-5
    alphaa = 1.85
    alphap = 16.9
    alphab = 217
    alphav = 281
    m.equation(prev*K*l == prev*(B/Q*xb+B/V*xd))
    m.equation(prev*Rmin == prev*(alphaa*xd*Sac/(60.05*l)/(k-alphaa)+alphap*xd*Spro/(74*l)/(k-alphap)+alphab*xdxd*Sbu/(88.1*l)/(k-alphap))+alphab*xdxd*Sva/(109*l)/(k-alphav))
    #acetate in bottoms 
    K2 = 1.3e-5
    m.equation(prev2*K2*l == prev2*(B/Q*xb+B/V*xd))
    m.equation(prev2*Rmin == prev2*(alphap*xd*Spro/(74*l)/(k-alphap)+alphab*xdxd*Sbu/(88.1*l)/(k-alphap))+alphab*xdxd*Sva/(109*l)/(k-alphav))

    #propionic in bottoms
    K3= 2.0e-5
    m.equation(prev2*K3*l == prev2*(B/Q*xb+B/V*xd))
    m.equation(prev3*Rmin == prev3*(alphab*xdxd*Sbu/(88.1*l)/(k-alphap))+alphab*xdxd*Sva/(109*l)/(k-alphav))
    m.equation(prev+prev2+prev3 == 1)
    qcond = (D*(Rmin+1)*(394.3*xd*Sac/(60.05*l)+774*xd*Spro/(74*l)+570*xd*Spro/(88.1*l)+531*xd*Spro/(74*l)))
    Vb = (D*(Rmin+1)-V)/B
    qreb = B*Vb*(394.3*xd*Sac/(60.05*l)+774*xd*Spro/(74*l)+570*xd*Spro/(88.1*l)+531*xd*Spro/(74*l))
    Q1 = Vb
else:
    Qvfa = Q
    Q1 = Qvfa 


#Nitrogen Adsorption
#######################################################################################
#Solve at steady state, require 3 days uptime 

#Co is influent nitrgoen is SIN from AnMBR 
kth = 0.02 #from lab work 
qe = 2.51 #from lab work *******Changed for sensitivity analysis***********
time = 3
Ct = m.Var(SINo,0)
M = m.Var(0) #mass in column
P_nitro = m.Var(0)
SINrem = m.Var(SIN,0)
backwash = m.Var(0)
NH4Fert = m.Var(0)
#Qin = Qvfa 
t = m.Var(3) #set min days to 3
m.Equation(t.dt()==1)
m.Equation(Ct.dt()*Q1/3 == 1/(2.42**((kth*qe*M)/Q1)-kth*SIN*t))
#assume 3 columns split flow equally 
m.Equation(backwash == M*3.41e2) # lab data for backwash for clinoptilolite
m.Equation(NH4Fert == (SIN-SINrem)*Q1)
m.Equation(P_nitro == 3*Q1*1000*9.81*23912.5352*1.8*30/6/(3.6*10**6))

#Phosphorus removal
#####################################################################
#Just mass balance here

CaO = m.Var(0)
P_phos = m.Var(0)
ratio = 1/2
maxPhosphorus = Spho #Change if setting effluent limit - not touched in CWM1
Sph = m. Var(Spho,0)

m.Equation(Sph == Spho*Q1 - CaO*ratio)
Pfert = CaO*ratio*Q1/1000/1000
CaO_mass =CaO*Q1/1000/1000
m.Equation(P_phos == Q1*1000*9.81*23912.5352*1.8*30/6/(3.6*10**6)+1482)


#CWM1
##########################################################################################
#note - for bacteria the influent isn't the max value, we're not filtering and they're growing
So_wet = m.Var(0)
Sf_wet = m.Var(0)
Sa_wet = m.Var(0)
SI_wet = m.Var(0)
Snh_wet= m.Var(0)
Sno_wet = m.Var(0)
Sso4_wet = m.Var(0)
Sh2s_wet = m.Var(0)
Xs_wet =m.Var(0)
XI_wet = m.Var(0)
Xh_wet = m.Var(0)
Xa_wet = m.Var(0)
Xfb_wet = m.Var(0)
Xamb_wet = m.Var(0)
Xasrb_wet = m.Var(0)
Xsob_wet = m.Var(0)
Sch4_wet = m.Var(0)
t= m.Var(0,0)
P_constructed    = m.Var(0)

#Not sure why ADM1 was set upt like this 
conc = {'So_wet':So_wet,
        'Sf_wet':Sf_wet,
        'Sa_wet':Sa_wet,
        'SI_wet':SI_wet,
        'Snh_wet':Snh_wet,
        'Sno_wet':Sno_wet,
        'Sso4_wet':Sso4_wet,
        'Sh2s_wet':Sh2s_wet,
        'Xs_wet':Xs_wet,
        'XI_wet':XI_wet,
        'Xh_wet':Xh_wet,
        'Xa_wet':Xa_wet,
        'Xfb_wet':Xfb_wet,
        'Xamb_wet':Xamb_wet,
        'Xasrb_wet':Xasrb_wet,
        'Xsob_wet':Xsob_wet} 

So_weto = 2000
Sf_weto = Ssu
Sa_weto = Sac+Saa+Spro+Sva+Sc4
SI_weto = SIC
Snh_weto = SIN
Sno_weto = 0
Sso4_weto = 28
Sh2s_weto =0 
Xs_weto =0
XI_weto = 0
Xh_weto = 0
Xa_weto = 0
Xfb_weto = 0
Xamb_weto = 0
Xasrb_weto = 0
Xsob_weto = 0


#Clean mess into matrices


f_hyd_si = 0 
f_bm_sf = 0.05
f_bm_xi = 0.1
Yh = 0.63
Ya =0.24
Yfb = 0.053
Yamb = 0.032
Yasrb = 0.05
Ysob = 0.12
kh = 3
Kx = 0.1
mu_h = 0.1
Ksf = 2
Koh = 0.2
Knhh = 0.05
Knh = 0.05
Kh2sh = 140
mu_g = 0.2 #I'm 95% sure this is just incorrectly written in the published CWM1, somteims they mention an ng and it also appears to be
#a correction factor like mu_h is, I think it's just an error in prtinging. I'm treating mu_g as n_g
Knoh = 0.5
Ksa = 4
bh = 0.4
mu_a = 1
Koa =1 
Knha = 0.5
Kh2sa = 140
ba = 0.5
mu_fb = 3 #Another error in CWM1
Ksfb = 28
Kh2sfb = 140
Kofb = 0.2
Knofb = 0.5
Knhfb = 0.01
bfb = 0.02
mu_amb = 0.085
Ksamb = 56
Kh2samb = 140
Koamb =0.0002
Knoamb = 0.0005
Knhamb =  0.01
bamb = 0.008
mu_asrb = 0.18
Ksasrb = 24
Ksosarb = 19
Kh2asrb = 140
Koasrb = 0.0002
Knoasrb = 0.0005
Knhasrb =0.01
basrb = 0.012
mu_sob = 5.28
Kssob = 0.24
Kosob = 0.2
Knhsob = 0.05
eta_sob = 1 #This datat wasn't in the model. Using correction factor = 1 so this is the same as the other bacteria w/o correction factors
Knosob = 0.5
bsob = 0.15

i_n_sf = 0.03
i_n_si = 0.01
i_n_xs = 0.04
i_n_xi = 0.03
i_n_bm = 0.07

rho_wet = {'w1':kh*(Xs_wet/(Xh_wet+Xfb_wet)/(Kx+(Xs_wet/(Xh_wet+Xfb_wet))))*(Xh_wet+mu_h*Xfb_wet),
    'w2':mu_h*Sf_wet/(Ksf+Sf_wet)*(Sf_wet/(Sf_wet*Sa_wet))*So_wet/(Koh+So_wet)*Snh_wet/(Knhh+Snh_wet)*Kh2sh/(Kh2sh+Sh2s_wet)*Xh_wet,
    'w3': mu_g*mu_h*Sf_wet/(Ksf+Sf_wet)*Sf_wet/(Sf_wet+Sa_wet)*Koh/(Koh+So_wet)*Sno_wet/(Knoh+Sno_wet)*Kh2sh/(Kh2sh+Sh2s_wet)*Xh_wet,
    'w4':mu_h*Sa_wet/(Ksa+Sa_wet)*Sa_wet/(Sf_wet+Sa_wet)*So_wet/(Koh+So_wet)*Snh_wet/(Knh+Snh_wet)*Kh2sh/(Kh2sh+Sh2s_wet)*Xh_wet,
    'w5':mu_g*mu_h*Sa_wet/(Ksa+Sa_wet)*Sa_wet/(Sf_wet+Sa_wet)*Koh/(Koh+So_wet)*Sno_wet/(Knoh+Sno_wet)*Snh_wet/(Knh+Snh_wet)*Kh2sh/(Kh2sh+Sh2s_wet)*Xh_wet,
    'w6':bh*Xamb_wet,
    'w7':mu_a*Snh_wet/(Knha+Snh_wet)*So_wet/(Koa+So_wet)*Kh2sa/(Kh2sa+Sh2s_wet)*Xa_wet,
    'w8':ba*Xa_wet,
    'w9':mu_fb*Sf_wet/(Ksfb+Sf_wet)*Kh2sfb/(Kh2sfb+Sh2s_wet)*Kofb/(Kofb+So_wet)*Knofb/(Knofb+Sno_wet)*Snh_wet/(Knhfb+Snh_wet)*Xfb_wet,
    'w10':bfb*Xfb_wet,
    'w11':mu_amb*Sa_wet/(Ksamb+Sa_wet)+Kh2samb/(Kh2sfb+Sh2s_wet)+Koamb/(Koamb+So_wet)*Knoamb/(Knoamb+Sno_wet)*Snh_wet/(Knhamb+Snh_wet)*Xamb_wet,
    'w12':bamb*Xamb_wet,
    'w13':mu_asrb*Sa_wet/(Ksasrb +Sa_wet)*Sso4_wet/(Koasrb+Sso4_wet)*Kh2asrb/(Kh2asrb+Sh2s_wet)*Koasrb/(Koasrb+So_wet)*Knoasrb/(Knoasrb+Sno_wet)*Snh_wet/(Knhasrb+Snh_wet)*Xasrb_wet,
    'w14':basrb*Xasrb_wet,
    'w15':mu_sob * Sh2s_wet/(Kssob+Sh2s_wet)*So_wet/(Kosob+So_wet)*Snh_wet/(Knhsob+Snh_wet)*Xsob_wet,
    'w16':mu_sob*eta_sob*Sh2s_wet/(Kssob+Sh2s_wet)*Sno_wet/(Knosob + Sno_wet)*Kosob/(Kosob+So_wet)*Snh_wet/(Knhsob+Snh_wet)*Xsob_wet,
    'w17':bsob*Xsob_wet
    }

nu_wet = {'n1':i_n_xs - (1-f_hyd_si*i_n_sf)-f_hyd_si*i_n_si,
        'n2': i_n_sf/Yh-i_n_bm,
        'n3':i_n_sf/Yh-i_n_bm,
        'n4': -i_n_bm,
        'n5':-i_n_bm,
        'n6':i_n_bm-f_bm_sf*i_n_sf -(1- f_bm_sf-f_bm_xi)*i_n_xs-f_bm_xi*i_n_xi,
        'n7':-i_n_bm-1/Ya,
        'n8':i_n_bm-f_bm_sf*i_n_sf -(1- f_bm_sf-f_bm_xi)*i_n_xs-f_bm_xi*i_n_xi,
        'n9':i_n_sf/Yfb-i_n_bm,
        'n10':i_n_bm-f_bm_sf*i_n_sf -(1- f_bm_sf-f_bm_xi)*i_n_xs-f_bm_xi*i_n_xi,
        'n11':-i_n_bm,
        'n12':i_n_bm-f_bm_sf*i_n_sf -(1- f_bm_sf-f_bm_xi)*i_n_xs-f_bm_xi*i_n_xi,
        'n13':-i_n_bm,
        'n14':i_n_bm-f_bm_sf*i_n_sf -(1- f_bm_sf-f_bm_xi)*i_n_xs-f_bm_xi*i_n_xi,
        'n15':-i_n_bm,
        'n16':-i_n_bm,
        'n17':i_n_bm-f_bm_sf*i_n_sf -(1- f_bm_sf-f_bm_xi)*i_n_xs-f_bm_xi*i_n_xi
}

m.Equation(So_wet.dt() == So_weto + rho_wet.get('w2')*(1-Yh)+ rho_wet.get('w4')*(1-Yh)-rho_wet.get('w7')*(4.57-Ya)/Ya-rho_wet.get('w15')*(2-Ysob)/Ysob)
m.Equation(Sf_wet.dt() == Sf_weto + rho_wet.get('w1')*(1-f_hyd_si)-rho_wet.get('w2')*(1/Yh)-rho_wet.get('w3')*(1/Yh)+rho_wet.get('w6')*f_bm_sf+rho_wet.get('w8')*f_bm_sf-rho_wet.get('w9')*(1/Yfb)+rho_wet.get('w10')*f_bm_sf+rho_wet.get('w12')*f_bm_sf+rho_wet.get('w14')*f_bm_sf+rho_wet.get('w17')*f_bm_sf)
m.Equation(Sa_wet.dt() == Sa_weto - rho_wet.get('w4')/Yh- rho_wet.get('w5')/Yh+rho_wet.get('w9')*(1-Yfb)/Yfb- rho_wet.get('w11')/Yamb)- rho_wet.get('w13')/Yasrb
m.Equation(SI_wet.dt() == SI_weto + f_hyd_si*rho_wet.get('w1'))
m.Equation(Snh_wet.dt() == Snh_weto + rho_wet.get('w1')*nu_wet.get('n1')+ rho_wet.get('w2')*nu_wet.get('n2')+ rho_wet.get('w3')*nu_wet.get('n3')+ rho_wet.get('w4')*nu_wet.get('n4')+ rho_wet.get('w5')*nu_wet.get('n5')+rho_wet.get('w6')*nu_wet.get('n6')+ rho_wet.get('w7')*nu_wet.get('n7')+ 
    rho_wet.get('w8')*nu_wet.get('n8')+ rho_wet.get('w9')*nu_wet.get('n9')+ rho_wet.get('w10')*nu_wet.get('n10')+ rho_wet.get('w11')*nu_wet.get('n11')+ rho_wet.get('w12')*nu_wet.get('n12')+ rho_wet.get('w13')*nu_wet.get('n13')+ rho_wet.get('w14')*nu_wet.get('n14')+ rho_wet.get('w15')*nu_wet.get('n15')+ rho_wet.get('w16')*nu_wet.get('n16')+ rho_wet.get('w17')*nu_wet.get('n17'))
m.Equation(Sno_wet.dt() == Sno_weto- rho_wet.get('w3')*(1-Yh)/(2.86*Yh)- rho_wet.get('w5')*(1-Yh)/(2.86*Yh)+rho_wet.get('w7')/Ya-rho_wet.get('w16')*(1-Ysob)/(0.875*Ysob))
m.Equation(Sso4_wet.dt() == Sso4_weto-rho_wet.get('w13')*(1-Yasrb)/(2*Yasrb)+rho_wet.get('w15')/Ysob+rho_wet.get('w16')/Ysob)
m.Equation(Sh2s_wet.dt() == Sh2s_weto-rho_wet.get('w13')*(1-Yasrb)/(2*Yasrb)+rho_wet.get('w15')/Ysob+rho_wet.get('w16')/Ysob)
m.Equation(Xs_wet.dt() == Xs_weto -rho_wet.get('w1')+rho_wet.get('w6')*(1-f_bm_sf-f_bm_xi)+rho_wet.get('w8')*(1-f_bm_sf-f_bm_xi)+rho_wet.get('w10')*(1-f_bm_sf-f_bm_xi)+rho_wet.get('w12')*(1-f_bm_sf-f_bm_xi)+rho_wet.get('w14')*(1-f_bm_sf-f_bm_xi)+rho_wet.get('w17')*(1-f_bm_sf-f_bm_xi))
m.Equation(XI_wet.dt() == XI_weto +rho_wet.get('w6')*f_bm_xi+rho_wet.get('w8')*f_bm_xi+rho_wet.get('w10')*f_bm_xi+rho_wet.get('w12')*f_bm_xi+rho_wet.get('w14')*f_bm_xi+rho_wet.get('w17')*f_bm_xi)
m.Equation(Xh_wet.dt() == Xh_weto + rho_wet.get('w2')+ rho_wet.get('w2')+ rho_wet.get('w3')+ rho_wet.get('w4')+ rho_wet.get('w5')-rho_wet.get('w6'))
m.Equation(Xa_wet.dt() == Xa_weto + rho_wet.get('w7')-rho_wet.get('w8'))
m.Equation(Xfb_wet.dt() == Xfb_weto +rho_wet.get('w9')-rho_wet.get('w10'))
m.Equation(Xamb_wet.dt() == Xamb_weto+rho_wet.get('w11')-rho_wet.get('w12'))
m.Equation(Xasrb_wet.dt() == Xasrb_weto + rho_wet.get('w13')-rho_wet.get('w14'))
m.Equation(Xsob_wet.dt() == Xsob_weto + rho_wet.get('w15')+rho_wet.get('w16')-rho_wet.get('w17'))
m.Equation(Sch4_wet.dt() == 141*Q1+rho_wet.get('w16')-rho_wet.get('w17'))
m.Equation(1 == t.dt())
m.Equation(P_constructed== Q1*1000*9.81*23912.5352*1.8/(3.6*10**6))




###########################################################################################
#LCA
###########################################################################################



LCI = {'electricity':0.054, #kg CO2/ kWh
    'lime':0.94, #kg/kg
    'citric':5.24-8, #kg/kg,
    'bleach':3.81e-5, #kg/kg dry,
    'clinoptilolite':0.000207, #kg/kg as zeolite,
    'brine':5.09e-6, #kg/kg dry,
    'methane':29.8, #kg/kg,
    'nox':3,#mg/m2day,
    'anhydrousAmmonia':10.59e-5, #kg/kg,
    'phosphorus':8.81e-5, #kg/kg,
    'heating':2e-6,#kg/MJ)
    }

Q_gal = Q*0.264172
in_pump_elec = 1600*Q_gal**0.4631

AnMBR_emissions = m.Var(0)
elec_degasser = 0

if scenario == 1:
    elec_dissolved_methane = Sc4gas*Q/1000/1000*0.9*55.5*0.8*0.4375
    AnMBR_meth = 0.1*Sc4gas*Q/1000/1000 #biogenic so subtract one mole co2 for one mole methane
    AnMBR_CO2 = 0 #biogenic
    elec_degasser = 0.8*Q*3*R*temp/2*(5)/1000000
elif scenario == 2:
    elec_dissolved_methane = 0
    AnMBR_meth = 0.1*Sc4gas*Q/1000/1000 #90% flare 
    AnMBR_CO2 = 0 #biogenic
else:
    elec_dissolved_methane = 0
    AnMBR_meth = Sc4gas*Q 
    AnMBR_CO2 = 0 #biogenic

heat_balance = m.Var(0)
heat_AnMBR = m.Var(0)
m.Equation(heat_balance == ((T-temp)*4.18*Q*1000)/1000000-Sc4gas*55.5*0.8)
heat_AnMBR = m.if2(heat_balance,((T-temp)*4.18*Q*1000)/1000000-Sc4gas*55.5*0.8,0)
elec_prod_AnMBR = m.if2(heat_balance,0,(Sc4gas*55.5*0.8-((T-temp)*4.18*Q*1000)/1000000)*0.4375) 
elec_con_AnMBR = in_pump_elec*3+qgas*0.264172*1600**0.4631+ 1/30*bfreq*0.26300*Q_gal**0.4631
elec_AnMBR = m.Var(0)
m.Equation(elec_AnMBR == elec_con_AnMBR-(elec_prod_AnMBR+elec_dissolved_methane))

m.Equation(AnMBR_emissions == heat_AnMBR*LCI.get('heating')+elec_AnMBR*LCI.get('electricity')/3.6*elec_degasser*LCI.get('electricity')/3.6+AnMBR_meth*(LCI.get('methane')-1)+citric_acid*LCI.get('citric')+bleach*LCI.get('bleach'))

VFA_emissions = m.Var(0)

m.Equation(VFA_emissions == qcond*LCI.get('heating')+qreb*LCI.get('heating')+in_pump_elec*LCI.get('electricity'))

nitrogen_emissions = m.Var(0)

m.Equation(nitrogen_emissions == P_nitro*LCI.get('electricity')-NH4Fert*LCI.get('anhydrousAmmonia')+M*LCI.get('clinoptilolite'))

p_emissions = m.Var(0)

m.Equation(p_emissions == P_phos*LCI.get('electricity')-Pfert*LCI.get('phosphorus')+CaO_mass*LCI.get('lime'))

CW_emissions = m.Var(0,lb=0)
m.Equation(CW_emissions == P_constructed*LCI.get('electricity')+Sch4_wet*Q1*LCI.get('methane'))





#################################################
#TEA

Rs = 0.3
Pe  = 0.073
Mc = 80
MA = 15
Rc = 474
CVc = 150
Pc = 98.5
perm = 1578
sludge = 518
VDtank =158000
VDpump =4200
CWc = 154800*1000000/100000
AA = 153100*3
FS = 15894
D = 0.133681*5
Caoc = 1.12
clinopc = 0.28
citricc = 14.8
bleachc = 3
heatc = 0.18/648000
fertc = 1300
Pfertc =600
years  = 20

hanmbr = m.Var(0,lb=0,ub=40)
RManmbr = m.Var(0,lb=0)
PWanmbr = m.Var(0,lb=0)
CVanmbr = m.Var(0,lb=0)
Danmbr = m.Var(0,lb=0)
Rmpr = m.Var(0,lb=0)
hpr =m.Var(0,lb=0,ub= 40)
A = m.Var(0,lb = 0)
AM= m.Var(0,lb=0)
CVpr = m.Var(0,lb=0)
Wpr = m.Var(0,lb=0)
Pwpr = m.Var(0,lb=0)
Capital_cost = m.Var(0,lb=0)
CWA = m.Var(0, lb = 0)
Annual_op = m.Var(0,lb =0)
NPV = m.Var(0,lb=0)


Vanmbr = 0.3*Vliq*(1-0.3)
m.Equation(Danmbr == 2*(Vanmbr/3.14/hanmbr))
m.Equation(RManmbr == hanmbr/Rs*2*3.14*Danmbr/2+(Danmbr/2*3.14)**(1/2)/scenario)
m.Equation(CVanmbr == 0.275*(hanmbr+4.5)*Danmbr*(7.4+0.5*hanmbr)+0.0675*Danmbr**2*(8+0.24*15.5*(hanmbr))+Danmbr*3.14)
m.Equation(PWanmbr==0.275*(hanmbr+4.5)*(Danmbr*(7.5+0.5*hanmbr))+0.0675*Danmbr**2*(8+0.24*15.5*hanmbr)+Danmbr*hanmbr)
m.Equation(AM == Q/MA)
m.Equation(Rmpr == hpr/Rs*A+A/Rs)
m.Equation(CVpr == 0.275*(hpr+4.5)*Wpr*(7.5+0.5*hpr))
m.Equation(Pwpr == hpr +2.17*Wpr+36+0.5*(2.24*Wpr+70))
m.Equation(CWA ==Q*t/D)
m.Equation(Annual_op == CaO_mass*Caoc+M*clinopc+citric_acid*citricc+bleach*bleachc +elec_AnMBR*Pe+heat_AnMBR*heatc +(P_nitro+P_phos+P_constructed)*Pe-NH4Fert*fertc-Pfert*Pc)

m.Equation(Capital_cost ==FS+ RManmbr*Rc+CVanmbr*CVc+PWanmbr*Pc+Rmpr*Rc+12*perm+3*sludge+MA*Mc+AA+CVpr*CVc+Pwpr*Pc+Rmpr*Rc+VDtank+VDpump +CWA*CWc)
m.Equation(NPV == -Capital_cost+Annual_op/(1+0.1)**1+Annual_op/(1+0.1)**2+Annual_op/(1+0.1)**3+Annual_op/(1+0.1)**4+Annual_op/(1+0.1)**5-MA*Mc/(1+0.1)**5+Annual_op/(1+0.1)**6+Annual_op/(1+0.1)**7+Annual_op/(1+0.1)**8+Annual_op/(1+0.1)**9+Annual_op/(1+0.1)**10-MA*Mc/(1+0.1)**10+Annual_op/(1+0.1)**11+Annual_op/(1+0.1)**12+Annual_op/(1+0.1)**13+Annual_op/(1+0.1)**14+Annual_op/(1+0.1)**15-MA*Mc/(1+0.1)**15+Annual_op/(1+0.1)**16+Annual_op/(1+0.1)**17+Annual_op/(1+0.1)**18+Annual_op/(1+0.1)**19+Annual_op/(1+0.1)**20)


m.time = [1,5,10,25,50,100,150,200]      # time points if doing non-steady state (these are in hours, adding more time points will make model run slower)
m.options.IMODE = 2 # simulation mode 
#For simulation mode if trying to do long term steady state 2 else use 5
m.options.NODES = 4 
m.options.CSV_WRITE=2  
m.options.MAX_ITER = 100  #Increase if not converging 
                        #   with internal nodes

m.Minimize(NPV)
#m.Minimize(AnMBR_emissions + VFA_emissions  +nitrogen_emissions +p_emissions +CW_emissions  )
                            

m.solve(disp=True)     # solve


