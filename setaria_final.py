import scobra
import csv
import pickle
from cobra import Reaction
import numpy as np
import pandas as pd
import os

def Datapath(Path='.\Data'):
    return Path


def SplitRevRxn(model):
    #modify.convert_to_irreversible(model)
    reactions_to_add = []
    for reaction in model.reactions:
        if reaction.reversibility == True:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = min(0, reaction.upper_bound) * -1
            reverse_reaction.upper_bound = reaction.lower_bound * -1
            reaction.lower_bound = 0
            reaction.upper_bound = max(0, reaction.upper_bound)
            #Make the directions aware of each other
            reaction.reflection = reverse_reaction
            reverse_reaction.reflection = reaction
            reaction_dict = dict([(k, v*-1) for k, v in reaction._metabolites.items()])
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction._gene_reaction_rule = reaction._gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    return model


def MakeInvGEWtDict(m, DWt='Mean',
                GEbase="GSMx2rxn_baseGE_N0_edited.txt",
                GElwr="GSMx2rxn_lwrGE_N0_edited.txt",
                GEmid="GSMx2rxn_midGE_N0_edited.txt",
                GEtip="GSMx2rxn_tipGE_N0_edited.txt"):
    '''m= Scobra Model File,
    GEfiles= give GE files path/ name
    DWt='Mean' or '10Mean' or 0.1Mean or 'Any Other Number'
    '''
    Path=Datapath()
    GEbase=os.path.join(Path,GEbase)
    GElwr=os.path.join(Path,GElwr)
    GEmid=os.path.join(Path,GEmid)
    GEtip=os.path.join(Path,GEtip)


    M=m.copy()
    #######SPLITTING REVERSIBLE REACTIONS
    M_split=M.copy()
    M_split=SplitRevRxn(M_split)

    ##making dict of GE values
    b1=open(GEbase,'r')
    b2=[x.strip() for x in b1.readlines()]
    b3=[y.split('\t') for y in b2]
    baseGE = {(b[0]+'_base'):(1/float(b[1])) for b in b3}
    b1.close()
    #print(baseGE)

    l1=open(GElwr,'r')
    l2=[x.strip() for x in l1.readlines()]
    l3=[y.split('\t') for y in l2]
    lwrGE = {(l[0]+'_lwr'):(1/float(l[1])) for l in l3}
    l1.close()
    #print(lwrGE)

    z1=open(GEmid,'r')
    z2=[x.strip() for x in z1.readlines()]
    z3=[y.split('\t') for y in z2]
    midGE = {(z[0]+'_mid'):(1/float(z[1])) for z in z3}
    z1.close()
    #print(midGE)

    t1=open(GEtip,'r')
    t2=[x.strip() for x in t1.readlines()]
    t3=[y.split('\t') for y in t2]
    tipGE = {(t[0]+'_tip'):(1/float(t[1])) for t in t3}
    t1.close()
    #print(tipGE)

    allGE={}
    allGE.update(baseGE)
    allGE.update(lwrGE)
    allGE.update(midGE)
    allGE.update(tipGE)

    #adding _reverse reactions in GE dict
    revrxn=[r for r in M_split.Reactions() if "_reverse" in r]
    revtemp=[x.replace("_reverse","") for x in revrxn]
    revGE={v+"_reverse":allGE[v] for v in revtemp if v in allGE.keys()}
    allGE.update(revGE)

    rest=[r for r in M_split.Reactions() if r not in allGE.keys()]
    NOGE=[]
    for rx in rest:
        if '_CP_BS' in rx or '_MPBS' in rx:
            NOGE.append(rx)
        if '_tx' in rx or '_biomass' in rx:
            NOGE.append(rx)
        if 'Lipids' in rx or 'AminoAcid' in rx or 'NuclicAcid' in rx or 'CarbosLignin' in rx:
            NOGE.append(rx)

    y=[1/Val for Val in allGE.values()]
    if DWt == 'Mean':
        YV=1/np.mean(y)
        restGE={x:YV for x in rest if x not in NOGE}
    elif DWt == '10Mean':
        YV=10/np.mean(y)
        restGE={x:YV for x in rest if x not in NOGE}
    elif DWt == '0.1Mean':
        YV=0.1/np.mean(y)
        restGE={x:YV for x in rest if x not in NOGE}
    elif type(DWt) == int or type(DWt) == float:
        YV=DWt
        restGE={x:YV for x in rest if x not in NOGE}
    else:
        YV=eval(DWt)
        restGE={x:YV for x in rest if x not in NOGE}
    allGE.update(restGE)

    return allGE


def MakeInvGEWtDictNOTRANSPORT(m,
                GEbase="GSMx2rxn_baseGE_N0_edited.txt",
                GElwr="GSMx2rxn_lwrGE_N0_edited.txt",
                GEmid="GSMx2rxn_midGE_N0_edited.txt",
                GEtip="GSMx2rxn_tipGE_N0_edited.txt"):
    '''m= Scobra Model File,
    GEfiles= give GE files path/ name
    '''
    Path=Datapath()
    GEbase=os.path.join(Path,GEbase)
    GElwr=os.path.join(Path,GElwr)
    GEmid=os.path.join(Path,GEmid)
    GEtip=os.path.join(Path,GEtip)



    M=m.copy()
    #######SPLITTING REVERSIBLE REACTIONS
    M_split=M.copy()
    M_split=SplitRevRxn(M_split)

    ##making dict of GE values
    b1=open(GEbase,'r')
    b2=[x.strip() for x in b1.readlines()]
    b3=[y.split('\t') for y in b2]
    baseGE = {(b[0]+'_base'):(1/float(b[1])) for b in b3}
    b1.close()
    #print(baseGE)

    l1=open(GElwr,'r')
    l2=[x.strip() for x in l1.readlines()]
    l3=[y.split('\t') for y in l2]
    lwrGE = {(l[0]+'_lwr'):(1/float(l[1])) for l in l3}
    l1.close()
    #print(lwrGE)

    z1=open(GEmid,'r')
    z2=[x.strip() for x in z1.readlines()]
    z3=[y.split('\t') for y in z2]
    midGE = {(z[0]+'_mid'):(1/float(z[1])) for z in z3}
    z1.close()
    #print(midGE)

    t1=open(GEtip,'r')
    t2=[x.strip() for x in t1.readlines()]
    t3=[y.split('\t') for y in t2]
    tipGE = {(t[0]+'_tip'):(1/float(t[1])) for t in t3}
    t1.close()
    #print(tipGE)

    allGE={}
    allGE.update(baseGE)
    allGE.update(lwrGE)
    allGE.update(midGE)
    allGE.update(tipGE)

    #adding _reverse reactions in GE dict
    revrxn=[r for r in M_split.Reactions() if "_reverse" in r]
    revtemp=[x.replace("_reverse","") for x in revrxn]
    revGE={v+"_reverse":allGE[v] for v in revtemp if v in allGE.keys()}
    allGE.update(revGE)

    rest=[r for r in M_split.Reactions() if r not in allGE.keys()]
    NOGE=[]
    for rx in rest:
        if '_CP_BS' in rx or '_MPBS' in rx:
            NOGE.append(rx)
        if '_tx' in rx or '_biomass' in rx:
            NOGE.append(rx)
        if 'Lipids' in rx or 'AminoAcid' in rx or 'NuclicAcid' in rx or 'CarbosLignin' in rx:
            NOGE.append(rx)
        if '_pc' in rx or '_mc' in rx or '_vc' in rx or '_xc' in rx:
            NOGE.append(rx)

    y=[1/Val for Val in allGE.values()]
    YV=1/np.mean(y)
    restGE={x:YV for x in rest if x not in NOGE}
    allGE.update(restGE)

    return allGE

def SetPhoton(m, ll):
    '''m=scobra model (to set light upper limit)
        ll= lifgt upper limit float value'''
    ll=ll
    ##Light_fix
    P={i:(0,ll) for i in m.Reactions() if 'Photon_tx' in i}
    m.SetConstraints(P)
    #return m

def SetOpenBiomass(m):
    biomass = {x:(-100.0, -0.0) for x in m.Reactions() if '_biomass' in x}
    #Set BIOMASS constraint
    m.SetConstraints(biomass)
    #return m

def SetMaintenance(m):
    mentainanceATP={i:(-8.5,-8.5) for i in m.Reactions() if 'ATPase_tx' in i}
    mentainanceNADP={i:(-0.944,-0.944) for i in m.Reactions() if 'NADPHoxc_tx' in i or 'NADPHoxm_tx' in i or 'NADPHoxp_tx' in i }
    NADPHoxx={i:(0,0) for i in m.Reactions() if 'NADPHoxx_tx' in i}

    #Set Mentainance constraint
    m.SetConstraints(mentainanceATP)
    m.SetConstraints(mentainanceNADP)
    #BLOCK NADPHoxx_tx
    m.SetConstraints(NADPHoxx)
    #return m

def SetEQUALPhoton(m):
    P={i:1 for i in m.Reactions() if 'Photon_tx' in i}
    m.SetReacsFixedRatio(P)
    #return m

def SetMPBSRatio(m):
    m.SetReacsFixedRatio({'Cell_biomass_MP_tip':1,'Cell_biomass_BS_tip':1})
    m.SetReacsFixedRatio({'Cell_biomass_MP_mid':1,'Cell_biomass_BS_mid':1})
    m.SetReacsFixedRatio({'Cell_biomass_MP_lwr':1,'Cell_biomass_BS_lwr':1})
    m.SetReacsFixedRatio({'Cell_biomass_MP_base':1,'Cell_biomass_BS_base':1})

def SetMPRubiscoRatio(m):
    m.SetReacsFixedRatio({'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN_p_MP_base':3,'RXN-961_p_MP_base':1})
    m.SetReacsFixedRatio({'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN_p_MP_lwr':3,'RXN-961_p_MP_lwr':1})
    m.SetReacsFixedRatio({'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN_p_MP_mid':3,'RXN-961_p_MP_mid':1})
    m.SetReacsFixedRatio({'RIBULOSE-BISPHOSPHATE-CARBOXYLASE-RXN_p_MP_tip':3,'RXN-961_p_MP_tip':1})

def SetBIOMASSRatio(m, r_tmlb=[]):
    '''m=scobra model to set biomas ratio,
    r_tmlb = LIST of ratio values, tip , mid, lwr, base, IN THIS ORDER'''
    D=r_tmlb
    T,M,L,B=D[0],D[1],D[2],D[3]
    m.SetReacsFixedRatio({'Cell_biomass_BS_tip':T,'Cell_biomass_BS_mid':M,'Cell_biomass_BS_lwr':L,'Cell_biomass_BS_base':B})


def BlockExtraRxns(m, rxnlist=[]):
    blklist={i:(0,0) for i in rxnlist}
    m.SetConstraints(blklist)



def TwoStepBiomassGradientSolve(m, P=200.0, r_tmlb=[], EqualPhoton=True, DWt='Mean',
                                Path='.\Data', SaveModel= False, ReturnSols= False,
                                SaveSols= False,
                                GEbase="GSMx2rxn_baseGE_N0_edited.txt",
                                GElwr="GSMx2rxn_lwrGE_N0_edited.txt",
                                GEmid="GSMx2rxn_midGE_N0_edited.txt",
                                GEtip="GSMx2rxn_tipGE_N0_edited.txt"):
    '''m=scobra model to run,
        P= photon amount in float
        DWt='Mean' or '10Mean' or 'Any Other Number'
        r_tmlb= LIST of ratio values, tip , mid, lwr, base, IN THIS ORDER,
        EqualPhoton= True/ False (all cells will use equal photon or not),
        GEfiles= give GE files path/ name'''

    Path=Datapath(Path)
    GEbase=os.path.join(Path,GEbase)
    GElwr=os.path.join(Path,GElwr)
    GEmid=os.path.join(Path,GEmid)
    GEtip=os.path.join(Path,GEtip)

    if SaveModel or SaveSols:
        os.mkdir('.\\Results')


    m4=m.copy()

    SetPhoton(m4, P)
    SetOpenBiomass(m4)
    SetMaintenance(m4)

    ##BLOCK specific rxns
    rxns2blk=[x for x in m4.Reactions() if 'Plastoquinol_Oxidase_p' in x or 'NADPH_Dehydrogenase_p' in x or 'MALIC-NADP-RXN_c' in x or 'PEPCARBOXYKIN-RXN_c' in x]
    BlockExtraRxns(m4, rxns2blk)
    adlrxns2blk=[x for x in m4.Reactions() if 'CO2_tx_BS' in x or 'O2_tx_BS' in x or '4.1.1.32-RXN_c_BS' in x or '1.1.1.39-RXN_m_BS' in x]
    BlockExtraRxns(m4, adlrxns2blk)

    M_split=m4.copy()
    M_split=SplitRevRxn(M_split)
    M_split1=M_split.copy()

    if EqualPhoton:
        SetEQUALPhoton(M_split1)
    SetBIOMASSRatio(M_split1, r_tmlb=r_tmlb)
    SetMPBSRatio(M_split1)
    SetMPRubiscoRatio(M_split1)
    #SetBiomassMaintenanceRatio(M_split1)

    #SET BIOMASS MAXIMIZATION OBJECTIVE
    Rxnbiomass=[x for x in M_split1.Reactions() if '_biomass' in x]
    M_split1.SetObjective(Rxnbiomass)
    M_split1.SetObjDirec(direc='Min') #AS BIOMASS FLUX IS -ve, MINIMIZING THE FLUX is MAXIMIZING BIOMASS

    #solve
    M_split1.MinFluxSolve()
    SolwoGEBioMax=M_split1.GetSol()

    print ("\n\ntx\n")
    M_split1.PrintSol('_tx')
    print "\n\nbiomass\n"
    M_split1.PrintSol('_biomass')
    print("\n\nCPs\n")
    M_split1.PrintSol('_CP')
    print("\n")

    if SaveModel:
        M_split1.WriteModel('.\\Results\\After_Run_BioMax_model_run6.xlsx')
        with open('.\\Results\\BioMax_Model_run6.pkl', 'wb') as outf:
            pickle.dump(M_split1, outf)
    if SaveSols:
        with open('.\\Results\\BioMax_Sol_run6.pkl', 'wb') as outf2:
            pickle.dump(SolwoGEBioMax, outf2)
        with open('.\\Results\\BioMax_Sol_run6.csv','w') as outf3:
            writer = csv.writer(outf3, delimiter='\t')
            for key, value in SolwoGEBioMax.items():
                writer.writerow([key, value])

    ### SET Biomass values for GE Integration
    BMax=SolwoGEBioMax.Filter('Cell_biomass_BS_tip')
    Bvalue=BMax['Cell_biomass_BS_tip']


    M_split2=M_split.copy()

    SetPhoton(M_split2, 1000.0)
    if EqualPhoton:
        SetEQUALPhoton(M_split2)
    SetBIOMASSRatio(M_split2, r_tmlb=r_tmlb)
    SetMPBSRatio(M_split2)
    SetMPRubiscoRatio(M_split2)
    #SetBiomassMaintenanceRatio(M_split2)

    BStip_biomass = {x:(Bvalue, Bvalue) for x in M_split2.Reactions() if 'Cell_biomass_BS_tip' in x}
    M_split2.SetConstraints(BStip_biomass)

    ### Addition of GE values
    ## making dict of GE values
    allGE=MakeInvGEWtDict(M_split2, DWt=DWt, GEbase=GEbase, GElwr=GElwr, GEmid=GEmid, GEtip=GEtip)
    print("GEDict length: "+str(len(allGE))+"\n")

    ## applying weight of GE values
    M_split2.SetObjective(allGE)
    M_split2.SetObjDirec(direc='Min')

    M_split2.MinFluxSolve()
    SolwGE2=M_split2.GetSol()

    print ("\nSolution with GE (GIVEN GE ONLY + Grad "+str(r_tmlb)+" )\n")
    if EqualPhoton:
        print("EqualPhoton: YES\n")
    else:
        print("EqualPhoton: NO\n")
    print ("\n\ntx\n")
    M_split2.PrintSol('_tx')
    print "\n\nbiomass\n"
    M_split2.PrintSol('_biomass')
    print("\n\nCPs\n")
    M_split2.PrintSol('_CP')

    if SaveModel:
        M_split2.WriteModel('.\\Results\\After_Run_wGE_model_run6.xlsx')
        with open('.\\Results\\wGE_Model_run6.pkl', 'wb') as outff:
            pickle.dump(M_split2, outff)
    if SaveSols:
        with open('.\\Results\\wGE_Sol_run6.pkl', 'wb') as outff2:
            pickle.dump(SolwGE2, outff2)
        with open('.\\Results\\wGE_Sol_run6.csv','w') as outff3:
            writer = csv.writer(outff3, delimiter='\t')
            for key, value in SolwGE2.items():
                writer.writerow([key, value])



    if ReturnSols:
        return SolwoGEBioMax, SolwGE2
    else:
        return M_split2, SolwGE2
