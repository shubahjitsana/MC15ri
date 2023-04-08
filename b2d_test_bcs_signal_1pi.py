#!/usr/bin/env python3
import sys
import basf2 as b2
import modularAnalysis as ma
# import stdV0s
import stdPi0s
import flavorTagger as ft
from variables import variables as vm  # shorthand for VariableManager
import variables.collections as vc
import variables.utils as vu
import vertex as vx
import variables.MCGenTopo as vmc   # module( interface of basf2 to TopoAna)
                                    # to extract out Topological information
# get input file number from the command line
file_name = sys.argv[1]
if "sub00/" in file_name :
    file_number = file_name.split("sub00/")
else :
    file_number = file_name.split("sub01/")
input_filename = file_number[1]

# get output folder name from the command line
output_foldername = sys.argv[2]
output_filename = input_filename

# weight file
weightfile = sys.argv[3]
#################################################################################

# set analysis global tag (needed for flavor tagging) [S10]
b2.conditions.prepend_globaltag(ma.getAnalysisGlobaltag())  # [E10]

# Perform analysis.
main = b2.create_path()

ma.inputMdstList(environmentType="default",filelist=f"{file_name}",path=main)


# fill final state particle lists
ma.fillParticleList("pi+:myPions2", "abs(d0)<0.2 and abs(z0)<1 and cosTheta >= -0.6", path=main)
ma.fillParticleList("K+", "abs(d0)<0.2 and abs(z0)<1 and cosTheta >= -0.6", path=main)

#D0 reconstruction
ma.reconstructDecay("anti-D0 -> K+ pi-:myPions2",cut="1.84 < M < 1.89",path=main)
ma.variablesToExtraInfo("anti-D0", variables={'M': 'M_before_fit'}, path=main)
vx.kFit("anti-D0", conf_level=0.0, fit_type='massvertex', path=main)
ma.matchMCTruth("anti-D0", path=main)

#B reconstruction
ma.reconstructDecay("B+ -> anti-D0 pi+:myPions2 pi-:myPions2 pi+:myPions2",cut="5.22 < Mbc < 5.3 and abs(deltaE) < 0.17",path=main)
ma.matchMCTruth("B+", path=main)

ma.buildRestOfEvent(target_list_name="B+", path=main)
cleanMask = ("cleanMask", "nCDCHits > 0 and useCMSFrame(p)<=3.2", "p >= 0.05 and useCMSFrame(p)<=3.2")
ma.appendROEMasks(list_name="B+", mask_tuples=[cleanMask], path=main)
ma.buildContinuumSuppression(list_name="B+", roe_mask="cleanMask", path=main) 


####################################################################################################################
####################################################################################################################
# ############### All veto (TOGETHER) starts here #########################
# roe_path = b2.create_path()                                                                                                                  
# deadEndPath = b2.create_path()                             
# # ma.signalSideParticleFilter('B+', '', roe_path, deadEndPath)  #??????????????????????

# ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path)     
# ma.fillSignalSideParticleList(outputListName='pi+:1st', decayString='B+ -> anti-D0 ^pi+ pi- pi+',path=roe_path)     
# ma.fillSignalSideParticleList(outputListName='pi-:2nd', decayString='B+ -> anti-D0 pi+ ^pi- pi+',path=roe_path)     
# ma.fillSignalSideParticleList(outputListName='pi+:3rd', decayString='B+ -> anti-D0 pi+ pi- ^pi+',path=roe_path)     

# ma.reconstructDecay(decayString='D*-:veto ->  anti-D0:sig pi-:2nd', cut='',path=roe_path)          
# ma.reconstructDecay(decayString='anti-D_10:veto ->  anti-D0:sig pi-:2nd pi+:3rd', cut='',path=roe_path)
# ma.reconstructDecay(decayString='D_s+:veto ->  pi+:1st pi-:2nd pi+:3rd', cut='',path=roe_path)
# ma.reconstructDecay(decayString='K_S0:veto ->  pi-:2nd pi+:3rd', cut='',path=roe_path)


# ma.rankByLowest(particleList='D*-:veto', variable='massDifference(0)', outputVariable='Dstrminus2D_md_rank', path=roe_path)
# ma.rankByLowest(particleList='anti-D_10:veto', variable='massDifference(0)', outputVariable='D_10D_md_rank', path=roe_path)
# ma.rankByLowest(particleList='D_s+:veto', variable='formula((InvM - 1.9683) **2 )', outputVariable='D_s_InvM_dif_rank', path=roe_path)
# ma.rankByLowest(particleList='K_S0:veto', variable='formula((InvM - 0.49761) **2 )', outputVariable='K_S0_23_InvM_dif_rank', path=roe_path)
# vm.addAlias('Dstrminus2D_md_rank', 'extraInfo(Dstrminus2D_md_rank)')
# vm.addAlias('D_10D_md_rank', 'extraInfo(D_10D_md_rank)')
# vm.addAlias('D_s_InvM_dif_rank', 'extraInfo(D_s_InvM_dif_rank)')
# vm.addAlias('K_S0_23_InvM_dif_rank', 'extraInfo(K_S0_23_InvM_dif_rank)')

# ma.variableToSignalSideExtraInfo(particleList='D*-:veto', varToExtraInfo={'massDifference(0)': 'Dstrminus2D_md_veto'}, path=roe_path)
# ma.variableToSignalSideExtraInfo(particleList='anti-D_10:veto', varToExtraInfo={'massDifference(0)': 'D_10D_md_veto'}, path=roe_path)
# ma.variableToSignalSideExtraInfo(particleList='D_s+:veto', varToExtraInfo={'InvM': 'D_s_InvM_veto'}, path=roe_path)
# ma.variableToSignalSideExtraInfo(particleList='K_S0:veto', varToExtraInfo={'InvM': 'K_S0_23_InvM_veto'}, path=roe_path)

# # execute roe_path for each RestOfEvent in the event                  
# # main.for_each('RestOfEvent', 'RestOfEvents', roe_path)
# vm.addAlias('DstrminusD_md', 'extraInfo(Dstrminus2D_md_veto)')
# vm.addAlias('D_10D_md', 'extraInfo(D_10D_md_veto)')
# vm.addAlias('D_s_InvM', 'extraInfo(D_s_InvM_veto)')
# vm.addAlias('K_S0_InvM_23', 'extraInfo(K_S0_23_InvM_veto)')

# veto_var = ['Dstrminus2D_md_rank', 'D_10D_md_rank', 'D_s_InvM_dif_rank', 'K_S0_23_InvM_dif_rank', 'DstrminusD_md', 'D_10D_md', 'D_s_InvM', 'K_S0_InvM_23']
# ################ All veto (TOGETHER) ends here ############################

############### Seperate veto starts here #########################
roe_path_Dstrminus2 = b2.create_path()                                                                                                                  
# deadEndPath_Dstr = b2.create_path()                             
# ma.signalSideParticleFilter('B+', '', roe_path_Dstrminus2, deadEndPath_Dstr)  #??????????????????????
ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path_Dstrminus2)
ma.fillSignalSideParticleList(outputListName='pi-:2nd', decayString='B+ -> anti-D0 pi+ ^pi- pi+',path=roe_path_Dstrminus2)     
ma.reconstructDecay(decayString='D*-:veto ->  anti-D0:sig pi-:2nd', cut='',path=roe_path_Dstrminus2)   
ma.rankByLowest(particleList='D*-:veto', variable='massDifference(0)', outputVariable='Dstrminus2D_md_rank', path=roe_path_Dstrminus2)
vm.addAlias('Dstrminus2D_md_rank', 'extraInfo(Dstrminus2D_md_rank)')
ma.variableToSignalSideExtraInfo(particleList='D*-:veto', varToExtraInfo={'massDifference(0)': 'Dstrminus2D_md_veto'}, path=roe_path_Dstrminus2)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_Dstrminus2)
vm.addAlias('Dstrminus2D_md', 'extraInfo(Dstrminus2D_md_veto)')

# roe_path_Dstrminus1 = b2.create_path()                                                                                                                  
# # deadEndPath_Dstr = b2.create_path()                             
# # ma.signalSideParticleFilter('B+', '', roe_path_Dstrminus1, deadEndPath_Dstr)  #??????????????????????
# ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path_Dstrminus1)
# ma.fillSignalSideParticleList(outputListName='pi+:1st', decayString='B+ -> anti-D0 ^pi+ pi- pi+',path=roe_path_Dstrminus1)     
# ma.reconstructDecay(decayString='D*-:veto ->  anti-D0:sig pi+:1st', cut='',path=roe_path_Dstrminus1)   
# ma.rankByLowest(particleList='D*-:veto', variable='massDifference(0)', outputVariable='Dstrminus1D_md_rank', path=roe_path_Dstrminus1)
# vm.addAlias('Dstrminus1D_md_rank', 'extraInfo(Dstrminus1D_md_rank)')
# ma.variableToSignalSideExtraInfo(particleList='D*-:veto', varToExtraInfo={'massDifference(0)': 'Dstrminus1D_md_veto'}, path=roe_path_Dstrminus1)
# main.for_each('RestOfEvent', 'RestOfEvents', roe_path_Dstrminus1)
# vm.addAlias('Dstrminus1D_md', 'extraInfo(Dstrminus1D_md_veto)')

# roe_path_Dstrminus3 = b2.create_path()                                                                                                                  
# # deadEndPath_Dstr = b2.create_path()                             
# # ma.signalSideParticleFilter('B+', '', roe_path_Dstrminus3, deadEndPath_Dstr)  #??????????????????????
# ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path_Dstrminus3)
# ma.fillSignalSideParticleList(outputListName='pi+:3rd', decayString='B+ -> anti-D0 pi+ pi- ^pi+',path=roe_path_Dstrminus3)     
# ma.reconstructDecay(decayString='D*-:veto ->  anti-D0:sig pi+:3rd', cut='',path=roe_path_Dstrminus3)   
# ma.rankByLowest(particleList='D*-:veto', variable='massDifference(0)', outputVariable='Dstrminus3D_md_rank', path=roe_path_Dstrminus3)
# vm.addAlias('Dstrminus3D_md_rank', 'extraInfo(Dstrminus3D_md_rank)')
# ma.variableToSignalSideExtraInfo(particleList='D*-:veto', varToExtraInfo={'massDifference(0)': 'Dstrminus3D_md_veto'}, path=roe_path_Dstrminus3)
# main.for_each('RestOfEvent', 'RestOfEvents', roe_path_Dstrminus3)
# vm.addAlias('Dstrminus3D_md', 'extraInfo(Dstrminus3D_md_veto)')

roe_path_D_10 = b2.create_path()                                                                                                                  
# deadEndPath_D_10 = b2.create_path()                             
# ma.signalSideParticleFilter('B+', '', roe_path_D_10, deadEndPath_D_10)  #??????????????????????
ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path_D_10)
# ma.fillSignalSideParticleList(outputListName='pi+:1st', decayString='B+ -> anti-D0 ^pi+ pi- pi+',path=roe_path)     # need to include this pion for veto also
ma.fillSignalSideParticleList(outputListName='pi-:2nd', decayString='B+ -> anti-D0 pi+ ^pi- pi+',path=roe_path_D_10)     
ma.fillSignalSideParticleList(outputListName='pi+:3rd', decayString='B+ -> anti-D0 pi+ pi- ^pi+',path=roe_path_D_10)  
ma.reconstructDecay(decayString='anti-D_10:veto ->  anti-D0:sig pi-:2nd pi+:3rd', cut='',path=roe_path_D_10)
ma.rankByLowest(particleList='anti-D_10:veto', variable='massDifference(0)', outputVariable='D_10D_md_rank', path=roe_path_D_10)
vm.addAlias('D_10D_md_rank', 'extraInfo(D_10D_md_rank)')
ma.variableToSignalSideExtraInfo(particleList='anti-D_10:veto', varToExtraInfo={'massDifference(0)': 'D_10D_md_veto'}, path=roe_path_D_10)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_D_10)
vm.addAlias('D_10D_md', 'extraInfo(D_10D_md_veto)')

roe_path_D_s = b2.create_path()                                                                                                                  
# deadEndPath_D_s = b2.create_path()   
# ma.signalSideParticleFilter('B+', '', roe_path_D_s, deadEndPath_D_s)  #??????????????????????
ma.fillSignalSideParticleList(outputListName='pi+:1st', decayString='B+ -> anti-D0 ^pi+ pi- pi+',path=roe_path_D_s)     
ma.fillSignalSideParticleList(outputListName='pi-:2nd', decayString='B+ -> anti-D0 pi+ ^pi- pi+',path=roe_path_D_s)   
ma.fillSignalSideParticleList(outputListName='pi+:3rd', decayString='B+ -> anti-D0 pi+ pi- ^pi+',path=roe_path_D_s)  
ma.reconstructDecay(decayString='D_s+:veto ->  pi+:1st pi-:2nd pi+:3rd', cut='',path=roe_path_D_s)
ma.rankByLowest(particleList='D_s+:veto', variable='formula((InvM - 1.9683) **2 )', outputVariable='D_s_InvM_dif_rank', path=roe_path_D_s)
vm.addAlias('D_s_InvM_dif_rank', 'extraInfo(D_s_InvM_dif_rank)')
ma.variableToSignalSideExtraInfo(particleList='D_s+:veto', varToExtraInfo={'InvM': 'D_s_InvM_veto'}, path=roe_path_D_s)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_D_s)
vm.addAlias('D_s_InvM', 'extraInfo(D_s_InvM_veto)')

roe_path_K_S0_23 = b2.create_path()                                                                                                                  
# deadEndPath_K_S0 = b2.create_path()                             
# ma.signalSideParticleFilter('B+', '', roe_path_K_S0_23, deadEndPath_K_S0)  #??????????????????????
# ma.fillSignalSideParticleList(outputListName='pi+:1st', decayString='B+ -> anti-D0 ^pi+ pi- pi+',path=roe_path_K_S0_23)     
ma.fillSignalSideParticleList(outputListName='pi-:2nd', decayString='B+ -> anti-D0 pi+ ^pi- pi+',path=roe_path_K_S0_23)     
ma.fillSignalSideParticleList(outputListName='pi+:3rd', decayString='B+ -> anti-D0 pi+ pi- ^pi+',path=roe_path_K_S0_23)    
ma.reconstructDecay(decayString='K_S0:veto ->  pi-:2nd pi+:3rd', cut='',path=roe_path_K_S0_23)
ma.rankByLowest(particleList='K_S0:veto', variable='formula((InvM - 0.49761) **2 )', outputVariable='K_S0_23_InvM_dif_rank', path=roe_path_K_S0_23)
vm.addAlias('K_S0_23_InvM_dif_rank', 'extraInfo(K_S0_23_InvM_dif_rank)')
ma.variableToSignalSideExtraInfo(particleList='K_S0:veto', varToExtraInfo={'InvM': 'K_S0_23_InvM_veto'}, path=roe_path_K_S0_23)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_K_S0_23)
vm.addAlias('K_S0_23_InvM', 'extraInfo(K_S0_23_InvM_veto)')

roe_path_K_S0_12 = b2.create_path()                                                                                                                  
# deadEndPath_K_S0 = b2.create_path()                             
# ma.signalSideParticleFilter('B+', '', roe_path_K_S0_12, deadEndPath_K_S0)  #??????????????????????
# ma.fillSignalSideParticleList(outputListName='pi+:1st', decayString='B+ -> anti-D0 ^pi+ pi- pi+',path=roe_path_K_S0_12)     
ma.fillSignalSideParticleList(outputListName='pi-:2nd', decayString='B+ -> anti-D0 pi+ ^pi- pi+',path=roe_path_K_S0_12)     
ma.fillSignalSideParticleList(outputListName='pi+:1st', decayString='B+ -> anti-D0 ^pi+ pi- pi+',path=roe_path_K_S0_12)    
ma.reconstructDecay(decayString='K_S0:veto -> pi+:1st pi-:2nd', cut='',path=roe_path_K_S0_12)
ma.rankByLowest(particleList='K_S0:veto', variable='formula((InvM - 0.49761) **2 )', outputVariable='K_S0_12_InvM_dif_rank', path=roe_path_K_S0_12)
vm.addAlias('K_S0_12_InvM_dif_rank', 'extraInfo(K_S0_12_InvM_dif_rank)')
ma.variableToSignalSideExtraInfo(particleList='K_S0:veto', varToExtraInfo={'InvM': 'K_S0_12_InvM_veto'}, path=roe_path_K_S0_12)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_K_S0_12)
vm.addAlias('K_S0_12_InvM', 'extraInfo(K_S0_12_InvM_veto)')

############### anti-D0* (gamma) veto starts here #########################
roe_path_Dstr0_gamma = b2.create_path()                                                                                                                  
deadEndPath_Dstr_gamma = b2.create_path()                             
ma.signalSideParticleFilter('B+', '', roe_path_Dstr0_gamma, deadEndPath_Dstr_gamma)
vm.addAlias("goodFWDGamma", "passesCut(clusterReg == 1 and clusterE > 0.075)")
vm.addAlias("goodBRLGamma", "passesCut(clusterReg == 2 and clusterE > 0.05)")
vm.addAlias("goodBWDGamma", "passesCut(clusterReg == 3 and clusterE > 0.1)")
vm.addAlias("goodGamma", "passesCut(goodFWDGamma or goodBRLGamma or goodBWDGamma)")
ma.fillParticleList(decayString="gamma:roe", cut="isInRestOfEvent == 1 and goodGamma", path=roe_path_Dstr0_gamma)
ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path_Dstr0_gamma)     
ma.reconstructDecay(decayString='anti-D*0:veto ->  anti-D0:sig gamma:roe', cut='',path=roe_path_Dstr0_gamma)          
ma.rankByLowest(particleList='anti-D*0:veto', variable='massDifference(0)',numBest=1, path=roe_path_Dstr0_gamma)
# vm.addAlias('D_str0_gammaD_md_rank', 'extraInfo(D_str0_gammaD_md_rank)')
ma.variableToSignalSideExtraInfo(particleList='anti-D*0:veto', varToExtraInfo={'massDifference(0)': 'D_str0_gammaD_md'}, path=roe_path_Dstr0_gamma)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_Dstr0_gamma)
vm.addAlias('D_str0_gammaD_md', 'extraInfo(D_str0_gammaD_md)')
#ma.printVariableValues('B+', ['dm'], path=main)
#ma.applyCuts('B+','dm > 0.15',path=main) #############COMPLIMENT
################ anti-D0* (gamma) veto ends here ############################

############### anti-D0* (pi0) veto starts here #########################
roe_path_Dstr0_pi0 = b2.create_path()                                                                                                                  
deadEndPath_Dstr0_pi0 = b2.create_path()                             
ma.signalSideParticleFilter('B+', '', roe_path_Dstr0_pi0, deadEndPath_Dstr0_pi0)
stdPi0s.stdPi0s(path=roe_path_Dstr0_pi0)  # ???????????cut="isInRestOfEvent == 1 ???????????
ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path_Dstr0_pi0)     
ma.reconstructDecay(decayString='anti-D*0:veto ->  anti-D0:sig pi0:eff60_May2020', cut='',path=roe_path_Dstr0_pi0)          
ma.rankByLowest(particleList='anti-D*0:veto', variable='massDifference(0)',numBest=1, path=roe_path_Dstr0_pi0)
# vm.addAlias('D_str0_pi0D_md_rank', 'extraInfo(D_str0_pi0D_md_rank)')
ma.variableToSignalSideExtraInfo(particleList='anti-D*0:veto', varToExtraInfo={'massDifference(0)': 'D_str0_pi0D_md'}, path=roe_path_Dstr0_pi0)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_Dstr0_pi0)
vm.addAlias('D_str0_pi0D_md', 'extraInfo(D_str0_pi0D_md)')
#ma.printVariableValues('B+', ['dm'], path=main)
#ma.applyCuts('B+','dm > 0.15',path=main) #############COMPLIMENT
################ anti-D0* (pi0) veto ends here ############################

############### D*- veto starts here #########################
roe_path_Dstrminus2_roe = b2.create_path()                                                                                                                  
deadEndPath_Dstrminus_roe = b2.create_path()                             
ma.signalSideParticleFilter('B+', '', roe_path_Dstrminus2_roe, deadEndPath_Dstrminus_roe)
ma.fillParticleList(decayString="pi-:roe", cut="isInRestOfEvent == 1 and  abs(d0)<0.2 and abs(z0)<1 and cosTheta >= -0.6", path=roe_path_Dstrminus2_roe)
ma.fillSignalSideParticleList(outputListName='anti-D0:sig', decayString='B+ -> ^anti-D0 pi+ pi- pi+',path=roe_path_Dstrminus2_roe)     

ma.reconstructDecay(decayString='D*-:veto ->  anti-D0:sig pi-:roe', cut='',path=roe_path_Dstrminus2_roe)
ma.rankByLowest(particleList='D*-:veto', variable='massDifference(0)',numBest=1, path=roe_path_Dstrminus2_roe)
# vm.addAlias('DstrminusroeD_md_rank', 'extraInfo(DstrminusroeD_md_rank)')
ma.variableToSignalSideExtraInfo(particleList='D*-:veto', varToExtraInfo={'massDifference(0)': 'DstrminusroeD_md'}, path=roe_path_Dstrminus2_roe)
main.for_each('RestOfEvent', 'RestOfEvents', roe_path_Dstrminus2_roe)
vm.addAlias('DstrminusroeD_md', 'extraInfo(DstrminusroeD_md)')
#ma.printVariableValues('B+', ['dm'], path=main)
#ma.applyCuts('B+','dm > 0.15',path=main) #############COMPLIMENT
################ D*- veto ends here ############################




veto_var = ['Dstrminus2D_md_rank', 'D_10D_md_rank', 'D_s_InvM_dif_rank',
            # 'Dstrminus1D_md_rank', 'Dstrminus3D_md_rank', 
            'K_S0_23_InvM_dif_rank','K_S0_12_InvM_dif_rank',
            'Dstrminus2D_md', 'D_10D_md', 'D_s_InvM', 'K_S0_23_InvM', 'K_S0_12_InvM',
            # 'Dstrminus1D_md', 'Dstrminus3D_md', 
            'D_str0_gammaD_md','D_str0_pi0D_md', 'DstrminusroeD_md']
################ Seperate veto ends here ############################
####################################################################################################################
####################################################################################################################

vx.kFit("B+", conf_level=0.0, fit_type='vertex', path=main)
vx.TagV("B+", 'breco', path=main)
ft.flavorTagger(particleLists=["B+"],weightFiles='B2nunubarBGx1',combinerMethods=['TMVA-FBDT'], path=main)
vm.addAlias("abs_qr","abs(qrOutput(FBDT))")

ma.buildEventShape(inputListNames=["B+"], path=main)

main.add_module("MVAExpert", listNames=["B+"], extraInfoName="SignalProbability", identifier=f"{weightfile}")  # name of the weightfile used

############################### Best candidate selection Start#########################
vm.addAlias('BCS3', 'formula(((Mbc - 5.27933)/0.0026) ** 2  + ((daughter(0,InvM) - 1.86483)/0.0040276) ** 2)')
ma.rankByLowest(particleList='B+', variable='BCS3', outputVariable='BCS3_rank', path=main)
vm.addAlias('BCS3_rank', 'extraInfo(BCS3_rank)')
# ma.applyCuts('B+','BCS3_rank==1',path=main)
############################### Best candidate selection Stop#########################

##############################calculation of invariant mass##########################
vm.addAlias('InvMD0and1stpi', 'daughterInvM(0, 1)')
vm.addAlias('InvMD0and2ndpi', 'daughterInvM(0, 2)')
vm.addAlias('InvMD0and3rdpi', 'daughterInvM(0, 3)')
vm.addAlias('InvM1stand2ndpi', 'daughterInvM(1, 2)')
vm.addAlias('InvM1stand3rdpi', 'daughterInvM(1, 3)')
vm.addAlias('InvM2ndand3rdpi', 'daughterInvM(2, 3)')
vm.addAlias('InvM1st2nd3rdpi', 'daughterInvM(1, 2, 3)')
vm.addAlias('InvMD0and1st2ndpi', 'daughterInvM(0, 1, 2)')
vm.addAlias('InvMD0and2nd3rdpi', 'daughterInvM(0, 2, 3)')
vm.addAlias('InvMD0and1st3rdpi', 'daughterInvM(0, 1, 3)')

invM_var = ['InvMD0and1stpi', 'InvMD0and2ndpi', 'InvMD0and3rdpi', 'InvM1stand2ndpi', 'InvM1stand3rdpi',
    'InvM2ndand3rdpi', 'InvM1st2nd3rdpi', 'InvMD0and1st2ndpi', 'InvMD0and2nd3rdpi', 'InvMD0and1st3rdpi']
##############################calculation of invariant mass##########################

############################## Calculation of angle between daughteres#############
vm.addAlias('AngelD0and1stpi', 'daughterAngle(0, 1)')
vm.addAlias('AngelD0and2ndpi', 'daughterAngle(0, 2)')
vm.addAlias('AngelD0and3rdpi', 'daughterAngle(0, 3)')
vm.addAlias('Angel1stand2ndpi', 'daughterAngle(1, 2)')
vm.addAlias('Angel1stand3rdpi', 'daughterAngle(1, 3)')
vm.addAlias('Angel2ndand3rdpi', 'daughterAngle(2, 3)')

vm.addAlias('AngelKpi', 'daughterAngle(0:0, 0:1)')

Angel_var = ['AngelD0and1stpi', 'AngelD0and2ndpi', 'AngelD0and3rdpi', 'Angel1stand2ndpi', 'Angel1stand3rdpi',
    'Angel2ndand3rdpi', 'AngelKpi']
############################## Calculation of angle between daughteres#############

# Create list of variables to save into the output file
vm.addAlias('ethoo0', 'formula((KSFWVariables(et) - KSFWVariables(hoo0))/(KSFWVariables(et) + KSFWVariables(hoo0)))')
vm.addAlias('hso00mm2', 'formula((KSFWVariables(hso00) - KSFWVariables(mm2))/(KSFWVariables(hso00) + KSFWVariables(mm2)))')
vm.addAlias('hoo3hoo1', 'formula((KSFWVariables(hoo3) + KSFWVariables(hoo1)))')
vm.addAlias('hso01hso03', 'formula((KSFWVariables(hso01) + KSFWVariables(hso03)))')
vm.addAlias('ethoo0TBz', 'formula((ethoo0 - cosTBz)/(ethoo0 + cosTBz))')
vm.addAlias('ethoo0hso00mm2', 'formula((ethoo0 - hso00mm2)/(ethoo0 + hso00mm2))')
vm.addAlias('hso12hso02hoo2', 'formula((KSFWVariables(hso12) + KSFWVariables(hso02) + KSFWVariables(hoo2)))')
vm.addAlias('hso12hso02hoo2thrustOm_lc', 'formula((hso12hso02hoo2 - thrustOm)/(hso12hso02hoo2 + thrustOm))')
vm.addAlias('hso12hso02hoo2thrustOm_p', 'formula((hso12hso02hoo2 + thrustOm))')

###################################R2
vm.addAlias('R2cosTBTO_lc', 'formula((R2 - cosTBTO)/(R2 + cosTBTO))')
vm.addAlias('R2cosTBTO_p', 'formula((R2 + cosTBTO))')
vm.addAlias('R2thrustBm_lc', 'formula((R2 - thrustBm)/(R2 + thrustBm))')
vm.addAlias('R2thrustBm_p', 'formula((R2 + thrustBm))')
vm.addAlias('TBTOhso12', 'formula((cosTBTO - KSFWVariables(hso12))/(cosTBTO + KSFWVariables(hso12)))')
vm.addAlias('thrustOmhoo2', 'formula((thrustOm - KSFWVariables(hoo2))/(thrustOm + KSFWVariables(hoo2)))')

vm.addAlias('R2thrustBmTBTOhso12', 'formula((R2thrustBm_lc - TBTOhso12)/(R2thrustBm_lc + TBTOhso12))')
vm.addAlias('R2thrustBmthrustOmhoo2', 'formula((R2thrustBm_lc - thrustOmhoo2)/(R2thrustBm_lc + thrustOmhoo2))')
vm.addAlias('TBTOhso12thrustOmhoo2', 'formula((TBTOhso12 - thrustOmhoo2)/(TBTOhso12 + thrustOmhoo2))')

vm.addAlias('R2thrustBmTBTOhso12thrustOmhoo2', 'formula((R2thrustBmTBTOhso12 - thrustOmhoo2)/(R2thrustBmTBTOhso12 + thrustOmhoo2))')

##########################R2

vm.addAlias('CMS_cosTheta', 'useCMSFrame(cosTheta)')
simpleCSVariables = [
    "ethoo0",
    "hso00mm2",
    "ethoo0hso00mm2",
    "hoo3hoo1",
    "hso01hso03",
    "ethoo0TBz",
    "hso12hso02hoo2",
    "hso12hso02hoo2thrustOm_lc",
    "hso12hso02hoo2thrustOm_p",
    "R2cosTBTO_lc",
    "R2cosTBTO_p",
    "R2thrustBm_lc",
    "R2thrustBm_p",
    "TBTOhso12",
    "thrustOmhoo2",
    "R2thrustBmTBTOhso12",
    "R2thrustBmthrustOmhoo2",
    "TBTOhso12thrustOmhoo2",
    "R2thrustBmTBTOhso12thrustOmhoo2",
    "abs_qr",       # calculated using ft.flavorTagger fn
    "sphericity",   # calculated using ma.buildEventShape fn
    "DeltaZ",       # calculated using vx.TagV fn
    "R2",           # rest of the variables are calculated using ma.buildContinuumSuppression fn
    "thrustBm",  #must    #significance
    "thrustOm",  #must    #significance
    "cosTBTO",  #must    #significance
    "cosTBz",  #must    #significance
    "CMS_cosTheta",  #must    #significance
    "KSFWVariables(et)",    #significance
    "KSFWVariables(mm2)",
    "KSFWVariables(hso00)",
    "KSFWVariables(hso01)",
    "KSFWVariables(hso02)",  #must    #significance     #############nearly_same##########
    # "KSFWVariables(hso03)",
    "KSFWVariables(hso04)",    #significance
    "KSFWVariables(hso10)",    #significance
    "KSFWVariables(hso12)",  #must    #significance     #############nearly_same##########
    "KSFWVariables(hso14)",    #significance
    # "KSFWVariables(hso20)",    #significance
    # "KSFWVariables(hso22)",
    # "KSFWVariables(hso24)",
    "KSFWVariables(hoo0)",  #must    #significance
    # "KSFWVariables(hoo1)",
    "KSFWVariables(hoo2)",    #significance
    # "KSFWVariables(hoo3)",
    # "KSFWVariables(hoo4)",    #significance
    # "CleoConeCS(1)",
    # "CleoConeCS(2)",
    # "CleoConeCS(3)",
    # "CleoConeCS(4)",
    # "CleoConeCS(5)",
    # "CleoConeCS(6)",
    # "CleoConeCS(7)",
    # "CleoConeCS(8)",
    # "CleoConeCS(9)",
]
vm.addAlias("D0M_BF", "extraInfo(M_before_fit)")
vm.addAlias("SigProb", "extraInfo(SignalProbability)")
vm.addAlias('PID_bin_kaon', 'ifNANgiveX(pidPairProbabilityExpert(321, 211, ALL), 0.5)')
vm.addAlias('chiSqrd', 'extraInfo(chiSquared)')
vm.addAlias('fitNdf', 'extraInfo(ndf)')


b_vars = []
bcs_var = ['BCS3_rank']
b_vars  += bcs_var
b_vars  += invM_var
b_vars  += Angel_var
b_vars += veto_var

other_var = ['p', 'E', 'isSignal', 'M', 'InvM', 'Mbc', 'deltaE', 'chiProb']
b_vars += other_var
b_vars += simpleCSVariables
track_variables = ['PID_bin_kaon', 'pionID', 'kaonID']
b_vars += vu.create_aliases_for_selected(
    track_variables,
    "B+ -> [anti-D0 -> ^K+ pi-] pi+ pi- pi+",
    prefix=["Kp"],
)

# D_var = vc.inv_mass + vc.deltae_mbc + vc.mc_truth + vc.vertex + vc.mc_vertex + ['chiSqrd', 'fitNdf']
D_var = ['M', 'D0M_BF', 'InvM', 'Mbc', 'deltaE', 'isSignal', 'dr', 'dz',
        'mcDecayVertexX', 'mcDecayVertexY', 'mcDecayVertexZ', 'x', 'y', 'z',
        'chiProb', 'chiSqrd', 'fitNdf']
b_vars += vu.create_aliases_for_selected(
    D_var,
    "B+ -> [^anti-D0 -> K+ pi-] pi+ pi- pi+",
    prefix=["D0_bar"],
)


ma.variablesToNtuple(
    decayString="B+",
    variables=b_vars + vmc.mc_gen_topo(200) +["SigProb", "isContinuumEvent"],  # + vmc.mc_gen_topo(200) 
    filename=f"{output_foldername}/{output_filename}",
    treename="tree",
    path=main,
)

b2.process(main)

# print out the summary
print(b2.statistics)