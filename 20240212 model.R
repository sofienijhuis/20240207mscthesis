#Purpose: PBK Model carbendazim, built with in vitro derived parameter values
#Species: Human (gender mixed) and rat (male, Sprague-Dawley)
#Compiled by: Sofie Nijhuis
#Organisation: Wageningen University


#Note: to switch from human to rat by commenting/uncommenting all parts that start with #HUMAN or #RAT
#-------------------------------------------------------------------------------------------------------------------------------------------
library(patchwork)
library(rxode2)
library(tidyverse)
library(ggplot2)

#### Absorption calculations #####
#Absorption
TPSACBZ = 67.01      #in Å^2 ref: Chemicalize (Chemaxon), 01/2024
TPSAHBC = 87.24      #in Å^2 ref: Chemicalize (Chemaxon), 01/2024
TPSAAB = 54.7        #in Å^2 ref: Chemicalize (Chemaxon), 01/2024

LogPappCBZ = -4.36 - 0.01*TPSACBZ            #in cm/s
LogPappHBC = -4.36 - 0.01*TPSAHBC
LogPappAB = -4.36 - 0.01*TPSAAB

LogPappCBZ106 = log10((10^LogPappCBZ) * 10^6) #in 10^-6 cm/s


LogPeffCBZ = 0.4926*LogPappCBZ106-0.1454   #in 10^-4 cm/s
LogPeffHBC = 0.4926*LogPappHBC*10^-6-0.1454
LogPeffAB = 0.4926*LogPappAB*10^-6-0.1454

PeffCBZ = 10^LogPeffCBZ                      #in 10-4 cm/s
PeffHBC = 10^LogPeffHBC
PeffAB = 10^LogPeffAB
# #HUMAN
R = 1.75 # intestinal radius in cm ref: Yu et al. (1999)
Tsi = 3.32 #intestinal transit time in h ref: Punt et al. (2022)

#RAT
# R = 0.18 # ref: Punt et al. (2022)
# Tsi = 88 / 60 #intestinal transit time in h ref: Grandoni et al. (2019), supplementary material

kaCBZ = (PeffCBZ*10^(-4))*2/R*3600
kaHBC =  PeffHBC*10^-4*2/R*3600
kaAB  = PeffAB*10^-4*2/R*3600




faCBZ = 0.825 #Sharma et al. (2022)

##### Initialization #####

# INITIALIZATION - Initialize values and set to 0
InitModel <- function(){
  inits <<- c(
    Ast = 0,
    AFCBZ = 0,
    AFHBC = 0,
    AFAB = 0,
    ALCBZ = 0,
    ALHBC = 0,
    ALAB = 0,
    AMHBCL = 0,
    AMABL = 0,
    ABCBZ = 0,
    ABHBC = 0,
    ABAB = 0,
    ARCBZ = 0,
    ARHBC = 0,
    ARAB = 0,
    ASCBZ = 0,
    ASHBC = 0,
    ASAB = 0,
    ABlCBZ = 0,
    ABlHBC = 0,
    ABlAB = 0,
    AKCBZ = 0,
    AKHBC = 0,
    AKAB = 0
    
  )
#####  #####
##### Model  #####
  
Model <<- RxODE({
  ##### Physiological  #####
  #HUMAN
  #Physiological parameters

  #Tissues reference: Brown et al, (1997)
  BW=83.3;              #body weight Caucasian in kg
  VBc = 0.02            #fraction of brain tissue unitless
  VFc=0.214             #fraction of fat tissue unitless
  VLc=0.026             #fraction of liver tissue unitless
  VKc = 0.0044	        #fraction of kidney tissue unitless
  VBlc=0.079            #fraction of blood unitless
  VRc=0.032             #fraction of richly perfused tissue  (GI-tract, Heart, Lungs,spleen) unitless
  VSc=0.56              #fraction of slowly perfused tissue (Bone,skin,Muscle) unitless


  #Cardiac parameters reference: Brown et al. (1997)
  COf = 80 / 1000 * 60  #Cardiac output in L/h/kg ref: Walton et al. (2004)
  CO = COf * BW         #cardiac output in L/h
  QBc = 0.114		        #fraction of blood flow to brain
  QFc = 0.052		        #fraction of blood flow to fat	unitless
  QLc = 0.227		        #fraction of blood flow to liver	unitless
  QKc = 0.175		        #fraction of blood flow to kidney unitless
  QRc = 0.43		        #fraction of blood flow to richly perfused tissue (1-rest) unitless
  QSc = 0.291 		      #fraction of blood flow to slowly perfused tissue (muscle, skin, bone) unitless

  #Renal excretion
  fGFR = 1.8 / 1000 * 60 #Glomerular filtration rate in L/h/kg ref: Walton et al. (2004)

# #RAT
# #Tissues reference: Brown et al, (1997)
# BW=0.185               #body weight kg
# VBc = 0.0057          #fraction of brain tissue unitless
# VFc=0.07              #fraction of fat tissue unitless
# VLc=0.037             #fraction of liver tissue unitless
# VKc = 0.8859	        #fraction of kidney tissue unitless
# VBlc=0.074            #fraction of blood unitless
# VRc=0.037             #fraction of richly perfused tissue  (GI-tract, Heart, Lungs,spleen) unitless
# VSc=0.6546            #fraction of slowly perfused tissue (Bone,skin,Muscle) unitless
# 
# 
# #Cardiac parameters reference: Brown et al. (1997)
# COf = 300 / 1000 * 60 #Cardiac output in L/h/kg ref: Walton et al. (2004)
# CO = BW * COf         #cardiac output in L/h
# QBc = 0.02		        #fraction of blood flow to brain
# QFc = 0.07		        #fraction of blood flow to fat	unitless
# QLc = 0.183		        #fraction of blood flow to liver	unitless
# QKc = 0.141		        #fraction of blood flow to kidney unitless
# QRc = 0.128		        #fraction of blood flow to richly perfused tissue (1-rest) unitless
# QSc = 0.458 		      #fraction of blood flow to slowly perfused tissue (muscle, skin, bone) unitless
# 
#   #Renal excretion
#   fGFR = 5.2 / 1000 * 60 #Glomerular filtration rate in L/h/kg ref: Walton et al. (2004)
# 
  # #Calculatd physiological parameters
  #
  VB = VBc*BW		        #volume of  brain tissue (calculated)
  VF = VFc*BW	          #volume of fat tissue (calculated)
  VL = VLc*BW		        #volume of liver tissue (calculated)
  VK = VKc*BW		        #volume of kidney tissue (calculated)
  VBl = VBlc*BW         #volume of blood (calculated)
  VR = VRc*BW	 	        #volume of  richly perfused tissue (calculated)
  VS = VSc*BW		        #volume of  slowly perfused tissue (calculated)

  QB = QBc*CO			      #blood flow to brain tissue (calculated) in L/h
  QF = QFc*CO			      #blood flow to fat tissue (calculated) in L/h
  QL = QLc*CO		        #blood flow to liver tissue (calculated) in L/h
  QK = QKc*CO			      #blood flow to kidneys (calculated) in L/h
  QS = QSc*CO 			    #blood flow to slowly perfused tissue (calculated) in L/h
  QR = QRc*CO 			    #blood flow to richly perfused tissue (calculated) in L/h


  #-------------------------------------------------------------------------------------------------------------------------------------------
  
  ##### Physicochemical  #####
  #Physicochemical parameters
  #HUMAN
  #Partition coefficients carbendazim (tissue:plasma)
  #Zwitterionic, pKa(acid) = 10, pKa(base) = 4.2, log P (pH 7 & 9) = 1.5, Mw = 191.21 g/mol
  PBCBZ = 0.81		      #brain/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PFCBZ = 0.09		      #fat/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PLCBZ = 0.56		      #liver/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PKCBZ =	0.58		      #kidney/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PRCBZ = 0.74  		    #richly perfused tissue/plasma partition coefficient of carbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  PSCBZ = 0.49  		    #slowly perfused tissue/plasma partition coefficient of carbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  FUPCBZ = 0.098		    #fraction unbound plasma calculated using QPPR of Lobell & Sivarajah (2003)
  BPRCBZ = 1

  #Partition coefficients 5-hydroxycarbendazim (tissue:plasma)
  #Zwitterionic, pKa(acid) = 8.5, pKa(base) = 5.5, log P = 1.45, Mw = 207.19 g/mol
  PBHBC = 0.73		      #brain/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PFHBC = 0.09   		    #fat/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PLHBC = 0.52		      #liver/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PKHBC = 0.55		      #kidney/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PRHBC = 0.68		      #richly perfused tissue/plasma partition coefficient of 5-hydroxycarbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  PSHBC = 0.46   		    #richly perfused tissue/plasma partition coefficient of 5-hydroxycarbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  FUPHBC = 0.402		    #fraction unbound plasma calculated using QPPR of Lobell & Sivarajah (2003)
  BPRHBC = 1

  #partition coefficients 2-aminobenzimidazole (tissue:plasma)
  #Diprotic acid, pKa1(acid) = 5.7, pKa2(acid) = 6.9, log P = 1.04, Mw = 133.15 g/mol

  PBAB = 0.08		        #brain/plasma partition coefficient of 2-aminobenzimidazole calculated using QPPR of Rodgers & Rowland (2006)
  PFAB = 0.06		        #fat/plasma partition coefficient of 2-aminobenzimidazole calculated using QPPR of Rodgers & Rowland (2006)
  PLAB = 0.11		        #liver/plasma partition coefficient of 2-aminobenzimidazole calculated using QPPR of Rodgers & Rowland (2006)
  PKAB = 0.16		        #kidney/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  PRAB = 0.19 		      #richly perfused tissue/plasma partition coefficient of 2-aminobenzimidazole (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  PSAB = 0.08   		    #slowly perfused tissue/plasma partition coefficient of 2-aminobenzimidazole (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  FUPAB = 0.139		      #fraction unbound plasma calculated using QPPR of Lobell & Sivarajah (2003)
  BPRAB = 1

  
  # #RAT
  # 
  # #Partition coefficients carbendazim (tissue:plasma)
  # #Zwitterionic, pKa(acid) = 10, pKa(base) = 4.2, log P (pH 7 & 9) = 1.5, Mw = 191.21 g/mol
  # PBCBZ = 0.81		      #brain/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PFCBZ = 0.09		      #fat/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PLCBZ = 0.56		      #liver/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PKCBZ =	0.58		      #kidney/plasma partition coefficient of carbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PRCBZ = 0.80  		    #richly perfused tissue/plasma partition coefficient of carbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  # PSCBZ = 0.67  		    #slowly perfused tissue/plasma partition coefficient of carbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  # FUPCBZ = 0.098		    #fraction unbound plasma calculated using QPPR of Lobell & Sivarajah (2003)
  # BPRCBZ = 1
  # 
  # #Partition coefficients 5-hydroxycarbendazim (tissue:plasma)
  # #70% glucuronated, 30% sulfonated
  # #Zwitterionic, pKa(acid) = 8.5, pKa(base) = 5.5, log P = 1.45, Mw = 207.19 g/mol
  # PBHBC = 0.73		      #brain/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PFHBC = 0.083   		    #fat/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PLHBC = 0.19		      #liver/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PKHBC = 0.25		      #kidney/plasma partition coefficient of 5-hydroxycarbendazim calculated using QPPR of Rodgers & Rowland (2006)
  # PRHBC = 0.27		      #richly perfused tissue/plasma partition coefficient of 5-hydroxycarbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  # PSHBC = 0.23   		    #slowly perfused tissue/plasma partition coefficient of 5-hydroxycarbendazim (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  # FUPHBC = 0.38		    #fraction unbound plasma calculated using QPPR of Lobell & Sivarajah (2003)
  # BPRHBC = 0.55
  # 
  # #partition coefficients 2-aminobenzimidazole (tissue:plasma)
  # #Diprotic acid, pKa1(acid) = 5.7, pKa2(acid) = 6.9, log P = 1.04, Mw = 133.15 g/mol
  # 
  # PBAB = 0.08		        #brain/plasma partition coefficient of 2-aminobenzimidazole calculated using QPPR of Rodgers & Rowland (2006)
  # PFAB = 0.06		        #fat/plasma partition coefficient of 2-aminobenzimidazole calculated using QPPR of Rodgers & Rowland (2006)
  # PLAB = 0.11		        #liver/plasma partition coefficient of 2-aminobenzimidazole calculated using QPPR of Rodgers & Rowland (2006)
  # PKAB = 0.16		        #kidney/plasma partition coefficient of 2-aminobenzimidazole calculated using QPPR of Rodgers & Rowland (2006)
  # PRAB = 0.187		      #richly perfused tissue/plasma partition coefficient of 2-aminobenzimidazole (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  # PSAB = 0.0882  		    #richly perfused tissue/plasma partition coefficient of 2-aminobenzimidazole (weighted average) calculated using QPPR of Rodgers & Rowland (2006)
  # FUPAB = 0.139		      #fraction unbound plasma calculated using QPPR of Lobell & Sivarajah (2003)
  # BPRAB = 1



##### Kinetic  #####
  #Kinetic parameters 
  
  #Renal excretion
  GFR = BW * fGFR #Glomerular filtration rate
  

  #HUMAN
  #Metabolic parameters CBZ -> HBC
  MPL = 32                                                          #scaling factor of human liver microsome in mg microsomal protein /g liver, Al-Malahmeh, A. J et al.,(2017)
  ClecHBC= 4.695                                                    #in uL/min/mg microsomal protein, data derived from experiment
  CleHBC = ClecHBC *60 * MPL * 1000 * VL /1000 /1000              #in L/h,/1000 to correct for pmol/umol, 60 -> min to hour, MPL -> mg microsomal protein /g liver,1000-> g to kg = L, VL = volume liver
  CleAB = 0
  # #RAT
  # MPL = 35                                                         #scaling factor of rat liver microsome in mg microsomal protein /g liver, Al-Malahmeh, A. J et al.,(2017)
  # VMax1cHBC= 181.3                                 	               #in pmol/min/mg microsomal protein, data derived from experiment
  # VMax1HBC = VMax1cHBC/1000000*60*MPL*VL*1000                      #in µmol/hr/liver, 1000000 -> pmol to umol, 60 -> min to hour, MPL -> mg microsomal protein /g liver, 1000 -> g to kg (liver), VL -> volume of liver tissue (VLc*BW)
  # Km1HBC = 8.87	                                                   #in µmol/L data derived from  experiment
  # #Metabolic parameters CBZ -> AB
  # VMax2cAB = 0	                                                   #in pmol/min/mg microsomal protein, data derived from experiment
  # VMax2AB = VMax2cAB/1000000*60*MPL*VL*1000                        #in µmol/hr/liver, 1000000 -> pmol to umol, 60 -> min to hour, MPL -> mg microsomal protein /g liver, 1000 -> g to kg (liver), VL -> volume of liver tissue (VLc*BW)
  # Km2AB = 0                                                        #in µmol/L data derived from  experiment

##### Compartments #####
  

  #fat compartment
  #Concentrations of CBZ, HBC and AB in fat and venous blood from fat
  CFCBZ = AFCBZ/VF
  CVFCBZ = (CFCBZ/PFCBZ)*BPRCBZ
    
  CFHBC = AFHBC/VF
  CVFHBC = (CFHBC/PFHBC)*BPRHBC
  
  CFAB = AFAB/VF
  CVFAB = (CFAB/PFAB)*BPRAB
  
  #liver compartment
  #Concentrations of CBZ, HBC and AB in liver and venous blood from liver
  CLCBZ = ALCBZ/VL
  CVLCBZ = (CLCBZ/PLCBZ)*BPRCBZ
  
  CLHBC = ALHBC/VL
  CVLHBC = (CLHBC/PLHBC)*BPRHBC
  
  CLAB = ALAB/VL
  CVLAB = (CLAB/PLAB)*BPRAB
  
  #Metabolism in liver
  
  #HUMAN
  #Amount metabolized CBZ -> 5HBC
  AMHBCL =  CleHBC * CVLCBZ
  AMABL = CleAB * CVLCBZ
  # #RAT
  # #Amount metabolized CBZ -> 5HBC
  # AMHBCL = (VMax1HBC*CVLCBZ)/(Km1HBC + CVLCBZ)
  # 
  # #Amount metabolized CBZ -> 2AB
  # AMABL = (VMax2AB*CVLCBZ)/(Km2AB + CVLCBZ)
  # #Metabolic parameters CBZ -> AB in blood
  # VMax3cAB = 0
  # VMax3AB = VMax3cAB/1000000*60*MPL*VL*1000
  # Km3AB = 0
  #
  #brain compartment
  #Concentrations of CBZ, HBC and AB in brain and venous blood from brain
  CBCBZ = ABCBZ/VB
  CVBCBZ = (CBCBZ/PBCBZ)*BPRCBZ
  
  CBHBC = ABHBC/VB
  CVBHBC = (CBHBC/PBHBC)*BPRHBC
  
  CBAB = ABAB/VB
  CVBAB = (CBAB/PBAB)*BPRAB
  
  #tissue compartment richly perfused tissue
  #Concentrations of CBZ, HBC and AB in richly perfused tissues and venous blood from richly perfused tissues
  CRCBZ = ARCBZ/VR
  CVRCBZ = (CRCBZ/PRCBZ)*BPRCBZ
  
  CRHBC = ARHBC/VR
  CVRHBC = (CRHBC/PRHBC)*BPRHBC
  
  CRAB = ARAB/VR
  CVRAB = (CRAB/PRAB)*BPRAB
  
  #tissue compartment slowly perfused tissue
  #Concentrations of CBZ, HBC and AB in slowly perfused tissues and venous blood from slowly perfused tissues
  CSCBZ = ASCBZ/VS
  
  CVSCBZ = (CSCBZ/PSCBZ)*BPRCBZ
  
  CSHBC = ASHBC/VS
  CVSHBC = (CSHBC/PSHBC)*BPRHBC
  
  CSAB = ASAB/VS
  CVSAB = (CSAB/PSAB)*BPRAB
  
  #blood compartment
  #Concentrations of CBZ, HBC and AB in blood
  CBlCBZ = ABlCBZ/VBl
  CBlHBC = ABlHBC/VBl
  CBlAB = ABlAB/VBl
  
  #Concentration in plasma
  CPlCBZ = CBlCBZ/BPRCBZ
  CPlHBC = CBlHBC/BPRHBC
  CPlAB = CBlAB/BPRAB
  
  #unbound fractions in plasma
  CPfCBZ = CPlCBZ*FUPCBZ
  CPfHBC = CPlHBC*FUPHBC
  CPfAB = CPlAB*FUPAB
  
  #kidney compartment
  #Concentrations of CBZ, HBC and AB in kidneys and venous blood from kidneys
  CKCBZ = AKCBZ/VK
  CVKCBZ = (CKCBZ/PKCBZ)*BPRCBZ
  
  CKHBC = AKHBC/VK
  CVKHBC = (CKHBC/PKHBC)*BPRHBC
  
  CKAB = AKAB/VK
  CVKAB = (CKAB/PKAB)*BPRAB
  
  #Renal excretion
  ARECBZ = CKCBZ * FUPCBZ * GFR
  AREHBC = CKHBC * FUPHBC * GFR
  AREAB = CKAB * FUPAB * GFR
  
  
#####Differential equations  #####
  #Stomach
  #HUMAN
  kaCBZ = 0.8844
  #RAT
  #kaCBZ = 8.598
  
  d/dt(Ast) = -kaCBZ * Ast
  
  #Cumulative metabolism in blood
  #d/dt(SAMABPl) = AMABPl
  AMABPl = 0
  #Amount in blood
  #Sum of the amount of compound in venous blood - amount already in blood
  d/dt(ABlCBZ) = QF * CVFCBZ + QL * CVLCBZ + QK * CVKCBZ + QB * CVBCBZ + QS * CVSCBZ + QR * CVRCBZ - (QF+QL+QK+QB+QS+QR) * CBlCBZ - AMABPl
  d/dt(ABlHBC) = QF * CVFHBC + QL * CVLHBC + QK * CVKHBC + QB * CVBHBC + QS * CVSHBC + QR * CVRHBC - (QF+QL+QK+QB+QS+QR) * CBlHBC - AMABPl
  d/dt(ABlAB) = QF * CVFAB + QL * CVLAB + QK * CVKAB + QB * CVBAB + QS * CVSAB + QR * CVRAB - (QF+QL+QK+QB+QS+QR) * CBlAB - AMABPl
  
  #Metabolism in plasma
  #AMABPl = (VMax3AB*CPfCBZ)/(Km3AB+CPfCBZ)
  
  #Fat
  d/dt(AFCBZ) = QF * (CBlCBZ - CVFCBZ)
  d/dt(AFHBC) = QF * (CBlHBC - CVFHBC)
  d/dt(AFAB) = QF * (CBlAB - CVFAB) 
  
  #Cumulative metabolism liver
  d/dt(SAMHBCL) = AMHBCL
  d/dt(SAMABL) = AMABL
  
  #Liver
  d/dt(ALCBZ) = QL * (CBlCBZ - CVLCBZ) + kaCBZ * Ast - AMHBCL- AMABL
  d/dt(ALHBC) = QL * (CBlHBC - CVLHBC) + AMHBCL
  d/dt(ALAB) = QL * (CBlAB - CVLAB) + AMABL
  
  #Brain
  d/dt(ABCBZ) = QB * (CBlCBZ - CVBCBZ)
  d/dt(ABHBC) = QB * (CBlHBC - CVBHBC)
  d/dt(ABAB) = QB * (CBlAB - CVBAB)
  
  #Richly perfused
  d/dt(ARCBZ) = QR * (CBlCBZ - CVRCBZ)
  d/dt(ARHBC) = QR * (CBlHBC - CVRHBC)
  d/dt(ARAB) = QR * (CBlAB - CVRAB)
  
  #Slowly perfused
  d/dt(ASCBZ) = QS * (CBlCBZ - CVSCBZ)
  d/dt(ASHBC) = QS * (CBlHBC - CVSHBC)
  d/dt(ASAB) = QS * (CBlAB - CVSAB)
  
  #Renal excretion
  d/dt(SARECBZ) <- ARECBZ;
  d/dt(SAREHBC) <- AREHBC;
  d/dt(SAREAB) <- AREAB;
  
  #Kidney
  d/dt(AKCBZ) = QK * (CBlCBZ - CVKCBZ) - ARECBZ
  d/dt(AKHBC) = QK * (CBlHBC - CVKHBC) - AREHBC
  d/dt(AKAB) = QK * (CBlAB - CVKAB) - AREAB
  
})
}
#### DOSING ####
#Molecular weight
MWCBZ = 191.21		    #Molecular weight carbendazim
MWHBC = 207.19		    #Molecular weight 5-hydroxycarbendazim
MWAB = 	133.15		    #Molecular weight 2-aminobenzimidazole

dosemg = 12 * 83.3                          #dose of CBZ in mg (dose in mg/kg * BW)
doseumol = dosemg / 1000 / MWCBZ * 10^6   #dose of CBZ in umol
TDumol = faCBZ * doseumol                 #absorbed dose taking into account the fractional absorption


ev <<- eventTable(amount.units = "umol", time.units = "h") %>% 
  et(seq(0,72, by = 0.1))
#Oral dosing
ev$add.dosing(
  dose = TDumol,
  dosing.to = "Ast"
)
#IV dosing
# ev$add.dosing(
#   dose = doseumol,
#   dosing.to = "ABlCBZ"
# )

#### Solving ####

#ModelOutputrativ12  <<- solve(Model, NULL, events = ev, inits = inits, cores = 4)
#ModelOutputratoral12  <<- solve(Model, NULL, events = ev, inits = inits, cores = 4)
#ModelOutputratoral1000  <<- solve(Model, NULL, events = ev, inits = inits, cores = 4)

ModelOutputhumanoral12  <<- solve(Model, NULL, events = ev, inits = inits, cores = 4)

##### Import in vivo data  #####
#Get data from paths
ivblood <- "C:/Users/sofie/OneDrive - Wageningen University & Research/MSc Thesis tox/20240112 Krechniak 1986 blood.csv"
ivliver <- "C:/Users/sofie/OneDrive - Wageningen University & Research/MSc Thesis tox/20240112 Krechniak 1986 liver.csv" 
ivkidney <- "C:/Users/sofie/OneDrive - Wageningen University & Research/MSc Thesis tox/20240113Krechniak1986kidney.csv"
ivexcr <- "C:/Users/sofie/OneDrive - Wageningen University & Research/MSc Thesis tox/20240112 Krechniak 1986 excretion oral.csv"
oralexcr <- "C:/Users/sofie/OneDrive - Wageningen University & Research/MSc Thesis tox/20240112 Krechniak 1986 excretion iv.csv"
oralserum <- "C:/Users/sofie/OneDrive - Wageningen University & Research/MSc Thesis tox/20240114 jia 2003 serum CBZ.csv"

# Read the CSV file into a dataframe
#Rat, dose = 12 mg/kg, BW = 0.22 
ivblood_df <- read.csv(ivblood,header=TRUE,sep=",",nrows=8)
ivliver_df <- read.csv(ivliver,header=TRUE,sep=",",nrows=8)
ivkidney_df <- read.csv(ivkidney,header=TRUE,sep=",",nrows=8)
#Rat, dose = 12 mg/kg, BW = 0.22 
oralexcr_df <- read.csv(oralexcr,header=TRUE,sep=",",nrows=6)
ivexcr_df <- read.csv(ivexcr,header=TRUE,sep=",",nrows=6)
#Rat, dose = 1000 mg/kg, BW = 
oralserum_df <- read.csv(oralserum, header = TRUE,sep=",", nrows=7)

#prevent axis labeling in scientific notation
options(scipen=100000)


#### read qivivetools.wur.nl data ####
qivivetoolsrat12 <- "C:/Users/sofie/OneDrive - Wageningen University & Research/MSc Thesis tox/202040207 qiviveresult rat liver.csv"
ivratqt_df <- read.csv(qivivetoolsrat12, header = TRUE, sep=",")
view(ivratqt_df)


##### IV dosing plots  #####
#Plot liver iv dosing
RatCBZliver2 = ggplot() +
  geom_point(data=ivliver_df,
             aes(x = ivliver_df[,1], y = ivliver_df[,3], shape= "CBZ")
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = CLCBZ, linetype = "CBZ"),
            
  ) +
  geom_point(data = ivliver_df,
             aes(x = ivliver_df[,4], y = ivliver_df[,6],shape = "5-HBC"),
             
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = CLHBC,linetype = "5-HBC"),
            
  ) +
  geom_point(data  = ivliver_df,
             aes(x = ivliver_df[,7], y = ivliver_df[,9],shape = "2-AB")
             
  ) +
  geom_line(data = ivratqt_df,
            aes(x=time, y=Cliver,linetype = "CBZ GM"),
            )+
  labs(title = "Liver",
       x="Time (h)",
       y = "Concentration (µM)",
       shape = "In vivo",
       linetype = "In silico") +
  theme_bw() +
  theme(legend.position = "none")+
  xlim(0,24) + 
  scale_y_continuous(trans='log10',limits = c(0.01,50))


print(RatCBZliver2)

#Plot blood iv dosing
RatCBZblood = ggplot() + 
  geom_point(data = ivblood_df,
             aes(x = ivblood_df[,1], y = ivblood_df[,3], shape = "CBZ"),
             
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = CBlCBZ,linetype = "CBZ")
            
  ) +
  geom_point(data = ivblood_df,
             aes(x = ivblood_df[,4], y = ivblood_df[,6], shape = "5-HBC")
             
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = CBlHBC, linetype = "5-HBC")
            
  ) +
  geom_point(data  = ivblood_df,
             aes(x = ivblood_df[,7], y = ivblood_df[,9], shape = "2-AB")
             
  ) +
  geom_line(data = ivratqt_df,
            aes(x=time, y=Cplasmavenous, linetype = "CBZ GM")
  )+
  labs(title = "Blood",
       x = "Time (h)",
       y = "Concentration (µM)",
       shape = "In vivo",
       linetype = "In silico") +
  theme_bw() +
  #theme(legend.position = "none") +
  xlim(0,24)+ 
  scale_y_continuous(trans='log10',limits = c(0.01,1000))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

print(RatCBZblood)
                        
#Plot kidney iv dosing
RatCBZkidney2 = ggplot() + 
  geom_point(data = ivkidney_df,
             aes(x = ivkidney_df[,1], y = ivkidney_df[,3], shape = "CBZ"),
             
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = CKCBZ, linetype = "CBZ")
            
  ) +
  geom_point(data = ivkidney_df,
             aes(x = ivkidney_df[,4], y = ivkidney_df[,6], shape = "5-HBC")
             
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = CKHBC, linetype = "5-HBC"),
            
  ) +
  geom_point(data  = ivkidney_df,
             aes(x = ivkidney_df[,7], y = ivkidney_df[,9],shape = "2-AB")
             ) +
  geom_line(data = ivratqt_df,
            aes(x=time, y=Ckidney, linetype = "CBZ GM")
            )+
  labs(title = "Kidney",
       x = "Time (h)",
       y = "Concentration (µM)",
       shape = "In vivo",
       linetype = "In silico") +
  theme_bw() +
  xlim(0,24) +
  scale_y_continuous(trans='log10',limits = c(0.01,100))

print(RatCBZkidney2)


#Plot urinal excretion iv dosing
RatCBZivexcr = ggplot() + 
  geom_point(data = ivexcr_df,
             aes(x = ivexcr_df[,1], y = ivexcr_df[,3], shape = "CBZ")
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = ARECBZ, linetype = "CBZ")
  ) +
  geom_point(data = ivexcr_df,
             aes(x = ivexcr_df[,4], y = ivexcr_df[,6], shape = "5-HBC")
  ) +
  geom_line(data = ModelOutputrativ12,
            aes(x = time, y = AREHBC, linetype = "5-HBC")
  ) +
  geom_point(data = ivexcr_df,
             aes(x = ivexcr_df[,7], y = ivexcr_df[,9], shape = "2-AB")
  ) +
  labs(title = "Urinal excretion",
       x = "Time (h)",
       y = "Excretion rate (µmol/h)",
       shape = "In vivo",
       linetype = "In silico") +
  theme_bw() +
  #theme(legend.position = "none") +
  xlim(0,75)+ 
  scale_y_continuous(trans='log10',limits = c(0.00001,1))

print(RatCBZivexcr)

#### Both oral/iv #####

RativCBZexcrcumulative = ggplot() +
  geom_line(data=ModelOutput,
            aes(x = time, y = SARECBZ),
            linetype = 1
  ) +
  geom_line(data=ModelOutput,
            aes(x = time, y = SAREHBC),
            linetype = 2
  ) +
  labs(x = "Time (h)",
       y = "Amount (µmol)") +
  theme_bw() +
  xlim(0,75)

print(RativCBZexcrcumulative)

RatoralCBZexcrcumulative = ggplot() +
  geom_line(data=ModelOutput,
            aes(x = time, y = SARECBZ),
            linetype = 1
  ) +
  geom_line(data=ModelOutput,
            aes(x = time, y = SAREHBC),
            linetype = 2
  ) +
  labs(x = "Time (h)",
       y = "Amount (µmol)") +
  theme_bw() +
  xlim(0,75)

print(RatoralCBZexcrcumulative)

HumanoralCBZexcrcumulative = ggplot() +
  geom_line(data=ModelOutput,
            aes(x = time, y = SARECBZ),
            linetype = 1
  ) +
  geom_line(data=ModelOutput,
            aes(x = time, y = SAREHBC),
            linetype = 2
  ) +
  labs(x = "Time (h)",
       y = "Amount (µmol)") +
  theme_bw() +
  xlim(0,75)

print(HumanoralCBZexcrcumulative)


  

  
##### Oral dosing plots  #####
# Plot serum, dose = 1000*BW, BW = 0.185
RatCBZserum = ggplot() + 
  geom_point(data = oralserum_df,
             aes(x = oralserum_df[,1], y = oralserum_df[,3], shape = "CBZ")
             
  ) +
  geom_line(data = ModelOutputratoral1000,
            aes(x = time, y = CPlCBZ, linetype = "CBZ"),
            
  ) +
  geom_line(data=ModelOutputratoral1000,
            aes(x=time, y=CPlHBC, linetype = "5-HBC"),
            
  ) +
  labs(title = "Serum, dose = 1000 mg/kg BW",
       x = "Time (h)",
       y = "Concentration (µM)",
       shape = "In vivo",
       linetype = "In silico"
  ) +
  theme_bw() +
  xlim(0,75)+ 
  scale_y_continuous(trans='log10',limits = c(0.01,5000))
  
print(RatCBZserum)

#Plot urinal excretion oral dosing
RatCBZoralexcr = ggplot() + 
  geom_line(data = ModelOutputratoral12,
            aes(x = time, y = ARECBZ, linetype = "CBZ")
            
  ) +
  geom_point(data = oralexcr_df,
             aes(x = oralexcr_df[,1], y = oralexcr_df[,3],shape = "CBZ")
             
  ) +
  geom_line(data = ModelOutputratoral12,
            aes(x = time, y = AREHBC, linetype = "5-HBC")
            
  ) +
  geom_point(data = oralexcr_df,
             aes(x = oralexcr_df[,4], y = oralexcr_df[,6], shape = "5-HBC")
            
  ) +
  geom_point(data = oralexcr_df,
             aes(x = oralexcr_df[,7], y = oralexcr_df[,9], shape = "2-AB")
             
  )+
  labs(title = "Urinal excretion, dose = 12 mg/kg BW",
       x = "Time (h)",
       y = "Excretion rate (µmol/h)",
       shape = "In vivo",
       linetype = "In silico") +
  theme_bw() +
  xlim(0,75)+ 
  scale_y_continuous(trans='log10',limits = c(0.00001,1))


print(RatCBZoralexcr)

#### Predicting plots oral, no in vivo data ####

RatCBZbrainpred = ggplot() +
  geom_line(data = ModelOutputratoral12,
            aes(x = time, y = CBCBZ),
            linetype = 1
  ) +
  geom_line(data = ModelOutputratoral12,
            aes(x = time, y = CBHBC),
            linetype = 2
  ) +
  labs(title = "Brain",
       x = "Time (h)",
       y = "Concentration (µM)") +
  theme_bw() +
  xlim(0,24)

print(RatCBZbrainpred)

HumanCBZbrainpred = ggplot() +
  geom_line(data = ModelOutputhumanoral12,
            aes(x = time, y = CBCBZ, linetype = "CBZ")
  ) +
  geom_line(data = ModelOutputhumanoral12,
            aes(x = time, y = CBHBC, linetype = "5-HBC")
  ) +
  labs(title = "Brain",
       x = "Time (h)",
       y = "Concentration (µM)",
       linetype = "In silico") +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(0,24)+ 
  scale_y_continuous(trans='log10',limits = c(0.001,100))

print(HumanCBZbrainpred)

HumanCBZbloodpred = ggplot() + 
  geom_line(data = ModelOutputhumanoral12,
            aes(x = time, y = CBlCBZ, linetype = "CBZ")
  ) +
  geom_line(data = ModelOutputhumanoral12,
            aes(x = time, y = CBlHBC, linetype = "HBC")
  ) +
  labs(title = "Blood",
       x = "Time (h)",
       y = "Concentration (µM)",
       linetype = "In silico") +
  theme_bw() +
  theme(legend.position = "none") +
  xlim(0,24)+ 
  scale_y_continuous(trans='log10',limits = c(0.001,100))

print(HumanCBZbloodpred)

HumanCBZexcrpred = ggplot() +
  geom_line(data = ModelOutputhumanoral12,
            aes(x = time, y = ARECBZ, linetype = "CBZ")
  ) +
  geom_line(data = ModelOutputhumanoral12,
            aes(x = time, y = AREHBC, linetype = "HBC")
  )+
  labs(title = "Urinal excretion",
       x = "Time (h)",
       y = "Excretion rate (µmol/h)",
       linetype = "In silico") +
  theme_bw() +
  #theme(legend.position = "none") +
  xlim(0,75)+ 
  scale_y_continuous(trans='log10',limits = c(0.0001,1000))

print(HumanCBZexcrpred)

#### Plots for presentation####
Plotsrativ12pp = (RatCBZivexcr|RatCBZblood)
print(Plotsrativ12pp)

#### Plots for report ####
Plotsrativ12 = ((RatCBZliver2|RatCBZkidney2)/(RatCBZblood|RatCBZivexcr))
Plotsratoral = (RatCBZserum|RatCBZoralexcr)
Plotshumanoral12 = (HumanCBZbrainpred|HumanCBZbloodpred|HumanCBZexcrpred)

                
#Ratio for exporting = 1000:750
Plotsrativ12 + plot_annotation(tag_levels = 'A')
#Ratio = 1000:375
Plotsratoral + plot_annotation(tag_levels = 'A')
Plotshumanoral12 + plot_annotation(tag_levels = 'A')

#### MASSBALANCE ####
colsCBZ <-c("Ast", "AFCBZ","ALCBZ","ABCBZ","ARCBZ","ASCBZ","ABlCBZ","AKCBZ","SARECBZ","SAMHBCL", "SAMABL")
colsHBC <-c("AFHBC","ALHBC","ABHBC","ARHBC","ASHBC","ABlHBC","AKHBC","SAREHBC")
colsAB <-c("AFAB","ALAB","AMABL","ABAB","ARAB","ASAB","ABlAB","AKAB","SAREAB")
colsAll <- c("Ast", "AFCBZ","ALCBZ","ABCBZ","ARCBZ","ASCBZ","ABlCBZ","AKCBZ","SARECBZ",
             "AFHBC","ALHBC","ABHBC","ARHBC","ASHBC","ABlHBC","AKHBC","SAREHBC",
             "AFAB","ALAB","AMABL","ABAB","ARAB","ASAB","ABlAB","AKAB","SAREAB")


CalcBodAll <- rowSums(ModelOutput[, colsAll])
CalcBodCBZ <- rowSums(ModelOutput[, colsCBZ])
CalcBodHBC <- rowSums(ModelOutput[, colsHBC])
CalcBodAB <- rowSums(ModelOutput[, colsAB])
# Calculate what goes in. This is the administered amount of compound, and for the metabolite submodels, this is the amount metabolized.

TotBodAll <- TDumol # Amount administered (in umol)
TotBodCBZ <- TDumol
TotBodHBC <- ModelOutput[, "SAMHBCL"] # Amount metabolized from CBZ
TotBodAB   <- 0#rowSums(ModelOutput[, c("SAMABL")])
# Calculate mass balance
# Total administered or produced amount - minus what is present in the model + 1 (=100%). Done for each timepoint
MassBall    <- TotBodAll - CalcBodAll + 1
MassBalCBZ  <- TotBodCBZ - CalcBodCBZ + 1
MassBalHBC  <- TotBodHBC - CalcBodHBC + 1
MassBalAB  <- TotBodAB - CalcBodAB + 1
# Calculate mass error (% of mass lost)
# Difference between total amount administered or produced and the amount in the (sub)model. 
# This is divided by the total amount administered or produced and multiplied by 100 to give a percentage
ErrorAll    <- (TotBodAll - CalcBodAll) / (TotBodAll + 10^(-30)) * 100
ErrorCBZ    <- (TotBodCBZ - CalcBodCBZ) / (TotBodCBZ + 10^(-30)) * 100
ErrorHBC    <- (TotBodHBC - CalcBodHBC) / (TotBodHBC + 10^(-30)) * 100
ErrorAB    <- (TotBodAB - CalcBodAB) / (TotBodAB + 10^(-30)) * 100
# Compile Mass Balance data into a dataframe
MassBal  <- cbind("time" = ModelOutput[,"time"], MassBall, ErrorAll, MassBalCBZ, ErrorCBZ, MassBalHBC, 
                   ErrorHBC, MassBalAB, ErrorAB)
MassBal <- as.data.frame(MassBal)

##### Mass balance plots  #####

ggplot(data = MassBal) +
  geom_line(aes(x = time,
                y = MassBall)) +
  geom_line(aes(x = time,
                y = ErrorAll)) +
  labs(title = "MB All",
       x = "Time (h)",
       y = "Concentration (umol/L)")

 ggplot(data = MassBal) +
  geom_line(aes(x = time,
                y = MassBalCBZ)) +
  geom_line(aes(x = time,
                y = ErrorCBZ)) +
  labs(title = "MB CBZ",
       x = "Time (h)",
       y = "Concentration (umol/L)")

ggplot(data = MassBal) +
  geom_line(aes(x = time,
                y = MassBalHBC)) +
  geom_line(aes(x = time,
                y = ErrorHBC)) +
  labs(title = "MB HBC",
       x = "Time (h)",
       y = "Concentration (umol/L)")

