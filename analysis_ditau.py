###########################################
#########  H_stastau Analysis   ###########
###########################################

import ROOT 
import numpy as np 
import math
import pickle
import os
import operator
from array import array
ROOT.gROOT.SetBatch(True)
import logging
# Initialize the xAOD infrastructure:
if(not ROOT.xAOD.Init().isSuccess()): print "Failed xAOD.Init()"
ROOT.gInterpreter.ProcessLine('#include "/cluster/home/bataju/compare/mt2/CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"')
ROOT.gSystem.Load("/cluster/home/bataju/compare/mt2/CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h")

#################################
########### Functions ###########
#################################

def find_child(p):
	#takes a particles container and retuns a list of the child particle container
	for i in range(p.nChildren()):
		if type(p.child(i).pdgId()) == []:
			continue
		elif p.child(i).pdgId() != p.pdgId():
			p_child = [p.child(i) for i in range(p.nChildren())]
			return p_child
			#return [p.child(i) for i in range(p.nChildren())]
		else:  
			return find_child(p.child(i))

def find_parent(p):
	#takes a particles container and retruns a list of the child particle container
	for i in range(p.nParents()):
		if p.parent(i).pdgId() == p.pdgId():
			#notice this is a list containing p.child(i) elements
			return [p.child(i) for i in range(p.nParents())]
		else:
			return find_parent(p.parent(i))

def pdg(n):
	#takes a list of particles container and returns the pdgId in a list

	return [i.pdgId() for i in n if i is not None]
def MASS(n):
	#takes a  particles container and returns the MASS in a list
	return [i.p4().M() for i in n] 

def Apdg(n):
	#takes a list of particles container and returns the pdgId in a list
	
	return [i.absPdgId() for i in n if i is not None]

def get_eta(n):
	#takes a list of particles container and returns the Eta in a list
	
	return [i.eta() for i in n]

def get_phi(n):
	#takes a list of particles container and returns the Phi in a list
	
	return [i.phi() for i in n]

def PT(n):
	#takes a  particles container and returns the pt in a list
	return [i.p4().Pt() for i in n if i is not None]

def output_True_from_list(test_list):
	#shows information about list 
	return list(filter(lambda i: test_list[i], range(len(test_list))))

def DeltaR(v1,v2):
	#calculates deltaR from two .p4()
	p1=v1
	p2=v2
	return p1.DeltaR(p2)
def DeltaPhi(v1,v2):
	#calculates deltaR from two .p4()
	p1=v1
	p2=v2
	return p1.DeltaPhi(p2)
def mT(v1,v2):
	# calcaulate transverse mass for p4 and p4
	p1=v1
	p2=v2
	return (p1+p2).Mt()

def mT_(v1,v2):
	# calcaulate transverse mass for p4 and met
	p1=v1
	p2=v2
	return (p1+p2).Mt()

def find_vis(p):
	invisible = [12,14,16]
	vis = []
	pp = find_child(p)
	a = ROOT.TLorentzVector(0,0,0,0)
	for i in pp:
		if i.absPdgId() not in invisible:
			vis.append(i)
			#print i
	assert len(vis), " No children found for {}".format(pdg(p))
	for i in vis:
		a +=i.p4()	
	return a

#################################
######### Histograms ############
#################################


#################################
######### Opening File ##########
#################################

# sig_file = ['/cluster/home/bataju/375_800/lxplux_aod/10k/Merge._10k.pool.root']	#800
# sig_file.append('/cluster/home/bataju/375_800/lxplux_aod/20k/Merge._20k.pool.root')	#800

# sig_file =['/cluster/home/bataju/lxplus_aod/500/Merge.aod.pool.root ']	#500
# sig_file.append('/cluster/home/bataju/lxplus_aod/500/Merge.10k.pool.root ')	#500
# sig_file.append('/cluster/home/bataju/filter/500_aod/merge._20k.pool.root')	#500

# sig_file = ['/cluster/home/bataju/BENCHMARK/1000/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 1000
# sig_file = ['/cluster/home/bataju/BENCHMARK/1500/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 1500
# sig_file = ['/cluster/home/bataju/BENCHMARK/2000/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 2000
# sig_file = ['/cluster/home/bataju/BENCHMARK/2500/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 2500
# sig_file = ['/cluster/home/bataju/BENCHMARK/700/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 700
# sig_file = ['/cluster/home/bataju/BENCHMARK/800/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 800

# sig_file = ['/cluster/home/bataju/compare/eventgeneration/different_stau_mass/1500/TRUTH1/DAOD_TRUTH1.test.pool.root']	#1500
# sig_file = ['/cluster/home/bataju/compare/eventgeneration/different_stau_mass/1000/TRUTH1/DAOD_TRUTH1.test.pool.root']	#1000
# sig_file = ['/cluster/home/bataju/compare/eventgeneration/different_stau_mass/2000/TRUTH1/DAOD_TRUTH1.test.pool.root']	#2000

# for filename in os.listdir(sig_300):
# 	if filename.startswith("Merge"):
# 		print filename
# 		sig_file.append(sig_300+filename)

#bkg 
# ttbar   =   ['/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.12458455._000337.pool.root.1']
# ttbar.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.12458455._004964.pool.root.1')
# ttbar.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.14761208._000530.pool.root.1')
# ttbar.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.14761208._004126.pool.root.1')

# ztautau =     ['/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182601._001043.pool.root.1']
# ztautau.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182601._001055.pool.root.1')
# ztautau.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182673._001851.pool.root.1')
# ztautau.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_la13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182673._003111.pool.root.1')
tautau = []
f1 = [0,1,2,4,5]
tautau = [['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341881.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH600_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4482_e5984_a875_r10724_r10726_p3759/DAOD_HIGG4D4.17073421._00000{}.pool.root.1'.format(i+1) for i in xrange(6)]] 	#600
# tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341885.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH1000_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4298_e5984_a875_r10201_r10210_p3749/DAOD_HIGG4D4.17259964._00000{}.pool.root.1'.format(i+1) for i in xrange(6)]) #1000
tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.345288.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH2000_yb2_tautauhh.deriv.DAOD_HIGG4D4.e5686_e5984_a875_r10724_r10726_p3759/DAOD_HIGG4D4.17076570._00000{}.pool.root.1'.format(i+1) for i in xrange(6)])	#2000
tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.345292.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH2500_yb2_tautauhh.deriv.DAOD_HIGG4D4.e5686_e5984_a875_r10724_r10726_p3759/DAOD_HIGG4D4.17073341._00000{}.pool.root.1'.format(i+1) for i in xrange(6)])	#2500
# tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341874.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH125_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4482_a875_r10724_p3759/DAOD_HIGG4D4.17073198._00000{}.pool.root.1'.format(i+1) for i in xrange(6)])	#125
tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341875.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH200_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4482_a875_r10724_p3759/DAOD_HIGG4D4.17075704._00000{}.pool.root.1'.format(i+1) for i in xrange(6)])	#200
# tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341877.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH300_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4482_a875_r10724_p3759/DAOD_HIGG4D4.17075878._00000{}.pool.root.1'.format(i+1) for i in xrange(6)])	#300
tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341879.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH400_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4298_a875_r10724_p3759/DAOD_HIGG4D4.17073318._00000{}.pool.root.1'.format(i+1) for i in [0,1,2,6,5]])	#400
tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341920.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH1500_yb2_tautauhh.deriv.DAOD_HIGG4D4.e5314_a875_r10724_p3759/DAOD_HIGG4D4.17075106._0000{}.pool.root.1'.format(i) for i in [11,12,13,14,15]])	#1500
# tautau.append(['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341879.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH400_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4298_a875_r10724_p3759/DAOD_HIGG4D4.17073318._00000{}.pool.root.1'.format(i+1)  for i in xrange(6)])
temp = [600,2000,2500,200,400,1500]
# temp = [400]
# f1 = [3,4,5,7,11]
f1 = [0,1,2,4,5]
# tautau = [['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341920.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH1500_yb2_tautauhh.deriv.DAOD_HIGG4D4.e5314_a875_r10724_p3759/DAOD_HIGG4D4.17075106._00000{}.pool.root.1'.format(i) for i in f1]]	#1500
# tautau = [['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341879.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH400_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4298_a875_r10724_p3759/DAOD_HIGG4D4.17073318._00000{}.pool.root.1'.format(i+1)  for i in f1]]	#400

# temp = [400]
# temp = [1500]


# file_list = sig_file[:]
file_list = tautau[:]

# print len(file_list) - len(sig_file), "the number of bkg files"
# print len(sig_file), "the number of sig files"
print file_list
print "number of files", len(file_list)

#################################
############ LIST ###############
#################################

count 			= 0
skiped 			= 0


########### Analysis Code Starts ############
for m in xrange(len(temp)):
	print temp[m]
	outF = ROOT.TFile.Open("H_2tau_{}.root".format(temp[m]), "RECREATE")
	outTree = ROOT.TTree('T','Test TTree')
	
	# creating vectors 
	Higgs 	= ROOT.vector(ROOT.TLorentzVector)()
	mets  	= ROOT.vector(ROOT.TLorentzVector)()
	tau1   	= ROOT.vector(ROOT.TLorentzVector)()
	tau2   	= ROOT.vector(ROOT.TLorentzVector)()
	vistau1 = ROOT.vector(ROOT.TLorentzVector)()
	vistau2 = ROOT.vector(ROOT.TLorentzVector)()

	# create branches

	outTree.Branch("higgs",Higgs )
	outTree.Branch("met",mets )
	outTree.Branch("tau1",tau1 )
	outTree.Branch("tau2",tau2 )
	outTree.Branch("vistau1",vistau1 )
	outTree.Branch("vistau2",vistau2 )
	Higgs.clear()
	mets.clear()
	tau1.clear()
	tau2.clear()
	vistau1.clear()
	vistau2.clear()
	
	for k in xrange(len(file_list[m])):
		f = ROOT.TFile.Open(file_list[m][k], "READONLY")
		t = ROOT.xAOD.MakeTransientTree(f, "CollectionTree")

		for entry in xrange(t.GetEntries()):
			if entry == 2000: break
			#========================= Standard Prints
			# print "Number of input events:", t.GetEntries()
			# print "Working on ", file_list[m][k]
			print  entry*100/t.GetEntries(), "% complete." , "{}/{} files".format(k+1,len(file_list[m]))
			print "This is entry number: ", entry
			#========================= Standard Prints
			
			t.GetEntry(entry)
			Met = t.MET_Truth
			all_particles = t.TruthParticles
			tau_particles = t.TruthTaus
			metvector1 = ROOT.TLorentzVector(0,0,0,0)
			metvector1.SetPtEtaPhiM(Met.get(1).met(), 0, Met.get(1).phi(), 0)
   			METp = ROOT.TLorentzVector(0,0,0,0)
			METp.SetPxPyPzE(Met.get(1).mpx(),Met.get(1).mpy(),0,Met.get(1).sumet())
		
			higgs = [i for i in all_particles if i.absPdgId() in [25]]
			print pdg(higgs)
			print "__________________________________________"
			higgs_with_children = []
			if len(higgs) < 1: continue
			for i  in higgs:
				if i.nChildren() == 0: continue
				elif i.child(0).pdgId() == i.pdgId():	continue
				elif i.child(0).absPdgId() in [15]: choosen_higgs = [i]
				# elif i.child(0).absPdgId() in [1000015,2000015]: choosen_higgs = [i.child(p) for p in xrange(i.nChildren())]
				else:	higgs_with_children.append(i)
			assert len(choosen_higgs) == 1
			
			
			higgs_last = find_child(choosen_higgs[0])
			

			tau_h = [i for i in higgs_last if i.absPdgId() == 15]
			if len(tau_h) < 2:	continue
			# assert len(tau_h) >= 2, "Less than one tau!!!"
			# print "Di-stau Pt of all the stau", (tau_h[0].p4() +tau_h[1].p4()).Pt()
			# print "parent of stau [0] \t" , tau_h[0].parent().pdgId() , "barcode \t" , tau_h[0].parent().barcode()
			# print "parent of stau [1] \t" , tau_h[1].parent().pdgId() , "barcode \t"  , tau_h[1].parent().barcode()
			t_h1 = [] #hadronic taus
			t_h2 = [] #hadronic taus
			tau_vis = []
			if len(tau_h) < 2: continue
			tau_1 = tau_h[0]
			tau_2 = tau_h[1]
			checkelmu = []

			
			for i in Apdg(find_child(tau_1)):
				checkelmu.append(i in [11,13])
			if sum(checkelmu) < 1:
				t_h1 = tau_1
				tau_vis.append(tau_1)
			
			checkelmu = []
			for i in Apdg(find_child(tau_2)):
				checkelmu.append(i in [11,13])
			if sum(checkelmu) < 1:	
				t_h2 = tau_2
				tau_vis.append(tau_2)	

			if len(tau_vis) < 2: continue	

			Higgs.push_back(choosen_higgs[0].p4())
			tau1.push_back(tau_h[0].p4())
			tau2.push_back(tau_h[1].p4())
			mets.push_back(metvector1)
			vistau1.push_back(find_vis(tau_h[1]))
			vistau2.push_back(find_vis(tau_h[1]))

	outTree.Fill()
	

	outF.Write()
	outF.Close()	
	print "Finished."
	