###########################################
#########  H_stastau Analysis   ###########
###########################################


import os
# print os.environ
import ROOT 
import numpy as np 
import math
import pickle
from array import array
ROOT.gROOT.SetBatch(True)
import logging
# Initialize the xAOD infrastructure:
if(not ROOT.xAOD.Init().isSuccess()): print "Failed xAOD.Init()"


#################################
########### Functions ###########
#################################

def find_child(p):
	#takes a particles container and retuns a list of the child particle container
	for i in range(p.nChildren()):
		if type(p.child(i).pdgId()) == []:
			continue
		elif p.child(i).pdgId() != p.pdgId():
			#notice this is a list containing p.child(i) elements
			#new
			p_child = [p.child(i) for i in range(p.nChildren())]
			return p_child
			#return [p.child(i) for i in range(p.nChildren())]
		else:
			return find_child(p.child(i))

def find_parent(p):
	#takes a particles container and retruns a list of the child particle container
	for i in range(p.nParents()):
		if p.parent(i).pdgId() != p.pdgId():
			#notice this is a list containing p.child(i) elements
			return [p.child(i) for i in range(p.nParents())]
		else:

			return find_parent(p.parent(i))

def pdg(n):
	#takes a list of particles container and returns the pdgId in a list

	return [i.pdgId() for i in n]

def Apdg(n):
	#takes a list of particles container and returns the pdgId in a list
	if isinstance(n[0],ROOT.xAOD.TruthParticle_v1):
		return [i.absPdgId() for i in n]
	else:
		pass
def get_eta(n):
	#takes a list of particles container and returns the Eta in a list
	
	return [i.eta() for i in n]

def get_phi(n):
	#takes a list of particles container and returns the Phi in a list
	
	return [i.phi() for i in n]

def PT(n):
	#takes a  particles container and returns the pt in a list
	if isinstance(n[0],ROOT.xAOD.TruthParticle_v1):
		return [i.p4().Pt() for i in n]
	elif isinstance(n[0],ROOT.TLorentzVector):
		return [i.Pt() for i in n]

def MASS(n):
	#takes a  particles container and returns the MASS in a list
	return [i.p4().M() for i in n]

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

# temp = [700,730,800,830,930,1000,1130,1280,1500,1700,2000,2200,2500]
temp = [1700]
sig_file = ['/cluster/home/bataju/BENCHMARK/'+str(i)+'/truth/DAOD_TRUTH1.test.pool.root' for i in temp]

file_list = sig_file[:]

print len(file_list) - len(sig_file), "the number of bkg files"
print len(sig_file), "the number of sig files"
print file_list


#################################
############ LIST ###############
#################################


count 			= 0
skiped 			= 0
outF = ROOT.TFile.Open("H_2stau_2tau_old_git_1700.root", "RECREATE")
outTree = ROOT.TTree('T','Test TTree')

# creating vectors 
Higgs = ROOT.vector(ROOT.TLorentzVector)()
stau1 = ROOT.vector(ROOT.TLorentzVector)()
stau2 = ROOT.vector(ROOT.TLorentzVector)()
outTree.Branch("higgs",Higgs )	
outTree.Branch("stau1",stau1 )
outTree.Branch("stau2",stau2 )

########### Analysis Code Starts ############

for k in xrange(len(file_list)):
	print "############"
	print file_list[k]

	f = ROOT.TFile.Open(file_list[k], "READONLY")
	t = ROOT.xAOD.MakeTransientTree(f, "CollectionTree")
	decay_chain = {}
	sorted_chain ={}
	for entry in xrange(t.GetEntries()):
		
		t.GetEntry(entry)
		Met = t.MET_Truth
		all_particles = t.TruthParticles
		metvector1 = ROOT.TLorentzVector(0,0,0,0)
		metvector1.SetPtEtaPhiM(Met.get(1).met(), 0, Met.get(1).phi(), 0)
		higgs = [i for i in all_particles if i.absPdgId() in [36,35]]
		assert len(higgs) ,"NO HIGGS BOSONS FOUND!!!"
		
		higgs_last = find_child(higgs[-1])
  		assert higgs_last is not None, " This higgs has not children"
		decay_chain["higgs"]=[higgs[-1]]
		
		Higgs.push_back(higgs[-1].p4()) 
			
		stau = [i for i in higgs_last if i.absPdgId() in [1000015,2000015]]
		assert len(stau) == 2, "There should be two staus!!!"
		decay_chain["stau"]=stau
		
		tau_h	= [] # tau list
		n1 = [i for i in find_child(stau[0]) if i.absPdgId() ==1000022]
		n2 = [i for i in find_child(stau[1]) if i.absPdgId() ==1000022]
		decay_chain["nue1"]=[n1[0],n2[0]]
		tau_1 = [i for i in find_child(stau[0]) if i.absPdgId() ==15]
		tau_2 = [i for i in find_child(stau[1]) if i.absPdgId() ==15]
		
		dPhi_tau1_tau2	= array('f',[])
		t_h1 = [] #hadronic taus
		t_h2 = [] #hadronic taus
		checkelmu = [] 
		for i in Apdg(find_child(tau_1[0])):
			checkelmu.append(i in [11,13])
		if sum(checkelmu) < 1:
			t_h1 = tau_1[:]
			tau_h.append(tau_1[0])
		
		checkelmu = []

		for i in Apdg(find_child(tau_2[0])):
			checkelmu.append(i in [11,13])
		if sum(checkelmu) < 1:	
			t_h2 = tau_2[:]
			tau_h.append(tau_2[0])
		
		#####################################################
		#t_lep t_lep
		if len(tau_h) < 1:	continue
		#####################################################
		#t_had t_lep
		elif len(tau_h) < 2:	continue
		#####################################################
		#t_had t_had 
		else:	
			tau_h_1 = find_vis(tau_h[0])
			tau_h_2 = find_vis(tau_h[1])
			decay_chain['tau']=[tau_h_1,tau_h_2]
			
			if decay_chain['tau'][0].Pt() > decay_chain['tau'][1].Pt():
				sorted_chain = dict(decay_chain)
			elif decay_chain['tau'][0].Pt() < decay_chain['tau'][1].Pt():
				sorted_chain['tau'] = list(reversed(decay_chain['tau']))
				sorted_chain['nue1'] = list(reversed(decay_chain['nue1']))
				sorted_chain['stau'] = list(reversed(decay_chain['stau']))
				sorted_chain['higgs'] = list(reversed(decay_chain["higgs"]))
			
			stau1.push_back(stau[0].p4())
			stau2.push_back(stau[1].p4())

	outTree.Fill()
	Higgs.clear()
	stau1.clear()
	stau2.clear()
outF.Write()
outF.Close()	
print "Finished."