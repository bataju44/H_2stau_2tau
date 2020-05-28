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
import logging
ROOT.gROOT.SetBatch(True)

# Initialize the xAOD infrastructure:
if(not ROOT.xAOD.Init().isSuccess()): print "Failed xAOD.Init()"
# for Mt2
ROOT.gInterpreter.ProcessLine('#include "/cluster/home/bataju/compare/mt2/CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h"')
ROOT.gSystem.Load("/cluster/home/bataju/compare/mt2/CalcGenericMT2/CalcGenericMT2/MT2_ROOT.h")

#################################
########### Functions ###########
#################################

def find_child(p):
	#doesnot return the immidate parent
	for i in range(p.nChildren()):
		if not p.child(i):
			continue
		elif p.child(i).pdgId() != p.pdgId():
			p_child = [p.child(i) for i in range(p.nChildren())]
			return p_child
			#return [p.child(i) for i in range(p.nChildren())]
		else:
			return find_child(p.child(i))

def find_parent(p):
	for i in range(p.nParents()):
		if p.parent(i).pdgId() != p.pdgId():
			return [p.child(i) for i in range(p.nParents())]
		else:

			return find_parent(p.parent(i))

def pdg(n):
	#takes a list of particles container and returns the pdgId in a list
	return [i.pdgId() for i in n if i is not None]

def Apdg(n):
	if not n:	
		return []
	elif isinstance(n[0],ROOT.xAOD.TruthParticle_v1):
		return [i.absPdgId() for i in n if i is not None]
	else:
		pass
def get_eta(n):
	#takes a list of particles container and returns the Eta in a list
	
	return [i.eta() for i in n if i is not None]

def get_phi(n):
	#takes a list of particles container and returns the Phi in a list
	
	return [i.phi() for i in n]

def PT(n):
	if not n:	
		return []
	elif isinstance(n[0],ROOT.xAOD.TruthParticle_v1):
		return [i.p4().Pt() for i in n if i is not None] 
	elif isinstance(n[0],ROOT.TLorentzVector):
		return [i.Pt() for i in n if i is not None]

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

def fill_hist(hist,array):
	[hist.Fill(_) for _ in array]

def scale_hist(hist_bkg_list,norm=1):
	if hist_bkg_list.Integral() != 0 :
		hist_bkg_list.Scale((norm/(hist_bkg_list.Integral())))
	# scale = norm/hist_bkg_list.Integral()
	# hist_bkg_list.Scale(scale)

#################################
######### Histograms ############
#################################


#################################
######### Opening File ##########
#################################

#
temp = [700,730,800,830,930,1000,1130,1280,1500,1700,2000,2200,2500]
# temp = [1500,1700,2000,2200,2500]
# temp = [730,800,1700]
sig_file = ['/cluster/home/bataju/BENCHMARK/'+str(i)+'/truth/DAOD_TRUTH1.test.pool.root' for i in temp]

file_list = sig_file[:]

print len(file_list) - len(sig_file), "the number of bkg files"
print len(sig_file), "the number of sig files"
print file_list


#################################
############ LIST ###############
#################################
h =ROOT.TH1F("h","h",30,0,1000E3)
h1 =ROOT.TH1F("h2","h1",30,0,1000E3)
h2 = ROOT.TH1F("h2",'h2',30,0,800E3)
count 			= 0
skiped 			= 0

mcweight_dict  ={} 
########### Analysis Code Starts ############

for k in xrange(len(file_list)):
	print "############"
	print file_list[k]

	f = ROOT.TFile.Open(file_list[k], "READONLY")
	t = ROOT.xAOD.MakeTransientTree(f, "CollectionTree")
	outF = ROOT.TFile.Open("H_2stauhad_{}.root".format(temp[k]), "RECREATE")
	outTree = ROOT.TTree('T','Test TTree')

	# creating vectors 
	Higgs = ROOT.vector(ROOT.TLorentzVector)()
	stau1 = ROOT.vector(ROOT.TLorentzVector)()
	stau2 = ROOT.vector(ROOT.TLorentzVector)()
	neu1  = ROOT.vector(ROOT.TLorentzVector)()
	neu2  = ROOT.vector(ROOT.TLorentzVector)()
	mets  = ROOT.vector(ROOT.TLorentzVector)()
	tau1  = ROOT.vector(ROOT.TLorentzVector)()
	tau2  = ROOT.vector(ROOT.TLorentzVector)()
	vistau1 = ROOT.vector(ROOT.TLorentzVector)()
	vistau2 = ROOT.vector(ROOT.TLorentzVector)()

	# create branches

	outTree.Branch("higgs",Higgs )
	outTree.Branch("stau1",stau1 )
	outTree.Branch("stau2",stau2 )
	outTree.Branch("neu1",neu1 )
	outTree.Branch("neu2",neu2 )
	outTree.Branch("met",mets )
	outTree.Branch("tau1",tau1 )
	outTree.Branch("tau2",tau2 )
	outTree.Branch("vistau1",vistau1 )
	outTree.Branch("vistau2",vistau2 )
	Higgs.clear()
	stau1.clear()
	stau2.clear()
	neu1.clear()
	neu2.clear()
	mets.clear()
	tau1.clear()
	tau2.clear()
	vistau1.clear()
	vistau2.clear()

	mt2_c = array('f',[])
	mt2_s = array('f',[])

	decay_chain = {}
	sorted_chain = {}

	# h = ROOT.TH1F('h1','h1',10,-8,8)
	mcweight = []
	sum_of_higgs_last = 0
	nevent = 0

	for entry in xrange(t.GetEntries()):

		# if entry == 20: break
		# print '========================= Standard Prints'
		# print "Number of input events:", t.GetEntries()
		# print "Working on ", file_list[k]
		print  entry*100/t.GetEntries(), "%% complete. %s/%s files"%(k+1,len(file_list))
		print "This is entry number: ", entry
		#========================= Standard Prints
		# logging.info("%s %% complete.  %s/%s files"%(entry*100/t.GetEntries(),k+1,len(file_list)))
    
		t.GetEntry(entry)
		Met = t.MET_Truth
		all_particles = t.TruthParticles
		metvector1 = ROOT.TLorentzVector(0,0,0,0)
		metvector1.SetPtEtaPhiM(Met.get(1).met(), 0, Met.get(1).phi(), 0)
		METp = ROOT.TLorentzVector(0,0,0,0)
		METp.SetPxPyPzE(Met.get(1).mpx(),Met.get(1).mpy(),0,Met.get(1).sumet())


		higgs = [i for i in all_particles if i.absPdgId() in [36,35]]
		higgs_with_childern_as_higgs = []
		higgs_with_children = []
		choosen_higgs = []
		for i  in higgs:
			if i.child(0).pdgId() == i.pdgId():	higgs_with_childern_as_higgs.append(i)
			elif i.child(0).absPdgId() in [1000015,2000015]:  choosen_higgs.append(i)
			else:	higgs_with_children.append(i)
		assert len(choosen_higgs) == 1,  " More than two higgs with immediate stau children."
	
		h.Fill(choosen_higgs[0].child(0).p4().Pt())	
		h1.Fill(choosen_higgs[0].child(1).p4().Pt())	
		h2.Fill(choosen_higgs[0].p4().Pt())

		for i in choosen_higgs:
			for j in xrange(i.nChildren()):
				print " Choosen higgs child: ", j ,"\t", i.child(j).pdgId()
		stau = [ i.child(j) for i in choosen_higgs for j in xrange(i.nChildren())]
		
		assert len(stau) == 2, " There should be two staus."
	
		#check for if di-stau pt is same as higgs pt
		if stau[0].parent().barcode() != choosen_higgs[0].barcode(): continue 
		tau_h	= [] # tau list
		
		n1 = [i for i in find_child(stau[0]) if i.absPdgId() ==1000022]
		n2 = [i for i in find_child(stau[1]) if i.absPdgId() ==1000022]

		tau_1 = [i for i in find_child(stau[0]) if i.absPdgId() ==15]
		tau_2 = [i for i in find_child(stau[1]) if i.absPdgId() ==15]
		


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
			# continue
			tau_h_1 = find_vis(tau_h[0])
			tau_h_2 = find_vis(tau_h[1])
			Higgs.push_back(choosen_higgs[0].p4()) 	
			stau1.push_back(stau[0].p4())
			stau2.push_back(stau[1].p4())
			neu1.push_back(n1[0].p4())
			neu2.push_back(n2[0].p4())
			tau1.push_back(tau_1[0].p4())
			tau2.push_back(tau_2[0].p4())
			mets.push_back(metvector1)
			
			vistau1.push_back(tau_h_1[0])
			vistau2.push_back(tau_h_2[0])
			
			print "PT of higgs that is being picked \t", PT(choosen_higgs) , '\t', 'barcode: ', '\t',choosen_higgs[0].barcode()
			print "PT of higgs that is not picked \t", [i.p4().Pt() for i in higgs_with_children ] , '\t', 'barcode: ', '\t',[ i.barcode() for i in higgs_with_children ]
			print "Immidate children of higgs ", [ (i.pdgId(),i.child(j).pdgId()) for  i in higgs for j in xrange(i.nChildren())]
			print "PT of all the higgs present \t", PT(higgs), '\t', 'barcode: ', '\t',[i.barcode() for i in higgs]
			# print PT(higgs), ' \n Pt higgs '
			print "Pdg higgs \t", pdg(higgs) 
			print "PT of child of all higgs \t" , [ PT(find_child(i)) for i in higgs if i is not None]
			print "Pdgid of child of all higgs \t" , [ pdg(find_child(i)) for i in higgs if i is not None]
			print "barcode  of all the Childen from higgs \t" , [ j.barcode() for i in higgs for j in find_child(i)  if i is not None]
			

			##  MT2
			# combined = float(520E3)
			# single = float(170E3)
			
			# Visap = sorted_chain["tau"][0]
			# Visbp = sorted_chain["tau"][1]
			
			# # print Visap.E() , 'visap'
			# # print METp.E(), " met E"
			# mt2c = ROOT.ComputeMT2(Visap,Visbp,METp,combined,combined)
			# mt2s = ROOT.ComputeMT2(Visap,Visbp,METp,single,single)

			# # print mt2c.Compute() , "mt2 com"
			# # print mt2s.Compute() , 'mt2 single'
			# mt2_c.append(mt2c.Compute())
			# mt2_s.append(mt2s.Compute())

	outTree.Fill()
	outF.Write()
	outF.Close()	

# can = ROOT.TCanvas('','',600,800)
# h.Draw("HIST")
# can.Print("stau1.png")
# h1.Draw("HIST")
# can.Print("stau2.png")
# h2.Draw("HIST")
# can.Print("pthiggs.png")