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
		if p.child(i).pdgId() != p.pdgId():
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
	#takes a list of particles container and returns the pdgId in a list
	if isinstance(n[0],ROOT.xAOD.TruthParticle_v1):
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
	#takes a  particles container and returns the pt in a list
	if isinstance(n[0],ROOT.xAOD.TruthParticle_v1):
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

# sig_file =['/cluster/home/bataju/lxplus_aod/500/Merge.aod.pool.root ']	#500
# sig_file.append('/cluster/home/bataju/lxplus_aod/500/Merge.10k.pool.root ')	#500
# sig_file.append('/cluster/home/bataju/filter/500_aod/merge._20k.pool.root')	#500

# sig_file = ['/cluster/home/bataju/375_800/lxplux_aod/10k/Merge._10k.pool.root']	#800
# sig_file.append('/cluster/home/bataju/375_800/lxplux_aod/20k/Merge._20k.pool.root')	#800

# sig_file = ['/cluster/home/bataju/compare/eventgeneration/different_stau_mass/1500/TRUTH1/DAOD_TRUTH1.test.pool.root']	#1500
# sig_file = ['/cluster/home/bataju/compare/eventgeneration/different_stau_mass/1000/TRUTH1/DAOD_TRUTH1.test.pool.root']	#1000
# sig_file = ['/cluster/home/bataju/compare/eventgeneration/different_stau_mass/2000/TRUTH1/DAOD_TRUTH1.test.pool.root']	#2000
# temp = [1500,1000,2000]
# sig_file = ['/cluster/home/bataju/compare/eventgeneration/different_stau_mass/'+str(i)+'/TRUTH1/DAOD_TRUTH1.test.pool.root' for i in temp]


# sig_file = ['/cluster/home/bataju/BENCHMARK/700/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 700
# sig_file = ['/cluster/home/bataju/BENCHMARK/800/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 800
# sig_file = ['/cluster/home/bataju/BENCHMARK/1000/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 1000
# sig_file = ['/cluster/home/bataju/BENCHMARK/1500/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 1500
# sig_file = ['/cluster/home/bataju/BENCHMARK/2000/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 2000
# sig_file = ['/cluster/home/bataju/BENCHMARK/2500/truth/DAOD_TRUTH1.test.pool.root']	#BenchMark 2500

# temp = [700,730,800,830,930,1000,1130,1280,1500,1700,2000,2200,2500]
# temp = [1500,1700,2000,2200,2500]
temp = [1000]
sig_file = ['/cluster/home/bataju/BENCHMARK/'+str(i)+'/truth/DAOD_TRUTH1.test.pool.root' for i in temp]

# temp = [250,300,500,650,850,1050,1450,2250]
# sig_file = ['/cluster/home/bataju/PAPER/'+str(i)+'/truth/DAOD_TRUTH1.test.pool.root' for i in temp]

# tautau = ['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341881.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH600_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4482_e5984_a875_r10724_r10726_p3759/DAOD_HIGG4D4.17073421._00000{}.pool.root.1'.format(i+1) for i in xrange(6)] 	#600
# tautau = ['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.341885.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH1000_yb2_tautauhh.deriv.DAOD_HIGG4D4.e4298_e5984_a875_r10201_r10210_p3749/DAOD_HIGG4D4.17259964._00000{}.pool.root.1'.format(i+1) for i in xrange(6)] #1000
# tautau = ['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.345288.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH2000_yb2_tautauhh.deriv.DAOD_HIGG4D4.e5686_e5984_a875_r10724_r10726_p3759/DAOD_HIGG4D4.17076570._00000{}.pool.root.1'.format(i+1) for i in xrange(6)]	#2000
# tautau = ['/cluster/home/bataju/tautau/DAOD/mc16_13TeV.345292.aMcAtNloPythia8EvtGen_A14NNPDF23LO_bbH2500_yb2_tautauhh.deriv.DAOD_HIGG4D4.e5686_e5984_a875_r10724_r10726_p3759/DAOD_HIGG4D4.17073341._00000{}.pool.root.1'.format(i+1) for i in xrange(6)]	#2500

#bkg 
# ttbar   =   ['/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.12458455._000337.pool.root.1']
# ttbar.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.12458455._004964.pool.root.1')
# ttbar.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.14761208._000530.pool.root.1')
# ttbar.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.410470.PhPy8EG_A14_ttbar_hdamp258p75_nonallhad.recon.AOD.e6337_s3126_r9364/AOD.14761208._004126.pool.root.1')

# ztautau =     ['/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182601._001043.pool.root.1']
# ztautau.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182601._001055.pool.root.1')
# ztautau.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182673._001851.pool.root.1')
# ztautau.append('/cluster/home/bataju/stau_analaysis/aod_files/mc16_13TeV.361108.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Ztautau.recon.AOD.e3601_s3126_r9364/AOD.11182673._003111.pool.root.1')

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
	outF = ROOT.TFile.Open("H_2stau_2tau_O_{}.root".format(temp[k]), "RECREATE")
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

	mt2_c = array('f',[])
	mt2_s = array('f',[])

	decay_chain = {}
	sorted_chain = {}

	# h = ROOT.TH1F('h1','h1',10,-8,8)
	mcweight = []
	sum_of_higgs_last = 0
	nevent = 0

	for entry in xrange(t.GetEntries()):
		logging.basicConfig(filename='{}B.log'.format(temp[k]), filemode='wb',level=logging.DEBUG, force =True)
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
		assert len(higgs) ,"NO HIGGS BOSONS FOUND!!!"
		# print '		Heavy Higgs Found: 				%s'% pdg(higgs)
		# print '		PT of the HIGGS:				%s'%PT(higgs)
		# print '		Mass of the HIGGS:				%s'%MASS(higgs)
		# print '		Choosing the last higgs			'
		# print [find_child(i) for i in higgs]
		# print reversed([i.p4().Pt() for i in higgs])
		# for i in higgs:
		# 	print i.child(0).pdgId()
		# 	if i.child(0).pdgId() == i.pdgId():
		# 		print "\t" ,i.child(0).child(0).pdgId()
		print "__________________________________________"
		higgs_with_children = []
		print higgs
		
		for i  in higgs:
			if i.child(0).pdgId() == i.pdgId():	continue
			elif i.child(0).absPdgId() in [1000015,2000015]: choosen_higgs = [i]
			else:	higgs_with_children.append(i)
		assert len(choosen_higgs) == 1

		
		print "PT of higgs that is being picked \t", PT(choosen_higgs) , '\t', 'barcode: ', '\t',choosen_higgs[0].barcode()
		print "PT of higgs that is not picked \t", [i.p4().Pt() for i in higgs_with_children ] , '\t', 'barcode: ', '\t',[ i.barcode() for i in higgs_with_children ]
		print "Immidate children of higgs ", [ (i.pdgId(),i.child(j).pdgId()) for  i in higgs for j in xrange(i.nChildren())]
		print "PT of all the higgs present \t", PT(higgs), '\t', 'barcode: ', '\t',[i.barcode() for i in higgs]
		# print PT(higgs), ' \n Pt higgs '
		print "Pdg higgs \t", pdg(higgs) 
		print "PT of child of all higgs \t" , [ PT(find_child(i)) for i in higgs if i is not None]
		print "Pdgid of child of all higgs \t" , [ pdg(find_child(i)) for i in higgs if i is not None]
		print "barcode  of all the Childen from higgs \t" , [ j.barcode() for i in higgs for j in find_child(i)  if i is not None]

		# Higgs.push_back(choosen_higgs[0].p4()) 
		# this a where to pick my Higgs, and i think that was a mistake 

		decay_chain["higgs"]=[choosen_higgs]
		sum_of_higgs_last +=1
		# print '		Child of the choosen_higgs:			%s' %pdg(choosen_higgs)
		# print '		Mass of the child of choosen_higgs:	%s' %MASS(choosen_higgs)
		# print '		PT of the child of choosen_higgs:	%s' %PT(choosen_higgs)
		# print PT(choosen_higgs)
		# print pdg(choosen_higgs)
		print "Choosen HIggs:", choosen_higgs
		# stau = [i for i in choosen_higgs if i.absPdgId() in [1000015,2000015]]
		print "child from the choosen higgs", choosen_higgs[0].child(0).pdgId(), choosen_higgs[0].child(1).pdgId()
		h.Fill(choosen_higgs[0].child(0).p4().Pt())	
		h1.Fill(choosen_higgs[0].child(1).p4().Pt())	
		h2.Fill(choosen_higgs[0].p4().Pt())

		for i in choosen_higgs:
			for j in xrange(i.nChildren()):
				print " Choose higgs child", i.child(j).pdgId()
		stau = [ i.child(j) for i in choosen_higgs for j in xrange(i.nChildren())]
		
		print pdg(stau)
		print "Di-stau Pt of all the stau", (stau[0].p4() +stau[1].p4()).Pt()
		print "Di-stau mass ", (stau[0].p4() +stau[1].p4()).M()

		print "parent of stau [0] \t" , stau[0].parent().pdgId() , "barcode \t" , stau[0].parent().barcode()
		print "parent of stau [1] \t" , stau[1].parent().pdgId() , "barcode \t"  , stau[1].parent().barcode()
		print [i.barcode() == stau[0].parent().barcode() for i in higgs ]
			
		
		# assert False
		assert len(stau) == 2, "There should be two staus!!!"
		decay_chain["stau"]=stau
		
		# assert False
		if stau[0].parent().barcode() != choosen_higgs[0].barcode(): continue #check for if di-stau pt is same as higgs pt
		tau_h	= [] # tau list

		Higgs.push_back(choosen_higgs[0].p4()) 
			
		stau1.push_back(stau[0].p4())
		stau2.push_back(stau[1].p4())
		
		n1 = [i for i in find_child(stau[0]) if i.absPdgId() ==1000022]
		n2 = [i for i in find_child(stau[1]) if i.absPdgId() ==1000022]
		
		decay_chain["nue1"]=[n1[0],n2[0]]
		tau_1 = [i for i in find_child(stau[0]) if i.absPdgId() ==15]
		tau_2 = [i for i in find_child(stau[1]) if i.absPdgId() ==15]
		
		dPhi_tau1_tau2	= array('f',[])
		t_h1 = [] #hadronic taus
		t_h2 = [] #hadronic taus
		# print '		Child of tau1: %s ,	Child tau2: %s'% (pdg(find_child(tau_1[0])),pdg(find_child(tau_2[0])))
		### checking for hadronic taus
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
			# tau_h.sort(reverse=True, key=lambda x: x.p4().Pt()) 
			print "Found a hadronic decay h h ####"
			logging.info("		This is entry number: %s"%entry)
			logging.info("		Working on %s"% file_list[k])
			logging.info('		***Found a hadronic decay***		')
			logging.info('		Higgs Found: 					%s' %pdg(higgs))
			logging.info('		PT of the HIGGS:				%s' %PT(higgs))
			logging.info('		Mass of the HIGGS:				%s' %MASS(higgs))
			logging.info('		Child of the choosen_higgs:			%s'%pdg(choosen_higgs))
			logging.info('		PT of the child of choosen_higgs:	%s'%PT(choosen_higgs))
			logging.info('		Mass of the child of choosen_higgs Stau mass:	%s'%MASS(choosen_higgs))
   			logging.info('		***Found a hadronic decay***	')
			logging.info('		Child of stau[0]: %s'%pdg(find_child(stau[0])))
			logging.info('		Child of stau[1]: %s'%pdg(find_child(stau[1])))

			logging.info('		PT of Child of first stau[0]:-	n1: %s	tau pt: %s'%((PT(n1),PT(tau_1))) )
			logging.info('		Mass of Child of first stau[0]:-	n1: %s	tau mass: %s'%(MASS(n1),MASS(tau_1)) )
			logging.info('		PT of Child of second stau[1]:-	n2: %s	tau pt: %s'%(PT(n2),PT(tau_2)) )
			logging.info('		Mass of Child of second stau[1]:-	n2: %s	tau mass: %s'%(MASS(n2),MASS(tau_2) ))
			logging.info('		Child of tau1: %s ,	Child tau2: %s'%(pdg(find_child(tau_1[0])),pdg(find_child(tau_2[0]))))
			
			tau_h_1 = find_vis(tau_h[0])
			tau_h_2 = find_vis(tau_h[1])
			
			decay_chain['tau']=[tau_h_1,tau_h_2]
			
			# print decay_chain.keys()
			# print [Apdg(i) for i in decay_chain.values()], " pdgId visible protion of tau so no pdgId"			
			# print [PT(i) for i in decay_chain.values()], " before sorting "			
			
			if decay_chain['tau'][0].Pt() > decay_chain['tau'][1].Pt():
				sorted_chain = dict(decay_chain)
			elif decay_chain['tau'][0].Pt() < decay_chain['tau'][1].Pt():
				sorted_chain['tau'] = list(reversed(decay_chain['tau']))
				sorted_chain['nue1'] = list(reversed(decay_chain['nue1']))
				sorted_chain['stau'] = list(reversed(decay_chain['stau']))
				sorted_chain['higgs'] = list(reversed(decay_chain["higgs"]))
			# print [PT(i) for i in sorted_chain.values()], " after sorting "
			# Higgs.push_back(choosen_higgs[0].p4()) 
			
			# stau1.push_back(sorted_chain['stau'][0].p4())
			# stau2.push_back(sorted_chain['stau'][1].p4())

			neu1.push_back(sorted_chain['nue1'][0].p4())
			neu2.push_back(sorted_chain['nue1'][1].p4())
			nevent +=1

			tau1.push_back(tau_h[0].p4())
			tau2.push_back(tau_h[1].p4())

			vistau1.push_back(sorted_chain['tau'][0])
			vistau2.push_back(sorted_chain['tau'][1])
			
			mets.push_back(metvector1)
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
		
	print "total number of events : \t", sum_of_higgs_last
	print "number of passed events: \t", nevent
	outF.Write()
	outF.Close()	
	print "Finished."

can = ROOT.TCanvas('','',600,800)
h.Draw("HIST")
can.Print("stau1.png")
h1.Draw("HIST")
can.Print("stau2.png")
h2.Draw("HIST")
can.Print("pthiggs.png")