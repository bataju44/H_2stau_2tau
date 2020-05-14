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
temp = [1700]
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


count 			= 0
skiped 			= 0
outF = ROOT.TFile.Open("H_2stau_2tau_old_1700.root", "RECREATE")
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
	pt_tau1 		= array('f',[])
	pt_tau2 		= array('f',[])
	mt_tau1_tau2	= array('f',[])
	mt_tau1_met		= array('f',[])
	mt_tau2_met		= array('f',[])
	met_Tr 			= array('f',[])
	mt_tot			= array('f',[])
	phi_tau1		= array('f',[])
	phi_tau2		= array('f',[])
	phi_met			= array('f',[])
	eta_tau1		= array('f',[])
	eta_tau2		= array('f',[])
	dR_tau1_tau2	= array('f',[])
	dPhi_tau1_tau2	= array('f',[])
	M_tau1_tau2 	= array('f',[])
	num_of_tau  	= array('f',[])
	num_of_t_th 	= array('f',[])
	MET_phi = array('f',[])
	px_t1 = array('f',[])
	py_t1 = array('f',[])
	px_t2 = array('f',[])
	py_t2 = array('f',[])
	mt_t1 = array('f',[])
	mt_t2 = array('f',[])
	met_px = array('f',[])
	met_py = array('f',[])
	met_tau0_dphi = array('f',[])
	met_tau1_dphi = array('f',[])
	Mt_2 = array('f',[])
	higgs_mass = array('f',[])
	stau_mass = array('f',[])
	leading_tau = array('f',[])
	sleading_tau = array('f',[])
	dR_stau1_tau1 = array('f',[])
	dR_stau2_tau2 = array('f',[])
	mass_higgs = array('f',[])
	mass_stau1 = array('f',[])
	mass_stau2 = array('f',[])
	pt_n1 = array('f',[])
	pt_n2 = array('f',[])
	pt_stau1 = array('f',[])
	pt_stau2 = array('f',[])
	pt_higgs = array('f',[])
	eta_stau1 = array('f',[])
	eta_stau2 = array('f',[])
	leading_neu1 = array('f',[])
	sleading_neu1 = array('f',[])
	dR_ltau_neu1 = array('f',[])
	dR_stau_neu1 = array('f',[])
	leading_stau = array('f',[])
	sleading_stau = array('f',[])
	tau1_px = array('f',[])
	tau2_py = array('f',[])
	tau1_py = array('f',[])
	tau2_px = array('f',[])
	mass_tau1 = array('f',[])
	mass_tau2 = array('f',[])
	total_nue_pt = array('f',[])
	total_nue_px = array('f',[])
	total_nue_py = array('f',[])
	decay_chain = {}
	sorted_chain = {}

	
	for entry in xrange(t.GetEntries()):
		logging.basicConfig(filename='{}B.log'.format(temp[k]), filemode='wb',level=logging.DEBUG, force =True)
		# if entry == 50: break
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
		higgs = [i for i in all_particles if i.absPdgId() in [36,35]]
		assert len(higgs) ,"NO HIGGS BOSONS FOUND!!!"
		print '		Heavy Higgs Found: 				%s'% pdg(higgs)
		print '		PT of the HIGGS:				%s'%PT(higgs)
		print '		Mass of the HIGGS:				%s'%MASS(higgs)
		print '		Choosing the last higgs			'
		
		higgs_last = find_child(higgs[-1])
  		assert higgs_last is not None, " This higgs has not children"
		decay_chain["higgs"]=[higgs[-1]]
		print '		Child of the higgs[-1]:			%s' %pdg(higgs_last)
		print '		Mass of the child of higgs[-1]:	%s' %MASS(higgs_last)
		print '		PT of the child of higgs[-1]:	%s' %PT(higgs_last)
		
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
		print '		Child of tau1: %s ,	Child tau2: %s'% (pdg(find_child(tau_1[0])),pdg(find_child(tau_2[0])))
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
			logging.info('		Child of the higgs[-1]:			%s'%pdg(higgs_last))
			logging.info('		PT of the child of higgs[-1]:	%s'%PT(higgs_last))
			logging.info('		Mass of the child of higgs[-1] Stau mass:	%s'%MASS(higgs_last))
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
			print decay_chain.keys()
			print [Apdg(i) for i in decay_chain.values()], " pdgId visible protion of tau so no pdgId"			
			print [PT(i) for i in decay_chain.values()], " before sorting "			
			if decay_chain['tau'][0].Pt() > decay_chain['tau'][1].Pt():
				sorted_chain = dict(decay_chain)
			elif decay_chain['tau'][0].Pt() < decay_chain['tau'][1].Pt():
				sorted_chain['tau'] = list(reversed(decay_chain['tau']))
				sorted_chain['nue1'] = list(reversed(decay_chain['nue1']))
				sorted_chain['stau'] = list(reversed(decay_chain['stau']))
				sorted_chain['higgs'] = list(reversed(decay_chain["higgs"]))
			print [PT(i) for i in sorted_chain.values()], " after sorting "

			stau1.push_back(stau[0].p4())
			stau2.push_back(stau[1].p4())

			pt_higgs.append(PT(higgs_last)[-1]) # PT of higgs 
			pt_stau1.append(stau[0].p4().Pt()) # PT of stau1
			pt_stau2.append(stau[1].p4().Pt())	# PT of stau2
			pt_n1.append(n1[0].p4().Pt())	# PT of n1
			pt_n2.append(n2[0].p4().Pt())	# PT of n2

			logging.info('		Visible part of  tau1: %s	tau2: %s'%(tau_h_1.Pt(),tau_h_2.Pt()))
			M_tau1tau2 = (tau_h_1+tau_h_2).M()

			tau1_px.append(tau_h[0].p4().Px())
			tau1_py.append(tau_h[1].p4().Py())
			tau2_px.append(tau_h_2.Px())
			tau2_py.append(tau_h_2.Py())
			met_px.append(Met.get(1).mpx())
			met_py.append(Met.get(1).mpy())
			mass_tau1.append(tau_h[0].p4().M())
			mass_tau2.append(tau_h[1].p4().M())

			print metvector1.Pt() , 'met pt'
			logging.info('		MET[1]:		%s'%metvector1.Pt())
			logging.info('		End of Decay')
			mt_tau1andtau2 = mT(tau_h_1,tau_h_2)
			mt_tau2andmet  = mT_(tau_h_2,metvector1)
			mt_tau1andmet  = mT_(tau_h_1,metvector1)
			total_mass_rec = math.sqrt( mt_tau1andtau2**2 + mt_tau2andmet**2 + mt_tau1andmet**2 )
			pt_tau1.append(tau_h_1.Pt())
			pt_tau2.append(tau_h_2.Pt())
			M_tau1_tau2.append((tau_h_1+tau_h_2).M())
			mt_tau1_tau2.append(mt_tau1andtau2)
			mt_tau1_met.append(mt_tau1andmet)
			mt_tau2_met.append(mt_tau2andmet)
			mt_tot.append(total_mass_rec)
			met_Tr.append(metvector1.Pt())
			phi_tau1.append(tau_h_1.Phi())
			phi_tau2.append(tau_h_2.Phi())
			eta_tau1.append(tau_h_1.Eta())
			eta_tau2.append(tau_h_2.Eta())
			dR_tau1_tau2.append(DeltaR(tau_h_1,tau_h_2))
			dPhi_tau1_tau2.append(DeltaPhi(tau_h_1,tau_h_2))
			dR_stau1_tau1.append(DeltaR(tau_h_1,stau[0].p4()))
			dR_stau2_tau2.append(DeltaR(tau_h_2,stau[0].p4()))
			mass_higgs.append(higgs[-1].p4().M())
			mass_stau1.append(stau[0].p4().M())
			mass_stau2.append(stau[1].p4().M())
			# tau_h.sort(reverse=True, key=lambda x: x.p4().Pt())  #sorting
			# l_tau = find_vis(tau_h[0])
			# s_tau = find_vis(tau_h[1])
			# leading_tau.append(l_tau.Pt())
			# sleading_tau.append(s_tau.Pt())
			leading_tau.append(sorted_chain["tau"][0].Pt())
			sleading_tau.append(sorted_chain["tau"][1].Pt())
			leading_neu1.append(sorted_chain["nue1"][0].p4().Pt())
			sleading_neu1.append(sorted_chain["nue1"][1].p4().Pt())
			dR_ltau_neu1.append(sorted_chain["tau"][0].DeltaR(sorted_chain["nue1"][0].p4()))
			dR_stau_neu1.append(sorted_chain["tau"][1].DeltaR(sorted_chain["nue1"][1].p4()))
			leading_stau.append(sorted_chain["stau"][0].p4().Pt())
			sleading_stau.append(sorted_chain["stau"][1].p4().Pt())
		
			t_n = (sorted_chain["nue1"][0].p4()+sorted_chain["nue1"][1].p4())
			total_nue_pt.append(t_n.Pt())
			total_nue_px.append(t_n.Px())
			total_nue_py.append(t_n.Py())
			tau1_px.append(tau_h[0].p4().Px())
			tau1_py.append(tau_h[1].p4().Py())
			tau2_px.append(tau_h_2.Px())
			tau2_py.append(tau_h_2.Py())
			met_px.append(Met.get(1).mpx())
			met_py.append(Met.get(1).mpy())
			mass_tau1.append(tau_h[0].p4().M())
			mass_tau2.append(tau_h[1].p4().M())
	
	outTree.Fill()
	Higgs.clear()
	stau1.clear()
	stau2.clear()
outF.Write()
outF.Close()	
print "Finished."
	# df = [leading_tau,sleading_tau,leading_neu1,sleading_neu1,dR_ltau_neu1,dR_stau_neu1,leading_stau,sleading_stau,met_Tr,dR_tau1_tau2,tau1_px,tau1_py,tau2_px,tau2_py,met_px,met_py,mass_tau1,mass_tau2,total_nue_pt,total_nue_px,total_nue_py]
	# var = ['leading_tau','sleading_tau','leading_neu1','sleading_neu1','dR_ltau_neu1','dR_stau_neu1','leading_stau','sleading_stau','met_Tr','dR_tau1_tau2','tau1_px','tau1_py','tau2_px','tau2_py','met_px','met_py','mass_tau1','mass_tau2','total_nue_pt','total_nue_px','total_nue_py']
	# data_ = {}
	# for i,j in zip(var,df):
	# 	data_[i]=j
	# with open('{}Bmt2.pkl'.format(temp[k]),'wb') as f:
	# 	pickle.dump(data_, f)


	# dt = [pt_tau1,pt_tau2,met_Tr,mt_tot,leading_tau,sleading_tau,pt_higgs,pt_stau1,pt_stau2,pt_n1,pt_n2,mass_higgs,mass_stau2,mass_stau1,dR_tau1_tau2,dPhi_tau1_tau2,eta_tau1,eta_tau2]
	# variable = ['pt_tau1','pt_tau2','met_Tr','mt_tot','leading_tau','sleading_tau','pt_higgs','pt_stau1','pt_stau2','pt_n1','pt_n2','mass_higgs','mass_stau2','mass_stau1','dR_tau1_tau2','dPhi_tau1_tau2','eta_tau1','eta_tau2']
# 	outfile = ROOT.TFile("H_2stau_Benchmark_{}.root".format(temp[k]),"RECREATE")
# 	ttree = ROOT.TTree("H_2stau_Benchmark_{}".format(temp[k]),"H_2stau_Benchmark_{}".format(temp[k]))
# 	for i in xrange(len(dt)):
# 		ttree.Branch('{}'.format(variable[i]),dt[i],'{}[{}]/F'.format(variable[i],len(dt[i])))
# 	for i in xrange(len(dt[0])):
# 		ttree.Fill()

# # print [len(i) for i in dt]
# 	ttree.Write()
# 	outfile.Write()
# 	outfile.Close()
	# data_dict = {}
	# for i in range(len(variable)):
	# 	data_dict[variable[i]]=dt[i]
	# with open('{}Bphi.pkl'.format(temp[k]),'wb') as f:
	# 	pickle.dump(data_dict, f)
	
