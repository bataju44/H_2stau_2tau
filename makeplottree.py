
import ROOT
import numpy as np
# import pandas as pd
import pickle
import logging
import matplotlib as plt
from array import array
ROOT.gROOT.SetBatch(True)

ROOT.gROOT.LoadMacro("/mnt/c/Users/Bond007/Downloads/atlasstyle/AtlasStyle.C")
ROOT.gROOT.LoadMacro("/mnt/c/Users/Bond007/Downloads/atlasstyle/AtlasUtils.C")
ROOT.gROOT.LoadMacro("/mnt/c/Users/Bond007/Downloads/atlasstyle/AtlasLabels.C")
ROOT.SetAtlasStyle()
# ROOT.gStyle.SetOptTitle(1)
# ROOT.gStyle.SetOptStat()
# Initialize the xAOD infrastructure:
# if(not ROOT.xAOD.Init().isSuccess()): print "Failed xAOD.Init()"


ROOT.gStyle.SetOptStat(1111)
ROOT.gStyle.SetPalette(57)
def open_pickle(f,name):
    with open(f) as file:
        name = pickle.load(file)
    return name

def fill_hist(hist,array):
	[hist.Fill(_) for _ in array]

def scale_hist(hist_bkg_list,norm=1):
	if hist_bkg_list.Integral() != 0 :
		hist_bkg_list.Scale((norm/(hist_bkg_list.Integral())))
	# scale = norm/hist_bkg_list.Integral()
	# hist_bkg_list.Scale(scale)

def H_mk(name,bin,minm,mixm):
	return ROOT.TH1F(name, "",bin,minm,mixm)

#bkg
with open('Ztautau_NONSORT', 'rb') as f:
    ztau = pickle.load(f)
with open('ttbar_NONSORT', 'rb') as f:
    ttbar = pickle.load(f)

ttbar_dict = {}
ztau_dict = {}
# for i in xrange(len(variable)):
#     ttbar_dict[variable[i]] = ttbar[i]
#     ztau_dict[variable[i]] = ztau[i]
# T = [800,2200]
#Benchmark
# Benchmark_list = [700,730,800,830,930,1000,1130,1500,2000,2500,1280,1700,2200]
# Benchmark_list = [730,800,830,930,1000,1130,1500,2000,2500,1700,2200]
# Benchmark_list=[T[0]]
# Benchmark_list=[2200,2500,1700,1500]
# Benchmark_list=[1000,2000,2500,1500]
# Benchmark_list = ["taulike1000"]
# tau_list       = Benchmark_list[:]
Benchmark_list = [800]
#PAPER
# Paper_list = [250,300,500,650,850,1050,1450,2250]
Paper_list = []
# temp1 = [125,200,300,400,600,1500,1000,2000,2500]
# temp1 = [200,300,400,600,1000,2000,2500]
# tau_list = [1000,2000,2500,1500]
# tau_list = [600]
tau_list = [400]
sig_file = ['/mnt/c/Users/Bond007/Downloads/H_2stau_2tau_'+str(i)+'.root' for i in Benchmark_list]
# sig_file = ['/mnt/c/Users/Bond007/Downloads/'+str(i)+'Bmt2eta.pkl' for i in Benchmark_list]
# sig_file = ['/mnt/c/Users/Bond007/Downloads/DAOD_TRUTH1.test.pool.root']
# tau_file = ['/mnt/c/Users/Bond007/Downloads/tauwhole'+str(i)+'.pkl' for i in tau_list]
tau_file = ['/mnt/c/Users/Bond007/Downloads/H_2tau_'+str(i)+'.root' for i in tau_list]

Benchmark = ['%s GeV Benchmark'%i for i in Benchmark_list]
tau = ['bbH %s GeV (A/H-taus)'%i for i in tau_list]

n = {}
# nt={}
for j in xrange(len(sig_file)):
    n[sig_file[j]] = Benchmark[j]
for j in xrange(len(tau_file)):
    with open(tau_file[j], 'rb') as f:
        # n[tau[j]] = pickle.load(f)
        n[tau_file[j]] = tau[j]



can = ROOT.TCanvas("can","can",800,600)
s3 =  ROOT.THStack("","")
l2 = ROOT.TLegend(0.65,0.65,0.90,0.90)
l2.SetBorderSize(0)
l2.SetFillStyle(0)
l2.SetTextSize(0.035)
h2 = ROOT.TH1F('#Delta_{#Phi}(stau)','h1',30,-6,6)
h9 = ROOT.TH1F('#Delte_{#Phi}(leading stau,di-stau)','h1',30,-6,6)
h10 = ROOT.TH1F('#Delte_{#Phi}(sub leading stau,di-stau)','h1',30,-6,6)
h3 = ROOT.TH1F('#Delta_{R}(stau)','h1',30,-2,6)
h1 = ROOT.TH1F('#Delta_{R}(distau,higgs) ','h1',30,-2,6)
h7 = ROOT.TH1F('leading_tau_distau_dR','h1',30,-2,6)
h8 = ROOT.TH1F('subleading_tau_distau_dR','h1',30,-2,6)
h4 = ROOT.TH1F('Higgs Pt','h1',30,0,        1600E3)
h5 = ROOT.TH1F('sub leading Stau','h1',30,0,1600E3)
h6 = ROOT.TH1F('Di-stau pt','h1',30,0,      1600E3)
h11 = ROOT.TH1F('Leading Stau Energy','h1',30,0,      1600E3)
h12 = ROOT.TH1F('Subleading Stau Energy','h1',30,0,   1600E3)
h0 = ROOT.TH1F('leading tau','h1',30,0, 1000E3)
h = ROOT.TH1F('leading stau','h1',30,0, 1000E3)
l2.AddEntry(h,str(Benchmark_list[0])+"GeV Benchmark",'lf')
l2.AddEntry(h0,str(tau_list[0])+"GeV bbHtautau",'lf')
h13 = ROOT.TH1F('(higgs.at(z).M())**2  + (higgs.at(z).Pt())**2 - (stau1.at(z).Pt())**2 - (stau2.at(z).Pt())**2)','h1',30,0,   1600E3)

print n
for files, masspoint in n.items():
    print masspoint
    if 'Benchmark'  in masspoint:
        f = ROOT.TFile.Open(files,"READONLY")

        for var in  f.Get("T"):
            # print len(var.stau1)
            for z in xrange(len(var.stau1)):
                print masspoint
                stau1 = var.stau1
                stau2 = var.stau2
                higgs = var.higgs
                neu1 =  var.neu1
                neu2 =  var.neu2
                met = var.met
                # print len(higgs)
                # print higgs[10].Pt()
                # print [i.Pt() for  i in higgs]
                h4.Fill(higgs.at(z).Pt())
                # h.Fill(stau1.at(z).Pt())
                print "\n______ Event no. {} _________\n".format(z)
                print higgs.at(z).M(), "\t Higgs mass"
                print higgs.at(z).Pt() , "\t Higgs Pt"
                print higgs.at(z).Phi() , "\t Higgs Phi"
                print higgs.at(z).Pz(), "\t Higgs Pz\n"

                h2.Fill(stau1.at(z).DeltaPhi(stau2.at(z)))
                h3.Fill(stau1.at(z).DeltaR(stau2.at(z)))
                # h.Fill(stau1.at(z).DeltaR(stau2.at(z)))

                print (stau1.at(z) + stau2.at(z)).M(), "\t di stau mass"
                print (stau1.at(z) + stau2.at(z)).Pt(), "\t di stau Pt"
                print (stau1.at(z) + stau2.at(z)).Phi(), "\t di stau Phi"
                print (stau1.at(z) + stau2.at(z)).Pz(), "\t di stau Pz\n"

                print neu1.at(z).M(), "\t neutralino1  mass"
                print neu1.at(z).Pt() , "\t neutralino1  Pt"
                print neu1.at(z).Phi() , "\t neutralino1  Phi"
                print neu1.at(z).Pz(), "\t neutralino1  Pz\n"

                print neu2.at(z).M(), "\t neutralinos2  mass"
                print neu2.at(z).Pt() , "\t neutralinos2  Pt"
                print neu2.at(z).Phi() , "\t neutralinos2  Phi"
                print neu2.at(z).Pz(), "\t neutralinos2  Pz\n"

                h6.Fill((stau1.at(z) + stau2.at(z)).Pt())
                # h.Fill((stau1.at(z) + stau2.at(z)).Pt())
                h1.Fill(higgs.at(z).DeltaR((stau1.at(z) + stau2.at(z))))
                h7.Fill(stau1.at(z).DeltaR((stau1.at(z) + stau2.at(z))))
                h8.Fill(stau2.at(z).DeltaR((stau1.at(z) + stau2.at(z))))
                # h9.Fill(stau1.at(z).DeltaPhi(higgs.at(z)))
                # h10.Fill(stau2.at(z).DeltaPhi(higgs.at(z)))
                h9.Fill(stau1.at(z).Phi())
                h10.Fill(stau1.at(z).Phi())
                # print stau1.at(z).M(), "\t stau1 mass"
                # print stau2.at(z).M(), "\t stau2 mass"
                # print higgs.at(z).M() ,"\t higgs mass"
                if stau1.at(z).Pt() > stau2.at(z).Pt():
                    h11.Fill(stau1.at(z).Energy())
                    print stau1.at(z).M() ,"\t Leading stau mass"
                    print stau1.at(z).Pt() ,"\t Leading stau Pt"
                    print stau1.at(z).Phi() ,"\t Leading stau Phi"
                    print stau1.at(z).Pz() ,"\t Leading stau Pz\n"
                    h.Fill(stau2.at(z).Pt())
                    h5.Fill(stau2.at(z).Pt())
                    h12.Fill(stau1.at(z).Energy())
                    print stau2.at(z).M() ,"\t sub Leading stau M"
                    print stau2.at(z).Pt() ,"\t sub Leading stau Pt"
                    print stau2.at(z).Phi() ,"\t sub Leading stau Phi"
                    print stau2.at(z).Pz() ,"\t sub Leading stau Pz\n"
                else:
                    h11.Fill(stau1.at(z).Energy())
                    print stau2.at(z).M(),"\t Leading stau mass"
                    print stau2.at(z).Pt(),"\t Leading stau Pt"
                    print stau2.at(z).Phi(),"\t Leading stau Phi"
                    print stau2.at(z).Pz(),"\t Leading stau Pz\n"
                    h.Fill(stau1.at(z).Pt())
                    h5.Fill(stau1.at(z).Pt())
                    h12.Fill(stau1.at(z).Energy())
                    print stau1.at(z).M() ,"\t sub Leading stau M"
                    print stau1.at(z).Pt() ,"\t sub Leading stau Pt"
                    print stau1.at(z).Phi() ,"\t sub Leading stau Phi"
                    print stau1.at(z).Pz() ,"\t sub Leading stau Pz\n"


                # print np.sqrt( (higgs.at(z).M()/2)**2 - (stau1.at(z).Pt())**2 - (higgs.at(z).Pt())**2 ),  "\t ( (mass of higgs /2)^2 - (pt of subleading tua)^2 - (pt of higgs)^2 ) ^ (1/2) \n"
                print np.sqrt( (higgs.at(z).M())**2  + (higgs.at(z).Pt())**2 - (stau1.at(z).Pt())**2 - (stau2.at(z).Pt())**2),  "\t ( (mass of higgs)^2 + (pt of higgs)^2 - (pt of subleading tua)^2 -(pt of leading tua)^2  ) ^ (1/2) \n"
                print stau1.at(z).DeltaPhi(stau2.at(z)), "\t DeltaPhi Stau \n"
                h13.Fill(np.sqrt( (higgs.at(z).M())**2  + (higgs.at(z).Pt())**2 - (stau1.at(z).Pt())**2 - (stau2.at(z).Pt())**2))
                
        scale_hist(h,norm=1)
        s3.Add(h)
        # print var
        # print files.split('/')[-1].split('.')[0]

    elif 'bbH' in masspoint:
        # continue
        print masspoint
        f = ROOT.TFile.Open(files,"READONLY")

        for var in  f.Get("T"):
            # print len(var.stau1)

            for z in xrange(len(var.tau1)):
                    higgs = var.higgs
                    mets = var.met
                    tau1 = var.tau1
                    tau2 = var.tau2
                    vistau1 = var.vistau1
                    vistau2 = var.vistau2
                    # h0.Fill(tau1.at(z).Pt())
                    # h0.Fill(tau2.at(z).Pt())
                    # h0.Fill((tau1.at(z) + tau2.at(z)).Pt())
                    # h0.Fill(tau1.at(z).DeltaR(tau2.at(z)))
                    print masspoint
                    print "\n______ Event no. {} _________\n".format(z)
                    print higgs.at(z).M(), "\t Higgs mass"
                    print higgs.at(z).Pt(), "\t Higgs Pt"
                    print higgs.at(z).Phi(), "\t Higgs Phi"
                    print higgs.at(z).Pz(), "\t Higgs Pz"

                    print (tau1.at(z) + tau2.at(z)).M(), "\t di tau mass"
                    print (tau1.at(z) + tau2.at(z)).Pt(), "\t di tau Pt"
                    print (tau1.at(z) + tau2.at(z)).Phi(), "\t di tau Phi"
                    print (tau1.at(z) + tau2.at(z)).Pz(), "\t di tau Pz"

                    if tau1.at(z).Pt() > tau2.at(z).Pt():
                        print tau1.at(z).M() ,"\t Leading tau mass"
                        print tau1.at(z).Pt() ,"\t Leading tau Pt"
                        print tau1.at(z).Phi() ,"\t Leading tau Phi"
                        print tau1.at(z).Pz() ,"\t Leading tau Pz"
                        # h11.Fill(tau1.at(z).Energy())
                        # print tau1.at(z).Pt() ,"   Leading tau Pt"
                        # print tau1.at(z).Phi() ,"  Leading stau Phi"
                    else:
                        print tau2.at(z).M() ,"\t Leading tau mass"
                        print tau2.at(z).Pt() ,"\t Leading tau Pt"
                        print tau2.at(z).Phi() ,"\t Leading tau Phi"
                        print tau2.at(z).Pz() ,"\t Leading tau Pz"
                        # h11.Fill(tau1.at(z).Energy())
                        # print tau2.at(z).Pt(),"    Leading tau Pt"
                        # print tau2.at(z).Phi(),"   Leading tau Phi"

                    if tau1.at(z).Pt() < tau2.at(z).Pt():
                        print tau1.at(z).M() ,"\t   sub Leading tau mass"
                        print tau1.at(z).Pt() ,"\t  sub Leading tau Pt"
                        print tau1.at(z).Phi() ,"\t  sub Leading tau Phi"
                        print tau1.at(z).Pz() ,"\t  sub Leading tau Pz"
                        h0.Fill(tau1.at(z).Pt())

                        # h5.Fill(tau1.at(z).Pt())
                        # h12.Fill(tau1.at(z).Energy())
                        # print tau1.at(z).Pt() ,"   sub Leading tau Pt"
                        # print tau1.at(z).Phi() ,"  sub Leading tau Phi"
                    else:
                        h0.Fill(tau2.at(z).Pt())
                        print tau2.at(z).M() ,"\t   sub Leading tau mass"
                        print tau2.at(z).Pt() ,"\t  sub Leading tau Pt"
                        print tau2.at(z).Phi() ,"\t  sub Leading tau Phi"
                        print tau2.at(z).Pz() ,"\t  sub Leading tau Pz"
                        # h5.Fill(stau2.at(z).Pt())
                    print tau1.at(z).DeltaPhi(tau2.at(z)), "\t DeltaPhi tau \n"

        scale_hist(h0,norm=1)
        s3.Add(h0)
        print var
        print files.split('/')[-1].split('.')[0]




# s3.SetMaximum(0.25)
h.Draw("HIST")
# h0.Draw("HIST SAME")
h.GetXaxis().SetTitle("#font[52]{P_{T} leading #tilde{#tau}/leading #tau}")
can.Print("leading_tau{}.png".format(Benchmark_list[0]))
h13.Draw("HIST")
h13.GetXaxis().SetTitle("#font[52](higgs.at(z).M())**2  + (higgs.at(z).Pt())**2 - (stau1.at(z).Pt())**2 - (stau2.at(z).Pt())**2)^1/2 }")
can.Print("Calc_{}.png".format(Benchmark_list[0]))


h11.Draw("HIST")
h11.GetXaxis().SetTitle("#font[52]{Energy leading #tilde{#tau}}")
can.Print("leading_tau_energy_{}.png".format(Benchmark_list[0]))

h12.Draw("HIST")
h12.GetXaxis().SetTitle("#font[52]{Energy sub leading #tilde{#tau}}")
can.Print("subleading_tau_energy_{}.png".format(Benchmark_list[0]))

h7.Draw("HIST")
h7.GetXaxis().SetTitle("#font[52]{#Delte_{R}(leading stau,di-stau)}")
can.Print("leading_tau_distau_dR{}.png".format(Benchmark_list[0]))
h8.Draw("HIST")
h8.GetXaxis().SetTitle("#font[52]{#Delte_{R}(subleading stau,di-stau)}")
can.Print("subleading_tau_distau_dR{}.png".format(Benchmark_list[0]))
h9.Draw("HIST")
h9.GetXaxis().SetTitle("#font[52]{#Phi leading stau)}")
# h9.GetXaxis().SetTitle("#font[52]{#Delte_{#Phi}(leading stau,higgs)}")
can.Print("leading_tau_phi{}.png".format(Benchmark_list[0]))
h10.Draw("HIST")
h10.GetXaxis().SetTitle("#font[52]{#Phi sun leading stau)}")
# h10.GetXaxis().SetTitle("#font[52]{#Delte_{#Phi}(subleading stau,higgs)}")
can.Print("subleading_Phi{}.png".format(Benchmark_list[0]))
h5.Draw("HIST")
h5.GetXaxis().SetTitle("#font[52]{P_{T} subleading #tilde{#tau}/subleading #tau}")
can.Print("subleading_tau{}.png".format(Benchmark_list[0]))
h4.Draw("HIST")
h4.GetXaxis().SetTitle("#font[52]{P_{T} H/A} [GeV]")
can.Print("pt_higgs{}.png".format(Benchmark_list[0]))
h6.Draw("HIST")
h6.GetXaxis().SetTitle("#font[52]{P_{T} di-stau } [GeV]")
can.Print("distau{}.png".format(Benchmark_list[0]))
h2.Draw("HIST")
h2.GetXaxis().SetTitle("#font[52]{#Delta_{R}(stau) } [GeV]")
can.Print("stau_staudphi{}.png".format(Benchmark_list[0]))
h1.Draw("HIST")
h1.GetXaxis().SetTitle("#font[52]{#Delta_{R}(distau,higgs) } [GeV]")
can.Print("distau_dR{}.png".format(Benchmark_list[0]))
h3.Draw("HIST")
h3.GetXaxis().SetTitle("#font[52]{#Delta_{#Phi}(stau) } [GeV]")
can.Print("stau_staudR{}.png".format(Benchmark_list[0]))

s3.Draw("PLC HIST nostack")
s3.GetYaxis().SetTitle("Events")
s3.GetXaxis().SetLabelSize(0.03)
s3.GetYaxis().SetLabelSize(0.02)
# s3.GetXaxis().SetTitle("#font[52]{#Delta_{R} taus and staus} ")
# s3.GetXaxis().SetTitle("#font[52]{leading #tilde{#chi}_{1}^{0}} [GeV]")
# s3.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
# s3.GetXaxis().SetTitle("#font[52]{MT2}")
# s3.GetXaxis().SetTitle("#font[52]{di-tau/ di-stau Pt} [GeV]")
# s3.GetXaxis().SetTitle("#font[52]{P_{T} H/A} [GeV]")
# s3.GetXaxis().SetTitle("#font[52]{P_{T} subleading #tilde{#tau}/subleading #tau}")
# s3.GetXaxis().SetTitle("#font[52]{P_{T} leading #tilde{#tau}/leading #tau}")
s3.GetXaxis().SetTitle("#font[52]{P_{T}  #tilde{#tau}/ #tau}")
# s3.GetXaxis().SetTitle("#font[52]{P_{T} leading #tau}")
l2.Draw()
# can.SetLogy()
can.Print("SLeading_stau_tau{}.png".format(Benchmark_list[0]))

assert False

# for j in xrange(len(sig_file_p)):
#     with open(sig_file_p[j], 'rb') as f:
#         n[Paper[j]] = pickle.load(f)
# n["Ztautatu"]=ztau
# n["ttbar"]=ttbar
# print len(ztau)
# print len(n)
# print len(Benchmark[0])
logging.basicConfig(filename='{}Bhisto.csv'.format(Benchmark_list[0]), filemode='wb',level=logging.DEBUG, force =True)
color = { "500 GeV (A/H-staus)":4,"850 GeV (A/H-staus)":2,"1050 GeV (A/H-staus)": 6,"1450 GeV (A/H-staus)":7, "2250 GeV (A/H-staus)":45,"650 GeV (A/H-staus)":30,"250 GeV (A/H-staus)":11,"300 GeV (A/H-staus)":39,
          "2500 GeV Benchmark":3,"2000 GeV Benchmark":5,"1500 GeV Benchmark":1,"800 GeV Benchmark":4,"1000 GeV Benchmark":46,"2200 GeV Benchmark":400,"700 GeV Benchmark":49,"730 GeV Benchmark":619 ,"830 GeV Benchmark":4,"930 GeV Benchmark":810,"1130 GeV Benchmark":41,"1280 GeV Benchmark":24,"1700 GeV Benchmark":30,
         "bbH 600 GeV (A/H-taus)":46,"bbH 1000 GeV (A/H-taus)":6,"bbH 2000 GeV (A/H-taus)":4,"bbH 2500 GeV (A/H-taus)":7,"bbH 125 GeV (A/H-taus)":880, "bbH 200 GeV (A/H-taus)":2,"bbH 300 GeV (A/H-taus)":840,"bbH 400 GeV (A/H-taus)":39,"bbH 1500 GeV (A/H-taus)":800,
         "Ztautatu":37,"ttbar":926,
         }


# with open('mcweight.pkl',"rb") as f:
# 	mc = pickle.load(f)
# print mc.keys()
# can = ROOT.TCanvas("can","can",800,600)
# s1 =  ROOT.THStack("","")
# l2 = ROOT.TLegend(0.65,0.65,0.90,0.90)
# l2.SetBorderSize(0)
# l2.SetFillStyle(0)
# l2.SetTextSize(0.03)
# print type(mc)
# for files, var in mc.items():
#     h = ROOT.TH1F('h1','h1',11,-2.5,8.5)
#     fill_hist(h,var)
#     scale_hist(h)
#     h.SetFillColor(color['%s GeV Benchmark'%files])
#     s1.Add(h)
#     l2.AddEntry(h,'%s GeV Benchmark'%files,'lf')

# s1.SetMaximum(10)
# s1.Draw("hist ")
# s1.GetYaxis().SetTitle("Events")
# s1.GetXaxis().SetLabelSize(0.03)
# s1.GetYaxis().SetLabelSize(0.02)
# s1.GetXaxis().SetTitle("#font[52]{mcEventWeight} ")
# s1.GetXaxis().SetTitle("#font[52]{#Delta_{R}} ")
# s1.GetXaxis().SetTitle("#font[52]{leading #tilde{#chi}_{1}^{0}} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{MT2} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{P_{T} H/A} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{sub leading #tilde{#tau} #eta/sub leading #tau #eta}")
# s1.GetXaxis().SetTitle("#font[52]{P_{T} sub leading #tau}")
# l2.Draw()
# can.SetLogy()
# can.Print("mcweight.png")

# assert False


# print n['2500 GeV Benchmark'].keys()
# print n['bbH 1000 GeV (A/H-taus)'].keys()
# print n['Ztautatu'].keys()
# assert False
# print len(n['Ztautatu'][0])


for files, var in n.items():
    can = ROOT.TCanvas("can","can",800,600)
    s3 =  ROOT.THStack("","")
    s1 =  ROOT.THStack("","")
    l2 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    l2.SetBorderSize(0)
    l2.SetFillStyle(0)
    l2.SetTextSize(0.03)
    l1 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    l1.SetBorderSize(0)
    l1.SetFillStyle(0)
    l1.SetTextSize(0.03)
    print files
    # sorted(var.keys())
    # if files == "Ztautatu":    continue
    met_150 = []
    met_100 = []
    pt_150 = []
    pt_100 = []
    dr_100 = []
    dr_100_l35 = []
    dr_100_g35 = []
    pt_100_l35 = []
    pt_100_g35 = []
    dr_100_l85 = []
    dr_100_g85 = []
    pt_100_l85 = []
    pt_100_g85 = []
    met_100_g85 = []
    met_100_l85 = []
    print len(var["met_Tr"])
    # print len(var["pt_tau2"])
    for i in xrange(len(var["sleading_tau"])):
        if var["met_Tr"][i] < 100E3:
            met_150.append(var["met_Tr"][i])
            pt_150.append(var['sleading_tau'][i])
        if var["met_Tr"][i] < 150E3:
            met_100.append(var["met_Tr"][i])
            pt_100.append(var['sleading_tau'][i])
            dr_100.append(var['dR_tau1_tau2'][i])
        # dr_100.append(var['dR_tau1_tau2'][i])
        if var['sleading_tau'][i] < 35E3:
            pt_100_l35.append(var['sleading_tau'][i])
            dr_100_l35.append(var['dR_tau1_tau2'][i])
        if var['sleading_tau'][i] > 35E3:
            pt_100_g35.append(var['sleading_tau'][i])
            dr_100_g35.append(var['dR_tau1_tau2'][i])
        if var['leading_tau'][i] < 85E3:
            pt_100_l85.append(var['leading_tau'][i])
            dr_100_l85.append(var['dR_tau1_tau2'][i])
            if var["met_Tr"][i] < 100E3:
                met_100_l85.append(var["met_Tr"][i])
        if var['leading_tau'][i] > 85E3:
            pt_100_g85.append(var['leading_tau'][i])
            dr_100_g85.append(var['dR_tau1_tau2'][i])
            if var["met_Tr"][i] < 100E3:
                met_100_g85.append(var["met_Tr"][i])
    h3 = ROOT.TH1F('h1','h1',20,-1,6)
    fill_hist(h3,var['dR_tau1_tau2'])

    # pt_100_l35 = [i for i in pt_100 if i<35E3 ]
    # pt_100_g35 = [i for i in pt_100 if i>35E3 ]
    TT= [dr_100_l35,dr_100_g35]
    h = ROOT.TH1F('h1','h1',20,-1,6)
    fill_hist(h,dr_100_g35)
    # scale_hist(h)
    # h.SetLineColor()
    h1 = ROOT.TH1F('h1','h1',20,-1,6)
    fill_hist(h1,dr_100_l35)
    # scale_hist(h1)
    # h1.SetLineColor(color[files])
    # h.SetMaximum(10)
    s1.Add(h1)
    s1.Add(h)

    # s1.Add(h3)
    l2.AddEntry(h,"dR, tau>35",'lf')
    l2.AddEntry(h1,"dR, tau<35",'lf')
    # l2.AddEntry(h3,"dR all",'lf')

    logging.info("For mA, total, >85 , <85 >85>100, <85>100 ,{}, {}, {}, {}, {}, {} ".format(files[:4],len(var['leading_tau']),len(pt_100_l85),len(pt_100_g85),len(met_100_l85),len(met_100_g85)))
    s1.Draw("PLC hist ")
    s1.GetYaxis().SetTitle("Events")
    s1.GetXaxis().SetLabelSize(0.03)
    s1.GetYaxis().SetLabelSize(0.02)
    # s1.GetXaxis().SetTitle("#font[52]{mcEventWeight} ")
    s1.GetXaxis().SetTitle("#font[52]{#Delta_{R}} %s"%files[:4])
    # s1.GetXaxis().SetTitle("#font[52]{leading #tilde{#chi}_{1}^{0}} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{MT2} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{P_{T} H/A} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{sub leading #tilde{#tau} #eta/sub leading #tau #eta}")
    # s1.GetXaxis().SetTitle("#font[52]{P_{T} sub leading #tau}")
    l2.Draw()
    # can.SetLogy()
    # can.Print("B_lgplc_stackb_dr_{}.png".format(files[:4]))
# assert False
can = ROOT.TCanvas("can","can",800,600)
s3 =  ROOT.THStack("","")
l1 = ROOT.TLegend(0.65,0.65,0.90,0.90)
l1.SetBorderSize(0)
l1.SetFillStyle(0)
l1.SetTextSize(0.03)
TTn = [ "dR, tau<35", "dR, tau>35"]
j = 0
for i in TT:

    h = ROOT.TH1F('h1','h1',33,-2,8)
    fill_hist(h,i)
    scale_hist(h)
    h.SetLineColor(j+1)
    h.SetMaximum(10)
    s3.Add(h)
    l1.AddEntry(h,TTn[j],'lf')
    j+=1
s3.SetMaximum(10)
s3.Draw("hist")
s3.GetYaxis().SetTitle("Events")
s3.GetXaxis().SetLabelSize(0.03)
s3.GetYaxis().SetLabelSize(0.02)
# s1.GetXaxis().SetTitle("#font[52]{mcEventWeight} ")
s3.GetXaxis().SetTitle("#font[52]{#Delta_{R} %s }"%files[:4])
# s1.GetXaxis().SetTitle("#font[52]{leading #tilde{#chi}_{1}^{0}} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{MT2} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{P_{T} H/A} [GeV]")
# s1.GetXaxis().SetTitle("#font[52]{sub leading #tilde{#tau} #eta/sub leading #tau #eta}")
# s1.GetXaxis().SetTitle("#font[52]{P_{T} sub leading #tau}")
l1.Draw()
can.SetLogy()
# can.Print("B_met_dr_all{}.png".format(Benchmark_list[0]))


# assert False
for files, var in n.items():
    can = ROOT.TCanvas("can","can",800,600)
    s3 =  ROOT.THStack("","")
    s1 =  ROOT.THStack("","")
    l2 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    l2.SetBorderSize(0)
    l2.SetFillStyle(0)
    l2.SetTextSize(0.03)
    l1 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    l1.SetBorderSize(0)
    l1.SetFillStyle(0)
    l1.SetTextSize(0.03)
    ptlead_85 = []
    for i in xrange(len(var["leading_tau"])):
        if var["leading_tau"][i] < 85E3:
            ptlead_85.append(var["leading_tau"][i])
    ptslead_35 = []
    ptmet_35 = []
    ptdr_35 = []
    ptmet_l35 = []
    ptmet_l150 = []
    ptmet_g35 = []
    ptmet_g150 = []
    for i in xrange(len(var["sleading_tau"])):
        if var["sleading_tau"][i] > 35E3:
            ptslead_35.append(var["sleading_tau"][i])
            ptmet_35.append(var["met_Tr"][i])
            ptdr_35.append(var["dR_tau1_tau2"][i])
            if var["met_Tr"][i] < 100E3:
                ptmet_l35.append(var["met_Tr"][i])
            if var["met_Tr"][i] > 100E3:
                ptmet_g35.append(var["met_Tr"][i])
            if var["met_Tr"][i] < 150E3:
                ptmet_l150.append(var["met_Tr"][i])
            if var["met_Tr"][i] > 150E3:
                ptmet_g150.append(var["met_Tr"][i])
    logging.info("For mA,total, tau>35, met<100, met>100,{},{}, {}, {},{}, {}, {}  ".format(files[:4],len(var['sleading_tau']),len(ptslead_35),len(ptmet_l35),len(ptmet_g35),len(ptmet_l150),len(ptmet_g150)))

    ptslead_ = []
    ptmet_ = []
    ptdr_ = []
    for i in xrange(len(var["sleading_tau"])):
        if var["sleading_tau"][i] > 35E3:
            ptslead_.append(var["sleading_tau"][i])
            ptmet_.append(var["met_Tr"][i])
            ptdr_.append(var["dR_tau1_tau2"][i])
    h1 = ROOT.TH1F('h1','h1',30,-2,8)
    fill_hist(h1,ptdr_ )
    # scale_hist(h1)
    # h1.SetLineColor(2)
    # h1.SetMaximum(10)
    s1.Add(h1)
    l1.AddEntry(h1,"tau>35GeV {}".format(len(ptdr_)),'lf')
    # s1.SetMaximum(10)
    h = ROOT.TH1F('h1','h1',30,-2,8)
    fill_hist(h,ptdr_35)
    # scale_hist(h)
    # h.SetLineColor()
    # h.SetMaximum(10)
    s1.Add(h)
    l1.AddEntry(h,"tau<35GeV {}".format(len(ptdr_35)),'lf')

    s1.Draw("plc nostack")
    s1.GetYaxis().SetTitle("Events")
    s1.GetXaxis().SetLabelSize(0.03)
    s1.GetYaxis().SetLabelSize(0.02)
    s1.GetXaxis().SetTitle("#font[52]{#Delta_{R} } %s"%files[:4])
    # s1.GetXaxis().SetTitle("#font[52]{leading #tilde{#chi}_{1}^{0}} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{MT2} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{P_{T} H/A} [GeV]")
    # s1.GetXaxis().SetTitle("#font[52]{sub leading #tilde{#tau} #eta/sub leading #tau #eta}")
    # s1.GetXaxis().SetTitle("#font[52]{P_{T} sub leading #tau}")
    l1.Draw()
    # can.SetLogy()
    # can.Print("B_staudr{}.png".format(files[:4]))
# assert False
can = ROOT.TCanvas("can","can",800,600)
s3 =  ROOT.THStack("","")
l2 = ROOT.TLegend(0.65,0.65,0.90,0.90)
l2.SetBorderSize(0)
l2.SetFillStyle(0)
l2.SetTextSize(0.03)
for files, var in n.items():

    # logging.info("For mA, ltau, 85 ,{}, {}, {} ".format(files,len(var['sleading_tau']),len(pt_35)))

    # h = ROOT.TH1F('h1','h1',30,0,300E3)
    # fill_hist(h,met_100)
    # # scale_hist(h)
    # h.SetLineColor(color[files])
    # s3.Add(h)
    # l2.AddEntry(h,files+str(len(met_100)),'lf')
    # fill_hist(h,l_tau)
    # # scale_hist(h)
    # h.SetLineColor(color[files])
    # s3.Add(h)
    # l2.AddEntry(h,files+str(len(l_tau)),'lf')

    for var_name, var_list in var.items():
        # print var_name, type(var_list)

        # if '(A/H-taus)' in files and var_name== 'pt_higgs':
        #     fill_hist(h,var_list)
        #     # h.SetFillColorAlpha(color[files],0.35)
        #     h.SetLineColor(color[files])

        #     l2.AddEntry(h,files,'lf')
        #     scale_hist(h)
        #     s3.Add(h)
        # if var_name == 'mt2_c' or var_name == 'mt2_s':
        #     # print var_list
        #     fill_hist(h,var_list)
        #     l2.AddEntry(h,files,'lf')
        #     h.SetLineColor(color[files])
        #     scale_hist(h)
        #     s3.Add(h)
        # if files == "Ztautatu" or files == "ttbar":
        #         continue
        # if var_name == 'dR_tau1_tau2':
        #     # print var_list
        #     fill_hist(h,var_list)
        #     l2.AddEntry(h,files,'lf')
        #     h.SetLineColor(color[files])
        #     # scale_hist(h,norm=1)
        #     s3.Add(h)
        # if '(A/H-taus)'  not in files:   continue
        # if var_name == 'leading_tau':
        print var_name
        h = ROOT.TH1F('h1','h1',30,0,1000E3)
        if '(A/H-taus)' in files:
            print var["leading_whole_tau"] == var["sleading_whole_tau"]
        # if var_name == 'pt_higgs':
        if (var_name == 'leading_stau' and 'Benchmark' in files  ) or (var_name == 'leading_whole_tau' and '(A/H-taus)' in files):
            # print var_list
            fill_hist(h,var_list)
            h.SetLineColor(color[files])
            scale_hist(h,norm=1)
            if files == "Ztautatu" or files == "ttbar":
                continue
                # h.SetFillColorAlpha(color[files],0.45)
            # # h.SetMaximum(10)
            s3.Add(h)
            l2.AddEntry(h,files,'lf')
            # # h.SetMaximum(10)
    # s3.SetMaximum(10)

s3.Draw("HIST nostack")
s3.GetYaxis().SetTitle("Events")
s3.GetXaxis().SetLabelSize(0.03)
s3.GetYaxis().SetLabelSize(0.02)
# s3.GetXaxis().SetTitle("#font[52]{#Delta_{R}} ")
# s3.GetXaxis().SetTitle("#font[52]{leading #tilde{#chi}_{1}^{0}} [GeV]")
# s3.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
# s3.GetXaxis().SetTitle("#font[52]{MT2}")
# s3.GetXaxis().SetTitle("#font[52]{P_{T} H/A} [GeV]")
# s3.GetXaxis().SetTitle("#font[52]{P_{T} subleading #tilde{#tau}/subleading #tau}")
s3.GetXaxis().SetTitle("#font[52]{P_{T} leading #tilde{#tau}/leading #tau}")
# s3.GetXaxis().SetTitle("#font[52]{P_{T} leading #tau}")
l2.Draw()
# can.SetLogy()
can.Print("B_l_stautau{}.png".format(Benchmark_list[0]))
assert False

# for files, var in n.items():
#     j =0
#     c2 = ROOT.TCanvas('c2','c2',800,600)
#     h = ROOT.TH1F('h1','h1',100,0,max(var['sleading_tau']))
#     leg1 = ROOT.TLegend(0.65,0.65,0.90,0.90)
#     leg1.SetBorderSize(0)
#     leg1.AddEntry(h,files,'lf')
#     for i in var['sleading_tau']:
#         h.Fill(i)
#     scale_hist(h)
#     h.Draw("HIST")
#     h.GetYaxis().SetTitle("Events")
#     h.GetXaxis().SetLabelSize(0.03)
#     h.GetYaxis().SetLabelSize(0.02)
#     h.GetXaxis().SetTitle("#font[52]{ subleading #tau} mA %s [GeV]"%files[:4])
#     ymax = h.GetMaximum()
#     p 	= array('d',[0.15,0.25,0.5,0.75])
#     q 	= array('d',[0,0,0,0])
#     h.GetQuantiles(4,q,p)
#     # logging.info("For mA, 15%,25%,50%,75%,{}, {}, {}, {} ,{}, {}".format(files[:4],q[0],q[1],q[2],q[3],ymax))

#     l = ROOT.TLine(q[0],0,q[0],ymax)
#     l.SetLineColor(4)
#     l.SetLineWidth(2)
#     l.Draw("")
#     l1 = ROOT.TLine(q[1],0,q[1],ymax)
#     l1.SetLineColor(1)
#     l1.SetLineWidth(2)
#     l1.Draw("")
#     l2 = ROOT.TLine(q[2],0,q[2],ymax)
#     l2.SetLineColor(2)
#     l2.SetLineWidth(2)
#     l2.Draw("")
#     l3 = ROOT.TLine(q[3],0,q[3],ymax)
#     l3.SetLineColor(3)
#     l3.SetLineWidth(2)
#     l3.Draw("")
#     leg1.Draw()
#     t = ROOT.TText()
#     t.DrawText(q[0],ymax/4,'15%% %0.5s'%q[0]).SetTextAngle(90)
#     t.DrawText(q[1],ymax/2,'25%% %0.5s'%q[1]).SetTextAngle(90)
#     t.DrawText(q[2],ymax/2,'50%% %0.5s'%q[2]).SetTextAngle(90)
#     t.DrawText(q[3],ymax/2,'75%% %0.5s'%q[3]).SetTextAngle(90)
#     # c2.Print("B_slead_div_{}.png".format(files[:4]))
met_a = []
dr_a = []
c35 = ROOT.TCanvas("c","c",800,600)
h35 = ROOT.THStack("h1","")
l35 = ROOT.TLegend(0.65,0.65,0.90,0.90)
l35.SetBorderSize(0)
l35.SetFillStyle(0)
l35.SetTextSize(0.03)
# assert False
for files, var in n.items():
    i=0

    # for var_name, var_list in var.items():
    #     minm = 0		#changing the minimum of the graph
    #     maxm = 600E3
    #     c2 = ROOT.TCanvas('c2','c2',800,600)
    #     hs1 = ROOT.THStack("hs1","")
    #     leg1 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    #     leg1.SetBorderSize(0)
    #     if "dR_ltau_neu1" in var_name or 'dR_stau_neu1' in var_name:
    #         maxm = 6
    #         minm = -2

    #     # ROOT.gStyle.SetPalette(67)
    #     leg1.SetFillStyle(0)
    #     leg1.SetTextSize(0.02)
    #     print var_name
    #     o = ROOT.TH1F(var_name,"",100,minm,maxm)
    #     fill_hist(o,var_list)
    #     leg1.AddEntry(o,files,'lf')
    #     o.SetLineColor(color[var_name])
    #     i=+1
    #     scale_hist(o)
    #     hs1.Add(o)

    #     hs1.Draw("HIST,nostack")
    #     hs1.GetYaxis().SetTitle("Events")
    #     hs1.GetXaxis().SetLabelSize(0.02)
    #     hs1.GetYaxis().SetLabelSize(0.02)
    #     hs1.GetXaxis().SetTitle("#font[52]{%s} [GeV]"% var_name)
    #     c2.SetLogy()
    #     # c2.Print("B_{}_{}.png".format(files,var_name))
    # c3 = ROOT.TCanvas("c3",'c3',800,600)
    # h2 = ROOT.TGraph(len(var.values()[0]),var["leading_tau"],var["leading_neu1"])
    # h2.Draw("A*")
    # h2.GetXaxis().SetTitle("#font[52]{leading tau} [GeV]")
    # h2.GetYaxis().SetTitle("#font[52]{nue1 from leading tau} [GeV]")
    # h2.GetXaxis().SetLabelSize(0.03)
    # h2.GetYaxis().SetLabelSize(0.03)
    # text = ROOT.TText(max(var["leading_tau"])*.9,max(var["leading_neu1"])*0.9,files)
    # text.Draw()
    # # c3.Print("B_lead_tau_nue1_{}.png".format(files))

    # c3 = ROOT.TCanvas("c3",'c3',800,600)
    # h3 = ROOT.TGraph(len(var.values()[0]),var["sleading_tau"],var["leading_neu1"])
    # h3.Draw("A*")
    # h3.GetXaxis().SetTitle("#font[52]{sub leading #tau} [GeV]")
    # h3.GetYaxis().SetTitle("#font[52]{#tilde{#chi}_{1}^{0} from leading #tau} [GeV]")
    # h3.GetXaxis().SetLabelSize(0.03)
    # h3.GetYaxis().SetLabelSize(0.03)
    # text = ROOT.TText(max(var["sleading_tau"])*.9,max(var["sleading_neu1"])*0.9,files)
    # text.Draw()
    # c3.Print("B_slead_tau_lead_nue1_{}.png".format(files[:4]))

    c4 = ROOT.TCanvas("c4",'c4',800,600)
    ROOT.gStyle.SetPalette(1)
    h4 = ROOT.TH2F("h4","",50,0,900E3,50,0,1000E3)
    for i,j in zip(var["leading_tau"],var["met_Tr"]):
        h4.Fill(i,j)
    h4.Draw("COL Z CJUST")
    h4.GetXaxis().SetTitle("#font[52]{leading #tau p_{T}(%s)}"%files[:4])
    h4.GetYaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
    h4.GetXaxis().SetLabelSize(0.03)
    h4.GetYaxis().SetLabelSize(0.03)
    h4.GetZaxis().SetLabelSize(0.03)
    c4.SetWindowSize(800,600)
    c4.Print("B_2d_lead_tua_met_{}.png".format(files[:4]))

    c5 = ROOT.TCanvas("c5",'c5',800,600)
    h5 = ROOT.TH2F("h5","",50,0,800E3,50,0,1400E3)
    for i,j in zip(var["sleading_tau"],var["sleading_neu1"]):
        h5.Fill(i,j)
    h5.Draw("COL Z CJUST")
    ttext = ROOT.TText(600E3,1200E3,files)

    ROOT.gPad.SetRightMargin(0.13)
    h5.GetXaxis().SetTitle("#font[52]{p_{T} subleading#tau} [GeV]")
    h5.GetYaxis().SetTitle("#font[52]{p_{T} #tilde{#chi}_{1}^{0} } [GeV]")
    h5.GetXaxis().SetLabelSize(0.03)
    h5.GetYaxis().SetLabelSize(0.03)
    h5.GetZaxis().SetLabelSize(0.03)
    ttext.Draw()
    # c5.SetWindowSize(800,600)
    # c5.Print("B_2d_sleadin_stua_snue_{}.png".format(files[:4]))
assert False
for files, var in n.items():
    s_tau = [i for i in var["sleading_tau"] if i < 35E3]
    s_met = []
    s_dR = []

    s_tau25 = [i for i in var["sleading_tau"] if i < 25E3]
    s_met25 = []
    s_dR25 = []
    print files
    for i in xrange(len(var["met_Tr"])):
        if var["sleading_tau"][i] < 35E3:
            s_met.append(var["met_Tr"][i])
    for i in xrange(len(var["met_Tr"])):
        if var["sleading_tau"][i] < 25E3:
            s_met25.append(var["met_Tr"][i])
    for i in xrange(len(var["dR_tau1_tau2"])):
        if var["sleading_tau"][i] < 35E3:
            s_dR.append(var["dR_tau1_tau2"][i])
    for i in xrange(len(var["dR_tau1_tau2"])):
        if var["sleading_tau"][i] < 25E3:
            s_dR25.append(var["dR_tau1_tau2"][i])

    # logging.info("For mA, len(tau),tau<35GeV,tau<25GeV,{}, {}, {}, {} ".format(files,len(var['sleading_tau']),len(s_tau),len(s_tau25)))

    # c7 = ROOT.TCanvas("c7","c7",800,600)

    c8 = ROOT.TCanvas("c8","c8",800,600)
    h6 = ROOT.TH1F("","",30,0,700E3)
    fill_hist(h6,s_met)
    scale_hist(h6,norm=1)
    h6.SetLineColor(2)
    h1 = ROOT.TH1F("","",30,0,700E3)
    fill_hist(h1,var["met_Tr"])
    scale_hist(h1,norm=1)
    # h1.SetLineColor(color[files])
    h1.SetLineColor(1)

    met_a.append(h1)
    h8 = ROOT.TH1F("","",30,0,700E3)
    fill_hist(h8,s_met25)
    scale_hist(h8,norm=1)
    h8.SetLineColor(3)
    hs1 = ROOT.THStack("h1","")
    leg4 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    leg4.SetBorderSize(0)
    leg4.SetFillStyle(0)
    leg4.SetTextSize(0.03)
    leg4.AddEntry(h6,"tau<35GeV %s"%len(s_met),'lf')
    leg4.AddEntry(h8,"tau<25GeV %s"%len(s_met25),'lf')
    leg4.AddEntry(h1,"%s  %s"%(files,len(var['met_Tr'])),'lf')
    hs1.Add(h6)
    hs1.Add(h8)
    hs1.Add(h1)
    hs1.Draw("HIST")
    hs1.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} mA %s [GeV]"%files[:4])
    hs1.GetYaxis().SetTitle("#font[52]{Events}")
    hs1.GetXaxis().SetLabelSize(0.03)
    hs1.GetYaxis().SetLabelSize(0.03)
    leg4.Draw()
    c8.Print("B_met_25_25_%s.png"%files[:4])

    c9 = ROOT.TCanvas("c9","c9",800,600)

    h7 = ROOT.TH1F("","",30,-1,6)
    fill_hist(h7,s_dR)
    h7.SetLineColor(color[files])
    scale_hist(h7,norm=1)
    h35.Add(h7)
    l35.AddEntry(h7,files,'lf')

    h9 = ROOT.TH1F("","",30,-1,6)
    fill_hist(h9,s_dR25)
    scale_hist(h9,norm=1)
    h9.SetLineColor(3)
    h2 = ROOT.TH1F("","%s  %s"%(files,len(var['dR_tau1_tau2'])),30,-1,6)
    fill_hist(h2,var["dR_tau1_tau2"])
    scale_hist(h2,norm=1)
    # h2.SetLineColor(color[files])
    h2.SetLineColor(1)
    dr_a.append(h2)
    hs2 = ROOT.THStack("h2","")
    leg5 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    leg5.SetBorderSize(0)
    leg5.SetFillStyle(0)
    leg5.SetTextSize(0.03)
    leg5.AddEntry(h7,"tau<35GeV %s"%len(s_dR),'lf')
    leg5.AddEntry(h9,"tau<25GeV %s"%len(s_dR25),'lf')
    leg5.AddEntry(h2,"%s  %s"%(files,len(var['dR_tau1_tau2'])),'lf')
    hs2.Add(h7)
    hs2.Add(h2)
    hs2.Add(h9)
    hs2.Draw("HIST ")
    hs2.GetXaxis().SetTitle("#font[52]{#Delta_{R}(leading#tau,subleading#tau)} mA %s"%files[:4])
    hs2.GetYaxis().SetTitle("#font[52]{Events}")
    hs2.GetXaxis().SetLabelSize(0.03)
    hs2.GetYaxis().SetLabelSize(0.03)
    leg5.Draw()
    c9.Print("B_dr_25_25_{}.png".format(files[:4]))

h35.Draw("HIST ")
h35.GetXaxis().SetTitle("#font[52]{#Delta_{R}(leading#tau,subleading#tau)} ")
h35.GetYaxis().SetTitle("#font[52]{Events}")
h35.GetXaxis().SetLabelSize(0.03)
h35.GetYaxis().SetLabelSize(0.03)
l35.Draw()
c9.Print("B_dr_35_35_all.png")
assert False
c = ROOT.TCanvas("c","c",800,600)
s1 =  ROOT.THStack("","")
le = ROOT.TLegend(0.65,0.65,0.90,0.90)
le.SetBorderSize(0)
le.SetFillStyle(0)
le.SetTextSize(0.03)
for i in xrange(len(met_a)):
    s1.Add(met_a[i])
    le.AddEntry(met_a[i],n.keys()[i],'lf')
s1.Draw("HIST ")
le.Draw()
s1.GetXaxis().SetTitle("#font[52]{E_{T}^{miss}} [GeV]")
s1.GetYaxis().SetTitle("#font[52]{Events}")
s1.GetXaxis().SetLabelSize(0.03)
s1.GetYaxis().SetLabelSize(0.03)
c.Print("B_MET_ALL.png")
le.Clear()
c1 = ROOT.TCanvas("c1","c1",800,600)
s1 =  ROOT.THStack("","")
le = ROOT.TLegend(0.65,0.65,0.90,0.90)
le.SetBorderSize(0)
le.SetFillStyle(0)
le.SetTextSize(0.03)
for i in xrange(len(dr_a)):
    s1.Add(dr_a[i])
    le.AddEntry(dr_a[i],n.keys()[i],'lf')
s1.Draw("HIST NOSTACK")
le.Draw()
s1.GetXaxis().SetTitle("#font[52]{#Delta_{R}(leading#tau,sub-leading#tau)}")
s1.GetYaxis().SetTitle("#font[52]{Events}")
s1.GetXaxis().SetLabelSize(0.03)
s1.GetYaxis().SetLabelSize(0.03)
c1.Print("B_DR_ALL.png")
le.Clear()
# le1.Clear()
st1 = ROOT.THStack("","")
st2 = ROOT.THStack("","")
le = ROOT.TLegend(0.65,0.65,0.90,0.90)
le.SetBorderSize(0)
le.SetFillStyle(0)
le.SetTextSize(0.02)
le1 = ROOT.TLegend(0.65,0.65,0.90,0.90)
le1.SetBorderSize(0)
le1.SetFillStyle(0)
le1.SetTextSize(0.02)
for files, var in n.items():
    for var_name, var_list in var.items():
        h2 = ROOT.TH1F("","",35,0,600E3)
        fill_hist(h2,var["sleading_tau"])
        scale_hist(h2)
        h2.SetLineColor(color[files])
        h3 = ROOT.TH1F("","",35,0,600E3)
        fill_hist(h3,var["leading_tau"])
        scale_hist(h3)
        h3.SetLineColor(color[files])
    le1.AddEntry(h3,files,'lf')
    st2.Add(h3)
    le.AddEntry(h2,files,'lf')
    st1.Add(h2)
c8 = ROOT.TCanvas("c7","c7",800,600)
st1.Draw("HIST NOSTACK")
le.Draw()
st1.GetXaxis().SetTitle("#font[52]{P_{T}(sub-leading#tau)}")
st1.GetYaxis().SetTitle("#font[52]{Events}")
st1.GetXaxis().SetLabelSize(0.03)
st1.GetYaxis().SetLabelSize(0.03)
c8.Print("B_STAU_ALL.png")
c8 = ROOT.TCanvas("c7","c7",800,600)
st2.Draw("HIST NOSTACK")
le1.Draw()
st2.GetXaxis().SetTitle("#font[52]{P_{T}(leading#tau)}")
st2.GetYaxis().SetTitle("#font[52]{Events}")
st2.GetXaxis().SetLabelSize(0.03)
st2.GetYaxis().SetLabelSize(0.03)
c8.Print("B_LTAU_ALL.png")
    # c8 = ROOT.TCanvas("c7","c7",800,600)
    # h8 = ROOT.TH1F("","",50,0,27E3)
    # leg4 = ROOT.TLegend(0.65,0.65,0.90,0.90)
    # leg4.SetBorderSize(0)
    # leg4.SetFillStyle(0)
    # leg4.SetTextSize(0.03)
    # fill_hist(h8,s_tau)
    # leg4.AddEntry(h8,files,'lf')
    # scale_hist(h8)
    # h8.GetXaxis().SetTitle("#font[52]{subleading #tau}")
    # h8.GetYaxis().SetTitle("#font[52]{Events}")
    # h8.GetXaxis().SetLabelSize(0.03)
    # h8.GetYaxis().SetLabelSize(0.03)
    # h8.Draw("HIST")
    # leg4.Draw()
    # # c8.Print("H_pt_slead_{}.png".format(files))
can = ROOT.TCanvas("can","can",800,600)
s3 =  ROOT.THStack("","")
l2 = ROOT.TLegend(0.65,0.65,0.90,0.90)
l2.SetBorderSize(0)
l2.SetFillStyle(0)
l2.SetTextSize(0.03)
for files, var in n.items():
    for var_name, var_list in var.items():
        print files, var_name



assert False
# def fill_hist(hist,array):
# 	[hist.Fill(_) for _ in array]
# def scale_hist(hist_bkg_list,norm=1):
# 	# if hist_bkg_list.Integral() != 0 :
# 	# 	hist_bkg_list.Scale((1/(hist_bkg_list.Integral())))
# 	scale = norm/hist_bkg_list.Integral()
# 	hist_bkg_list.Scale(scale)
# def H_mk(name,bin,minm,mixm):
#     return ROOT.TH1F(name, "",bin,minm,mixm)

# def open_pickle(f):
#     with open(f) as file:
#         return pickle.load(file)

# var = ['pt_tau1','pt_tau2','met_Truth','mt_tot','dR_tau1_tau2']
# color = ['kBlue', 'kRed', 'kGreen', 'kPink', 'kOrange','kPurple']
# # file_list_T = ['1000T.pkl','1500T.pkl','2000T.pkl','500T.pkl','800T.pkl']
# # file_list_B = ['1000B.pkl','1500B.pkl','2000B.pkl']

# file_list_T = ['500T.pkl','800T.pkl','1000T.pkl','1500T.pkl','2000T.pkl']
# file_list_B = ['1000B.pkl','1500B.pkl','2000B.pkl','2500B.pkl']

# open_file = [open_pickle(i) for i in file_list_T]

# # data_dict = dict(zip(var,))
# data_dict = {}

# for i in range(len(var)):
#     hist_list = []
#     for j in xrange(len(file_list_T)):
#         hist =  H_mk(var[i]+' '+file_list_T[j],30,min(open_file[j][i]),max(open_file[j][i]))
#         fill_hist(hist,open_file[j][i])
#         scale_hist(hist)
#         hist_list.append(hist)
#         # hist.Draw("HIST")
#         # print hist.GetMaximum()
#         # print hist.GetXaxis().GetXmax(), "###"
#         print file_list_T[j]
#         print var[i]
#     # hist_list = hist_list.sort(reverse=True, key=lambda x: x.GetMaximum())
#     print hist_list
#     data_dict[var[i]] = hist_list
# # print data_dict['pt_tau2']

# for v in var:
#     c = ROOT.TCanvas()
#     hs = ROOT.THStack(v,v)
#     # print data_dict[v].GetMaximum()
#     # mem = data_dict[v].sort(reverse=True, key=lambda x: x.GetMaximum())
#     for i in data_dict[v]:
#         hs.Add(i)
#     hs.Draw("HIST SAME C")
#     c.Print("{}_tfile.pdf".format(v))

# assert False



# histogram = [H_mk(var[i],30,min(hstau_800_truth_noTrig[i]),max(hstau_800_truth_noTrig[i])) for i in xrange(len(var))]
# for i in xrange(len(var)):
#     fill_hist(histogram[i],hstau_800_truth_noTrig[i])
#     scale_hist(histogram[i])
#     c1 = ROOT.TCanvas()
#     c1.cd()
#     print var[i]
#     histogram[i].Draw("HIST")
#     histogram[i].SetStats(0)
#     histogram[i].SetLineColor(ROOT.kRed)
#     histogram[i].SetTitle("{}_2000b".format(var[i]))
#     c1.Print("test_2000b_{}.pdf".format(var[i]))
