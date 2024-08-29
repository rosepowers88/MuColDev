from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TLorentzVector, TTree, TVector3
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TGraphErrors, TF1, TProfile, TProfile2D
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange, kViolet
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from math import *
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from array import array
import os
import logging
import itertools
import fnmatch
import csv
import pandas as pd
import mplhep as hep

parser= OptionParser()
parser.add_option("-m", "--inFile_0_50_10T", dest='inFile_0_50_10T',
                  default="histos_photon.root", help="Name of the ROOT file")
parser.add_option("-n", "--inFile_50_250_10T", dest='inFile_50_250_10T',
                  default="histos_photon.root", help="Name of the ROOT file")
parser.add_option("-o", "--inFile_250_1000_10T", dest='inFile_250_1000_10T',
                  default="histos_photon.root", help="Name of the ROOT file")

(options, args) = parser.parse_args()

#load files -- the slices
fFile_0_50 = TFile(options.inFile_0_50_10T, "READ")
fFile_50_250 = TFile(options.inFile_50_250_10T, "READ")
fFile_250_1000 = TFile(options.inFile_250_1000_10T, "READ")


#this is not really useful
######## ECAL dimensions #######
zlim = 230.7 #cm
R_inner_lim = 150.0 #cm
R_outer_lim = 185.7 #cm
X_0 = 8.897 #cm -- rad length for AL
###############################

####### ROOT gStyle ##########
gStyle.SetPadTopMargin(0.09)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetOptStat(0)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
################################

# Energy binning (EPJC style), angular binning
EBins = array('d', (0., 50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800., 850., 900., 950., 1000.))

#ThetaBins = array('d', (20.*TMath.Pi()/180, 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180., 90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., 160*TMath.Pi()/180))
ThetaBins = np.linspace(0.175,2.96,30)

# Declare response profiles and reso plots
h_response_profile_E = TProfile('E_profile', 'E_profile', len(EBins)-1, EBins, 0,1000)
h_response_profile_th = TProfile('Th_profile', 'Th_profile', len(ThetaBins)-1, ThetaBins, 0, 1000)
h_2d_response = TProfile2D('T-E_profile', 'T-E_profile', len(ThetaBins)-1, ThetaBins, len(EBins)-1, EBins)
h_2d_resp_corr = TProfile2D('TEcorr', 'TEcorr', len(ThetaBins)-1, ThetaBins, len(EBins)-1, EBins)
h_reso_E = TH1D('reso_E', 'reso_E', len(EBins)-1, EBins)
h_reso_theta = TH1D('reso_theta', 'reso_theta', len(ThetaBins)-1, ThetaBins)

h_E_v_theta = TH2D ('E_v_theta', 'E_v_theta', len(ThetaBins)-1, ThetaBins, len(EBins)-1, EBins)

#declare arrays
e_arr = array('d')
th_arr = array('d')
th_err_arr=array('d')
sigma_arr = array('d')
sigma_arr_th = array('d')
e_err_arr = array('d')
sigma_err_arr = array('d')
sigma_err_arr_th = array('d')
mu_arr = array('d')
mu_err_arr = array('d')
mu_arr_th = array('d')
mu_arr_th_err = array('d')

#none of the commented chunk below is used anymore
'''
#Define piecewise function for radlengths as fn of theta
def NX_0(theta):
    if 0.577 < theta and theta < 0.678:
        return (zlim/cos(theta)-R_inner_lim/sin(theta))/X_0
    elif 0.678 < theta and theta < TMath.Pi()-0.678:
        return ((R_outer_lim-R_inner_lim)/sin(theta))/X_0
    elif TMath.Pi()-0.678 < theta and theta < TMath.Pi()-0.577:
        return (zlim/cos(TMath.Pi()-theta)-R_inner_lim/sin(TMath.Pi()-theta))/X_0
    else:
        return 0

def corr_ratio(x,A):
    t = NX_0(x)
    return (2**t)*A

'''

#get files
files = [fFile_0_50, fFile_50_250, fFile_250_1000]


#####################################
# BINNED CALIBRATION                #
# double for loop, take the         #
# Etrue/Erec ratio for each theta/E #
# slice and apply to the 2D bins.   #
#####################################

corr_matrix = [] #for each theta bin, a list of avgs for each E bin
#uncomment the below to regenerate the calibration matrix
'''
for tbin in range(0, len(ThetaBins)-1):
    rmax = 0
    ThMin = ThetaBins[tbin]
    ThMax = ThetaBins[tbin+1]
    Elist = []
    for Ebin in range(0, len(EBins)-1):
        EMin = EBins[Ebin]
        EMax = EBins[Ebin+1]
        ratios = []
        for file in files:
            #tree = file.Get("photon_tree")
            tree = file.Get("jet_tree")
            for entry in tree:
                E_truth = entry.photon_E
                theta = entry.photon_theta_angle
                E_reco = entry.jet_energy
                theta_reco = entry.jet_theta
                #get rid of weird theta values
                if theta_reco < 0:
                    continue
                # get rid of negative energy values
                if E_reco < 0:
                    continue
                #get rid of lower tail, file-specific
                #if file == fFile_50_250 and E_reco < 35:
                #    continue
                #if file == fFile_250_1000 and E_reco < 100:
                #    continue
                #finally, dump the ratio into a list
                if EMin <= E_truth and E_truth < EMax:
                    if ThMin <= theta and theta < ThMax:
                        ratios.append(E_truth/E_reco)
        if len(ratios) > 0:
            avg_ratio = np.average(ratios)
        else:
            print("No entries in bin: ", Ebin, " Emin: ", EMin)
            avg_ratio = 0
        Elist.append(avg_ratio)
    print("Theta slice ", tbin, " done")
    print("Starting theta slice ",tbin+1)
    corr_matrix.append(Elist)


#write calibration map to a csv file
with open('calibMap_jets', 'w', newline = '') as csvfile:
    mapwriter = csv.writer(csvfile)
    mapwriter.writerows(corr_matrix)

#after initial calibrtion, read in the csv file and just pull values from there
# have block below uncommented if reading in file: corr_matrix will be defined differently
'''
with open('calibMap_jets.csv', 'r') as csvToRead:
    calibmap=csv.reader(csvToRead)
    corr_matrix = list(calibmap)


#plot the CalibMaps in both cases
#convert the matrix back into floats
corr_map = []
for thetaSlice in corr_matrix:
    corrfloats = []
    for corr in thetaSlice:
        corrfloats.append(float(corr))
    corr_map.append(corrfloats)



#plot with pcolormesh
plt.pcolormesh(EBins, ThetaBins, corr_map, vmin=np.min(corr_map), vmax=np.max(corr_map))
plt.colorbar()
plt.xlabel("True E [GeV]")
plt.ylabel("True Theta [rad]")
#use mplhep for labels
hep.cms.label(exp = "Muon Collider", data = False, rlabel = "MAIA Detector Concept", loc=0, italic=(1,0,0))
plt.savefig('calibmap_jets.pdf')
plt.close()


################### RESOLUTION STUDY #################


cx = TCanvas("", "", 800, 600)
gStyle.SetOptStat(1)

for Ebin in range(0, len(EBins)-1):
    EMin = EBins[Ebin]
    EMax = EBins[Ebin+1]
    proj_name = "E reso, "+str(EBins[Ebin])+"<E<"+str(EBins[Ebin+1])
    file_name = "Ereso"+str(Ebin)
    if  EMin < 750.:
        h_my_proj_2 = TH1D(proj_name, proj_name, 150, -0.1, 0.1)
    else:
        h_my_proj_2 = TH1D(proj_name, proj_name, 150, -.05, .05)


    for file in files:
        tree = file.Get("jet_tree")
        for entry in tree:
            E_truth = entry.photon_E
            theta = entry.photon_theta_angle
            E_reco = entry.jet_energy
            E_corr = entry.jet_energy
            theta_reco = entry.jet_theta
            #get rid of weird theta values
            if theta_reco < 0:
                continue
            # get rid of negative energy values
            if E_reco < 0:
                continue
            for tbin in range(0, len(ThetaBins)-1):
                ThMin = ThetaBins[tbin]
                ThMax = ThetaBins[tbin+1]
                #get corr factor from profile
                corr = float(corr_matrix[tbin][Ebin])
                #if entry is in the 2d bin, correct energy
                if EMin <= E_truth and E_truth < EMax:
                    if ThMin <= theta and theta < ThMax:
                        if corr > 0:
                            E_corr = E_reco * corr
                        #isolate the peak of the distribution so we can fit a Gaussian
                        neg_lim = -0.1
                        pos_lim = 0.1
                        if EMin < 50:
                            neg_lim = -0.08
                            pos_lim = 0.08
                        elif 50==EMin:
                            neg_lim = -0.05
                            pos_lim = 0.05
                        elif 100 <= EMin:
                            neg_lim = -0.02
                            pos_lim = 0.02
                        if (E_corr - E_truth)/E_truth > neg_lim and (E_corr - E_truth)/E_truth < pos_lim:
                            h_my_proj_2.Fill((E_corr - E_truth)/E_truth)
                            h_response_profile_E.Fill(E_truth, (E_truth)/E_corr)
                            h_response_profile_th.Fill(theta, (E_truth)/E_corr)
           

    #now start fitting the gaussians
    lim = 0.0
    if EMin<3000.:
        if EMin<250:
            gaussFit1 = TF1("gaussfit", "gaus", -0.5, 0.5)
        elif EMin>700:
            gaussFit1 = TF1("gaussfit", "gaus", -.02, 0.02)   
        else:
            gaussFit1 = TF1("gaussfit", "gaus", -0.1, 0.1)
        gaussFit1.SetLineColor(kRed)
        gaussFit1.SetParameter(1, 0.)
        gaussFit1.SetParameter(2, h_my_proj_2.GetRMS())
        h_my_proj_2.Fit(gaussFit1, "E")
        gStyle.SetOptFit(0o111);
        h_my_proj_2.Draw("HIST")
        gaussFit1.Draw("LSAME")
        cx.Update()
        cx.SaveAs("slices_corr_full"+file_name+".root")
        sigma = gaussFit1.GetParameter(2)
        sigma_err = gaussFit1.GetParError(2)
        bincenter = (EMax-EMin)/2
        e_arr.append(h_reso_E.GetBinCenter(Ebin+1))
        e_err_arr.append(bincenter)
        sigma_arr.append(sigma)
        mu_arr.append(gaussFit1.GetParameter(1))
        mu_err_arr.append(gaussFit1.GetParError(1))
        print("SIGMA=",sigma)
        #print(h1_reso_E.GetBinCenter(bin+1))
        sigma_err_arr.append(sigma_err)
        

#uncomment for theta resolution studies
'''
for bin in range(0,len(ThetaBins)-1):
    minTh = ThetaBins[bin]
    maxTh = ThetaBins[bin+1]

    proj_name = "E reso, "+str(round(ThetaBins[bin],3))+"<Theta<"+str(round(ThetaBins[bin+1],3))
    print(proj_name)
    #adjust bin range so we don't underbin
    binrange = 1.5
    #if ThetaBins[bin]< 10.0:
        #binrange = 1.5
    #elif ThetaBins[bin]<70.0:
        #binrange = 1.
   # else:
        #binrange = 0.5
    h_my_proj_2 = TH1D(proj_name, proj_name, 150, -binrange, binrange)

    for file in files:
        tree = file.Get("photon_tree")
        for entry in tree:
            E_truth = entry.E_truth
            E_reco = entry.E
            theta = entry.theta
            if correctEnergy: 
                E_reco = E_corr(E_reco, theta)
            #h_response_profile_E.Fill(E_truth, (E_truth-E_reco)/E_truth)
            #h_response_profile_th.Fill(theta, (E_truth-E_reco)/E_truth)
            if E_reco > 0:
                #if(entry.theta > 0.65 and  entry.theta < 1.01) or (entry.theta > 2.13 and entry.theta <2.53): #if entry.theta < 0.612:
                #if theta >1.01 and theta < 2.13:
                if (E_truth > minE) and (E_truth < maxE):
                    if (E_reco - E_truth)/E_truth > -0.9:
                        h_my_proj_2.Fill((E_reco - E_truth)/E_truth)


    #lim = 0.0
    #if minE<5000:
     #   if minE<50:
    lim = 0.5
      #  else:
       #     lim = 0.25
    gaussFit = TF1("gaussfit", "gaus", -1.0*lim, lim)
    gaussFit.SetLineColor(kBlue)
    gaussFit.SetParameter(1, 0.)
    gaussFit.SetParameter(2, h_my_proj_2.GetRMS())
    h_my_proj_2.Fit(gaussFit, "E")
    gStyle.SetOptFit(0o111);
    h_my_proj_2.Draw("HIST")
    gaussFit.Draw("LSAME")
    cx.SaveAs("slices_ph/"+proj_name+".pdf")
    sigma = gaussFit.GetParameter(2)
    sigma_err = gaussFit.GetParError(2)
    if bin > 0:
        bincenter = (maxTh-minTh)/2
        th_arr.append(h_reso_theta.GetBinCenter(bin+1))
        th_err_arr.append(bincenter)
        sigma_arr_th.append(sigma)
        mu_arr_th.append(gaussFit.GetParameter(1))
        mu_arr_th_err.append(gaussFit.GetParError(1))
        print("SIGMA=",sigma)
        #print(h1_reso_E.GetBinCenter(bin+1))
        sigma_err_arr_th.append(sigma_err)
    else:
        gaussFit = TF1("gaussfit", "gaus", -0.15, 0.15)
        gaussFit.SetLineColor(kBlue)
        gaussFit.SetParameter(1, 0.)
        gaussFit.SetParameter(2, h_my_proj_2.GetRMS())
        h_my_proj_2.Fit(gaussFit, "E")
        h_my_proj_2.Draw("HIST")
        gaussFit.Draw("LSAME")
        cx.SaveAs("slices_ph/"+proj_name+".pdf")
        sigma = gaussFit.GetParameter(2)
        sigma_err = gaussFit.GetParError(2)
        if bin > 0:
            th_arr.append(h_reso_theta.GetBinCenter(bin+1))
            th_err_arr.append(maxTh-minTh)
            sigma_arr_th.append(sigma)
            sigma_err_arr_th.append(sigma_err)
            mu_arr_th.append(gaussFit.GetParameter(1))
            mu_arr_th_err.append(gaussFit.GetParError(1))
            
'''     

#declare the graphs
gr_reso_E = TGraphErrors(len(e_arr), e_arr, sigma_arr, e_err_arr, sigma_err_arr)
#gr_reso_Th = TGraphErrors(len(th_arr), th_arr, sigma_arr_th,th_err_arr, sigma_err_arr_th)
#gr_response_E = TGraphErrors(len(e_arr), e_arr, mu_arr, e_err_arr, mu_err_arr)
#gr_response_Th = TGraphErrors(len(th_arr), th_arr, mu_arr_th, th_err_arr, mu_arr_th_err)



cE1=TCanvas("", "", 800, 600)
gr_reso_E.SetTitle(" ")
gr_reso_E.GetYaxis().SetTitle("Jet #sigma_{E}/E")
gr_reso_E.GetYaxis().SetTitleOffset(1.4)
gr_reso_E.GetXaxis().SetTitleOffset(1.2)
gr_reso_E.GetXaxis().SetTitle("True jet energy [GeV]")
gr_reso_E.SetLineColor(kBlue)
gr_reso_E.SetLineWidth(2)
gr_reso_E.Draw("AP")
cE1.Update()


#ROOT fit

resoFit = TF1(    "resofit", "sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 0., 1000., 3)
resoFit.SetLineColor(kBlue)
resoFit.SetParName(0, "Stochastic")
resoFit.SetParName(1, "Noise")
resoFit.SetParName(2, "Constant")
resoFit.SetParameter(0, 0.5)
resoFit.SetParameter(1, 0.)
resoFit.SetParameter(2, 0.)
gStyle.SetOptFit(0)
gr_reso_E.Fit(resoFit)
#resoFit10.SetLineColor(kBlue)
#resoFit.Write("ResoFit")
gr_reso_E.Draw("AP")
resoFit.Draw("LSAME")
#gr_reso_E.SaveAs("Ereso_gr_EE.root")
#resoFit.SaveAs("Ereso_fit_EE.root")
cE1.Update()
cE1.SaveAs("Ereso_matrixcorr_jets.root")


#MAIA format python plots (!!)
###############RESO PLOT IN PAPER GENERATED WITH THE BELOW####################

#redefine the fit function
def py_resofit(stoch, noise, const, E):
    return np.sqrt((stoch**2)/E+(noise/E)**2+const**2)

Epoints = np.linspace(25,1000,1000)
resopoints=[]
#grab the fit values from the ROOT fit
for energy in Epoints:
    resopoints.append(py_resofit(resoFit.GetParameter(0), resoFit.GetParameter(1), resoFit.GetParameter(2), energy))

plt.errorbar(e_arr, sigma_arr, yerr = sigma_err_arr, xerr = e_err_arr, fmt=".", color='blue')
plt.plot(Epoints, resopoints,'--', c='blue', label="Fit")
plt.xlabel("True Photon Energy [GeV]")
plt.ylabel("Photon Energy Resolution $\sigma/E$")
hep.cms.label(exp = "Muon Collider", data = False, rlabel = "MAIA Detector Design", loc=0, italic=(1,0,0))
plt.show()
plt.savefig("Ereso_EE_maia.pdf")
plt.close()



cE2=TCanvas("", "", 800, 600)
#gr_response_E.SetTitle(" ")
#gr_response_E.GetYaxis().SetTitle("Photon #mu_{E}/E")
#gr_response_E.GetYaxis().SetTitleOffset(1.4)
#gr_response_E.GetXaxis().SetTitleOffset(1.2)
#gr_response_E.GetXaxis().SetTitle("True jet energy [GeV]")
#gr_response_E.SetLineColor(kOrange)
#gr_response_E.SetLineWidth(2)
#gr_response_E.Draw("AP")
h_response_profile_E.SetTitle("")
h_response_profile_E.GetYaxis().SetTitle("Jet E_{true}/E_{reco}")
h_response_profile_E.GetYaxis().SetTitleOffset(1.4)
h_response_profile_E.GetXaxis().SetTitleOffset(1.2)
h_response_profile_E.GetXaxis().SetTitle("True photon energy [GeV]")
h_response_profile_E.SetLineColor(kRed)
h_response_profile_E.SetLineWidth(2)
h_response_profile_E.Draw("HIST PE")
cE2.Update()
cE2.SaveAs("Eresp_corr_jets.root")



'''
cTh1=TCanvas("", "", 800, 600)
gr_reso_Th.SetTitle(" ")
gr_reso_Th.GetYaxis().SetTitle("Photon #sigma_{E}/E")
gr_reso_Th.GetYaxis().SetTitleOffset(1.4)
gr_reso_Th.GetXaxis().SetTitleOffset(1.2)
gr_reso_Th.GetXaxis().SetTitle("Photon #theta [rad]")
gr_reso_Th.SetLineColor(kBlue)
gr_reso_Th.SetLineWidth(2)
gr_reso_Th.Draw("AP")
cTh1.Update()
cTh1.SaveAs("ThReso_corr.root")
'''

cTh2=TCanvas("", "", 800, 600)
#gr_response_Th.SetTitle(" ")
#gr_response_Th.GetYaxis().SetTitle("Photon #mu_{E}/E")
#gr_response_Th.GetYaxis().SetTitleOffset(1.4)
#gr_response_Th.GetXaxis().SetTitleOffset(1.2)
#gr_response_Th.GetXaxis().SetTitle("Photon #theta [rad]")
#gr_response_Th.SetLineColor(kOrange)
#gr_response_Th.SetLineWidth(2)
#gr_response_Th.Draw("AP")
h_response_profile_th.SetTitle("")
h_response_profile_th.GetYaxis().SetTitle("Jet E_{true}/E_{reco}")
h_response_profile_th.GetYaxis().SetTitleOffset(1.4)
h_response_profile_th.GetXaxis().SetTitleOffset(1.2)
h_response_profile_th.GetXaxis().SetTitle("True jet theta [rad]")
h_response_profile_th.SetLineColor(kRed)
h_response_profile_th.SetLineWidth(2)
h_response_profile_th.Draw("HIST PE")
cTh2.Update()
cTh2.SaveAs("ThResp_corr_jets.root")

