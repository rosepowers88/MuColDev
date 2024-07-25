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
ThetaBins = np.linspace(0.3,3.0,30)

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

#define corrected energy as a function of RECO energy and theta
def E_corr(E, theta, A):
    return E*corr_ratio(theta, A)

#plot what the ratio ~should~ look like
ths = np.linspace(0,TMath.Pi(),5000)
ratios=[]
ratios.clear()
for th in ths:
    ratios.append(corr_ratio(th,0.15))

#ratio_plot=plt.scatter(ths, ratios, s=4, c='k')
#plt.ylim(0,50)
#plt.xlabel(r'$\theta$')
#plt.ylabel("$E_{Truth}/E_{Reco}$")
#plt.savefig('ratio_plot.pdf')
#plt.show()
#plt.close(ratio_plot)
#corr_Es = array('d', (25, 75, 125, 175, 225, 275, 325, 375, 425, 475, 525, 575, 625, 675, 725, 775, 825, 875, 925))
#corr_factors = array('d',(1.43e-1, 1.42e-1, 1.39e-1, 1.23e-1, 1.02e-1, 9.62e-2, 9.31e-2, 8.56e-2, 8.21e-2, 1.50e-2, 1.50e-2, 1.49e-2, 1.48e-2, 1.94e-2, 1.99e-2, 2.73e-2, 2.61e-2, 4.09e-2, 4.55e-2))
#corr_plot = plt.scatter(corr_Es, corr_factors)
#plt.xlabel("E [GeV]")
#plt.ylabel("Correction factor")
#plt.savefig('corr_plot.pdf')
'''

#get files
files = [fFile_0_50, fFile_50_250, fFile_250_1000]


#####################################
# BINNED CALIBRATION ATTEMPT        #
# double for loop, simply take the  #
# Etrue/Erec ratio for each theta/E #
# slice and apply to the 2D bins.   #
#####################################

corr_matrix = [] #for each theta bin, a list of avgs for each E bin
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
            tree = file.Get("photon_tree")
            for entry in tree:
                E_truth = entry.E_truth
                theta = entry.theta
                E_reco = entry.E
                if EMin <= E_truth and E_truth < EMax:
                    if ThMin <= theta and theta < ThMax:
                        ratios.append(E_truth/E_reco)
        if len(ratios) > 0:
            if max(ratios) > rmax:
                rmax = max(ratios)
            #print(rmax)
            #plt.hist(ratios,bins=50, log=True)
            #plt.xlim(0,500)
            #plt.ylim(1,7000)
            #plt.xlabel('E_truth/E_reco')
            #plt.show()
            #plt.savefig('ratios_'+str(tbin)+'.pdf')
            #plt.close()
            #UNCOMMENT BELOW to remove outliers (exclusive correction matrix)
            #for ratio in ratios:
                #if ratio/np.average(ratios) > 10.0:
                    #ratios.remove(ratio)
            avg_ratio = np.average(ratios)
        else:
            print("No entries in bin: ", Ebin, " Emin: ", EMin)
            avg_ratio = 0
        Elist.append(avg_ratio)
    #plt.savefig('ratios_'+str(tbin)+'.pdf')
    #plt.close()
    print("Theta slice ", tbin, " done")
    print("Starting theta slice ",tbin+1)
    #print(Elist)
    corr_matrix.append(Elist)


#write calibration map to a csv file
with open('calibMap_unfiltered.csv', 'w', newline = '') as csvfile:
    mapwriter = csv.writer(csvfile)
    mapwriter.writerows(corr_matrix)

#after initial calibrtion, read in the csv file and just pull values from there
#uncomment below if reading in file: corr_matrix will be defined differently
'''
with open('calibMap_unfiltered.csv', 'r') as csvToRead:
    calibmap=csv.reader(csvToRead)
    corr_matrix = list(calibmap)




cx = TCanvas("", "", 800, 600)
gStyle.SetOptStat(1)

for Ebin in range(0, len(EBins)-1):
    EMin = EBins[Ebin]
    EMax = EBins[Ebin+1]
    proj_name = "E reso, "+str(EBins[Ebin])+"<E<"+str(EBins[Ebin+1])
    file_name = "Ereso"+str(Ebin)
    
    h_my_proj_2 = TH1D(proj_name, proj_name, 150, -1.5, 1.5)
    
    for file in files:
        tree = file.Get("photon_tree")
        for entry in tree:
            E_truth = entry.E_truth
            theta = entry.theta
            E_reco = entry.E
            E_corr = entry.E
            #restrict to barrel region
            #if theta < 1.01 or theta > 2.13:
                #continue
            for tbin in range(0, len(ThetaBins)-1):
                ThMin = ThetaBins[tbin]
                ThMax = ThetaBins[tbin+1]
                #print(ThMin, "< theta <", ThMax, ", ", EMin, "< E <", EMax)
                #get corr factor from profile
                corr = float(corr_matrix[tbin][Ebin])
                #print(corr)
                #if entry is in the 2d bin, correct energy
                if EMin <= E_truth and E_truth < EMax:
                    #print("in E Bin")
                    if ThMin <= theta and theta < ThMax:
                        #print("in theta bin")
                        if corr > 0:
                            E_corr = entry.E * corr
                            #print("Corrected")
                            #h_2d_resp_corr.Fill(theta, E_truth, E_truth/E_reco)
                        if (E_corr - E_truth)/E_truth > -0.2 and (E_corr - E_truth)/E_truth < 0.2:
                                #print("Passed last check")
                                #print((E_reco-E_truth)/E_truth)
                            h_my_proj_2.Fill((E_corr - E_truth)/E_truth)
                            h_response_profile_E.Fill(E_truth, (E_truth)/E_corr)
                            h_response_profile_th.Fill(theta, (E_truth)/E_corr)
           


    lim = 0.0
    if EMin<3000.:
        if EMin<100:
            gaussFit = TF1("gaussfit", "gaus", -1.5, 1.5)
        elif EMin<500 and EMin>550:
            gaussFit = TF1("gaussfit", "gaus", -.5, 0.5)
        else:
            gaussFit = TF1("gaussfit", "gaus", -0.2, 0.2)
        gaussFit.SetLineColor(kRed)
        gaussFit.SetParameter(1, 0.)
        gaussFit.SetParameter(2, h_my_proj_2.GetRMS())
        h_my_proj_2.Fit(gaussFit, "E")
        gStyle.SetOptFit(0o111);
        h_my_proj_2.Draw("HIST")
        gaussFit.Draw("LSAME")
        cx.Update()
        cx.SaveAs("slices_corr_full"+file_name+".root")
        sigma = gaussFit.GetParameter(2)
        sigma_err = gaussFit.GetParError(2)
        if Ebin > 0:
            bincenter = (EMax-EMin)/2
            e_arr.append(h_reso_E.GetBinCenter(Ebin+1))
            e_err_arr.append(bincenter)
            sigma_arr.append(sigma)
            mu_arr.append(gaussFit.GetParameter(1))
            mu_err_arr.append(gaussFit.GetParError(1))
            print("SIGMA=",sigma)
            #print(h1_reso_E.GetBinCenter(bin+1))
            sigma_err_arr.append(sigma_err)
        else:
            gaussFit = TF1("gaussfit", "gaus", -0.15, 0.15)
            gaussFit.SetLineColor(kRed)
            gaussFit.SetParameter(1, 0.)
            gaussFit.SetParameter(2, h_my_proj_2.GetRMS())
            h_my_proj_2.Fit(gaussFit, "E")
            h_my_proj_2.Draw("HIST")
            gaussFit.Draw("LSAME")
            cx.Update()
            cx.SaveAs("slices_corr_full"+file_name+".root")
            sigma = gaussFit.GetParameter(2)
            sigma_err = gaussFit.GetParError(2)
            if Ebin > 0:
                e_arr.append(h1_reso_E.GetBinCenter(Ebin+1))
                e_err_arr.append(EMax-EMin)
                sigma_arr.append(sigma)
                sigma_err_arr.append(sigma_err)
    
        

'''
#bool to turn energy calib on or off
correctEnergy = True

cx = TCanvas("", "", 800, 600)
gStyle.SetOptStat(1)
#start with energy
#correction_factors=[]
#const_factors=[]
for bin in range(0,len(EBins)-1):
    minE = EBins[bin]
    maxE = EBins[bin+1]
    
    proj_name = "E reso, "+str(EBins[bin])+"<E<"+str(EBins[bin+1])
    print(proj_name)
    #adjust bin range so we don't underbin
    binrange = 0
    if EBins[bin]< 10.0:
        binrange = 1.5
    elif EBins[bin]<70.0:
        binrange = 1.5
    else:
        binrange = 1.5
    h_my_proj = TH1D(proj_name, proj_name, 150, -binrange+1.5, binrange+1.5)
    h_resp_prof_binned=TProfile('Th_profile_b'+str(EBins[bin]), 'E_tru/E_rec, '+str(EBins[bin])+'<E<'+str(EBins[bin+1]), len(ThetaBins)-1, ThetaBins, 0, 1000)
    
    for file in files:
        tree = file.Get("photon_tree")
        for entry in tree:
            E_truth = entry.E_truth
            E_reco = entry.E
            theta = entry.theta
            #if correctEnergy: 
                #E_reco = E_reco*corr_ratio(theta, 0.14053)
            if theta >=0.8 and theta <= 2.4:
                if EBins[bin] < E_truth and E_truth > EBins[bin]+50.:
                    #h_response_profile_E.Fill(E_truth, E_truth/E_reco)
                    h_resp_prof_binned.Fill(theta, E_truth/E_reco)

    calib_fit=TF1("calibfit", "[0]*(2^(4.01/sin(x)))", 0.678, 2.46, 1)
    calib_lin=TF1("calibconstfit", "[0]+[1]*x", 0.7, 2.2, 1)
    calib_lin.SetLineColor(kRed)
    calib_lin.SetParName(0, "Constant")
    calib_lin.SetParName(1, "lin")
    #calib_lin.SetParameter(0,3.)
    calib_lin.FixParameter(1,0.)
    calib_fit.SetLineColor(kBlue)
    calib_fit.SetParName(0, "Constant")
    calib_fit.SetParameter(0,0.15)
    gStyle.SetOptFit(0)
    if EBins[bin] < 0:
        h_resp_prof_binned.Fit(calib_fit)
        print("Fit: ", calib_fit.GetParameter(0))
        correction_factor=calib_fit.GetParameter(0)
        h_resp_prof_binned.Draw("HIST PE")
        calib_fit.Draw("LSAME")
    else:
        bincounts=[]
        for tbin in range(0, len(ThetaBins)-1):
            if ThetaBins[tbin] >= 0.8 and ThetaBins[tbin] <= 2.4:
                content=h_resp_prof_binned.GetBinContent(tbin)
                bincounts.append(content)
        const_factor=np.average(bincounts)
        calib_lin.FixParameter(0,const_factor)
        h_resp_prof_binned.Fit(calib_lin)
        #print("Const: ", calib_lin.GetParameter(0))
        #print("Lin: ", calib_lin.GetParameter(1))
        h_resp_prof_binned.Draw("HIST PE")
        calib_lin.Draw("LSAME")
    h_resp_prof_binned.SaveAs("ThResp_"+str(EBins[bin])+".root")
    
    #if EBins[bin]<500:
        #h_my_proj_2 = TH1D(proj_name, proj_name, 150, 0,2.)
    #else:
        #h_my_proj_2 = TH1D(proj_name, proj_name, 150, -0.8, 0)


    h_my_proj_2 = TH1D(proj_name, proj_name, 150, -1., 0.2)
    for file in files:
        tree = file.Get("photon_tree")
        for entry in tree:
            E_truth = entry.E_truth
            E_reco = entry.E
            theta = entry.theta
            if E_reco > 0:
                #if(entry.theta > 0.65 and  entry.theta < 1.01) or (entry.theta > 2.13 and entry.theta <2.53): #if entry.theta < 0.612:
                if theta >0.678 and theta < 2.46:
                    if (E_truth > minE) and (E_truth < maxE):
                        if correctEnergy:
                            #if EBins[bin] >400:
                            E_reco = E_reco * const_factor
                                #print(const_factor, ", ", E_reco)
                            #else:
                                #E_reco = E_corr(E_reco, theta, correction_factor)
                                #for tbin in range(0, len(ThetaBins)-1):
                                    #if theta > ThetaBins[tbin] and theta <= ThetaBins[tbin+1]:
                                        #E_reco = E_reco / bincounts[tbin]
                        if (E_reco - E_truth)/E_truth > -0.9:
                            #print("Passed last check")
                            #print((E_reco-E_truth)/E_truth)
                            h_my_proj_2.Fill((E_reco - E_truth)/E_truth)
                            h_response_profile_E.Fill(E_truth, (E_truth-E_reco)/E_truth)
                            h_response_profile_th.Fill(theta, (E_truth-E_reco)/E_truth)
           


    lim = 0.0
    if minE<5000:
        #if minE<500:
            #gaussFit = TF1("gaussfit", "gaus", -0.5, 1.5)
        #else:
            #gaussFit = TF1("gaussfit", "gaus", -1., 0)
        gaussFit = TF1("gaussfit", "gaus", -1., -0.2)
        gaussFit.SetLineColor(kBlue)
        gaussFit.SetParameter(1, 0.)
        gaussFit.SetParameter(2, h_my_proj_2.GetRMS())
        h_my_proj_2.Fit(gaussFit, "E")
        gStyle.SetOptFit(0o111);
        h_my_proj_2.Draw("HIST")
        #gaussFit.Draw("LSAME")
        cx.SaveAs("slices_ph/"+proj_name+".pdf")
        sigma = gaussFit.GetParameter(2)
        sigma_err = gaussFit.GetParError(2)
        if bin > 0:
            bincenter = (minE-maxE)/2
            e_arr.append(h_reso_E.GetBinCenter(bin+1))
            e_err_arr.append(bincenter)
            sigma_arr.append(sigma)
            mu_arr.append(gaussFit.GetParameter(1))
            mu_err_arr.append(gaussFit.GetParError(1))
            print("SIGMA=",sigma)
            #print(h1_reso_E.GetBinCenter(bin+1))
            sigma_err_arr.append(sigma_err)
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
                e_arr.append(h1_reso_E.GetBinCenter(bin+1))
                e_err_arr.append(maxE-minE)
                sigma_arr.append(sigma)
                sigma_err_arr.append(sigma_err)
'''
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
gr_reso_E = TGraphErrors(len(e_arr), e_arr, sigma_arr, e_err_arr, sigma_err_arr)
#gr_reso_Th = TGraphErrors(len(th_arr), th_arr, sigma_arr_th,th_err_arr, sigma_err_arr_th)
#gr_response_E = TGraphErrors(len(e_arr), e_arr, mu_arr, e_err_arr, mu_err_arr)
#gr_response_Th = TGraphErrors(len(th_arr), th_arr, mu_arr_th, th_err_arr, mu_arr_th_err)




#add the fits in later





cE1=TCanvas("", "", 800, 600)
gr_reso_E.SetTitle(" ")
gr_reso_E.GetYaxis().SetTitle("Photon #sigma_{E}/E")
gr_reso_E.GetYaxis().SetTitleOffset(1.4)
gr_reso_E.GetXaxis().SetTitleOffset(1.2)
gr_reso_E.GetXaxis().SetTitle("True photon energy [GeV]")
gr_reso_E.SetLineColor(kBlue)
gr_reso_E.SetLineWidth(2)
gr_reso_E.Draw("AP")
cE1.Update()
#cE1.SaveAs("Ereso.root")



resoFit = TF1(
    "resofit", "sqrt([0]*[0]/x+[1]*[1]/(x*x)+[2]*[2])", 0., 1000., 3)
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

gr_reso_E.Draw("AP")
resoFit.Draw("LSAME")
cE1.Update()
cE1.SaveAs("Ereso_matrixcorr_full.root")


cE2=TCanvas("", "", 800, 600)
#gr_response_E.SetTitle(" ")
#gr_response_E.GetYaxis().SetTitle("Photon #mu_{E}/E")
#gr_response_E.GetYaxis().SetTitleOffset(1.4)
#gr_response_E.GetXaxis().SetTitleOffset(1.2)
#gr_response_E.GetXaxis().SetTitle("True photon energy [GeV]")
#gr_response_E.SetLineColor(kOrange)
#gr_response_E.SetLineWidth(2)
#gr_response_E.Draw("AP")
h_response_profile_E.SetTitle("")
h_response_profile_E.GetYaxis().SetTitle("Photon (E_{truth}/E_{reco}")
h_response_profile_E.GetYaxis().SetTitleOffset(1.4)
h_response_profile_E.GetXaxis().SetTitleOffset(1.2)
h_response_profile_E.GetXaxis().SetTitle("True photon energy [GeV]")
h_response_profile_E.SetLineColor(kViolet)
h_response_profile_E.SetLineWidth(2)
h_response_profile_E.Draw("HIST PE")
cE2.Update()
cE2.SaveAs("Eresp_corr.root")

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
h_response_profile_th.GetYaxis().SetTitle("Photon E_{true}/E_{reco}")
h_response_profile_th.GetYaxis().SetTitleOffset(1.4)
h_response_profile_th.GetXaxis().SetTitleOffset(1.2)
h_response_profile_th.GetXaxis().SetTitle("True photon energy [GeV]")
h_response_profile_th.SetLineColor(kViolet)
h_response_profile_th.SetLineWidth(2)
h_response_profile_th.Draw("HIST PE")
cTh2.Update()
cTh2.SaveAs("ThResp_corr.root")

'''
cThE=TCanvas("", "", 800, 600)
h_2d_response.SetTitle("")
h_2d_response.GetYaxis().SetTitle("Photon truth E [GeV]")
h_2d_response.GetYaxis().SetTitleOffset(1.4)
h_2d_response.GetXaxis().SetTitleOffset(1.2)
h_2d_response.GetXaxis().SetTitle("Photon #theta [rad]")
h_2d_response.Draw("COLZ")
cThE.Update()
cThE.SaveAs("ThE_2D_prof.root")
'''
'''
cThE2=TCanvas("", "", 800, 600)
h_2d_resp_corr.SetTitle("")
h_2d_resp_corr.GetYaxis().SetTitle("Photon truth E [GeV]")
h_2d_resp_corr.GetYaxis().SetTitleOffset(1.4)
h_2d_resp_corr.GetXaxis().SetTitleOffset(1.2)
h_2d_resp_corr.GetXaxis().SetTitle("Photon #theta [rad]")
h_2d_resp_corr.Draw("COLZ")
cThE2.Update()
cThE2.SaveAs("ThE_2D_prof_corr.root")
'''

#fit the middle range to get the constant
#then we can worry about edge effects later, but at least we can give Tova some decent Barrel plots
'''
calib_fit=TF1("calibfit", "[0]*(2^(4.01/sin(x)))", 0.678, 2.46, 1)
calib_fit.SetLineColor(kBlue)
calib_fit.SetParName(0, "Constant")
calib_fit.FixParameter(0,0.15)
gStyle.SetOptFit(0)
h_response_profile_th.Fit(calib_fit)
h_response_profile_th.Draw("HIST PE")
calib_fit.Draw("LSAME")
cTh2.Update()
cTh2.SaveAs("ThResp_fitted.root")

corr_factor=calib_fit.GetParameter(0)
print("Correction factor: ",corr_factor)

'''
