from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TFile, TLorentzVector, TMath, TTree, TVector3
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch

#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_photons.root',
                  type=str, default='ntup_photons.root')
(options, args) = parser.parse_args()

arrBins_theta = array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000., 2500., 5000.))

#Global Flag
do_calo_tree = False

# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_theta = TH1D('matched_theta', 'matched_theta', len(arrBins_theta)-1, arrBins_theta)
h_matched_E = TH1D('matched_E', 'matched_E', len(arrBins_E)-1, arrBins_E)

h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_E)-1, arrBins_E) #GeV: true generator pT


h_Npfo = TH1D('Npfo', "Npfo", 1000, 0, 1000)
h_EMpfo_E = TH1D('EMpfo_E', 'EMpfo_E', len(arrBins_E)-1, arrBins_E)
h_HADpfo_E = TH1D('HADpfo_E', 'HADpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedEMpfo_E = TH1D('matchedEMpfo_E', 'matchedEMpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedEMpfo_theta = TH1D('matchedEMpfo_theta', 'matchedEMpfo_theta',
                            len(arrBins_theta)-1, arrBins_theta)
h_matchedHADpfo_E = TH1D('matchedHADpfo_E', 'matchedHADpfo_E', len(arrBins_E)-1, arrBins_E)
h_matchedHADpfo_theta = TH1D('matchedHADpfo_theta', 'matchedHADpfo_theta',
                             len(arrBins_theta)-1, arrBins_theta)
h_pfo_type = TH1D('pfo_type', "pfo_type", 3000, 0, 3000)
h_jet_const_type = TH1D('jet_const_type', "jet_const_type", 3000, 0, 3000) #pfo object PDG number

h_deltaEM_E = TH1D('deltaEM_E', 'deltaEM_E', 250, -1000, 1000)
h_deltaHAD_E = TH1D('deltaHAD_E', 'deltaHAD_E', 250, -1000, 1000)
h_delta_E_sumE = TH1D('delta_E_sumE', 'delta_E_sumE', 250, -1000, 1000)

# Low-level sim hit distributions
h_ECAL_simhit_E = TH1D('ECAL_simhit_E', 'ECAL_simhit_E', 100, 0, 20)  # GeV
h_ECAL_simhit_layer = TH1D(
    'ECAL_simhit_layer', 'ECAL_simhit_layer', 100, 0, 100)
h_ECAL_simhit_layer_ele = TH1D(
    'ECAL_simhit_layer_ele', 'ECAL_simhit_layer_ele', 100, 0, 100)
h_ECAL_simhit_layer_gamma = TH1D(
    'ECAL_simhit_layer_gamma', 'ECAL_simhit_layer_gamma', 100, 0, 100)
h_ECAL_simhit_layer_other = TH1D(
    'ECAL_simhit_layer_other', 'ECAL_simhit_layer_other', 100, 0, 100)

h_HCAL_simhit_E = TH1D('HCAL_simhit_E', 'HCAL_simhit_E', 100, 0, 20)  # GeV
h_HCAL_simhit_layer = TH1D(
    'HCAL_simhit_layer', 'HCAL_simhit_layer', 100, 0, 100)

# Low-level digitised hit distributions
h_ECAL_hit_time = TH1D('ECAL_hit_time', 'ECAL_hit_time', 100, -10, 10)  # ns
h_ECAL_hit_E = TH1D('ECAL_hit_E', 'ECAL_hit_E', 100, 0, 20)  # GeV
h_ECAL_hit_R = TH1D('ECAL_hit_R', 'ECAL_hit_R', 100, 1700, 4000)  # m
h_ECAL_hit_layer = TH1D('ECAL_hit_layer', 'ECAL_hit_layer', 100, 0, 100)

h_HCAL_hit_time = TH1D('HCAL_hit_time', 'HCAL_hit_time', 100, -10, 10)  # ns
h_HCAL_hit_E = TH1D('HCAL_hit_E', 'HCAL_hit_E', 100, 0, 20)  # GeV
h_HCAL_hit_R = TH1D('HCAL_hit_R', 'HCAL_hit_R', 100, 1700, 4000)  # m
h_HCAL_hit_layer = TH1D('HCAL_hit_layer', 'HCAL_hit_layer', 100, 0, 100)

# Aggregated energy info
h_sumE = TH1D('sumE', 'sumE', 120, 0, 6000)  # GeV
h_ECAL_sumE = TH1D('ECAL_sumE', 'ECAL_sumE', 120, 0, 6000)  # GeV
h_HCAL_sumE = TH1D('HCAL_sumE', 'HCAL_sumE', 120, 0, 6000)  # GeV
h_EMfrac = TH1D('EMfrac', 'EMfrac', 100, 0, 1)  # GeV
h_EMfrac_PFO = TH1D('EMfrac_PFO', 'EMfrac_PFO', 100, 0, 1)  # GeV

### Jet Histos
h_jet_energy = TH1D('jet_energy', 'jet_energy', len(arrBins_E)-1, arrBins_E)
h_jet_dR = TH1D('jet_dR', 'jet_dR', 100, 0, TMath.Pi())
h_min_jet_dR = TH1D('min_jet_dR', 'min_jet_dR', 100, 0, TMath.Pi())
h_jet_multiplicity = TH1D('jet_multiplicity', 'jet_multiplicity', 10, 0, 10)


# Histo list for writing to outputs
histos_list = [h_truth_E, h_truth_theta, h_matched_theta,
               h_EMpfo_E, h_matched_E,
               h_HADpfo_E,
               h_matchedEMpfo_E, h_matchedEMpfo_theta,
               h_matchedHADpfo_E, h_matchedHADpfo_theta,
               h_deltaEM_E, h_deltaHAD_E, h_delta_E_sumE,
               h_Npfo, h_pfo_type,
               h_jet_const_type, h_truth_pT,
               h_ECAL_hit_time, h_ECAL_hit_E, h_ECAL_hit_R,
               h_HCAL_hit_time, h_HCAL_hit_E, h_HCAL_hit_R,
               h_ECAL_simhit_E, h_HCAL_simhit_E,
               h_ECAL_sumE, h_HCAL_sumE,
               h_EMfrac, h_EMfrac_PFO,
               h_ECAL_hit_layer, h_HCAL_hit_layer,
               h_ECAL_simhit_layer, h_ECAL_simhit_layer_ele, h_ECAL_simhit_layer_gamma, h_ECAL_simhit_layer_other,
               h_HCAL_simhit_layer, 
               h_jet_energy, h_jet_dR, h_min_jet_dR, h_jet_multiplicity
               ]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
photon_tree = TTree("photon_tree", "photon_tree")
E = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
E_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
pT = array('d', [0])
pT_truth = array('d', [0])
photon_tree.Branch("E",  E,  'var/D')
photon_tree.Branch("phi", phi, 'var/D')
photon_tree.Branch("theta", theta, 'var/D')
photon_tree.Branch("E_truth",  E_truth,  'var/D')
photon_tree.Branch("phi_truth", phi_truth, 'var/D')
photon_tree.Branch("theta_truth", theta_truth, 'var/D')

####################################
calo_tree = TTree("calo_tree", "calo_tree")
E_cell = array('d', [0])
#theta_cell = array('d', [0])
E_sim = array('d', [0])
#theta_sim = array('d', [0])
isECAL = array('i', [0])
calo_tree.Branch("E_cell",  E_cell,  'var/D')
#calo_tree.Branch("theta_cell", theta_cell, 'var/D')
calo_tree.Branch("E_sim",  E_sim,  'var/D')
#calo_tree.Branch("theta_sim", theta_sim, 'var/D')
calo_tree.Branch("isECAL",  isECAL,  'var/I')

##############jet tree
jet_tree = TTree("jet_tree", "jet_tree")
photon_E = array('d', [0]) #initial neutron energy
photon_pT = array('d', [0]) #initial neutron pT
photon_phi_angle = array('d', [0]) #initial neutron phi
photon_theta_angle = array('d', [0]) #initial neutron theta
#neutron_index = array('d', [0]) #initial neutron indexing number
jet_energy = array('d', [0]) #reconstructed jet energy
jet_pT = array('d', [0]) #reconstructed jet pT
jet_phi = array('d', [0]) #reconstructed jet phi
jet_theta = array('d', [0]) #reconstructed jet theta
jet_dR = array('d', [0]) #angular distance between jet and neutron

jet_tree.Branch("photon_E", photon_E, 'var/D')
jet_tree.Branch("photon_pT", photon_pT, 'var/D')
jet_tree.Branch("photon_phi_angle", photon_phi_angle, 'var/D')
jet_tree.Branch("photon_theta_angle", photon_theta_angle, 'var/D')
#jet_tree.Branch("neutron_index", neutron_index, 'var/I')
jet_tree.Branch("jet_energy", jet_energy, 'var/D')
jet_tree.Branch("jet_pT", jet_pT, 'var/D')
jet_tree.Branch("jet_phi", jet_phi, 'var/D')
jet_tree.Branch("jet_theta", jet_theta, 'var/D')
jet_tree.Branch("jet_dR", jet_dR, 'var/D')

to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

filenum=0

for file in to_process:
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    reader.open(file)
    filenum=filenum+1
    if filenum == 138:
        continue
    if filenum == 815:
        continue
    # loop over all events in the file
    for ievt, event in enumerate(reader):
        if filenum == 139 and ievt == 0:
            continue
        if ievt % 1 == 0:
            print(" ")
            print("File "+str(filenum))
            print("Processing event " + str(ievt))

        # Fill the truth-level histos, the first particle is always the gun
        mcpCollection = event.getCollection('MCParticle')
        h_truth_E.Fill(mcpCollection[0].getEnergy())
        dp3 = mcpCollection[0].getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcpCollection[0].getEnergy())
        h_truth_theta.Fill(tlv.Theta())

        E_truth[0] = mcpCollection[0].getEnergy()
        phi_truth[0] = tlv.Phi()
        theta_truth[0] = tlv.Theta()

        '''
        for mcp in mcpCollection:
            if mcp.getGeneratorStatus() == 1 and len(mcp.getParents()) == 0:
                dp3 = mcp.getMomentum()
                tlv = TLorentzVector()
                tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], mcp.getEnergy())
                h_truth_E.Fill(mcp.getEnergy())
                h_truth_theta.Fill(tlv.Theta())
        '''
        
        # Fill the reco-level histos
        pfoCollection = event.getCollection('PandoraPFOs')
        h_Npfo.Fill(len(pfoCollection))

        # Match true pfo with closest reco PFO in deltaR
        matchedEM_E = -1.
        matchedEM_theta = -1.
        matchedHAD_E = -1.
        matchedHAD_theta = -1.
        allEM_E = 0.
        allHAD_E = 0.

        minDREM = 999999.
        minDRHAD = 999999.

        E[0] = -1
        phi[0] = -4
        theta[0] = -1
        
        #noncentral = True
        
        for pfo in pfoCollection:
            h_pfo_type.Fill(abs(pfo.getType()))
            dp3 = pfo.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], pfo.getEnergy())

            if abs(pfo.getType()) == 22:
                allEM_E = allEM_E + pfo.getEnergy()
            elif abs(pfo.getType()) == 2112:
                allHAD_E = allHAD_E + pfo.getEnergy()

            dR = tlv_pfo.DeltaR(tlv)

            if pfo.getEnergy()/mcpCollection[0].getEnergy() > 0.: #change threshold here if needed
                if dR < minDREM and abs(pfo.getType()) == 22:
                    minDREM = dR
                    matchedEM_E = pfo.getEnergy()
                    matchedEM_theta = tlv_pfo.Theta()
                    E[0] = pfo.getEnergy()
                    phi[0] = tlv_pfo.Phi()
                    theta[0] = tlv_pfo.Theta()
                    #if abs(theta[0]) < TMath.Pi()/6:
                     #   print("central")
                      #  noncentral=False
                if dR < minDRHAD and abs(pfo.getType()) == 2112:
                    minDRHAD = dR
                    matchedHAD_E = pfo.getEnergy()
                    matchedHAD_theta = tlv_pfo.Theta()

       
        print(allEM_E, allHAD_E, matchedEM_E, matchedHAD_E)

        h_EMpfo_E.Fill(allEM_E)
        h_HADpfo_E.Fill(allHAD_E)

        if matchedEM_E > 0:
            h_matchedEMpfo_E.Fill(matchedEM_E)
            h_matchedEMpfo_theta.Fill(matchedEM_theta)
            h_deltaEM_E.Fill(matchedEM_E-mcpCollection[0].getEnergy())
        if matchedHAD_E > 0:
            h_matchedHADpfo_E.Fill(matchedHAD_E)
            h_matchedHADpfo_theta.Fill(matchedHAD_theta)
            h_deltaHAD_E.Fill(matchedHAD_E-mcpCollection[0].getEnergy())

        if allHAD_E+allEM_E > 0:
            h_EMfrac_PFO.Fill(allEM_E/(allHAD_E+allEM_E))


        jetCollection = event.getCollection('JetOut')
        minDRjet = 999999.

        E[0] = -1
        pT[0] = -1
        phi[0] = -4
        theta[0] = -1
        #final_minDRjet[0] = -TMath.Pi()

        jet_energy[0] = -1
        jet_dR[0] = -1

        #jet_multiplicity[0] = len(jetCollection)
        h_jet_multiplicity.Fill(len(jetCollection))
        #n_ind = 0
        for jet in jetCollection:
            #n_ind += 1
            dp3 = jet.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            dR = tlv_pfo.DeltaR(tlv)
            h_jet_energy.Fill(jet.getEnergy())
            h_jet_dR.Fill(dR)

            photon_E[0] = E_truth[0]
            photon_pT[0] = pT_truth[0]
            photon_phi_angle[0] = phi_truth[0]
            photon_theta_angle[0] = theta_truth[0]
            #neutron_index[0] = n_ind
            jet_energy[0] = jet.getEnergy()
            jet_pT[0] = tlv_pfo.Perp()
            jet_phi[0] = tlv_pfo.Phi()
            jet_theta[0] = tlv_pfo.Theta()
            jet_dR[0] = dR

            '''if jet_n > 0:
                #print("a jet is found!")
                if jet_energy[0] == -1:
                    jet_energy[0] = jet.getEnergy()
                else:
                    jet_energy.append(jet.getEnergy())
                jet_dR.append(dR)

            if jet_energy[0] == -1:
                jet_energy[0] = jet.getEnergy()
                #print("1st jet found: jet_energy = {:3e}".format(jet_energy[0]))
            #else:
                #jet_energy.append(jet.getEnergy())
                #print("2nd jet found: jet_energy = {:3e}".format(jet_energy[-1]))

            if jet_dR[0] == -1:
                jet_dR[0] = dR
            else:
                jet_dR.append(dR)'''

            if dR < minDRjet:
                minDRjet = dR
                E[0] =jet.getEnergy()
                pT[0] = tlv.Perp()
                phi[0] = tlv_pfo.Phi()
                theta[0] = tlv_pfo.Theta()

            for constituent in jet.getParticles():
                h_jet_const_type.Fill(abs(constituent.getType()))

            jet_tree.Fill()

        #final_minDRjet[0] = minDRjet
        h_min_jet_dR.Fill(minDRjet)
        if minDRjet<0.2:
            h_matched_theta.Fill(tlv.Theta())
            h_matched_E.Fill(tlv.Energy())

        clusterCollection = event.getCollection('PandoraClusters')
        minDRclus = 99999.
        for cluster in clusterCollection:

            px = cluster.getEnergy()*sin(cluster.getITheta())*cos(cluster.getIPhi())
            py = cluster.getEnergy()*sin(cluster.getITheta())*sin(cluster.getIPhi())
            pz = cluster.getEnergy()*cos(cluster.getITheta())

            tlv_clus = TLorentzVector()
            tlv_clus.SetPxPyPzE(px, py, pz, cluster.getEnergy())

            dR = tlv_clus.DeltaR(tlv)

            if dR < minDRclus:
                minDRclus = dR
                '''
                E[0] = cluster.getEnergy()
                phi[0] = cluster.getIPhi()
                theta[0] = cluster.getITheta()
                '''
        #if noncentral:
        photon_tree.Fill()

        for jet in jetCollection:
            dp3 = jet.getMomentum()
            tlv_pfo = TLorentzVector()
            tlv_pfo.SetPxPyPzE(dp3[0], dp3[1], dp3[2], jet.getEnergy())

            dR = tlv_pfo.DeltaR(tlv)

            if dR < minDRjet:
                minDRjet = dR
        
        if minDRjet<0.4:
            h_matched_theta.Fill(tlv.Theta())
            h_matched_E.Fill(tlv.Energy())
        
        # Fill the hit-level histos and aggregated energy
        ECAL_sumE = 0.
        ecal_coll = ['EcalBarrelCollectionSel','EcalEndcapCollectionSel']
        ecal_relcoll = ['EcalBarrelRelationsSimSel','EcalEndcapRelationsSimSel']
        
        for icoll, coll in enumerate(ecal_coll):

            try:
                ECALhitCollection = event.getCollection(coll)

                relationCollection = event.getCollection(ecal_relcoll[icoll])
                relation = UTIL.LCRelationNavigator(relationCollection)

                encoding = ECALhitCollection.getParameters(
                ).getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)

                for hit in ECALhitCollection:
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder["layer"].value()

                    h_ECAL_hit_time.Fill(hit.getTime())
                    h_ECAL_hit_E.Fill(hit.getEnergy())
                    h_ECAL_hit_layer.Fill(layer, hit.getEnergy())

                    ECAL_sumE = ECAL_sumE + hit.getEnergy()
                    pos = hit.getPosition()
                    h_ECAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))

                    E_cell[0] = hit.getEnergy()
                    sim_hit = relation.getRelatedToObjects(hit)[0]
                    E_sim[0] = sim_hit.getEnergy()
                    isECAL[0] = 1

                    calo_tree.Fill()

                h_ECAL_sumE.Fill(ECAL_sumE)
            except:
                print("No", coll, "found")

        HCAL_sumE = 0.
        hcal_coll = ['HcalBarrelsCollectionSel','HcalEndcapsCollectionSel']
        hcal_relcoll = ['HcalBarrelRelationsSimSel','HcalEndcapRelationsSimSel']

        for icoll, coll in enumerate(hcal_coll):

            try:
                HCALhitCollection = event.getCollection(coll)

                relationCollection = event.getCollection(hcal_relcoll[icoll])
                relation = UTIL.LCRelationNavigator(relationCollection)

                encoding = HCALhitCollection.getParameters(
                ).getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)

                for hit in HCALhitCollection:
                    cellID = int(hit.getCellID0())
                    decoder.setValue(cellID)
                    layer = decoder["layer"].value()

                    h_HCAL_hit_time.Fill(hit.getTime())
                    h_HCAL_hit_E.Fill(hit.getEnergy())
                    h_HCAL_hit_layer.Fill(layer, hit.getEnergy())
                    HCAL_sumE = HCAL_sumE + hit.getEnergy()
                    pos = hit.getPosition()
                    h_HCAL_hit_R.Fill(sqrt(pos[0]*pos[0]+pos[1]*pos[1]))

                    E_cell[0] = hit.getEnergy()
                    sim_hit = relation.getRelatedToObjects(hit)[0]
                    E_sim[0] = sim_hit.getEnergy()
                    isECAL[0] = 0

                    calo_tree.Fill()

                h_HCAL_sumE.Fill(HCAL_sumE)
            except:
                print("No", coll, "found")

        print(ECAL_sumE, HCAL_sumE)

        h_sumE.Fill(ECAL_sumE+HCAL_sumE)

        if ECAL_sumE+HCAL_sumE > 0:
            h_EMfrac.Fill(ECAL_sumE/(ECAL_sumE+HCAL_sumE))
        else:
            h_EMfrac.Fill(0)
        h_delta_E_sumE.Fill(ECAL_sumE+HCAL_sumE-mcpCollection[0].getEnergy())

    reader.close()

# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
photon_tree.Write()
jet_tree.Write()
if do_calo_tree:
    calo_tree.Write()
output_file.Close()
