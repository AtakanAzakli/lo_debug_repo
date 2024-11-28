#include "TChain.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TMath.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootClasses.h"

void processFile(const char* inputFile, const char* outputFile)
{
    // Load the LHEF events file into a TChain
    TChain chain("LHEF");
    chain.Add(inputFile);
 
    // Open a file to save the output
    TFile *fout = new TFile(outputFile, "RECREATE");

    // Create a TTree to store results with all the new features
    TTree *tree = new TTree("eventData", "Event Data with Features");

    // Declare variables
    int l0_id, l1_id, q0_id, q1_id;
    double l0_pt, l1_pt, q0_pt, q1_pt;
    double l0_phi, l1_phi, q0_phi, q1_phi;
    double l0_eta, l1_eta, q0_eta, q1_eta;
    double l0_m, l1_m, q0_m, q1_m;
    double l0_e, l1_e, q0_e, q1_e;
    double met_et, met_phi, m_ll, m_qq;
    double pt_ll, pt_qq, d_phi_ll, d_phi_qq;
    double d_eta_ll, d_eta_qq, d_y_ll, d_y_qq;
    double sqrtHT, MET_sig, m_l0q0, m_l0q1, m_l1q0, m_l1q1;

    // Create branches
    tree->Branch("l0_id", &l0_id);
    tree->Branch("l1_id", &l1_id);
    tree->Branch("q0_id", &q0_id);
    tree->Branch("q1_id", &q1_id);
    tree->Branch("l0_pt", &l0_pt);
    tree->Branch("l1_pt", &l1_pt);
    tree->Branch("q0_pt", &q0_pt);
    tree->Branch("q1_pt", &q1_pt);
    tree->Branch("l0_phi", &l0_phi);
    tree->Branch("l1_phi", &l1_phi);
    tree->Branch("q0_phi", &q0_phi);
    tree->Branch("q1_phi", &q1_phi);
    tree->Branch("l0_eta", &l0_eta);
    tree->Branch("l1_eta", &l1_eta);
    tree->Branch("q0_eta", &q0_eta);
    tree->Branch("q1_eta", &q1_eta);
    tree->Branch("l0_m", &l0_m);
    tree->Branch("l1_m", &l1_m);
    tree->Branch("q0_m", &q0_m);
    tree->Branch("q1_m", &q1_m);
    tree->Branch("l0_e", &l0_e);
    tree->Branch("l1_e", &l1_e);
    tree->Branch("q0_e", &q0_e);
    tree->Branch("q1_e", &q1_e);
    tree->Branch("met_et", &met_et);
    tree->Branch("met_phi", &met_phi);
    tree->Branch("m_ll", &m_ll);
    tree->Branch("m_qq", &m_qq);
    tree->Branch("pt_ll", &pt_ll);
    tree->Branch("pt_qq", &pt_qq);
    tree->Branch("d_phi_ll", &d_phi_ll);
    tree->Branch("d_phi_qq", &d_phi_qq);
    tree->Branch("d_eta_ll", &d_eta_ll);
    tree->Branch("d_eta_qq", &d_eta_qq);
    tree->Branch("d_y_ll", &d_y_ll);
    tree->Branch("d_y_qq", &d_y_qq);
    tree->Branch("sqrtHT", &sqrtHT);
    tree->Branch("MET_sig", &MET_sig);
    tree->Branch("m_l0q0", &m_l0q0);
    tree->Branch("m_l0q1", &m_l0q1);
    tree->Branch("m_l1q0", &m_l1q0);
    tree->Branch("m_l1q1", &m_l1q1);

    // Create ExRootTreeReader to process the events
    ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader->GetEntries();

    // Access the "Particle" branch
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");

    // Loop over all events in the file
    for(Long64_t entry = 0; entry < numberOfEntries; ++entry)
    {
        // Read event data
        treeReader->ReadEntry(entry);
        if (entry % 100000 == 0) cout << "Processing Event " << entry << endl;

        // Initialize variables to store particles and their IDs
        TLorentzVector electron, positron;
        int electron_id = 0, positron_id = 0;
        TLorentzVector upQuark, downQuark;
        int upQuark_id = 0, downQuark_id = 0;

        // Variables for MET
        double met_px = 0.0;
        double met_py = 0.0;

        // Loop over particles in the event
        for(Int_t part_i = 0; part_i < branchParticle->GetEntries(); ++part_i)
        {
            TRootLHEFParticle *particle = (TRootLHEFParticle*) branchParticle->At(part_i);

            // Electrons and positrons
            if(particle->PID == 11) // Electron
            {
                electron.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                electron_id = particle->PID;
            }
            else if(particle->PID == -11) // Positron
            {
                positron.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                positron_id = particle->PID;
            }

            // Up and down quarks
            else if(particle->PID == 2) // Up quark
            {
                upQuark.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                upQuark_id = particle->PID;
            }
            else if(particle->PID == 1) // Down quark
            {
                downQuark.SetPtEtaPhiE(particle->PT, particle->Eta, particle->Phi, particle->E);
                downQuark_id = particle->PID;
            }

            // Neutrinos contribute to MET
            else if(abs(particle->PID) == 12 || abs(particle->PID) == 14 || abs(particle->PID) == 16) // Neutrinos
            {
                met_px += particle->PT * cos(particle->Phi);
                met_py += particle->PT * sin(particle->Phi);
            }
        }

        // Calculate MET
        met_et = sqrt(met_px * met_px + met_py * met_py);
        met_phi = atan2(met_py, met_px);

        // Check if all particles were found
        if(electron.Pt() > 0 && positron.Pt() > 0 && upQuark.Pt() > 0 && downQuark.Pt() > 0)
        {
            // Determine the most energetic lepton
            TLorentzVector l0, l1;
            int l0_id_temp, l1_id_temp;

            if(electron.Pt() > positron.Pt())
            {
                l0 = electron;
                l0_id_temp = electron_id;
                l1 = positron;
                l1_id_temp = positron_id;
            }
            else
            {
                l0 = positron;
                l0_id_temp = positron_id;
                l1 = electron;
                l1_id_temp = electron_id;
            }

            // Determine the most energetic quark
            TLorentzVector q0, q1;
            int q0_id_temp, q1_id_temp;

            if(upQuark.Pt() > downQuark.Pt())
            {
                q0 = upQuark;
                q0_id_temp = upQuark_id;
                q1 = downQuark;
                q1_id_temp = downQuark_id;
            }
            else
            {
                q0 = downQuark;
                q0_id_temp = downQuark_id;
                q1 = upQuark;
                q1_id_temp = upQuark_id;
            }

            // Transverse momentum
            l0_pt = l0.Pt();
            l1_pt = l1.Pt();
            q0_pt = q0.Pt();
            q1_pt = q1.Pt();

            // Phi
            l0_phi = l0.Phi();
            l1_phi = l1.Phi();
            q0_phi = q0.Phi();
            q1_phi = q1.Phi();

            // Eta
            l0_eta = l0.Eta();
            l1_eta = l1.Eta();
            q0_eta = q0.Eta();
            q1_eta = q1.Eta();

            // Mass
            l0_m = l0.M();
            l1_m = l1.M();
            q0_m = q0.M();
            q1_m = q1.M();

            // Energy
            l0_e = l0.E();
            l1_e = l1.E();
            q0_e = q0.E();
            q1_e = q1.E();

            // Particle IDs
            l0_id = l0_id_temp;
            l1_id = l1_id_temp;
            q0_id = q0_id_temp;
            q1_id = q1_id_temp;

            // Center-of-mass energies
            m_ll = (l0 + l1).M();
            m_qq = (q0 + q1).M();

            // Transverse momentum of systems
            pt_ll = (l0 + l1).Pt();
            pt_qq = (q0 + q1).Pt();

            // Angle differences
            d_phi_ll = l0.DeltaPhi(l1);
            d_phi_qq = q0.DeltaPhi(q1);

            // Pseudorapidity differences
            d_eta_ll = l0.Eta() - l1.Eta();
            d_eta_qq = q0.Eta() - q1.Eta();

            // Rapidity differences
            d_y_ll = l0.Rapidity() - l1.Rapidity();
            d_y_qq = q0.Rapidity() - q1.Rapidity();

            // sqrtHT
            sqrtHT = sqrt(l0_pt + l1_pt + q0_pt + q1_pt + met_et);

            // MET significance
            MET_sig = met_et / sqrtHT;

            // Center-of-mass energies for lepton-quark pairs
            m_l0q0 = (l0 + q0).M();
            m_l0q1 = (l0 + q1).M();
            m_l1q0 = (l1 + q0).M();
            m_l1q1 = (l1 + q1).M();

            // Fill the TTree
            tree->Fill();
        }
    }

    // Save the results to the output file
    fout->cd();
    tree->Write();
    fout->Close();
}

void ana()
{
    const char* inputFiles[] = {
        "/ML-study-hww-EFT-cHW/events_cHW-0.100000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW-1.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW-2.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW-5.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW-10.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW0.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW0.100000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW1.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW2.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW5.000000/unweighted_events.root",
        "/ML-study-hww-EFT-cHW/events_cHW10.000000/unweighted_events.root"
    };

    const char* outputFiles[] = {
        "output_bsm_-01.root",
        "output_bsm_-1.root",
        "output_bsm_-2.root",
        "output_bsm_-5.root",
        "output_bsm_-10.root",
        "output_sm.root",
        "output_bsm_01.root",
        "output_bsm_1.root",
        "output_bsm_2.root",
        "output_bsm_5.root",
        "output_bsm_10.root"
    };

    for(int i = 0; i < 11; ++i) // Adjust loop size as necessary
    {
        processFile(inputFiles[i], outputFiles[i]);
    }
}