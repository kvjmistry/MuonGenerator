#include <random>
#include <vector>
#include "TH2.h"
#include "TRotation.h"
#include "TVector3.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include "TStyle.h"
#include <fstream>

/* 
This code will read in a csv file containing muon events in zenith and azimuth
and will sample the distribution to for generating events. A check is made to 
the existing GetRandom function which is used in root. Events are smeared based
on the bin widths via a Gaussian to avoid binning artifacts in the generation of
the events 
*/

double pi = 3.14159;

// Main
void RandomGenerator(){

    // Load in the histogram
    TFile f("files/MuonAnaAllRuns.root");
    // TFile f("SimulatedMuonFile.root");
    TH2F *hist;
    f.GetObject("za", hist);
    hist->SetDirectory(0);
    f.Close();

    double y_high = hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY()+1);
    double x_high = hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX()+1);
    TH2D* hist_root = new TH2D("hist_root", ";Azimuth; Zenith", hist->GetNbinsX()+1,  0, x_high ,hist->GetNbinsY()+1, 0, y_high);

    TH2D* histXY = new TH2D("histXY", ";X; Y", 100, -1, 1 , 50, -1, 0 );
    TH2D* histXZ = new TH2D("histXZ", ";X; Z", 100, -1, 1 , 100, -1, 1 );
    TH2D* histYZ = new TH2D("histYZ", ";Y; Z", 50, -1, 0 , 100, -1, 1 );

    // Rotate the azimuth and zenith
    TRotation *rPhi = new TRotation();
    rPhi->RotateY(-150*pi/180);
    
    TVector3 dir;

    for (int i = 0; i < 1e6; i++){
        double zenith  = 0.;
        double azimuth = 0.;
        hist->GetRandom2(azimuth, zenith);


        zenith  *= pi;
        azimuth *= pi;
        hist_root->Fill(azimuth, zenith);

        dir.SetX(sin(zenith) * sin(azimuth));
        dir.SetY(-cos(zenith));
        dir.SetZ(-sin(zenith) * cos(azimuth));

        dir *= *rPhi;

        histXY->Fill(dir.X(), dir.Y());
        histXZ->Fill(dir.X(), dir.Z());
        histYZ->Fill(dir.Y(), dir.Z());

    }

    TCanvas *c = new TCanvas();
    hist_root->Draw("colz");

    // TCanvas *cXY = new TCanvas();
    // histXY->Draw("colz");

    TCanvas *cXZ = new TCanvas();
    histXZ->Draw("colz");

    // TCanvas *cYZ = new TCanvas();
    // histYZ->Draw("colz");

    // TCanvas *c2 = new TCanvas();
    // hist->Draw("colz");
    

    // File pointer
    
    // Simulated muons from MCeQ
    // std::ifstream fin("files/SimulatedMuonsProposalMCEq.csv");
    std::ifstream fin("files/SimulatedMuonsProposalMCEqEnergy.csv");
    // std::ifstream fin("SimulatedMuonFile.csv"); 
    
    // Data
    // std::ifstream fin("files/MuonAnaAllRuns.csv");
    
    // std::string dist_name_ = "za";
    std::string dist_name_ = "zae";

    // Check if file has opened properly
    if (!fin.is_open())
        std::cout << "Error: Could not read in the input data CSV file" << std::endl;
  
    // Read the Data from the file as strings
    std::string s_header, s_flux, s_azimuth, s_zenith, s_energy;
    std::string s_azimuth_smear, s_zenith_smear, s_energy_smear;
    std::vector<double> flux_, azimuths_, zeniths_, energies_; ///< Values of flux, azimuth and zenith from file
    std::vector<double> azimuth_smear_; ///< List of Azimuth bin smear values
    std::vector<double> zenith_smear_;  ///< List of Zenith bin smear values
    std::vector<double> energy_smear_;  ///< List of Energy bin smear values

    std::vector<double> zenith_bins, azimuth_bins, energy_bins; 

    // Loop over the lines in the file and add the values to a vector
    while (fin.peek()!=EOF) {

        std::getline(fin, s_header, ',');

        // Angle input only
        if (s_header == "value" && dist_name_ == "za"){
            std::getline(fin, s_flux, ',');
            std::getline(fin, s_azimuth, ',');
            std::getline(fin, s_zenith, ',');
            std::getline(fin, s_azimuth_smear, ',');
            std::getline(fin, s_zenith_smear, '\n');

            flux_.push_back(stod(s_flux));
            azimuths_.push_back(stod(s_azimuth));
            zeniths_.push_back(stod(s_zenith));
            azimuth_smear_.push_back(stod(s_azimuth_smear));
            zenith_smear_.push_back(stod(s_zenith_smear));
        }
        // Angle + Energy input
        if (s_header == "value" && dist_name_ == "zae"){
            std::getline(fin, s_flux, ',');
            std::getline(fin, s_azimuth, ',');
            std::getline(fin, s_zenith, ',');
            std::getline(fin, s_energy, ',');
            std::getline(fin, s_azimuth_smear, ',');
            std::getline(fin, s_zenith_smear, ',');
            std::getline(fin, s_energy_smear, '\n');

            flux_.push_back(stod(s_flux));
            azimuths_.push_back(stod(s_azimuth));
            zeniths_.push_back(stod(s_zenith));
            energies_.push_back(stod(s_energy));
            azimuth_smear_.push_back(stod(s_azimuth_smear));
            zenith_smear_.push_back(stod(s_zenith_smear));
            energy_smear_.push_back(stod(s_energy_smear));
        }
        if (s_header == "zenith"){
            std::getline(fin, s_zenith, '\n');
            zenith_bins.push_back(stod(s_zenith));

        }
        if (s_header == "azimuth"){
            std::getline(fin, s_azimuth, '\n');
            azimuth_bins.push_back(stod(s_azimuth));

        }
        if (s_header == "energy"){
            std::getline(fin, s_energy, '\n');
            energy_bins.push_back(stod(s_energy));

        }

    }

    int const nbins_zeni = zenith_bins.size()-1;
    double* edges_zeni = &zenith_bins[0]; // Cast to an array
    
    int const nbins_azimuth = azimuth_bins.size()-1;
    double* edges_azimuth = &azimuth_bins[0]; // Cast to an array 


    TH2D* hist_cpp2   = new TH2D("hist_cpp2", ";Azimuth; Zenith", nbins_azimuth, edges_azimuth , nbins_zeni, edges_zeni);
    // TH2D* hist_cpp2   = new TH2D("hist_cpp2", ";Azimuth; Zenith", hist->GetNbinsX()+1,  0, x_high ,hist->GetNbinsY()+1, 0, y_high);
    TH2D* histXY_cpp2 = new TH2D("histXY_cpp2", ";X; Y", 100, -1, 1 , 50, -1, 0 );
    TH2D* histXZ_cpp2 = new TH2D("histXZ_cpp2", ";X; Z", 100, -1, 1 , 100, -1, 1 );
    TH2D* histYZ_cpp2 = new TH2D("histYZ_cpp2", ";Y; Z", 50, -1, 0 , 100, -1, 1 );

    TH1D* histEnergy = new TH1D("histEnergy", ";Energy;Entries", 20, 200, 2000 );

    gStyle->SetOptStat(0);

    // Discrete distribution with bin index and intensity to sample from
    std::discrete_distribution<int> discr_dist(std::begin(flux_), std::end(flux_));
    std::mt19937 RN_engine;  // Random Number Generator Engine
    RN_engine.seed(17392); // Seed with the current time

    // Create random number generator to smear bins by a small amount in azimuth and zenith
    std::mt19937 RN_engine_az;
    RN_engine_az.seed(17392+1);

    std::mt19937 RN_engine_zen;
    RN_engine_zen.seed(17392+2); // Extra factor 3e4 to keep the seeds unique

    std::mt19937 RN_engine_E;
    RN_engine_E.seed(17392+3); // Extra factor 3e4 to keep the seeds unique
    
    int N_samp = 1e6;
    std::vector<int> samples(N_samp);

    // Index to set a unique seed to the generator
    int index = 0;

    for (auto & RN_indx: samples){

        bool invalid_evt = true;
        double az_smear  = 0;
        double zen_smear = 0;
        double energy_smear = 0;

        double az_samp = 0.0;
        double zen_samp = 0.0;
        double energy_samp = 0.0;

        while(invalid_evt){

            // Generate random index weighted by the bin contents
            RN_indx = discr_dist(RN_engine);

            std::normal_distribution<double> dist_az(0,  azimuth_smear_.at(RN_indx));
            std::normal_distribution<double> dist_zen(0, zenith_smear_.at(RN_indx));
            
            // Get the Gaussian smear values
            az_smear  = dist_az(RN_engine_az);
            zen_smear = dist_zen(RN_engine_zen);

            az_samp = azimuths_.at(RN_indx)+az_smear;
            zen_samp = zeniths_.at(RN_indx)+zen_smear;

            if (zen_samp < 0.0){
                invalid_evt = true;
            }
            else
                invalid_evt = false;

            if (dist_name_ == "zae"){
                std::normal_distribution<double> dist_E(0,  energy_smear_.at(RN_indx));
                energy_smear = dist_E(RN_engine_E);

                energy_samp = energies_.at(RN_indx)+energy_smear;
                histEnergy->Fill(energy_samp);

            }
        }
        

        dir.SetX(sin(zen_samp) * sin(az_samp));
        dir.SetY(-cos(zen_samp));
        dir.SetZ(-sin(zen_samp) * cos(az_samp));

        dir *= *rPhi;

        // Fill histogram
        hist_cpp2->Fill(az_samp, zen_samp);
    
        histXY_cpp2->Fill(dir.X(), dir.Y());
        histXZ_cpp2->Fill(dir.X(), dir.Z());
        histYZ_cpp2->Fill(dir.Y(), dir.Z());

        index++;

    }

    TCanvas *c4 = new TCanvas();
    hist_cpp2->Draw("colz");

    TCanvas *cXY_cpp2 = new TCanvas();
    histXY_cpp2->Draw("colz");

    TCanvas *cXZ_cpp2 = new TCanvas();
    histXZ_cpp2->Draw("colz");

    TCanvas *cYZ_cpp2 = new TCanvas();
    histYZ_cpp2->Draw("colz");

    if (dist_name_ == "zae"){
        TCanvas *cE_cpp2 = new TCanvas();
        histEnergy->Draw("hist");
    }

}
