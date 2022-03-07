#include <random>
#include <vector>
#include "TH2.h"
#include "TRotation.h"
#include "TVector3.h"

/* 
This code will read in a csv file containing muon events in zenith and azimuth
and will sample the distribution to for generating events. A check is made to 
the existing GetRandom function which is used in root. Events are smeared based
on the bin widths via a Gaussian to avoid binning artifacts in the generation of
the events 
*/

double pi = 3.14159;

// SEEDING in NEXUS:
// std::cout << "Krish: The RN is:" << CLHEP::HepRandom::getTheSeed() << " " << event->GetEventID() << std::endl;

// Get the bin widths to smear the events by
std::vector<double> GetBinWidths(std::vector<double> bins){

    // Vector of Bin Widths
    std::vector<double> BW = {};

    for (int i = 0; i < bins.size()-1; i++){
        BW.push_back((bins.at(i+1) - bins[i]));
        // std::cout << bins[i] << std::endl;
    }

    return BW;

}

void GenerateRandom(std::vector<double> weights, std::vector<double> az, std::vector<double> zen, std::vector<double> az_edges, std::vector<double> zen_edges, TH2D* hist, TH2D* histXY, TH2D* histXZ, TH2D* histYZ, TH1D* histZ){

    gStyle->SetOptStat(0);

    // Discrete distribution with bin index and intensity to sample from
    std::discrete_distribution<int> discr_dist(std::begin(weights), std::end(weights));
    std::mt19937 RN_engine;  // Random Number Generator Engine
    RN_engine.seed(17392); // Seed with the current time

    // Create random number generator to smear bins by a small amount in azimuth and zenith
    std::mt19937 RN_engine_az;
    RN_engine_az.seed(17392+1);

    std::mt19937 RN_engine_zen;
    RN_engine_zen.seed(17392+2); // Extra factor 3e4 to keep the seeds unique
    
    int N_samp = 1e7;
    std::vector<int> samples(N_samp);

    // Rotate the azimuth and zenith
    TRotation *rPhi = new TRotation();
    rPhi->RotateY(-150* pi/180);

    
    std::vector<double> az_BW  = GetBinWidths(az_edges);
    std::vector<double> zen_BW   = GetBinWidths(zen_edges);

    double zen_BW_smear{1.0e6}; // Amount to smear the randomly sampled zenith values by
    double az_BW_smear{1.0e6};  // Amount to smear the randomly sampled azimuth values by

    TVector3 dir;

    // Index to set a unique seed to the generator
    int index = 0;

    for (auto & RN_indx: samples){

        bool invalid_evt = true;
        double az_smear = 0.0;
        double zen_smear = 0.0;
        double az_samp = 0.0;
        double zen_samp = 0.0;

        while(invalid_evt){

            // Generate random index weighted by the bin contents
            RN_indx = discr_dist(RN_engine);

            // Loop over the zenith values and find the corresponding bin width to smear
            for (int i = 0; i < zen_edges.size()-1; i++){

                // Catch very last bin
                if (zen_BW_smear == 1e6 && i == zen_edges.size()-2){
                    if (zen[RN_indx] >= zen_edges[i] 
                        && zen[RN_indx] <= zen_edges[i+1]){
                    
                        zen_BW_smear = zen_BW[i];
                    }
                }
                else {
                    if (zen[RN_indx] >= zen_edges[i] 
                        && zen[RN_indx] < zen_edges[i+1]){
                        
                        zen_BW_smear = zen_BW[i];

                    }
                }

            }

            // Loop over the azimuth values and find the corresponding bin width to smear
            for (int i = 0; i < az_edges.size()-1; i++){
                // Include last bin in check
                if (az_BW_smear == 1e6 && i == az_edges.size()-2){
                    if (az[RN_indx] >= az_edges[i] 
                        && az[RN_indx] <= az_edges[i+1]){
                    
                        az_BW_smear = az_BW[i];

                    }
                }
                else {

                    if (az[RN_indx] >= az_edges[i] 
                        && az[RN_indx] < az_edges[i+1]){
                        
                        az_BW_smear = az_BW[i];

                    }
                }

            }

            if (az_BW_smear == 1.0e6 || zen_BW_smear == 1.0e6 )
                std::cout << "Error BW not set" << std::endl;

            std::normal_distribution<double> dist_az(0, az_BW_smear);
            std::normal_distribution<double> dist_zen(0, zen_BW_smear);
            
            // Get the Gaussian smear values
            az_smear  = dist_az(RN_engine_az);
            zen_smear = dist_zen(RN_engine_zen);

            az_samp = az.at(RN_indx)+az_smear;
            zen_samp = zen.at(RN_indx)+zen_smear;

            if (zen_samp < 0.0){
                invalid_evt = true;
            }
            else
                invalid_evt = false;
        }
        
        // Fill histogram
        hist->Fill(az_samp, zen_samp);

        dir.SetX(sin(zen_samp) * sin(az_samp));
        dir.SetY(-cos(zen_samp));
        dir.SetZ(-sin(zen_samp) * cos(az_samp));

        dir *= *rPhi;

        histXY->Fill(dir.X(), dir.Y());
        histXZ->Fill(dir.X(), dir.Z());
        histYZ->Fill(dir.Y(), dir.Z());
        histZ->Fill(dir.Z());

        index++;

    }

}


// Main
void RandomGenerator(){

    // Load in the histogram
    TFile f("MuonAnaAllRuns.root");
    TH2F *hist;
    f.GetObject("za", hist);
    hist->SetDirectory(0);
    f.Close();

    TH2D* hist_root = new TH2D("hist_root", ";Azimuth; Zenith", hist->GetNbinsX()+1,  0, 2*pi ,hist->GetNbinsY()+1, 0, 0.5*pi );

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
    std::ifstream fin("SimulatedMuonsProposalMCEq.csv");
    
    // Note these are equivalent!!
    // std::ifstream fin("MeasuredMuonsFromData.csv");
    // std::ifstream fin("MuonAnaAllRuns.csv");
    
    // Check if file has opened properly
    if (!fin.is_open())
        std::cout << "Error: Could not read in the input data CSV file" << std::endl;
  
    // Read the Data from the file 
    std::string s_intensity, s_beta, s_alpha;
    std::vector<double> intensity, beta, alpha;
    std::vector<double> azimuth_bins = {}; // List of Azimuth bin edges
    std::vector<double> zenith_bins = {};  // List of Zenith bin edges
  
    // Loop over the lines in the file and add the values to a vector
    while (fin.peek()!=EOF) {
  
        std::getline(fin, s_intensity, ',');

        // Load in alpha/zenith bin edges
        if (s_intensity == "zenith"){
            std::getline(fin, s_alpha, '\n');
            zenith_bins.push_back(stod(s_alpha));
        }
        // Load in beta/azimuth bin edges
        else if (s_intensity == "azimuth"){
            std::getline(fin, s_beta, '\n');
            azimuth_bins.push_back(stod(s_beta));
        }
        // Load in the histogram values
        else {
            std::getline(fin, s_beta, ',');
            std::getline(fin, s_alpha, '\n');

            intensity.push_back(stod(s_intensity));
            beta.push_back(stod(s_beta));
            alpha.push_back(stod(s_alpha));
        }


    }

    int const nbins_zeni = zenith_bins.size()-1;
    double* edges_zeni = &zenith_bins[0]; // Cast to an array
    
    int const nbins_azimuth = azimuth_bins.size()-1;
    double* edges_azimuth = &azimuth_bins[0]; // Cast to an array 


    // TH2D* hist_cpp2   = new TH2D("hist_cpp2", ";Azimuth; Zenith", nbins_azimuth, edges_azimuth , nbins_zeni, edges_zeni);
    TH2D* hist_cpp2   = new TH2D("hist_cpp2", ";Azimuth; Zenith", 50, 0, 2*pi, 50, 0, pi/2);
    TH2D* histXY_cpp2 = new TH2D("histXY_cpp2", ";X; Y", 100, -1, 1 , 50, -1, 0 );
    TH2D* histXZ_cpp2 = new TH2D("histXZ_cpp2", ";X; Z", 100, -1, 1 , 100, -1, 1 );
    TH2D* histYZ_cpp2 = new TH2D("histYZ_cpp2", ";Y; Z", 50, -1, 0 , 100, -1, 1 );

    TH1D* histZ_cpp2 = new TH1D("histZ_cpp2", ";Z; Entries", 75, -1, 1 );

    GenerateRandom(intensity, beta, alpha, azimuth_bins, zenith_bins, hist_cpp2, histXY_cpp2, histXZ_cpp2, histYZ_cpp2, histZ_cpp2);

    TCanvas *c4 = new TCanvas();
    hist_cpp2->Draw("colz");

    // TCanvas *cXY_cpp2 = new TCanvas();
    // histXY_cpp2->Draw("colz");

    TCanvas *cXZ_cpp2 = new TCanvas();
    // gStyle->SetPalette(kViridis);
    // cXZ_cpp2->SetLogz();
    histXZ_cpp2->Draw("colz");

    // TCanvas *cYZ_cpp2 = new TCanvas();
    // histYZ_cpp2->Draw("colz");

    TCanvas *cZ_cpp2 = new TCanvas();
    histZ_cpp2->Draw("hist");

}
