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
std::vector<double> GetSmearBins(std::vector<double> bins){

    std::vector<double> smears = {};

    for (int i = 0; i < bins.size(); i++){
        smears.push_back((bins[i+1] - bins[i]));
    }

    return smears;

}

void GenerateRandom(std::vector<double> weights, std::vector<double> azimuth, std::vector<double> zenith, std::vector<double> azi_edges, std::vector<double> zeni_edges, TH2D* hist, TH2D* histXY, TH2D* histXZ, TH2D* histYZ){

    gStyle->SetOptStat(0);

    // Discrete distribution with bin index and intensity to sample from
    std::discrete_distribution<int> dist(std::begin(weights), std::end(weights));
    std::default_random_engine gen;  // Generator
    gen.seed(time(0)); // Seed with the current time for uniqueness
    int N_samp = 1e6;
    std::vector<int> samples(N_samp);

     // Rotate the azimuth and zenith
    TRotation *rPhi = new TRotation();
    rPhi->RotateY(-140);

    
    std::vector<double> azi_smears = GetSmearBins(azi_edges);
    std::vector<double> zeni_smears = GetSmearBins(zeni_edges);

    double smear_pct_z{0.00}; // Percentage to smear the zenith values by
    double smear_pct_a{0.00}; // Percentage to smear the azimuth values by

    TVector3 dir;

    // Index to set a unique seed to the generator
    int index = 0;

    for (auto & ran_idx: samples){

        // Generate random index weighted by the bin contents
        ran_idx = dist(gen);

        // Loop over the zenith values and find the appropriate smear value
        for (int i = 0; i < zeni_edges.size()-1; i++){
            if (zenith.at(ran_idx) >= zeni_edges.at(i) && zenith.at(ran_idx) < zeni_edges.at(i+1)){
                smear_pct_z = zeni_smears.at(i);
            }
        }

        // Loop over the azimuth values and find the appropriate smear value
        for (int i = 0; i < azi_edges.size()-1; i++){
            if (azimuth.at(ran_idx) >= azi_edges.at(i) && azimuth.at(ran_idx) < azi_edges.at(i+1)){
                smear_pct_a = azi_smears.at(i);
            }
        }


        // Create random number generator to smear bins by a small amount in azimuth and zenith
        std::default_random_engine gen_gauss_azi;
        gen_gauss_azi.seed(1723+index);
        std::normal_distribution<double> dist_azi(0, smear_pct_a);

        std::default_random_engine gen_gauss_zeni;
        gen_gauss_zeni.seed(1723+index+3.0e4); // Extra factor 3e4 to keep the seeds unique
        std::normal_distribution<double> dist_zeni(0, smear_pct_z);
        
        // Get the Gaussian smear values
        double smear_azi  = dist_azi(gen_gauss_azi);
        double smear_zeni = dist_zeni(gen_gauss_zeni);

        double a_ = azimuth.at(ran_idx)+smear_azi;
        double z_ = zenith.at(ran_idx)+smear_zeni;
        
        // Fill histogram
        hist->Fill(a_, z_);

        dir.SetX(sin(z_) * sin(a_));
        dir.SetY(-cos(z_));
        dir.SetZ(-sin(z_) * cos(a_));

        dir *= *rPhi;

        histXY->Fill(dir.X(), dir.Y());
        histXZ->Fill(dir.X(), dir.Z());
        histYZ->Fill(dir.Y(), dir.Z());

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
    rPhi->RotateY(-140);
    
    TVector3 dir;

    for (int i = 0; i < 1e7; i++){
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

    // TCanvas *c = new TCanvas();
    // hist_root->Draw("colz");

    TCanvas *cXY = new TCanvas();
    histXY->Draw("colz");

    TCanvas *cXZ = new TCanvas();
    histXZ->Draw("colz");

    TCanvas *cYZ = new TCanvas();
    histYZ->Draw("colz");

    TCanvas *c2 = new TCanvas();
    hist->Draw("colz");
    

    // File pointer
    std::ifstream fin("MeasuredMuonsFromData.csv");
    // std::ifstream fin("SimulatedMuonsFromProposal.csv");
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


    TH2D* hist_cpp2   = new TH2D("hist_cpp2", ";Azimuth; Zenith", nbins_azimuth, edges_azimuth , nbins_zeni, edges_zeni);
    TH2D* histXY_cpp2 = new TH2D("histXY_cpp2", ";X; Y", 100, -1, 1 , 50, -1, 0 );
    TH2D* histXZ_cpp2 = new TH2D("histXZ_cpp2", ";X; Z", 100, -1, 1 , 100, -1, 1 );
    TH2D* histYZ_cpp2 = new TH2D("histYZ_cpp2", ";Y; Z", 50, -1, 0 , 100, -1, 1 );

    GenerateRandom(intensity, beta, alpha, azimuth_bins, zenith_bins, hist_cpp2, histXY_cpp2, histXZ_cpp2, histYZ_cpp2);

    TCanvas *c4 = new TCanvas();
    hist_cpp2->Draw("colz");

    TCanvas *cXY_cpp2 = new TCanvas();
    histXY_cpp2->Draw("colz");

    TCanvas *cXZ_cpp2 = new TCanvas();
    histXZ_cpp2->Draw("colz");

    TCanvas *cYZ_cpp2 = new TCanvas();
    histYZ_cpp2->Draw("colz");


}
