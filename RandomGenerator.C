// Test code to generate random Numbers 
#include <random>
#include <vector>
#include "TH2.h"
#include "TRotation.h"
#include "TVector3.h"

void GenerateRandom(std::vector<double> weights, std::vector<double> azimuth, std::vector<double> zenith, TH2D* hist, TH2D* histXY, TH2D* histXZ, TH2D* histYZ){

    gStyle->SetOptStat(0);

    // Discrete distribution with bin index and intensity to sample from
    std::discrete_distribution<int> dist(std::begin(weights), std::end(weights));
    std::default_random_engine gen;  // Generator
    gen.seed(time(0)); // Seed with the current time for uniqueness
    int N_samp = 1e6;
    std::vector<int> samples(N_samp);


    std::vector<double> zeni_edges= {0.        , 0.20033484, 0.28379411, 0.34816602, 0.40271584, 0.45102681,
                                     0.49493413, 0.53552665, 0.5735131 , 0.60938531, 0.64350111, 0.67613051,
                                     0.70748321, 0.73772597, 0.76699401, 0.79539883, 0.82303369, 0.84997757,
                                     0.87629806, 0.90205362, 0.92729522, 0.95206764, 0.97641053, 1.00035922,
                                     1.02394538, 1.04719755, 1.07014161, 1.09280113, 1.11519765, 1.13735101,
                                     1.15927948};


    std::vector<double> zeni_smears= {0.20033484232311688, 0.0834592668852091, 0.06437191206463344, 0.05454982030770095, 0.04831097021560099, 0.0439073145446332, 0.04059252797349239, 0.03798645010870905, 0.03587220360769827, 0.034115800762489434, 0.03262940076737708, 0.031352702218681694, 0.030242756673905857, 0.029268039408617774, 0.028404822322277035, 0.027634861950832534, 0.026943873789830985, 0.026320495243533637, 0.02575556242418431, 0.025241594409087353, 0.02477241812103337, 0.024342890671188888, 0.023948690601140243, 0.023586158703977844, 0.023252175097645322, 0.022944063193710784, 0.022659513885635807, 0.02239652512312884, 0.022153353325937486, 0.021928474002397946};

    
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

        // Azimuth has appropriate smear values
        smear_pct_a = 0.10833078;

        // Create random number generator to smear bins by a small amount in azimuth and zenith
        std::default_random_engine gen_gauss_azi;
        gen_gauss_azi.seed(time(0)*index);
        std::normal_distribution<double> dist_azi(0, smear_pct_a);

        std::default_random_engine gen_gauss_zeni;
        gen_gauss_zeni.seed(time(0)*index+100); // Extra factor 100 to keep the seeds unique
        std::normal_distribution<double> dist_zeni(0, smear_pct_z);
        
        // Get the Gaussian smear values
        double smear_azi  = dist_azi(gen_gauss_azi);
        double smear_zeni = dist_zeni(gen_gauss_zeni);

        // if (zenith.at(ran_idx) == 0)
        //     std::cout << smear_zeni << std::endl;

        double a_ = azimuth.at(ran_idx)+smear_azi;
        double z_ = zenith.at(ran_idx)+smear_zeni;
        
        // Fill histogram
        hist->Fill(a_, z_);

        dir.SetX(sin(z_) * sin(a_));
        dir.SetY(-cos(z_));
        dir.SetZ(-sin(z_) * cos(a_));

        histXY->Fill(dir.X(), dir.Y());
        histXZ->Fill(dir.X(), dir.Z());
        histYZ->Fill(dir.Y(), dir.Z());

        index++;

    }

}



void RandomGenerator(){

    double pi = 3.14159;

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

        // dir *= *rPhi;

        histXY->Fill(dir.X(), dir.Y());
        histXZ->Fill(dir.X(), dir.Z());
        histYZ->Fill(dir.Y(), dir.Z());

    }

    TCanvas *c = new TCanvas();
    hist_root->Draw("colz");

    TCanvas *cXY = new TCanvas();
    histXY->Draw("colz");

    TCanvas *cXZ = new TCanvas();
    histXZ->Draw("colz");

    TCanvas *cYZ = new TCanvas();
    histYZ->Draw("colz");

    TCanvas *c2 = new TCanvas();
    hist->Draw("colz");
    

    // Try C++ implementation
    std::vector<double> weights = {}; // Intensity in a given azimuth/zenith Bin
    std::vector<double> azimuth = {}; // List of Azimuth values
    std::vector<double> zenith = {};  // List of Zenith values
    std::vector<double> azimuth_bins = {}; // List of Azimuth bin edges
    std::vector<double> zenith_bins = {};  // List of Zenith bin edges

    
    // Histogram for testing if random number generator is working
    TH2D* hist_cpp = new TH2D("hist_cpp", ";Azimuth; Zenith", hist->GetNbinsX(), 0, 2.0 ,hist->GetNbinsY(), 0, 0.5);

    // Get Bin Weights, azimuth and zenith values into a vector
    for (unsigned int row = 0; row < hist->GetNbinsX()+1; row++){
        for (unsigned int col = 0; col < hist->GetNbinsY()+1; col++){
            weights.push_back(hist->GetBinContent(row, col));
            azimuth.push_back(hist->GetXaxis()->GetBinCenter(row+1));
            zenith.push_back(hist->GetYaxis()->GetBinCenter(col+1));
        }
    }

    // GenerateRandom(weights, azimuth, zenith, hist_cpp);

    TCanvas *c3 = new TCanvas();
    hist_cpp->Draw("colz");


    // File pointer
    // std::ifstream fin("MeasuredMuonsFromData.csv");
    std::ifstream fin("SimulatedMuonsFromProposal.csv");
    
    // Check if file has opened properly
    if (!fin.is_open())
        std::cout << "Error: Could not read in the input data CSV file" << std::endl;
  
    // Read the Data from the file 
    std::string s_intensity, s_beta, s_alpha;
    std::vector<double> intensity, beta, alpha;
  
    // Loop over the lines in the file and add the values to a vector
    while (fin.peek()!=EOF) {
  
        std::getline(fin, s_intensity, ',');

        // Load in alpha/zenith bin edges
        if (s_intensity == "alpha" || s_intensity == "zenith"){
            std::getline(fin, s_alpha, '\n');
            zenith_bins.push_back(stod(s_alpha));
        }
        // Load in beta/azimuth bin edges
        else if (s_intensity == "beta" || s_intensity == "azimuth"){
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

    GenerateRandom(intensity, beta, alpha, hist_cpp2, histXY_cpp2, histXZ_cpp2, histYZ_cpp2);

    TCanvas *c4 = new TCanvas();
    hist_cpp2->Draw("colz");

    TCanvas *cXY_cpp2 = new TCanvas();
    histXY_cpp2->Draw("colz");

    TCanvas *cXZ_cpp2 = new TCanvas();
    histXZ_cpp2->Draw("colz");

    TCanvas *cYZ_cpp2 = new TCanvas();
    histYZ_cpp2->Draw("colz");


}
