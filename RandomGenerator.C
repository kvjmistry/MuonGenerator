// Test code to generate random Numbers 
#include <random>
#include <vector>
#include "TH2.h"
#include "TRotation.h"
#include "TVector3.h"

void GenerateRandom(std::vector<double> weights, std::vector<double> azimuth, std::vector<double> zenith, TH2D* hist, TH2D* histXY, TH2D* histXZ, TH2D* histYZ){


    // Discrete distribution with bin index and intensity to sample from
    std::discrete_distribution<int> dist(std::begin(weights), std::end(weights));
    std::default_random_engine gen;  // Generator
    gen.seed(time(0)); // Seed with the current time for uniqueness
    int N_samp = 10000000;
    std::vector<int> samples(N_samp);
    
    double smear_pct{0.02}; // Percentage to smear the azimuth and zenith values by

    // Create random number generator to smear bins by a small amount in azimuth and zenith
    std::default_random_engine gen_gauss_azi;
    gen_gauss_azi.seed(time(0));
    std::normal_distribution<double> dist_azi(1.0,0.02);

    std::default_random_engine gen_gauss_zeni;
    gen_gauss_zeni.seed(time(0));
    std::normal_distribution<double> dist_zeni(1.0,0.02);

    TVector3 dir;

    for (auto & ran_idx: samples){
        
        // Generate random index weighted by the bin contents
        ran_idx = dist(gen);
        
        // Get the Gaussian smear values
        double smear_azi  = dist_azi(gen_gauss_azi);
        double smear_zeni = dist_zeni(gen_gauss_zeni);
        
        // Fill histogram
        hist->Fill(azimuth.at(ran_idx)*smear_azi, zenith.at(ran_idx)*smear_zeni);

        double a_ = azimuth.at(ran_idx)*smear_azi;
        double z_ = zenith.at(ran_idx)*smear_zeni;

        dir.SetX(sin(z_) * sin(a_));
        dir.SetY(-cos(z_));
        dir.SetZ(-sin(z_) * cos(a_));

        histXY->Fill(dir.X(), dir.Y());
        histXZ->Fill(dir.X(), dir.Z());
        histYZ->Fill(dir.Y(), dir.Z());


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
    rPhi->RotateY(-150);
    
    TVector3 dir;

    for (int i = 0; i < 1000000; i++){
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
    std::vector<double> azimuth = {}; // List of Azimuth bin centers
    std::vector<double> zenith = {};  // List of Zenith bin centers
    
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


    // Histogram for testing if random number generator is working
    TH2D* hist_cpp2   = new TH2D("hist_cpp2", ";Azimuth; Zenith", 50, 0, 2.0*3.141 ,50, 0, 0.5*3.141);
    TH2D* histXY_cpp2 = new TH2D("histXY_cpp2", ";X; Y", 50, -1, 1 , 25, -1, 0 );
    TH2D* histXZ_cpp2 = new TH2D("histXZ_cpp2", ";X; Z", 50, -1, 1 , 50, -1, 1 );
    TH2D* histYZ_cpp2 = new TH2D("histYZ_cpp2", ";Y; Z", 25, -1, 0 , 50, -1, 1 );

    // File pointer
    std::ifstream fin("MeasuredMuonsFromData.csv");
    
    // Check if file has opened properly
    if (!fin.is_open())
        std::cout << "Error: Could not read in the input data CSV file" << std::endl;
  
    // Read the Data from the file 
    std::string s_intensity, s_beta, s_alpha;
    std::vector<double> intensity, beta, alpha;
  
    // Loop over the lines in the file and add the values to a vector
    while (fin.peek()!=EOF) {
  
        std::getline(fin, s_intensity, ',');
        std::getline(fin, s_beta, ',');
        std::getline(fin, s_alpha, '\n');

        intensity.push_back(stod(s_intensity));
        beta.push_back(stod(s_beta));
        alpha.push_back(stod(s_alpha));

    }

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
