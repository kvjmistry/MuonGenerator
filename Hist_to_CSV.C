#include "TH2.h"
#include <iostream>


// Convert root histogram to csv file
void Hist_to_CSV(){

    double pi = 3.14159;

    // Load in the histogram
    TFile f("MuonAnaAllRuns.root");
    TH2F *hist;
    f.GetObject("za", hist);
    hist->SetDirectory(0);
    f.Close();

    std::vector<double> weights = {}; // Intensity in a given azimuth/zenith Bin
    std::vector<double> azimuth = {}; // List of Azimuth values
    std::vector<double> zenith = {};  // List of Zenith values
    std::vector<double> azimuth_bins = {}; // List of Azimuth bin edges
    std::vector<double> zenith_bins = {};  // List of Zenith bin edges

    
    // Histogram for testing if random number generator is working
    TH2D* hist_cpp = new TH2D("hist_cpp", ";Azimuth; Zenith", hist->GetNbinsX(), 0, 2.0 ,hist->GetNbinsY(), 0, 0.5);

    // Get Bin Weights, azimuth and zenith values into a vector
    for (unsigned int row = 1; row < hist->GetNbinsX()+1; row++){
        for (unsigned int col = 1; col < hist->GetNbinsY()+1; col++){
            weights.push_back(hist->GetBinContent(row, col));
            azimuth.push_back(hist->GetXaxis()->GetBinCenter(row+1)*pi);
            zenith.push_back(hist->GetYaxis()->GetBinCenter(col+1)*pi);

            if (row == 1){
                zenith_bins.push_back(hist->GetYaxis()->GetBinLowEdge(col)*pi);

                // Add the final bin edge
                if (col == hist->GetNbinsY())
                    zenith_bins.push_back(hist->GetYaxis()->GetBinLowEdge(col+1)*pi);
                
            }

        }

        azimuth_bins.push_back(hist->GetXaxis()->GetBinLowEdge(row)*pi);

        // Add the final bin edge
        if (row == hist->GetNbinsX())
            azimuth_bins.push_back(hist->GetXaxis()->GetBinLowEdge(row+1)*pi);

    }

    std::ofstream myfile;
    myfile.open ("MuonAnaAllRuns.csv");

    for (int i = 0; i < weights.size(); i++){
        myfile << weights.at(i) << "," << azimuth.at(i) << ","<< zenith.at(i)<< "\n";
    }

    for (int i = 0; i < zenith_bins.size(); i++){
        myfile << "zenith" << "," << zenith_bins.at(i) << "\n";
    }

    for (int i = 0; i < azimuth_bins.size(); i++){
        myfile << "azimuth" << "," << azimuth_bins.at(i) << "\n";
    }

    myfile.close();

}