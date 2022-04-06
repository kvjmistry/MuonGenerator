#include "TH2.h"
#include "TArray.h"
#include <iostream>

double GetBinSmearValue(std::vector<double> bins, double sampVal, std::vector<double> BinWidths);


// Convert root histogram to csv file
void Hist_to_CSV(){

    double pi = 3.14159;

    std::string mode = "data";
    // std::string mode = "sim";

    double corr_factor = pi; // choose 1 or pi
    int y_cut_off = 1;
    TFile *f;

    if (mode == "data"){
        f =  (TFile*)TFile::Open("MuonAnaAllRuns.root");
    }
    else {
        f =  (TFile*)TFile::Open("SimulatedMuonFile.root");
        corr_factor = 1.0;
        y_cut_off = 1;
    }

    // Load in the histogram
    // TFile f("MuonAnaAllRuns.root");
    
    TH2F *hist;
    f->GetObject("za", hist);
    hist->SetDirectory(0);
    f->Close();

    std::vector<double> weights = {}; // Intensity in a given azimuth/zenith Bin
    std::vector<double> azimuth = {}; // List of Azimuth values
    std::vector<double> zenith = {};  // List of Zenith values
    std::vector<double> azimuth_bins = {}; // List of Azimuth bin edges
    std::vector<double> zenith_bins = {};  // List of Zenith bin edges

    std::vector<double> azimuth_BW = {}; // List of Azimuth bin widths
    std::vector<double> zenith_BW  = {};  // List of Zenith bin widths

    std::vector<double> azimuth_smear = {}; // List of Azimuth bin smear values
    std::vector<double> zenith_smear  = {};  // List of Zenith bin smear values

    

    TArrayD* Ax_bins = (TArrayD*) hist->GetXaxis()->GetXbins();
    double xbins[hist->GetNbinsX()];

    for (unsigned int i =0; i < Ax_bins->GetSize(); i++){
        xbins[i] = Ax_bins->GetAt(i);
    }

    TArrayD* Ay_bins = (TArrayD*) hist->GetYaxis()->GetXbins();
    double ybins[hist->GetNbinsX()];

    for (unsigned int i =0; i < Ay_bins->GetSize(); i++){
        ybins[i] = Ay_bins->GetAt(i);
    }
    
    // Histogram for testing if random number generator is working
    TH2D* hist_cpp;

    if (mode == "data"){
        hist_cpp = new TH2D("hist_cpp", ";Azimuth; Zenith", hist->GetNbinsX(), 0, 2.0*pi ,hist->GetNbinsY(), 0, 0.5*pi);
    }
    else {
        hist_cpp = new TH2D("hist_cpp", ";Azimuth; Zenith", hist->GetNbinsX(), xbins ,hist->GetNbinsY(), ybins);
    }
    
     
    // Get Bin Weights, azimuth and zenith values into a vector
    for (unsigned int row = 1; row < hist->GetNbinsX()+1; row++){
        for (unsigned int col = y_cut_off; col < hist->GetNbinsY()+1; col++){
            weights.push_back(hist->GetBinContent(row, col) / (hist->GetXaxis()->GetBinWidth(row) * hist->GetYaxis()->GetBinWidth(col) ));
            azimuth.push_back(hist->GetXaxis()->GetBinCenter(row)*corr_factor);
            zenith.push_back(hist->GetYaxis()->GetBinCenter(col)*corr_factor);

            if (row == 1){
                zenith_bins.push_back(hist->GetYaxis()->GetBinLowEdge(col)*corr_factor);
                // std::cout << hist->GetYaxis()->GetBinLowEdge(col+1)<< std::endl;

                // Add the final bin edge
                if (col == hist->GetNbinsY())
                    zenith_bins.push_back(hist->GetYaxis()->GetBinLowEdge(col+1)*corr_factor);
                
            }

        }

        azimuth_bins.push_back(hist->GetXaxis()->GetBinLowEdge(row)*corr_factor);

        // Add the final bin edge
        if (row == hist->GetNbinsX())
            azimuth_bins.push_back(hist->GetXaxis()->GetBinLowEdge(row+1)*corr_factor);

    }

    // Get the bin widths to smear by
    for (unsigned int bin = 1; bin < hist->GetNbinsX()+1; bin++){
        azimuth_BW.push_back(hist->GetXaxis()->GetBinWidth(bin));
    }

    for (unsigned int bin = y_cut_off; bin < hist->GetNbinsY()+1; bin++){
        zenith_BW.push_back(hist->GetYaxis()->GetBinWidth(bin));
    }

    // Azimuth smear values
    for (int i = 0; i < azimuth.size(); i++){
       azimuth_smear.push_back(GetBinSmearValue(azimuth_bins, azimuth.at(i), azimuth_BW));
    }

    // Zenith smear values
    for (int i = 0; i < zenith.size(); i++){
       zenith_smear.push_back(GetBinSmearValue(zenith_bins, zenith.at(i), zenith_BW));
    }

    std::ofstream myfile;

    if (mode == "data"){
         myfile.open ("MuonAnaAllRuns.csv");
    }
    else {
        myfile.open ("SimulatedMuonsProposalMCEq.csv");
    }

   

    for (int i = 0; i < weights.size(); i++){
        myfile << "value," << weights.at(i) << "," << azimuth.at(i) << ","<< zenith.at(i) << "," <<azimuth_smear.at(i) << "," << zenith_smear.at(i) <<  "\n";
    }

    for (int i = 0; i < zenith_bins.size(); i++){
        myfile << "zenith" << "," << zenith_bins.at(i) << "\n";
    }

    for (int i = 0; i < azimuth_bins.size(); i++){
        myfile << "azimuth" << "," << azimuth_bins.at(i) << "\n";
    }

    myfile.close();

    for (unsigned int i = 0; i < weights.size(); i++){
        hist_cpp->Fill(azimuth.at(i), zenith.at(i), weights.at(i));
    }

    gStyle->SetOptStat(0);
    TCanvas *c = new TCanvas();
    // hist_cpp->Divide(hist);
    hist_cpp->Draw("COLZ");


}


double GetBinSmearValue(std::vector<double> bins, double sampVal, std::vector<double> BinWidths){

  double smearVal{std::numeric_limits<double>::lowest()};

  // Loop over the azimuth values and find the corresponding bin width to smear
  for (int i = 0; i < bins.size()-1; i++){
      
      // Include last bin edge in check
      if (smearVal == std::numeric_limits<double>::lowest() && i == bins.size()-2){
        if (sampVal >= bins[i] && sampVal <= bins[i+1]){
          
            smearVal = BinWidths[i];
            break;

        }
      }
      else {

        if (sampVal >= bins[i] && sampVal < bins[i+1]){
            
            smearVal = BinWidths[i];
            break;

        }
      }

  }

  return smearVal;

}