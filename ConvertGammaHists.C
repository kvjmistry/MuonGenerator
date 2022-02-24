
#include "TFile.h"
#include "TH1D.h"
#include <iostream>

// Convert Gamma histograms to CSV files
void ConvertGammaHists(){

    // Strings for the file names
    TSystemDirectory dir("hists", "hists");
    
    TList *files = dir.GetListOfFiles();
    
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files); 
        
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            
            if (!file->IsDirectory() && fname.EndsWith(".root")) {
                std::cout << fname.Data() << std::endl;

                TFile *f = TFile::Open(Form("hists/%s", fname.Data())); 
            
                // Get the histogram from the file
                TH1D *h = (TH1D*) f->Get("GammaEnergy");

                std::ofstream myfile;
                    
                    fname.ReplaceAll(".root", "");
                    myfile.open (Form("csv/%s.csv", fname.Data()));

                // Loop over the bins
                for (unsigned int bin = 1; bin < h->GetNbinsX()+1;bin++){
                    
                    if (bin == 1)
                        myfile << "E [keV],Entries" << "\n";

                    myfile << h->GetBinCenter(bin) << "," << h->GetBinContent(bin) << "\n";

                }

                myfile.close();

                f->Close();
            
            }
        }
    }

    // Load in a csv file and take the ratio to see if they are all 1.0
    TFile  *f = TFile::Open("hists/GammaEnergy_Bi214.root");
    TH1D* h = (TH1D*) f->Get("GammaEnergy");

    std::ifstream fin("csv/GammaEnergy_Bi214.csv");
    
    // Check if file has opened properly
    if (!fin.is_open())
        std::cout << "Error: Could not read in the input data CSV file" << std::endl;
  
    // Read the Data from the file 
    std::vector<double> E, events;
    std::string s_E, s_events;

    TH1D* h_csv = new TH1D("GammaEnergy_", "GammaEnergy", 10000, 0, 4000.);
  
    // Loop over the lines in the file and add the values to a vector
    while (fin.peek()!=EOF) {
  
        std::getline(fin, s_E, '\n');

        // Load in alpha/zenith bin edges
        if (s_E == "E [keV],Entries"){
            continue;
        }
        
        std::getline(fin, s_E, ',');
        std::getline(fin, s_events, '\n');

        E.push_back(stod(s_E));
        events.push_back(stod(s_events));

        h_csv->Fill(stod(s_E), stod(s_events));
    }

    h_csv->Divide(h);

    TCanvas *c = new TCanvas();

    h_csv->SetMaximum(2);

    h_csv->Draw("hist");


}