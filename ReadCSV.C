#include<iostream>
#include<vector>
#include<fstream>
#include<string>

// Function to read in a CSV file the C++ way
void ReadCSV(){

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

        std::cout << s_intensity << " " << s_beta << " "<< s_alpha << std::endl;
    }

}