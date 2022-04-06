void CompileMacro(){

  // gSystem->Load("FluggNtuple/FluxNtuple_C.so");
  gSystem->AddIncludePath(" -I/opt/homebrew/Cellar/clhep/2.4.5.1/include");
  TString t = ".L RandomGenerator.C+";
  gROOT->ProcessLine(t);
  gSystem->Exit(0);

}