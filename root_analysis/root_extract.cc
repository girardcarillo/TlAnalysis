#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TTree.h>
#include <string>
#include <cmath>
#include "/home/girardcarillo/Workdir/SNPlot/SNPlot.cc"
#include "root_extract.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <TList.h>
#include <TFile.h>
#include <TRandom.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TLegend.h>

Double_t A_208Tl = 54*31.5;//2*31.5;
Double_t exposure = 17.5;
Double_t nb_simulated_events = 0;
Double_t real_efficiency = 0.00312;

using namespace std;

Double_t GetBeta(Double_t energy);
Double_t GetTheoreticalTime(Double_t length, Double_t beta);
Double_t GetRandomSigma(Double_t random, Double_t energy);
Double_t GetSigma_L(Double_t theoretical_time, Double_t energy);
Double_t GetInternal_Chi2(Double_t delta_t, Double_t sigma);
Double_t GetInternalProbability(Double_t internal_Chi2);
void config_histo1D(TH1F *histo, const char *histoTitle, const char *XaxisTitle, const char *YaxisTitle,int linewidth,Color_t color,const char * draw_type);
void config_histo2D(TH2F *histo, const char *histoTitle, const char * XaxisTitle, const char *YaxisTitle,int linewidth,Color_t color,const char * draw_type);

void root_extract(Double_t gaus_sigma, Double_t sigma_L, string isotope,int IC=0){

  string const final_rate("final_rate.txt");
  ofstream final_flux(final_rate.c_str());

  Double_t tau = 0.294; //ns

  Double_t slope = 1/tau;
  string file_name;
  if (isotope == "208Tl") {
    if (IC==100) {
      nb_simulated_events=1.e5;
      file_name = "$WORKDIR/Tl_analysis/208Tl/root_file_simus_IC_100/208Tl-";
    }
    else if (IC==0) {
      nb_simulated_events=1.e6;
      file_name = "$WORKDIR/Tl_analysis/208Tl/root_file_simus_IC_0/208Tl_IC0-";
    }
  }
  else if (isotope == "0nubb") {
    nb_simulated_events=1.e6;
    file_name = "$WORKDIR/Tl_analysis/0nubb/0nubb_0ns_";
  }

  string file_extension = ".root";
  string file_name_new;
  TList *list = new TList();
  string file_adress;
  string tree_adress;

  ///Loop on all .root output files
  for (int i = 0; i < 10; ++i) {
    stringstream ss;
    ss << i;
    string str = ss.str();
    file_name_new = file_name+str+file_extension;
    file_adress = string("f")+str;
    tree_adress = string("tree")+str;
    const char* c_file_name_new = file_name_new.c_str();
    TFile *file_adress = TFile::Open(c_file_name_new,"READ");
    TTree *tree_adress = (TTree*)file_adress->Get("calorimeter_hit");
    list->Add(tree_adress);
  }

  ///Merge all trees in "newtree"
  // TFile* outputfile = TFile::Open("output.root", "recreate");
  TTree *newtree = TTree::MergeTrees(list);
  // newtree->Write();
  // outputfile->Close();
  // delete outputfile;

  Double_t time_difference_E;
  Double_t time;
  Double_t time_Emin;
  Double_t time_Emax;
  Double_t probability;
  Double_t length_Emin;// = 1049.02;
  Double_t length_Emax;// = 905.924;
  Double_t minimal_energy;// = 0.62446;
  Double_t maximal_energy;// = 2.45767;
  Double_t sigma_time_Emin;// = 4.29912e-05;
  Double_t sigma_time_Emax;// = 2.16706e-05;
  Int_t event_counter;

  newtree->SetBranchAddress("time_difference_E",&time_difference_E);
  newtree->SetBranchAddress("time_Emin",&time_Emin);
  newtree->SetBranchAddress("time_Emax",&time_Emax);
  newtree->SetBranchAddress("probability",&probability);
  newtree->SetBranchAddress("length_Emin",&length_Emin);
  newtree->SetBranchAddress("length_Emax",&length_Emax);
  newtree->SetBranchAddress("minimal_energy",&minimal_energy);
  newtree->SetBranchAddress("maximal_energy",&maximal_energy);
  newtree->SetBranchAddress("sigma_time_Emin",&sigma_time_Emin);
  newtree->SetBranchAddress("sigma_time_Emax",&sigma_time_Emax);
  newtree->SetBranchAddress("event_counter",&event_counter);

  ///Create histograms to be filled with reconstructed data outputs
  Long64_t nentries = newtree->GetEntries();
  Double_t nbins = 0;
  if (isotope=="0nubb") {
    nbins = 100;
  }
  else if (isotope=="208Tl") {
    nbins = 30;
  }

  TH1F *hsigma_time_Emax = new TH1F("","",nbins,-2,2);
  TH1F *hsigma_time_Emin = new TH1F("","",nbins,-2,2);
  TH1F *hnoised_time_difference = new TH1F("smired_calo_times","",nbins,-6,6);
  TH1F *htime_Emin = new TH1F("time_Emin","Time of e-(Emin)",nbins,1,6);
  TH1F *htime_Emax = new TH1F("time_Emax","Time of e-(Emax)",nbins,1,6);
  TH1F *hEmin = new TH1F("minimal_energy","Minimal energy",nbins,0,3);
  TH1F *hEmax = new TH1F("maximal_energy","Maximal energy",nbins,0,3);
  TH1F *hinternal_probability = new TH1F("internal_probability","",nbins,0,1);
  TH1F *hintegrated_probability = new TH1F("integrated_probability","Integrated probability",nbins,0,1);
  TH1F *hcalculated_internal_probability = new TH1F("calculated_internal_probability","",nbins,0,1);
  TH1F *hinternal_Chi2 = new TH1F("Internal_Chi2","Calculated internal Chi2",nbins,0,1000);
  TH2F *h2Probability = new TH2F("Biplot","P_{exp}(P_{int})",nbins,0,1,nbins,0,1);
  TH1F *hsigma_L_Emax = new TH1F("","",nbins,0,0.01);
  TH1F *hsigma_L_Emin = new TH1F("","",nbins,0,0.1);
  TH1F *hsigma_tot = new TH1F("","",nbins,0.02,0.1);

  ///Fill trees with reconstructed data outputs
  TRandom *grandom1 = new TRandom(1234);
  TRandom *grandom2 = new TRandom(1234);

  Double_t noised_time_Emin = 0;
  Double_t noised_time_Emax = 0;
  Double_t random_Emax = 0;
  Double_t random_Emin = 0;
  Double_t theoretical_time_Emax = 0;
  Double_t theoretical_time_Emin = 0;
  Double_t noised_time_difference = 0;
  Double_t sigma_L_Emax = 0;
  Double_t sigma_L_Emin = 0;
  // Double_t tau_sc = 2.500;//ns
  // Double_t FWHM_TTS = 2.250;//ns
  Double_t sigma_tot = 0;
  Double_t internal_Chi2 = 0;
  int counter_test = 0;
  Double_t calculated_internal_probability = 0;
  Double_t f_integrated_probability = 0;
  double maximumX = 0;
  double evalY = 0;
  double integral1 = 0;
  double integral2 = 0;
  double integral = 0;
  double mean_lower_bound = 0;
  double mean_upper_bound = 0;
  double ecart_moyenne = 0;

  TF1 *f_expogauss = new TF1("f_expogauss","([0]/2)*exp(([0]/2)*(2*[2]+[0]*pow([1],2)-2*(x)))*erfc(([2]+[0]*pow([1],2)-(x))/(sqrt(2)*[1]))",-6,6);
  TF1 *f_pint = new TF1("pint","erfc(x)",-6,6);

  double nb_tot_events = 0;
  double nb_cut_events = 0;

  cout << "total #event = " << newtree->GetEntries() << endl;
  for (Long64_t i=0;i<newtree->GetEntries();i++) {
    newtree->GetEntry(i);

    theoretical_time_Emax = GetTheoreticalTime(length_Emax,GetBeta(maximal_energy));
    theoretical_time_Emin = GetTheoreticalTime(length_Emin,GetBeta(minimal_energy));

    random_Emax = GetRandomSigma(grandom1->Gaus(0,gaus_sigma),maximal_energy);
    random_Emin = GetRandomSigma(grandom2->Gaus(0,gaus_sigma),minimal_energy);

    noised_time_Emin = time_Emin + random_Emin;
    noised_time_Emax = time_Emax + random_Emax;

    noised_time_difference = (noised_time_Emax - theoretical_time_Emax) - (noised_time_Emin - theoretical_time_Emin);

    sigma_L_Emax = GetSigma_L(theoretical_time_Emax,maximal_energy);
    sigma_L_Emin = GetSigma_L(theoretical_time_Emin,minimal_energy);

    sigma_tot = sqrt(pow(random_Emax,2)+pow(random_Emin,2)+pow(sigma_L_Emax,2)+pow(sigma_L_Emin,2)+pow(sigma_L,2));
    internal_Chi2 = GetInternal_Chi2(noised_time_difference,sigma_tot);
    calculated_internal_probability = GetInternalProbability(internal_Chi2);

    f_expogauss->SetParameters(slope,sigma_tot,0);

    maximumX = f_expogauss->GetMaximumX(-6,6);
    evalY = f_expogauss->Eval(noised_time_difference);

    if (noised_time_difference <= maximumX) {
      integral1 = f_expogauss->Integral(-6,noised_time_difference);
      integral2 = f_expogauss->Integral(f_expogauss->GetX(evalY,maximumX,6),6);
    }
    else {
      integral1 = f_expogauss->Integral(noised_time_difference,6);
      integral2 = f_expogauss->Integral(-6,f_expogauss->GetX(evalY,-6,maximumX));
    }

    integral = integral1+integral2;

    //if (calculated_internal_probability>integral) {

    // if (noised_time_difference<0) {

    hintegrated_probability->Fill(integral);
    hcalculated_internal_probability->Fill(calculated_internal_probability);
    hnoised_time_difference->Fill(noised_time_difference);
    h2Probability->Fill(calculated_internal_probability,integral);
    htime_Emin->Fill(time_Emin);
    htime_Emax->Fill(time_Emax);
    hEmin->Fill(minimal_energy);
    hEmax->Fill(maximal_energy);
    hinternal_probability->Fill(probability);
    hinternal_Chi2->Fill(internal_Chi2);
    hsigma_time_Emin->Fill(sigma_time_Emin);
    hsigma_time_Emax->Fill(sigma_time_Emax);
    hsigma_L_Emin->Fill(sigma_L_Emin);
    hsigma_L_Emax->Fill(sigma_L_Emax);
    hsigma_tot->Fill(sigma_tot);

    nb_cut_events++;
    // }
    //}

    final_flux << event_counter << endl;
    nb_tot_events++;


  }

  cout << f_expogauss->GetMaximumX(-6,6) << endl;

  double efficiency = (nb_cut_events/1e5)*real_efficiency;

  double nbr_bdf = A_208Tl*efficiency*exposure ;

  cout << "#total events = " << nb_tot_events << endl;
  cout << "#cut events = " << nb_cut_events << endl;

  cout << "Efficiency of Pint/Pexp cut = " << nb_cut_events/nb_tot_events << endl;

  cout << "Nback = " << nbr_bdf << endl;

  final_flux.close();

  ///Normalization of probability histograms
  Double_t scale_calo = 1./hnoised_time_difference->Integral();
  hnoised_time_difference->Scale(scale_calo);

  ///Drawing

  TCanvas *c_all = new TCanvas("canvas","canvas");
  TGaxis::SetMaxDigits(3);
  c_all->Divide(2,2);

  SNPlot a_plot;
  a_plot.Set_SNgstyle();

  c_all->cd(1);
  config_histo1D(hcalculated_internal_probability," ","probability","#Event",3,1,"HIST");
  config_histo1D(hintegrated_probability," ","probability","#Event",3,2,"SAME");

  auto legend1 = new TLegend(0.2,0.7,0.55,0.9);
  legend1->AddEntry(hcalculated_internal_probability,"P_{int}","l");
  legend1->AddEntry(hintegrated_probability,"P_{exp}","l");
  legend1->Draw();

  c_all->cd(2);
  config_histo1D(hinternal_Chi2,"Internal #chi^{2}","#chi^{2}_{int}","#Event",3,1,"HIST");

  c_all->cd(3);
  config_histo1D(hnoised_time_difference,"#Delta t distribution","#Delta t (ns)","#Event",3,1,"HIST");

  c_all->cd(4);
  config_histo2D(h2Probability,"P_{int} vs P_{exp}","P_{int}","P_{exp}",3,1,"colz");

  c_all->Modified();
  c_all->SaveAs("../plots/all.pdf");

  TCanvas *c_biplot = new TCanvas("canvas","canvas");
  config_histo2D(h2Probability,"^{208}Tl: P_{int} vs P_{exp} with #sigma_{t}=400 ps @1 MeV","P_{int}","P_{exp}",3,1,"colz");
  c_biplot->SaveAs("../plots/PintVSPexp.pdf");
  //f_expogauss->Draw();


}


Double_t GetBeta(Double_t energy){
  Double_t beta;
  Double_t electron_mass = 0.511;
  beta = sqrt(energy * (energy + 2.*electron_mass)) / (energy + electron_mass);
  return beta;
}

Double_t GetTheoreticalTime(Double_t length, Double_t beta){
  Double_t theoretical_time;
  Double_t light_velocity = 2.99792458e+2;
  theoretical_time = length/(beta*light_velocity);
  return theoretical_time;
}

Double_t GetRandomSigma(Double_t random, Double_t energy){
  Double_t random_sigma;
  random_sigma = random/sqrt(energy);
  return random_sigma;
}

Double_t GetSigma_L(Double_t theoretical_time, Double_t energy){
  Double_t sigma_L;
  Double_t electron_mass = 0.511;
  sigma_L = (theoretical_time/(electron_mass*pow(energy/electron_mass+1,3)*pow(GetBeta(energy),3./2.)))*((0.08*sqrt(energy))/(2*sqrt(2*log(2))));
  return sigma_L;
}

Double_t GetInternal_Chi2(Double_t delta_t, Double_t sigma){
  Double_t valeur;
  valeur = pow(delta_t,2)/pow(sigma,2);
  return valeur;
}

Double_t GetInternalProbability(Double_t internal_Chi2){
  Double_t valeur;
  valeur = erfc(sqrt(internal_Chi2)/sqrt(2));
  return valeur;
}

void config_histo1D(TH1F *histo, const char *histoTitle, const char * XaxisTitle, const char *YaxisTitle,int linewidth,Color_t color,const char * draw_type)
{
  histo->Draw(draw_type);
  histo->GetYaxis()->SetTitleSize(0.048);
  histo->GetXaxis()->SetTitleSize(0.048);
  histo->SetLineWidth(linewidth);
  histo->SetLineColor(color);
  histo->SetTitle(histoTitle);
  histo->GetXaxis()->SetTitle(XaxisTitle);
  histo->GetYaxis()->SetTitle(YaxisTitle);
}


void config_histo2D(TH2F *histo, const char *histoTitle, const char * XaxisTitle, const char *YaxisTitle,int linewidth,Color_t color,const char * draw_type)
{
  histo->Draw(draw_type);
  histo->GetYaxis()->SetTitleSize(0.048);
  histo->GetXaxis()->SetTitleSize(0.048);
  histo->SetLineWidth(linewidth);
  histo->SetLineColor(color);
  histo->SetTitle(histoTitle);
  histo->GetXaxis()->SetTitle(XaxisTitle);
  histo->GetYaxis()->SetTitle(YaxisTitle);
}
