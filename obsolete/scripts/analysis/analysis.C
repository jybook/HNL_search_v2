#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TFile.h"
#include <math.h>
#include "HNLDecayEnergy.h"

double weight_lifetime(double tau_min, double tau_max, double lifetime, double lifetime_proper);
double get_weight(double weight, double mass, double hnlEnergy, double dMin, double dMax, double lifetime, double U_tau4_sq, int nFiles, string mode);
double gamma_pi0_nu(double G_F, double f, double m_pi0, double m_HNL);
double gamma_H_l(double G_F, double V, double f, double m_H, double m_l, double m_HNL);
double gamma_rho_l(double g_p, double G_F, double V, double m_rho, double m_l, double m_HNL);
double gamma_rho0_nu(double g_p, double G_F, double m_rho0, double m_HNL);
double gamma_nu1_l1_l2(double G_F, double m_l1, double m_l2, double m_HNL);
double gamma_nu_ll(double G_F, double sin2thetaW, double m_l, double m_HNL);
double FullWidth(double m_HNL);

void analysis(string file, int nFiles){

	// Open the file containing the tree.
	const char *file_char = file.c_str();
	TFile *myFile = new TFile(file_char,"READ");

	// Create a TTreeReader for the tree
	TTreeReader myReader("MasterTree", myFile);
	TTreeReaderValue<double> weight(myReader, "weight.value");
	TTreeReaderValue<double> hnl_mass(myReader, "hnl_mass.value");
	TTreeReaderValue<double> hnl_energy(myReader, "hnl_energy.value");
	TTreeReaderValue<double> distance_min(myReader, "distance_min.value");
	TTreeReaderValue<double> distance_max(myReader, "distance_max.value");
	TTreeReaderValue<double> lifetime(myReader, "lifetime.value");
	TTreeReaderValue<double> cascade0Taupede_energy(myReader, "cascade0Taupede_energy.value");

	// Declare histograms
	TH1I *h_npulse = new TH1I("h_npulse", ";Number of pulses;Events", 500, 0, 500);
	TH1F *h_cascade0Taupede_energy = new TH1F("h_cascade0Taupede_energy", ";First cascade best-fit energy (GeV);Events", 50, 0, 500);

	// Choose a value of the mixing parameter
	double U_tau4_sq = 1e-03;

	// Loop over all entries of the TTree or TChain.
	while (myReader.Next()) {
		// Compute final weight
		double weight_final = get_weight(*weight,*hnl_mass,*hnl_energy,*distance_min,*distance_max,*lifetime,U_tau4_sq,nFiles,"");
		h_cascade0Taupede_energy->Fill(*cascade0Taupede_energy,weight_final);
	}

	h_cascade0Taupede_energy->Draw();

}

double weight_lifetime(double tau_min, double tau_max, double lifetime, double lifetime_proper){
	double pdf_uniform = 1./(tau_max-tau_min);
	double pdf_inverse = (1./(log(tau_max)-log(tau_min)))*(1./lifetime);
	double pdf_exp1 = 1./lifetime_proper;
	double pdf_exp2 = exp(-lifetime/lifetime_proper);
	double pdf_exp = pdf_exp1*pdf_exp2;
	return pdf_exp/pdf_inverse;
}

double get_weight(double weight, double mass, double hnlEnergy, double dMin, double dMax, double lifetime, double U_tau4_sq, int nFiles, string mode){
	double c = 299792458.;
	double hbar = 6.582119569E-25;
	double gamma = (hnlEnergy+mass)/mass;
	double speed = c*sqrt(1-pow(1./gamma,2));
	double tau_min = dMin/(gamma*speed);
	double tau_max = dMax/(gamma*speed);
	lifetime = 1e-09*lifetime;
	double lifetime_proper = hbar/(FullWidth(mass)*U_tau4_sq);
	double weight_lifetime_value = weight_lifetime(tau_min,tau_max,lifetime,lifetime_proper);
	return (10000.0/nFiles)*U_tau4_sq*weight_lifetime(tau_min,tau_max,lifetime,lifetime_proper)*weight;
	//if mode == 'event': weights = 8*365*24*3600*(1./1000.)*weights
}

double gamma_pi0_nu(double G_F, double f, double m_pi0, double m_HNL){
    if(m_HNL < m_pi0) return 0;
    double mass_ratio = pow(m_pi0/m_HNL,2);
    double factor = pow(G_F,2)*pow(f,2)*pow(m_HNL,3)/(32*M_PI);
    return factor*pow(1-mass_ratio,2);
}

double gamma_H_l(double G_F, double V, double f, double m_H, double m_l, double m_HNL){
    if(m_HNL < (m_H+m_l)) return 0;
    double mass_ratio1 = pow(m_l/m_HNL,2);
    double mass_ratio2 = pow(m_H/m_HNL,2);
    double mass_diff_minus = pow((m_H-m_l)/m_HNL,2);
    double mass_diff_plus = pow((m_H+m_l)/m_HNL,2);
    double factor = pow(G_F*V*f,2)*pow(m_HNL,3)/(16*M_PI);
    return factor*(pow(1-mass_ratio1,2)-(mass_ratio2*(1+mass_ratio1)))*sqrt((1-mass_diff_minus)*(1-mass_diff_minus));
}

double gamma_rho_l(double g_p, double G_F, double V, double m_rho, double m_l, double m_HNL){
    if(m_HNL < (m_rho+m_l)) return 0;
    double mass_ratio1 = pow(m_l/m_HNL,2);
    double mass_ratio2 = pow(m_rho/m_HNL,2);
    double mass_diff_minus = pow((m_rho-m_l)/m_HNL,2);
    double mass_diff_plus = pow((m_rho+m_l)/m_HNL,2);
    double mass_diff_minus_2 = (pow(m_l,2)-2*pow(m_rho,2))/pow(m_HNL,2);
    double factor = pow(g_p/m_rho,2)*pow(G_F,2)*pow(V,2)*pow(m_HNL,3)/(8*M_PI);
    return factor*(pow(1-mass_ratio1,2)-mass_ratio2*(1+mass_diff_minus_2))*sqrt((1-mass_diff_minus)*(1-mass_diff_minus));
}

double gamma_rho0_nu(double g_p, double G_F, double m_rho0, double m_HNL){
    if(m_HNL < m_rho0) return 0;
    double mass_ratio = pow(m_rho0/m_HNL,2);
    double factor = pow(g_p/m_rho0,2)*pow(G_F,2)*pow(m_HNL,3)/(16*M_PI);
    return factor*(1+2*mass_ratio)*pow(1-mass_ratio,2);
}

double gamma_nu1_l1_l2(double G_F, double m_l1, double m_l2, double m_HNL){
    if(m_HNL < (m_l1+m_l2)) return 0;
    double factor = pow(G_F,2)*pow(m_HNL,5)/(192*pow(M_PI,3));
    double mass_ratio = std::max(m_l1,m_l2)/m_HNL;
    double x_factor = 1-(8*pow(mass_ratio,2))+(8*pow(mass_ratio,6))-pow(mass_ratio,8)-(12*pow(mass_ratio,4)*log(pow(mass_ratio,2)));
    return factor*x_factor;
}

double gamma_nu_ll(double G_F, double sin2thetaW, double m_l, double m_HNL){
    if(m_HNL < 2*m_l) return 0;
    double x = pow(m_l/m_HNL,2);
    double L = 0;
    double num = 0, den = 0;
    if(x<=0.25){
        num = 1-(3*x)-(1-x)*sqrt(1-4*x);
        den = x*(1+sqrt(1-4*x));
    }
    if(den!=0 && num/den>0) L = log(num/den);
    double C1 = 0.25*(1-4*sin2thetaW+8*pow(sin2thetaW,2));
    double C2 = 0.5*sin2thetaW*(2*sin2thetaW-1);
    double C3 = 0.25*(1+4*sin2thetaW+8*pow(sin2thetaW,2));
    double C4 = 0.5*sin2thetaW*(2*sin2thetaW+1);
    double factor = pow(G_F,2)*pow(m_HNL,5)/(192*pow(M_PI,3));
    double gamma_same_flavor = (C3*((1-14*x-2*pow(x,2)-12*pow(x,3))*sqrt(1-4*x)+12*pow(x,2)*(pow(x,2)-1)*L)+(4*C4*(x*(2+10*x-12*pow(x,2))*sqrt(1-4*x)+6*pow(x,2)*(1-2*x+2*pow(x,2))*L)));
    double gamma_opposite_flavor = (C1*((1-14*x-2*pow(x,2)-12*pow(x,3))*sqrt(1-4*x)+12*pow(x,2)*(pow(x,2)-1)*L)+(4*C2*(x*(2+10*x-12*pow(x,2))*sqrt(1-4*x)+6*pow(x,2)*(1-2*x+2*pow(x,2))*L)));
    return factor*gamma_opposite_flavor;
}

double FullWidth(double m_HNL){
    // Declare global physical constants
    double G_F = 1.17E-05; // Fermi constant in GeV^-2
    double sin2thetaW = 0.23122; // sin^2(weak angle)

    return gamma_nu_ll(G_F,sin2thetaW,511E-06,m_HNL)+          // nu_ee
           gamma_pi0_nu(G_F,0.130,0.135,m_HNL)+                // pi0_nu
           gamma_nu_ll(G_F,sin2thetaW,0.1056,m_HNL)+           // nu_mumu
           gamma_pi0_nu(G_F,0.156,0.548,m_HNL)+                // eta_nu
           gamma_rho0_nu(0.102,G_F,0.775,m_HNL)+               // rho0_nu
           gamma_pi0_nu(G_F,-0.0585,0.958,m_HNL)+              // etaprime_nu
           gamma_nu1_l1_l2(G_F,511E-06,1.776,m_HNL)+           // nue_e_tau
           gamma_nu1_l1_l2(G_F,0.1056,1.776,m_HNL)+            // numu_mu_tau
           gamma_H_l(G_F,0.97427,0.130,0.140,1.776,m_HNL)+     // tau_pi
           gamma_H_l(G_F,0.224,0.130,0.494,1.776,m_HNL)+       // tau_K
           gamma_rho_l(0.102,G_F,9.74E-01,0.775,1.776,m_HNL); // tau_rho
}
