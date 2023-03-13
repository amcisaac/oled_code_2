#include <iostream>
#include <cmath>
#include <vector>
#include "tclap/CmdLine.h"


using namespace std;
using namespace TCLAP;

enum Species {Vacancy, Charge, Exciton};
enum Rates {kr, kEEA, k_EEA, kECA, k_ECA, kRec, k_Rec, ki, k_i, K_Rec, K_EEA, K_ECA, K_i, kDC, k_DC, kDE,  k_DE, NUMBER_OF_RATES};


void getConc(vector<double> &OnSiteAvg, vector<double> BathAvg, vector<double> BathRatesAvg, vector<double> rates ) {
    
    int NN = 4;
    
    // calculate D*x_i from rate equations
    OnSiteAvg[Vacancy] = NN * rates[kRec] * rates[kECA] * BathAvg[Charge] * BathAvg[Charge] + (BathRatesAvg[Charge] + rates[kDC] * BathAvg[Vacancy]) * (rates[k_i] + rates[kr] + NN * rates[kEEA] * BathAvg[Exciton] + NN * (rates[k_Rec] + rates[kDE]) * BathAvg[Vacancy]) + (rates[k_i] * rates[kRec] + rates[kr] * rates[kRec] + NN * (rates[kEEA] * rates[kRec] * BathAvg[Exciton] + rates[kECA] * BathRatesAvg[Charge] + rates[kDC] * rates[kECA] * BathAvg[Vacancy] + rates[kDE] * rates[kRec] * BathAvg[Vacancy])) * BathAvg[Charge];

//    OnSiteAvg[Vacancy] = (rates[k_i] + NN * BathRatesAvg[Charge] + NN * rates[kDC] * BathAvg[Vacancy]) * (rates[kr] + NN * rates[k_Rec] * BathAvg[Vacancy] + NN * rates[kEEA] * BathAvg[Exciton] + NN * rates[kECA] * BathAvg[Charge] + NN * rates[kDE] * BathAvg[Vacancy]) + NN * rates[kRec] * BathAvg[Charge] * (rates[kr] + NN * rates[kEEA] * BathAvg[Exciton] + NN * rates[kECA] * BathAvg[Charge] + NN * rates[kDE] * BathAvg[Vacancy]);

    OnSiteAvg[Charge] = (rates[k_DC] * BathAvg[Charge] + BathRatesAvg[Exciton]) * (rates[k_i] + rates[kr] + NN * rates[kECA] * BathAvg[Charge] + NN * rates[kEEA] * BathAvg[Exciton]) + ((rates[k_DC] * BathAvg[Charge] + BathRatesAvg[Exciton]) * rates[kDE] * NN + rates[k_Rec] * (rates[ki] + NN * ((rates[k_DC] + rates[k_ECA]) * BathAvg[Charge] + (rates[k_DE] + rates[k_EEA]) * BathAvg[Exciton] + BathRatesAvg[Exciton]))) * BathAvg[Vacancy];
    
//    OnSiteAvg[Charge] = (rates[ki] + NN * BathRatesAvg[Exciton] + NN * rates[k_DC] * BathAvg[Charge]) *(rates[kr] + NN * rates[k_Rec] * BathAvg[Vacancy] + NN * rates[kEEA] * BathAvg[Exciton] + NN * rates[kECA] * BathAvg[Charge] + NN * rates[kDE] * BathAvg[Vacancy]) + (NN * rates[k_EEA] * BathAvg[Exciton] + NN * rates[k_ECA] * BathAvg[Charge] + NN * rates[k_DE] * BathAvg[Exciton]) * NN * rates[k_Rec] * BathAvg[Vacancy];

    OnSiteAvg[Exciton] = NN * (rates[k_DC] + rates[k_ECA]) * rates[kRec] * BathAvg[Charge] * BathAvg[Charge] + (rates[ki] + NN * (rates[k_DE] + rates[k_EEA]) * BathAvg[Exciton]) * (BathRatesAvg[Charge] + rates[kDC] * BathAvg[Vacancy]) + (rates[ki] * rates[kRec] + NN * (((rates[k_DE] + rates[k_EEA]) * BathAvg[Exciton] + BathRatesAvg[Exciton]) * rates[kRec] + rates[k_ECA] * (BathRatesAvg[Charge] + rates[kDC] * BathAvg[Vacancy]))) * BathAvg[Charge];
    
//    OnSiteAvg[Exciton] = NN * rates[kRec] * BathAvg[Charge] * (rates[ki] + NN * BathRatesAvg[Exciton] + NN * rates[k_DC] * BathAvg[Charge]) + (NN * rates[k_EEA] * BathAvg[Exciton] + NN * rates[k_ECA] * BathAvg[Charge] + NN * rates[k_DE] * BathAvg[Exciton]) * (rates[k_i] + NN * rates[kRec] * BathAvg[Charge] + NN * BathRatesAvg[Charge] + NN * rates[kDC] * BathAvg[Vacancy]);
    
    // calculate D
    double denom = OnSiteAvg[Vacancy] + OnSiteAvg[Charge] + OnSiteAvg[Exciton];

    // returns x_i at a given site
    OnSiteAvg[Vacancy] /= denom; 
    
    OnSiteAvg[Charge] /= denom; 

    OnSiteAvg[Exciton] /= denom; 
    
    
}
int main(int argc, char** argv) {


	double xAvg;
	double yAvg;
	int N = 50;
	double sigmaX;
	double sigmaY;
//	vector<double> VacI(2*N);
//	vector<double> ChargeI(2*N);
//	vector<double> ExcitonI(2*N);

	vector<double> OnSiteAvg(3);
	vector<double> BathAvg(3);
	vector<double> BathRatesAvg(3);
	vector<double> rates(NUMBER_OF_RATES);
	int disorderedRate;


	try {
                CmdLine cmd("Enter rates and stuff");

                ValueArg<double> krArg("r", "kr", "Exciton radiative rate", true,0.1,"double");
                cmd.add(krArg);
                ValueArg<double> kEEAArg("a", "eea", "Exciton exciton annhilation rate", true,0.1,"double");
                cmd.add(kEEAArg);
                ValueArg<double> kECAArg("d", "eca", "Exicton charge annihlation rate", true,0.1,"double");
                cmd.add(kECAArg);
                ValueArg<double> k_ECAArg("", "beca", "Exicton charge annihlation backward rate", true,0.1,"double"); // set to 0
                cmd.add(k_ECAArg);
                ValueArg<double> kRecArg("e", "rec", "Charge charge recombination rate", true,0.1,"double");
                cmd.add(kRecArg);
                ValueArg<double> k_RecArg("b", "brec", "Charge charge recombination backward rate", true,0.1,"double"); // constant
                cmd.add(k_RecArg);
                ValueArg<double> kiArg("", "ki", "Charge injection rate", true,0.1,"double"); // constant 
                cmd.add(kiArg);
                ValueArg<double> k_iArg("i", "kbi", "Charge injection backward rate", true,50,"double");
                cmd.add(k_iArg);
                ValueArg<double> meanXArg("", "meanX", "bimolecular reaction rate mean", true, -6 ,"double");
                cmd.add(meanXArg);
                ValueArg<double> sigmaXArg("", "sigmaX", "bi molecular reaction rate sigma", true,0,"double");
                cmd.add(sigmaXArg);
                ValueArg<double> meanYArg("", "meanY", "injection reaction rate mean", true, -6 ,"double"); // ki (controls current)
                cmd.add(meanYArg);
                ValueArg<double> sigmaYArg("", "sigmaY", "injection reaction rate sigma", true,0,"double"); // set to 0
                cmd.add(sigmaYArg);
                ValueArg<double> k_EEAArg("z", "beea", "Exciton exciton annhilation backward rate", true,0,"double"); // set to 0
                cmd.add(k_EEAArg);
                ValueArg<double> kDCArg("p", "dc", "Charge Diffusion rate", true,0,"double"); // constant
                cmd.add(kDCArg);
                ValueArg<double> k_DCArg("x", "bdc", "Charge Diffusion backward rate", true,0,"double"); // constant, = dc
                cmd.add(k_DCArg);
                ValueArg<double> kDEArg("y", "de", "Exciton Diffusion rate", true,0,"double"); // constant, = dc
                cmd.add(kDEArg);
                ValueArg<double> k_DEArg("q", "bde", "Exciton Diffusion backward rate", true,0,"double"); //constant, = dc
                cmd.add(k_DEArg);
                ValueArg<int> DisArg("", "disr", "Bimolecular rate to disorder", true,1,"int"); // index of the rate w disorder
                cmd.add(DisArg);

        
		ValueArg<double> initVac("v","vac", "initial Vacancy pop", false, 0.45,"double");
		cmd.add(initVac);
		ValueArg<double> initC("C","cpop", "initial charge pop", false, 0.45,"double");
		cmd.add(initC);

                cmd.parse(argc, argv);
//enum Rates {kr, kEEA, k_EEA, kECA, k_ECA, kRec, k_Rec, ki, k_i, K_Rec, K_EEA, K_ECA, K_i, kDC, k_DC, kDE,  k_DE};
		disorderedRate = DisArg.getValue();

        // rates = list of rate constants
		rates[kr] = krArg.getValue();
		rates[kEEA] = kEEAArg.getValue();
		rates[k_EEA] = k_EEAArg.getValue();
		rates[kECA] = kECAArg.getValue();
		rates[k_ECA] = k_ECAArg.getValue();
		rates[kRec] = kRecArg.getValue();
		rates[k_Rec] = k_RecArg.getValue();
		rates[ki] = kiArg.getValue();
		rates[k_i] = k_iArg.getValue();


		//equilibrium constants at end of rates vector
		if (rates[k_Rec] == 0) rates[k_Rec] = 0.00000001;
		if (rates[k_EEA] == 0) rates[k_EEA] = 0.00000001;
		if (rates[k_ECA] == 0) rates[k_ECA] = 0.00000001;
		if (rates[k_i] == 0)   rates[k_i] = 0.00000001;
		//if (rates[k_DC] == 0)   rates[k_DC] = 0.00000001;
		//if (rates[k_DE] == 0)   rates[k_DE] = 0.00000001;


		//equilibrium constants 
		rates[K_i]   = rates[ki] / rates[k_i];
		rates[K_EEA] = rates[kEEA] / rates[k_EEA];
		rates[K_ECA] = rates[kECA] / rates[k_ECA];
		rates[K_Rec] = rates[kRec] / rates[k_Rec];

		//diffusion constants
		rates[kDC] = kDCArg.getValue();
		rates[k_DC] = k_DCArg.getValue();
		rates[kDE] = kDEArg.getValue();
		rates[k_DE] = k_DEArg.getValue();
        
		xAvg = meanXArg.getValue();
                sigmaX = sigmaXArg.getValue();

		yAvg = meanYArg.getValue();
                sigmaY = sigmaYArg.getValue();

		BathAvg[Vacancy] = initVac.getValue();
		BathAvg[Charge] = initC.getValue();
		BathAvg[Exciton] = 1 - BathAvg[Charge] - BathAvg[Vacancy];

        } catch (ArgException &e) {
                cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
                return 1;
        }

	cout.precision(17);

//	for (double sigma = 100*exp(-15); sigma < 1; sigma++) {
//		VacB = 0.45;
//		HB = VacB;
			
	BathRatesAvg[Charge] = BathAvg[Charge] * rates[kRec];
	BathRatesAvg[Exciton] = BathAvg[Exciton] * rates[k_Rec];

	if (sigmaX == 0) sigmaX = 0.0000001;
	if (sigmaY == 0) sigmaY = 0.0000000000000001;
	double dx = 0.15 * sqrt(sigmaX);
	double dy = 0.15 * sqrt(sigmaY);
	double Px;
	double Py;
	

    // BEGIN SCF: res = 0 when converged
	double res = 1;
	while (res > pow(10,-8)) {
		vector<double> Accum(3, 0);
		vector<double> RatesAccum(3, 0);
		
		double Norm = 0;
		
		for (int j = -0; j <= 0; j++){
			double y = yAvg + j*dy;
			Py = 1;
			dy = 1;
			//Py = exp(-(y-yAvg)*(y-yAvg)/(2*sigmaY)) / sqrt(2*sigmaY*3.1415926536);
			// First disorderd rate
			rates[ki] = rates[k_i]*exp(-y);
//			rates[ki] = exp(-y);
		
            // integrates over distribution of disordered rate	
			for (int i = -N; i <= N; i++){
				double x = xAvg + i*dx;
				Px = exp(-(x-xAvg)*(x-xAvg)/(2*sigmaX)) / sqrt(2*sigmaX*3.1415926536);
				


				// disorderd rate here
				rates[disorderedRate] = exp(-x);
				if (disorderedRate == kRec) rates[k_Rec] = rates[kRec] * rates[K_Rec];
							

                // get concentrations of each species				
				vector<double> OnSite(3);
				getConc(OnSite, BathAvg, BathRatesAvg, rates);

                // weighted avg over dist. to get new averages (eventually)
				// if i = N or -N, only 2 nearest neighbors (???)
                if ((i == -N) || (i == N)) {
					Accum[Vacancy] += OnSite[Vacancy] * Px * Py * 0.5;
					Accum[Charge] += OnSite[Charge] * Px * Py * 0.5;
					Accum[Exciton] += OnSite[Exciton] * Px * Py * 0.5;

					RatesAccum[Vacancy] += OnSite[Vacancy] * Px * Py *  rates[k_Rec] * 0.5;
					RatesAccum[Charge] += OnSite[Charge] * Px * Py * rates[kRec] * 0.5;
					RatesAccum[Exciton] += OnSite[Exciton] * Px * Py * rates[k_Rec] * 0.5;
					Norm += Px * Py * 0.5;
				}
				else {
					Accum[Vacancy] += OnSite[Vacancy] * Px * Py;
					Accum[Charge] += OnSite[Charge] * Px * Py;
					Accum[Exciton] += OnSite[Exciton] * Px * Py;

					RatesAccum[Vacancy] += OnSite[Vacancy] * Px * Py * rates[k_Rec];
					RatesAccum[Charge] += OnSite[Charge] * Px * Py * rates[kRec];
					RatesAccum[Exciton] += OnSite[Exciton] * Px * Py * rates[k_Rec];
					Norm += Px * Py;
				}
	       // cout << ki << endl;
			}
		}
        // new averages!
		Accum[Vacancy] *= dx * dy;
		Accum[Charge]  *= dx * dy;
		Accum[Exciton] *= dx * dy; 

        // new average rates (avg(krec*x))
		RatesAccum[Vacancy] *= dx * dy;
		RatesAccum[Charge]  *= dx * dy;
		RatesAccum[Exciton] *= dx * dy; 

		Norm *= dx * dy;

        // calc difference for scf convergence
		res = abs(Accum[Vacancy] - BathAvg[Vacancy]) + abs(Accum[Charge] - BathAvg[Charge]) +  abs(Accum[Exciton] - BathAvg[Exciton]); 
		
        // update with new avg's
		BathAvg[Vacancy] = Accum[Vacancy];
		BathAvg[Charge] = Accum[Charge];
		BathAvg[Exciton] = Accum[Exciton];

		BathRatesAvg[Vacancy] = RatesAccum[Vacancy];
		BathRatesAvg[Charge] = RatesAccum[Charge];
		BathRatesAvg[Exciton] = RatesAccum[Exciton];
	}
//	for (int i= 0; i < 2*N + 1;i++) {
//		cout << xAvg + (i - N) * dx << "\t " << VacI[i] << "\t " << HI[i]  << "\t " << H2I[i]<< endl;
//	}

    // if converged, return avg phi, c, e
	cout << BathAvg[Vacancy] << " " << BathAvg[Charge] << " " << BathAvg[Exciton] << endl;
//	}
	return 0;
}
