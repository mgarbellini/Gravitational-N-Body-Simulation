// F. Borando - francesco.borando@studenti.unimi.it
// M. Garbellini - matteo.garbellini@studenti.unimi.it
// Università degli Studi di Milano
// Dept. of Physics
// January 2019
//
// GRAVITATIONAL N-BODY SIMULATION
// tool that produces a graphic of energy distribution from an input file produced by energy_data_fixer.cpp

#include <iostream>
#include <fstream>
#include <cstring>
#include "TApplication.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

int main(int argc, char* argv[])
{
	if(!argv[1])
    {
        std::cerr << "Please FILENAME as argv[1]!" << std::endl;
        std::exit (EXIT_FAILURE);
    }


	TApplication app("App",NULL, NULL);
	TCanvas *c = new TCanvas("c", "Energy Plot",800, 800);

	TGraph *graph = new TGraph(argv[1], "%lg %lg");
	
	//Titles
	graph->SetTitle("Energy Plot");
	graph->GetXaxis()->SetTitle("Iteration");
	graph->GetYaxis()->SetTitle("Energy");

	//Marker style and Draw
	//graph->SetMarkerStyle(kCross);
	graph->Draw("AP");
    
    std::string title = "ener_";
    title += argv[1];
    title.erase(title.end()-4,title.end());
    title += ".pdf";
	c->Print(title.c_str());
	//app.Run(kFALSE);
	return 0;

}
