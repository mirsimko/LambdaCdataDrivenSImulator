/* *********************************************************************
 *  ROOT macro - generates 2D slices in pT from 3D distributions of DCA
 *  
 *  Authors:
 *            **Miroslav Simko
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/


#include "TH3D.h"
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"

#include <iostream>
using namespace std;

TH3D* getHisto(TFile &fDca1, int iEta, int iVz, int iCent) ;

// ----------------------------------------------------------------------------
void create2Dhistos()
{
  // contstants
  // const Int_t nParticles = 3;
  const int nVzs = 4;
  const int nEtas = 5;
  const int nCentDca = 6;
  const Int_t nPtBins = 24;
  const double epsilon = 0.0001;

  // input file and output file
  TFile fDca1("IncludeProton/DCA/Proton/NoBinWidth_3D_Dca_VsPt_Centrality_Eta_Phi_Vz_Zdcx_proton.root");
  TFile *outHistF = new TFile("IncludeProton/DCA/Proton/2DProjection_simCent_NoBinWidth_3D_Dca_VsPt_Centrality_Eta_Phi_Vz_Zdcx_proton_trimmed.root", "RECREATE");

  char n;
  cout << "Press any key" << endl;
  cin >> n;

  // edges of the pt bins in 3D DCA histograms
  const Double_t ptEdge[nPtBins + 1] =
  { 0. , 0.2 , 0.3 , 0.4 , 0.5 , 0.6 , 0.7 , 0.8 , 0.9 ,
    1. , 1.2 , 1.4 , 1.6 , 1.8 , 2.  , 2.2 , 2.4 , 2.6 , 2.8 ,
    3. , 3.5, 4.  , 4.5 , 5., 12.0};

  // edges of the centrality in 3D DCA histograms
  const Int_t centEdge[nCentDca + 1] = 
  { 0, 1, 2, 3, 4, 6, 8};

  for (int iEta = 0; iEta < nEtas; ++iEta)
  {
    for (int iVz = 0; iVz < nVzs; ++iVz)
    {
      // index of centrality in 3D histograms
      int iCent3d = 0;

      for(int iCent = 0; iCent < nCentDca; ++iCent)
      {
	// get the 3D histo
	TH3D *hist3D = getHisto(fDca1, iEta, iVz, iCent3d);

	// getting parameters of the 3D histogram (nBins, min and max of the axes, and title)
	int nXyBins = hist3D->GetYaxis()->GetNbins();
	double xyMin = hist3D->GetYaxis()->GetXmin();
	double xyMax = hist3D->GetYaxis()->GetXmax();
	const char* xyTitle = hist3D->GetYaxis()->GetTitle();
	int nZBins = hist3D->GetZaxis()->GetNbins();
	double zMin = hist3D->GetZaxis()->GetXmin();
	double zMax = hist3D->GetZaxis()->GetXmax();
	const char* zTitle = hist3D->GetZaxis()->GetTitle();
	const char *hist3dTitle = hist3D->GetTitle();

	// index of pT bin in 3D histogram
	int iPt3d = 1;

	// slicing 3D histogram into 2D slices
	for (int iPt = 0; iPt < nPtBins; ++iPt)
	{
	  // creating the 2D histogram
	  const char *h2dName = Form("mh2DcaPtCentPartEtaVz_%i_%i_%i_%i_%i", 2, iEta, iVz, iCent, iPt);
	  double midPt = (ptEdge[iPt+1]+ptEdge[iPt])*0.5; // middle of the pT bin
	  const char *h2dTitle = Form("%s_pT-%.2f", hist3dTitle, midPt);
	  TH2D hist2D(h2dName, h2dTitle, nXyBins, xyMin, xyMax, nZBins, zMin, zMax);
	  hist2D.GetXaxis()->SetTitle(xyTitle);
	  hist2D.GetYaxis()->SetTitle(zTitle);

	  while ( iCent3d  < centEdge[iCent+1] )
	  {
	    // get 3D histo
	    hist3D = getHisto(fDca1, iEta, iVz, iCent3d);
	    cout << "Title:" << h2dTitle << ", centrality: " << iCent << ", pT: " << iPt << endl;

	    while ( (hist3D->GetXaxis()->GetBinUpEdge(iPt3d)) - epsilon < ptEdge[iPt+1])
	    {
	      // filling the 2D histogram
	      for (int iXy = 0; iXy < nXyBins+2; ++iXy)
	      {
		for (int iZ = 0; iZ < nZBins+2; ++iZ)
		{
		  // copy bin content
		  double last =  hist2D.GetBinContent(iXy, iZ) ;
		  hist2D.SetBinContent(iXy, iZ, hist3D->GetBinContent(iPt3d,iXy,iZ) + last );
		}
	      }
	      // raise index of pT bin of 3D hist
	      ++iPt3d;
	    }
	    // raise index of Centrality 3D hist
	    ++iCent3d;
	    // delete the 3D histo
	    // delete hist3D;
	  }

	  // testing if some particles were recorded
	  // if (hist2D.Integral() == 0)
	  // {
	  //   // this means something went wrong
	  //   cerr << "Number of particles recorded: " << hist2D.Integral() 
	  //     << " from pT histogram: " << iPt 
	  //     << " between " << ptEdge[iPt] << " and " << ptEdge[iPt+1] << endl;
	  //   cerr << "Centrality of 2D: " << iCent << ", Centrality of 3D: " << iCent3d << endl;
	  //   cerr << ", pT bin number: " << iPt3d 
	  //     << " between " << hist3D->GetXaxis()->GetBinLowEdge(iPt3d-1)  << " and " << hist3D->GetXaxis()->GetBinUpEdge(iPt3d-1) << endl;
	  //   cerr << "Title: " << h2dTitle << endl;
	  //   cerr << "Name: " << h2dName << endl;
	  //   throw;
	  // }

	  // Saving the 2D histogram
	  outHistF->cd();
	  hist2D.Write();  
	}
      }
    }
    cout<<"Finished writing histograms for Eta " << iEta << endl;
  }
  cout << endl;
  cout << "Done..." << endl;
}

// ----------------------------------------------------------------------------
// Getting 3D histogram
TH3D* getHisto(TFile &fDca1, int iEta, int iVz, int iCent) 
{
  TH3D* hist3D = 0;
  const char *h3dName = Form("mh3DcaXyZPtCentPartEtaVz_%i_%i_%i_%i", 0, iEta, iVz, iCent+2); 
  hist3D = (TH3D*)(fDca1.Get(h3dName));
  if (!hist3D)
  {
    cerr << "getHisto: histogram \"" << h3dName << "\" not found." << endl;
    throw;
  }
  return hist3D;
}
// ----------------------------------------------------------------------------
