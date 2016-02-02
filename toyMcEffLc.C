/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for Lambda_c decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *            **Miroslav Simko (simkomir@fjfi.cvut.cz)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/

#include <iostream>
#include <fstream>
#include <new>
#include <cmath>

#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TF1.h"
#include "TClonesArray.h"
#include "TPythia6.h"
#include "TPythia6Decayer.h"
#include "TRandom3.h"
#include "TParticle.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TGraph.h"
#include "TMath.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"
#include "dataDrivenFastSimulator.h"

using namespace std;

void setDecayChannels(int const mdme);
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters);
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TLorentzVector const& pMom, TVector3 v00, int LambdaCcharge);
void getKinematics(TLorentzVector& b, double const mass);
float resMass(TLorentzVector const &pMom, TLorentzVector const &kMom, TLorentzVector const &piMom, int decayMode);

void bookObjects();
void write();
TPythia6Decayer* pydecay;
TNtuple* nt;
TNtuple* ntTMVA;
TFile* result;

TF1* fWeightFunction = NULL;

string outFileName;
std::pair<int, int> const decayChannels(4277, 4354); // first and last Lc decay channel number
std::pair<float, float> const momentumRange(0, 12);

float const acceptanceRapidity = 1.0;
float const M_KS = 0.49767;
DecayMode mDecayMode;

bool const saveNt = true;

// centrality and p_T distributions
TH1D *nBinCent;
float const nBin[nCent] = {1012, 805, 577, 365, 221, 127, 66.8, 32.4, 15.};
float maxNbin = 3221.2;
TF1* fLevy;


float const Dyield = 1.39751; /*D0 yiield per event*/ 
float const LambdaDratio = 0.275;/*ratio between Lambda_c and D0, taken from arXiv:hep-ex/0508019 - Table 4*/ 

//============== main  program ==================
void toyMcEffLc(unsigned long nEvts = 1000, int startCent = 0, int endCent = 8, const char *outputFileName = "Lc.toyMC.test.root", int modeOfDecay = 3)
{
   cout << "Ntuples will be saved to \"" << outputFileName << "\"" << endl;
   // cout << "Mass = " << M_LAMBDA_C_PLUS << endl;

   outFileName = outputFileName;
   mDecayMode = (DecayMode)modeOfDecay;
   
   gRandom->SetSeed();

   pydecay = TPythia6Decayer::Instance();
   // Reading new decay channels
   pydecay-> SetDecayTableFile("AddedDecays.list");
   pydecay-> ReadDecayTable();

   pydecay->Init();

   // char input;
   // cout << "Pres any key" << endl;
   // cin >> input;
   // TPythia6::Instance()->Pylist(12); // this is for writing the Decay table to std_out

   loadAllDistributions(startCent, endCent);
   bookObjects();

   
   // selecting decay channels
   float branchingRatio = 0;
   switch(mDecayMode)
   {
     case kKstarProton:
       TPythia6::Instance()->SetMDME(617,1,1);
       TPythia6::Instance()->SetMDME(618,1,0);
       TPythia6::Instance()->SetMDME(619,1,0);

       setDecayChannels(4294);

       branchingRatio = 0.016;
       break;
     case kLambda1520Pion:
       TPythia6::Instance()->SetMDME(4276,1,1);

       setDecayChannels(4344);

       branchingRatio = 0.018;
       break;
     case kDeltaPPkaon:
       TPythia6::Instance()->SetMDME(1052,1,1);

       setDecayChannels(4291);

       branchingRatio = 0.0086;
       break;
     case kPionKaonProton:
       setDecayChannels(4343);

       branchingRatio = 0.05;
       break;
     default:
       cerr << "Unknown decay mode: exiting" << endl;
       return;
       break;       
   }

   unsigned long const npart = (unsigned long) floor( Dyield * LambdaDratio * branchingRatio * (double)nEvts * 0.5*totalNBin/maxNbin); // has to be divided by 2 because we are creating LC+,LC- pairs
   cout << "Number of produced Lambda_C: " << npart*2 << endl;

   TLorentzVector* b_d = new TLorentzVector;
   TClonesArray ptl("TParticle", 10);
   for (unsigned long ipart = 0; ipart < npart; ipart++)
   {
      if (!(ipart % 5000))
         cout << "____________ 2*ipart = " << 2*ipart << " ________________" << endl;

      getKinematics(*b_d, M_LAMBDA_C_PLUS);

      decayAndFill(4122, b_d, M_LAMBDA_C_PLUS, ptl);
      decayAndFill(-4122, b_d, M_LAMBDA_C_PLUS, ptl);

      if (ipart%1000 == 1) // save
      {
	result->cd();
	if(saveNt)
	  nt->AutoSave("SaveSelf");
	// nt->FlushBaskets();
	ntTMVA->AutoSave("SaveSelf");
	// ntTMVA->FlushBaskets();
      }
   }

   write();
   delete b_d;
}
// =======================================================

void setDecayChannels(int const mdme)
{
   for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
   TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

// =======================================================
void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters)
{
   pydecay->Decay(kf, b);
   pydecay->ImportParticles(&daughters);

   TLorentzVector kMom;
   TLorentzVector piMom;
   TLorentzVector pMom;
   TVector3 v00;

   int nTrk = daughters.GetEntriesFast();
   int LambdaCcharge = 0;
   // cout << "Number of daughters: " <<  nTrk << endl;
   for (int iTrk = 0; iTrk < nTrk; ++iTrk)
   {
      TParticle* ptl0 = (TParticle*)daughters.At(iTrk);
      // cout << "Daughter PDG number: " << ptl0->GetPdgCode() << endl;

      switch (abs(ptl0->GetPdgCode()))
      {
         case 321:
            ptl0->Momentum(kMom);
            v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
	    // cout << "Kaon" << endl;
            break;
         case 211:
            ptl0->Momentum(piMom);
	    if(signbit(ptl0->GetPdgCode()))
	      LambdaCcharge = -1;
	    else
	      LambdaCcharge = 1;
	    // cout << "Pion" << endl;
            break;
         case 2212:
            ptl0->Momentum(pMom);
	    // cout << "Proton" << endl;
            break;
         default:
            break;
      }
   }
   daughters.Clear();

   try
   {
     fill(kf, b, weight, kMom, piMom, pMom, v00, LambdaCcharge);
   }
   catch(std::bad_alloc &ba)
   {
     cerr << "bad_alloc in \"fill\" function: " << ba.what() << endl;
     throw ba;
   }
}

// =======================================================
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TLorentzVector const& pMom, TVector3 v00, int LambdaCcharge)
{
   int const centrality = floor(nBinCent->GetRandom());

   TVector3 const vertex = getVertex(centrality);
   // smear primary vertex
   // float const sigmaVertex = sigmaVertexCent[cent];
   // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

   v00 += vertex;

   // smear momentum
   // cout << "piPt = " << piMom.Perp() << ", kPt = " << kMom.Perp() << ", pPt = "<< pMom.Perp() << endl;
   TLorentzVector const piRMom = smearMom(0, piMom);
   TLorentzVector const kRMom = smearMom(1, kMom);
   TLorentzVector const pRMom = smearMom(2, pMom);
   // cout << "piRPt = " << piRMom.Perp() << ", kRPt = " << kRMom.Perp() << ", pRPt = "<< pRMom.Perp() << endl;

   // smear position
   // cout << "Smearing position ..." << endl;
   // cout << "Pion" << endl;
   TVector3 const piRPos = smearPosData(0, vertex.z(), centrality, piRMom, v00);
   // cout << "Kaon" << endl;
   TVector3 const kRPos = smearPosData(1, vertex.z(), centrality, kRMom, v00);
   // cout << "Proton" << endl;
   TVector3 const pRPos = smearPosData(2, vertex.z(), centrality, pRMom, v00);
   // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
   // TVector3 const pRPos = smearPos(pMom, pRMom, v00);

   // reconstruct
   // cout << "Reconstructing..." << endl;
   TLorentzVector const rMom = pRMom + piRMom + kRMom;
   float const piDca = dca(piMom.Vect(), v00, vertex);
   float const kDca = dca(kMom.Vect(), v00, vertex);
   float const pDca = dca(pMom.Vect(), v00, vertex);
   float const piRDca = dca(piRMom.Vect(), piRPos, vertex);
   float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
   float const pRDca = dca(pRMom.Vect(), pRPos, vertex);

   // cout << "Calculating smeared vertex..." << endl;
   // smeared decay vertex
   TVector3 v0;
   // vertices between pairs
   TVector3 v012;
   TVector3 v023;
   TVector3 v013;

   // cout << "Calculating DCAs..." << endl;
   float const dca12 = dca1To2(kRMom.Vect(), kRPos, piRMom.Vect(), piRPos, v012);
   float const dca23 = dca1To2(piRMom.Vect(), piRPos, pRMom.Vect(), pRPos, v023);
   float const dca13 = dca1To2(kRMom.Vect(), kRPos, pRMom.Vect(), pRPos, v013);

   v0 = (1./3.) * (v012 + v023 + v013) ;

   // distances of the pairs vertices
   float const vDist1 = (v012 - v023).Mag();
   float const vDist2 = (v023 - v013).Mag();
   float const vDist3 = (v012 - v013).Mag();
   float vDistMax = vDist1 > vDist2 ? vDist1 : vDist2;
   vDistMax = vDistMax > vDist3 ? vDistMax : vDist3;

   float const decayLength = (v0 - vertex).Mag();
   float const dcaToPv = dca(rMom.Vect(), v0, vertex);
   float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());

   // cout << "is HFT track?" << endl;
   bool const isPiHft = matchHft(0, vertex.z(), centrality, piRMom);
   bool const isKhft =  matchHft(1, vertex.z(), centrality, kRMom );
   bool const isPhft =  matchHft(2, vertex.z(), centrality, pRMom );

   // is Tpc track?
   bool const isPiTpc = tpcReconstructed(0,  LambdaCcharge, centrality, piRMom);
   bool const isKTpc =  tpcReconstructed(1, -LambdaCcharge, centrality, kRMom );
   bool const isPTpc =  tpcReconstructed(2,  LambdaCcharge, centrality, pRMom );

                       // save

   // cout << "saving..." << endl;
   // error block
   try
   {
     if (saveNt)
     {
       float arr[100];
       int iArr = 0;
       arr[iArr++] = centrality;
       arr[iArr++] = vertex.X();
       arr[iArr++] = vertex.Y();
       arr[iArr++] = vertex.Z();

       arr[iArr++] = kf;
       arr[iArr++] = b->M();
       arr[iArr++] = b->Perp();
       arr[iArr++] = b->PseudoRapidity();
       arr[iArr++] = b->Rapidity();
       arr[iArr++] = b->Phi();
       arr[iArr++] = v00.X();
       arr[iArr++] = v00.Y();
       arr[iArr++] = v00.Z();

       arr[iArr++] = rMom.M();
       arr[iArr++] = rMom.Perp();
       arr[iArr++] = rMom.PseudoRapidity();
       arr[iArr++] = rMom.Rapidity();
       arr[iArr++] = rMom.Phi();

       arr[iArr++] = dca12;
       arr[iArr++] = dca23;
       arr[iArr++] = dca13;
       arr[iArr++] = decayLength;
       arr[iArr++] = dcaToPv;
       arr[iArr++] = cosTheta;

       arr[iArr++] = kMom.M();
       arr[iArr++] = kMom.Perp();
       arr[iArr++] = kMom.PseudoRapidity();
       arr[iArr++] = kMom.Rapidity();
       arr[iArr++] = kMom.Phi();
       arr[iArr++] = kDca;

       arr[iArr++] = kRMom.M();
       arr[iArr++] = kRMom.Perp();
       arr[iArr++] = kRMom.PseudoRapidity();
       arr[iArr++] = kRMom.Rapidity();
       arr[iArr++] = kRMom.Phi();
       arr[iArr++] = kRPos.X();
       arr[iArr++] = kRPos.Y();
       arr[iArr++] = kRPos.Z();
       arr[iArr++] = kRDca;

       arr[iArr++] = piMom.M();
       arr[iArr++] = piMom.Perp();
       arr[iArr++] = piMom.PseudoRapidity();
       arr[iArr++] = piMom.Rapidity();
       arr[iArr++] = piMom.Phi();
       arr[iArr++] = piDca;

       arr[iArr++] = piRMom.M();
       arr[iArr++] = piRMom.Perp();
       arr[iArr++] = piRMom.PseudoRapidity();
       arr[iArr++] = piRMom.Rapidity();
       arr[iArr++] = piRMom.Phi();
       arr[iArr++] = piRPos.X();
       arr[iArr++] = piRPos.Y();
       arr[iArr++] = piRPos.Z();
       arr[iArr++] = piRDca;

       arr[iArr++] = pMom.M();
       arr[iArr++] = pMom.Perp();
       arr[iArr++] = pMom.PseudoRapidity();
       arr[iArr++] = pMom.Rapidity();
       arr[iArr++] = pMom.Phi();
       arr[iArr++] = pDca;

       arr[iArr++] = pRMom.M();
       arr[iArr++] = pRMom.Perp();
       arr[iArr++] = pRMom.PseudoRapidity();
       arr[iArr++] = pRMom.Rapidity();
       arr[iArr++] = pRMom.Phi();
       arr[iArr++] = pRPos.X();
       arr[iArr++] = pRPos.Y();
       arr[iArr++] = pRPos.Z();
       arr[iArr++] = pRDca;

       arr[iArr++] = isKhft;
       arr[iArr++] = isPiHft;
       arr[iArr++] = isPhft;

       arr[iArr++] = isKTpc;
       arr[iArr++] = isPiTpc;
       arr[iArr++] = isPTpc;

       arr[iArr++] = resMass(pMom, kMom, piMom, mDecayMode);
       arr[iArr++] = resMass(pRMom, kRMom, piRMom, mDecayMode);

       arr[iArr++] = vDistMax;
       arr[iArr++] = vDist1;
       arr[iArr++] = vDist2;
       arr[iArr++] = vDist3;

       nt->Fill(arr);
     } // if (saveNt)
   }
   catch(std::bad_alloc &ba)
   {
     cerr << "bad_alloc in saving \"nt\": " << ba.what() << endl;
     throw ba;
   } // end of error block
   
   // error block
   try	
   {
     // __________________________________________
     // using cuts
     // __________________________________________
     bool const PtCut = pRMom.Perp() > 0.3 && kRMom.Perp() > 0.3 && piRMom.Perp() > 0.3;
     bool const dcaCut = dca12 < 200 && dca23 < 200 && dca13 < 200;
     bool const dLengthCut = decayLength > 30;
     bool const cosThetaCut = cosTheta > 0.98;
     bool const HftCut = isPhft && isKhft && isPiHft;
     bool const EtaCut = TMath::Abs(kRMom.PseudoRapidity()) < 1. && TMath::Abs(pRMom.PseudoRapidity()) < 1. && TMath::Abs(piRMom.PseudoRapidity()) < 1.;
     bool const tpcCut = isKTpc && isPiTpc && isPTpc;

     if ( !( PtCut && dcaCut && dLengthCut && cosThetaCut && HftCut && EtaCut && tpcCut) )
     {
       // Which one did not go?
       // cout << "| " ;
       // if(!PtCut)
       //   cout << "PtCut == 0       | ";
       // else
       //   cout << "                 | ";
       // if(!dcaCut)
       //   cout << "dcaCut == 0      | ";
       // else
       //   cout << "                 | ";
       // if(!dLengthCut)
       //   cout << "dLengthCut == 0  | ";
       // else
       //   cout << "                 | ";
       // if(!cosThetaCut)
       //   cout << "cosThetaCut == 0 | ";
       // else
       //   cout << "                 | ";
       // if(!HftCut)
       //   cout << "HftCut == 0      | ";
       // else
       //   cout << "                 | ";
       // if(!EtaCut)
       //   cout << "EtaCut == 0      | ";
       // else
       //   cout << "                 | ";
       // if(!tpcCut)
       //   cout << "tpcCut == 0      | ";
       // else
       //   cout << "                 | ";
  
       // cout << endl;

       // if(HftCut)
       //   cout << "********************* HFT track *********************" << endl;

       return;
     }
     else
       cout << "************************************" << endl;
       cout << "Good Lambda_c" << endl;
       cout << "************************************" << endl;

     // __________________________________________
     // end of cuts
     // filling TMVA histograms
     const float umToCm = 0.0001;

     float TMVA[100];
     int iTMVA = 0;

     TMVA[iTMVA++] = rMom.M();
     TMVA[iTMVA++] = rMom.Perp();
     TMVA[iTMVA++] = 1;
     TMVA[iTMVA++] = rMom.Phi();
     TMVA[iTMVA++] = rMom.PseudoRapidity();

     TMVA[iTMVA++] = dca12*umToCm;
     TMVA[iTMVA++] = dca23*umToCm;
     TMVA[iTMVA++] = dca13*umToCm;

     TMVA[iTMVA++] = cosTheta;
     TMVA[iTMVA++] = decayLength*umToCm;

     TMVA[iTMVA++] = kRMom.Perp();
     TMVA[iTMVA++] = pRMom.Perp();
     TMVA[iTMVA++] = piRMom.Perp();

     TMVA[iTMVA++] = kRDca*umToCm;
     TMVA[iTMVA++] = pRDca*umToCm;
     TMVA[iTMVA++] = piRDca*umToCm;

     TMVA[iTMVA++] = vDistMax*umToCm;

     ntTMVA->Fill(TMVA);
   }
   catch(std::bad_alloc &ba)
   {
     cerr << "bad_alloc in saving \"ntTMVA\": " << ba.what() << endl;
     throw ba;
   } // end of error block
}

//___________
float resMass(TLorentzVector const &pMom, TLorentzVector const &kMom, TLorentzVector const &piMom, int decayMode)
{
  TLorentzVector resMom;
  switch (decayMode)
  {
    case kPionKaonProton:
      return 0;
      break;
    case kKstarProton:
      resMom = kMom + piMom;
      return resMom.M();
      break;
    case kLambda1520Pion:
      resMom = pMom + kMom;
      return resMom.M();
      break;
    case kDeltaPPkaon:
      resMom = pMom + piMom;
      return resMom.M();
      break;
    default:
      return -1;
      break;
  };
}
//
//___________
void getKinematics(TLorentzVector& b, double const mass)
{
   float const pt = fLevy->GetRandom();
   // float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
   //cout << "Momentum: " << pt << endl;
   float const y = gRandom->Uniform(-acceptanceRapidity, acceptanceRapidity);
   float const phi = TMath::TwoPi() * gRandom->Rndm();

   float const mT = sqrt(mass * mass + pt * pt);
   float const pz = mT * sinh(y);
   float const E = mT * cosh(y);

   b.SetPxPyPzE(pt * cos(phi), pt * sin(phi) , pz, E);
}

//___________
void bookObjects()
{
   result = new TFile(outFileName.c_str(), "recreate");
   result->cd();

   int BufSize = (int)pow(2., 16.);
   if(saveNt)
   {
     nt = new TNtuple("nt", "", "cent:vx:vy:vz:"
		      "pid:m:pt:eta:y:phi:v0x:v0y:v0z:"
		      "rM:rPt:rEta:rY:rPhi:"
		      "dca12:dca23:dca13:decayLength:dcaToPv:cosTheta:" // Rc pair
		      "kM:kPt:kEta:kY:kPhi:kDca:" // MC Kaon
		      "kRM:kRPt:kREta:kRY:kRPhi:kRVx:kRVy:kRVz:kRDca:" // Rc Kaon
		      "piM:piPt:piEta:piY:piPhi:piDca:" // MC Pion
		      "piRM:piRPt:piREta:piRY:piRPhi:piRVx:piRVy:piRVz:piRDca:" // Rc Pion
		      "pM:pPt:pEta:pY:pPhi:pDca:" // MC Proton
		      "pRM:pRPt:pREta:pRY:pRPhi:pRVx:pRVy:pRVz:pRDca:" // Rc Proton
		      "kHft:piHft:pHft:"
		      "kTpc:piTpc:pTpc:"
		      "MResonance:MRResonance:"
		      "maxVertexDist:vDist1:vDist2:vDist3",
		      BufSize); // distances of vertices of track pairs
   }

   ntTMVA = new TNtuple("ntTMVA", "", "m:pt:charges:phi:eta:" // basic properties of Lambda_c
				      "dcaDaugthers31:dcaDaugthers23:dcaDaugthers12:" // dca daughters (pi-K, pi-p, p-K)
				      "cosPntAngle:dLength:" // cosTheta and decay Length
				      "p1pt:p2pt:p3pt:"
				      "p1Dca:p2Dca:p3Dca:" // daughters (K, p, pi)
				      "maxVertexDist", BufSize);

   // setting the centrality dependence histogram
   nBinCent = new TH1D("nBinCent", "Number of bins vs centrality", endCentrality - beginCentrality + 1, beginCentrality, endCentrality);
   for (int i = 0; i < endCentrality - beginCentrality + 1; ++i)   
   {
      nBinCent->SetBinContent(i+1, nBin[beginCentrality + i]);
      totalNBin+=nBin[beginCentrality + i];
   }

   // p_T spectrum
   fLevy = new TF1("fLevy", "[0]*TMath::Exp(-[1]/(x-[3]))*TMath::Power(x-[3],1-[2])*x", momentumRange.first, momentumRange.second);
   // in latex: $A  \left( \frac{\mu}{x-\phi} \right)^{1-\alpha} \exp \left(- \frac{\mu}{x-\phi}\right)$
   // parameters: A = 8.17808e+06, \mu = 1.92432e+01, \alpha = 1.39339e+01, \phi = -9.04949e-01
   //
   // this was fitted from published D0 data
   // The additional p_T is from Jacobian
   //
   fLevy->SetParameters(8.17808e+06, 1.92432e+01, 1.39339e+01, -9.04949e-01);
}
//___________
void write()
{
   result->cd();
   if(saveNt)
     nt->Write();
   ntTMVA->Write();
   result->Close();
}
//___________
