/* *********************************************************************
 *  ROOT macro - Toy Monte Carlo Simulation for D0 decay
 *  Includes Momentum Resolution, DCA, hft ration, TPC efficiency ...
 *  Example for D0 --> Kpi
 *
 *  Authors:
 *            Guannan Xie (guannanxie@lbl.gov)
 *            **Mustafa Mustafa (mmustafa@lbl.gov)
 *            Hao Qiu (hqiu@lbl.gov)
 *
 *  ** Code Maintainer
 *
 * *********************************************************************
*/

#include <iostream>
#include <fstream>

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
void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TVector3 v00);
void getKinematics(TLorentzVector& b, double const mass);

void bookObjects();
void write();
TPythia6Decayer* pydecay;
TNtuple* nt;
TFile* result;

TF1* fWeightFunction = NULL;

string outFileName = "KS.toyMc.root";
std::pair<int, int> const decayChannels(613,614);
std::pair<float, float> const momentumRange(0, 12);

float const acceptanceRapidity = 1.0;
float const M_KS = 0.49767;

bool const saveNt = false;
//============== main  program ==================
void toyMcEffLc(int npart = 100)
{
   gRandom->SetSeed();
   loadAllDistributions();
   bookObjects();

   pydecay = TPythia6Decayer::Instance();
   pydecay->Init();

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

   TLorentzVector* b_d = new TLorentzVector;
   TClonesArray ptl("TParticle", 10);
   for (int ipart = 0; ipart < npart; ipart++)
   {
      if (!(ipart % 100000))
         cout << "____________ ipart = " << ipart << " ________________" << endl;

      getKinematics(*b_d, M_LAMBDA_C_PLUS);

      decayAndFill(4122, b_d, ptl, mDecayMode);
      decayAndFill(-4122, b_d, ptl, mDecayMode);

      if (ipart%1000 == 1) // save
      {
	nt->AutoSave("SaveSelf");
	ntTMVA->AutoSave("SaveSelf");
      }
   }

   write();
}

void setDecayChannels(int const mdme)
{
   for (int idc = decayChannels.first; idc < decayChannels.second + 1; idc++) TPythia6::Instance()->SetMDME(idc, 1, 0);
   TPythia6::Instance()->SetMDME(mdme, 1, 1);
}

void decayAndFill(int const kf, TLorentzVector* b, double const weight, TClonesArray& daughters)
{
   pydecay->Decay(kf, b);
   pydecay->ImportParticles(&daughters);

   TLorentzVector kMom;
   TLorentzVector piMom;
   TLorentzVector pMom;
   TVector3 v00;

   int nTrk = daughters.GetEntriesFast();
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

   fill(kf, b, weight, kMom, piMom, pMom, v00);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& kMom, TLorentzVector const& piMom, TLorentzVector const& pMom, TVector3 v00)
{
   int const centrality = floor(nCent * gRandom->Rndm());

   TVector3 const vertex = getVertex(centrality);
   // smear primary vertex
   // float const sigmaVertex = sigmaVertexCent[cent];
   // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

   v00 += vertex;

   // smear momentum
   TLorentzVector const piRMom = smearMom(0, piMom);
   TLorentzVector const kRMom = smearMom(1, kMom);
   TLorentzVector const pRMom = smearMom(2, pMom);

   // smear position
   TVector3 const piRPos = smearPosData(0, vertex.z(), centrality, piRMom, v00);
   TVector3 const kRPos = smearPosData(1, vertex.z(), centrality, kRMom, v00);
   TVector3 const pRPos = smearPosData(2, vertex.z(), centrality, pRMom, v00);
   // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
   // TVector3 const pRPos = smearPos(pMom, pRMom, v00);

   // reconstruct
   TLorentzVector const rMom = pRMom + piRMom + kRMom;
   float const piDca = dca(piMom.Vect(), v00, vertex);
   float const kDca = dca(kMom.Vect(), v00, vertex);
   float const pDca = dca(pMom.Vect(), v00, vertex);
   float const piRDca = dca(piRMom.Vect(), piRPos, vertex);
   float const kRDca = dca(kRMom.Vect(), kRPos, vertex);
   float const pRDca = dca(pRMom.Vect(), pRPos, vertex);

   // smeared decay vertex
   TVector3 v0;
   // vertices between pairs
   TVector3 v012;
   TVector3 v023;
   TVector3 v013;

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

   bool const isPhft = matchHft(cent, pRMom);
   bool const isKhft = matchHft(cent, kRMom);
   bool const isPiHft = matchHft(cent, piRMom);

                       // save
   if (saveNt)
   {
     float arr[100];
     int iArr = 0;
     arr[iArr++] = cent;
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

     arr[iArr++] = resMass(pMom, kMom, piMom, decayMode);
     arr[iArr++] = resMass(pRMom, kRMom, piRMom, decayMode);

     arr[iArr++] = vDistMax;
     arr[iArr++] = vDist1;
     arr[iArr++] = vDist2;
     arr[iArr++] = vDist3;

     nt->Fill(arr);
   } // if (saveNt)

   // __________________________________________
   // using cuts
   // __________________________________________
   bool const PtCut = pRMom.Perp() > 0.3 && kRMom.Perp() > 0.3 && piRMom.Perp() > 0.3;
   bool const dcaCut = dca12 < 200 && dca23 < 200 && dca13 < 200;
   bool const dLengthCut = decayLength > 30;
   bool const cosThetaCut = cosTheta > 0.98;
   bool const HftCut = isPhft && isKhft && isPiHft;
   bool const EtaCut = TMath::Abs(kRMom.PseudoRapidity()) < 1. && TMath::Abs(pRMom.PseudoRapidity()) < 1. && TMath::Abs(piRMom.PseudoRapidity()) < 1.;
   
   if ( !( PtCut && dcaCut && dLengthCut && cosThetaCut && HftCut && EtaCut ) )
     return;

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

void getKinematics(TLorentzVector& b, double const mass)
{
   float const pt = gRandom->Uniform(momentumRange.first, momentumRange.second);
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
   nt = new TNtuple("nt", "", "cent:vx:vy:vz:"
                    "pid:m:pt:eta:y:phi:v0x:v0y:v0z:" // MC Ks
                    "rM:rPt:rEta:rY:rPhi:rV0x:rV0y:rV0z:" // Rc Ks
                    "dca12:decayLength:dcaD0ToPv:cosTheta:angle12:cosThetaStar:" // Rc pair
                    "p1M:p1Pt:p1Eta:p1Y:p1Phi:p1Dca:" // MC Pi +
                    "p1RM:p1RPt:p1REta:p1RY:p1RPhi:p1RVx:p1RVy:p1RVz:p1RDca:p1Tpc:" // Rc Pi +
                    "p2M:p2Pt:p2Eta:p2Y:p2Phi:p2Dca:" // MC Pi -
                    "p2RM:p2RPt:p2REta:p2RY:p2RPhi:p2RVx:p2RVy:p2RVz:p2RDca:p2Tpc:" // Rc Pi -
                    "p1Hft:p2Hft",BufSize);

   ntTMVA = new TNtuple("ntTMVA", "", "m:pt:charges:phi:eta:" // basic properties of Lambda_c
				      "dcaDaugthers31:dcaDaugthers23:dcaDaugthers12:" // dca daughters (pi-K, pi-p, p-K)
				      "cosPntAngle:dLength:" // cosTheta and decay Length
				      "p1pt:p2pt:p3pt:"
				      "p1Dca:p2Dca:p3Dca:" // daughters (K, p, pi)
				      "maxVertexDist");

}
//___________
void write()
{
   result->cd();
   nt->Write();
   result->Close();
}
