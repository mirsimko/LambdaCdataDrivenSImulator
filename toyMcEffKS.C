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
//============== main  program ==================
void toyMcEffKS(int npart = 100)
{
   gRandom->SetSeed();
   loadAllDistributions();
   bookObjects();

   pydecay = TPythia6Decayer::Instance();
   pydecay->Init();

   setDecayChannels(613); // K_S -> pi pi 

   TLorentzVector* b_d = new TLorentzVector;
   TClonesArray ptl("TParticle", 10);
   for (int ipart = 0; ipart < npart; ipart++)
   {
      if (!(ipart % 100000))
         cout << "____________ ipart = " << ipart << " ________________" << endl;

      getKinematics(*b_d, M_KS);

      decayAndFill(310, b_d, 1.0, ptl);
      if (ipart%1000 == 1) nt->AutoSave("SaveSelf");
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

   TLorentzVector p1Mom;
   TLorentzVector p2Mom;
   TVector3 v00;

   int nTrk = daughters.GetEntriesFast();
   for (int iTrk = 0; iTrk < nTrk; ++iTrk)
   {
     TParticle* ptl0 = (TParticle*)daughters.At(iTrk);

     switch (ptl0->GetPdgCode())
     {
       //Will only have pions
       case 211:
         ptl0->Momentum(p1Mom);
         v00.SetXYZ(ptl0->Vx() * 1000., ptl0->Vy() * 1000., ptl0->Vz() * 1000.); // converted to μm
         break;
       case -211:
         ptl0->Momentum(p2Mom);
         break;
       default:
         break;
     }
   }
   daughters.Clear();

   fill(kf, b, weight, p1Mom, p2Mom, v00);
}

void fill(int const kf, TLorentzVector* b, double weight, TLorentzVector const& p1Mom, TLorentzVector const& p2Mom, TVector3 v00)
{
   int const centrality = floor(nCent * gRandom->Rndm());

   TVector3 const vertex = getVertex(centrality);
   // smear primary vertex
   // float const sigmaVertex = sigmaVertexCent[cent];
   // TVector3 const vertex(gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex), gRandom->Gaus(0, sigmaVertex));

   v00 += vertex;

   // smear momentum
   TLorentzVector const p1RMom = smearMom(0, p1Mom);
   TLorentzVector const p2RMom = smearMom(0, p2Mom);

   // smear position
   TVector3 const p1RPos = smearPosData(0, vertex.z(), centrality, p1RMom, v00);
   TVector3 const p2RPos = smearPosData(0, vertex.z(), centrality, p2RMom, v00);
   // TVector3 const kRPos = smearPos(kMom, kRMom, v00);
   // TVector3 const pRPos = smearPos(pMom, pRMom, v00);

   // reconstruct
   TLorentzVector const rMom = p1RMom + p2RMom;
   float const p1Dca = dca(p1Mom.Vect(), v00, vertex);
   float const p2Dca = dca(p2Mom.Vect(), v00, vertex);
   float const p1RDca = dca(p1RMom.Vect(), p1RPos, vertex);
   float const p2RDca = dca(p2RMom.Vect(), p2RPos, vertex);

   TVector3 v0;
   float const dca12 = dca1To2(p1RMom.Vect(), p1RPos, p2RMom.Vect(), p2RPos, v0);
   float const decayLength = (v0 - vertex).Mag();
   float const dcaD0ToPv = dca(rMom.Vect(), v0, vertex);
   float const cosTheta = (v0 - vertex).Unit().Dot(rMom.Vect().Unit());
   float const angle12 = p1RMom.Vect().Angle(p2RMom.Vect());

   TLorentzVector p1RMomRest = p1RMom;
   TVector3 beta;
   beta.SetMagThetaPhi(rMom.Beta(), rMom.Theta(), rMom.Phi());
   p1RMomRest.Boost(-beta);
   float const cosThetaStar = rMom.Vect().Unit().Dot(p1RMomRest.Vect().Unit());

   // save
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
   arr[iArr++] = v0.X();
   arr[iArr++] = v0.Y();
   arr[iArr++] = v0.Z();

   arr[iArr++] = dca12;
   arr[iArr++] = decayLength;
   arr[iArr++] = dcaD0ToPv;
   arr[iArr++] = cosTheta;
   arr[iArr++] = angle12;
   arr[iArr++] = cosThetaStar;

   arr[iArr++] = p1Mom.M();
   arr[iArr++] = p1Mom.Perp();
   arr[iArr++] = p1Mom.PseudoRapidity();
   arr[iArr++] = p1Mom.Rapidity();
   arr[iArr++] = p1Mom.Phi();
   arr[iArr++] = p1Dca;

   arr[iArr++] = p1RMom.M();
   arr[iArr++] = p1RMom.Perp();
   arr[iArr++] = p1RMom.PseudoRapidity();
   arr[iArr++] = p1RMom.Rapidity();
   arr[iArr++] = p1RMom.Phi();
   arr[iArr++] = p1RPos.X();
   arr[iArr++] = p1RPos.Y();
   arr[iArr++] = p1RPos.Z();
   arr[iArr++] = p1RDca;
   arr[iArr++] = tpcReconstructed(0,1,centrality,p1RMom);

   arr[iArr++] = p2Mom.M();
   arr[iArr++] = p2Mom.Perp();
   arr[iArr++] = p2Mom.PseudoRapidity();
   arr[iArr++] = p2Mom.Rapidity();
   arr[iArr++] = p2Mom.Phi();
   arr[iArr++] = p2Dca;

   arr[iArr++] = p2RMom.M();
   arr[iArr++] = p2RMom.Perp();
   arr[iArr++] = p2RMom.PseudoRapidity();
   arr[iArr++] = p2RMom.Rapidity();
   arr[iArr++] = p2RMom.Phi();
   arr[iArr++] = p2RPos.X();
   arr[iArr++] = p2RPos.Y();
   arr[iArr++] = p2RPos.Z();
   arr[iArr++] = p2RDca;
   arr[iArr++] = tpcReconstructed(0,-1,centrality,p2RMom);

   arr[iArr++] = matchHft(1, vertex.z(), centrality, p1RMom);
   arr[iArr++] = matchHft(0, vertex.z(), centrality, p2RMom);

   nt->Fill(arr);
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
}
//___________
void write()
{
   result->cd();
   nt->Write();
   result->Close();
}
