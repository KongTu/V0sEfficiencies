#include <iostream>
#include <sstream>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include <vector>
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include <string>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHist.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

double ks_YieldCal( TH1F* inputHist ){

//Starting to do RooFit:

    RooRealVar x("x","mass",0.44,0.56);
    RooDataHist data("data","dataset",x, inputHist );
    RooPlot* xframe = x.frame(240);

    data.plotOn(xframe, Name("data"));

    RooRealVar mean("mean","mean",0.50,0.49,0.51);
    RooRealVar sigma1("sigma1","sigma1",0.003,0.001,0.01);
    RooRealVar sigma2("sigma2","sigma2",0.003,0.001,0.01);
    RooRealVar sig1("sig1","signal1",10,0,10000000);
    RooRealVar sig2("sig2","signal2",10,0,10000000);
    RooRealVar a("a","a",0,-100000,100000);
    RooRealVar b("b","b",0,-100000,100000);
    RooRealVar cp("cp","cp",0,-100000,100000);
    RooRealVar d("d","d",0,-100000,100000);

    RooRealVar f("f","f",0,-100000,100000);
    RooRealVar g("g","g",0,-100000,100000);
    RooRealVar h("h","h",0,-100000,100000);
    RooRealVar k("k","k",0,-100000,100000);

    RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
    RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
    RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
    RooRealVar polysig("polysig","polysig",10,0,10000000);
    RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

    x.setRange("cut",0.45,0.54);

    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));

    sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
    sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));

    xframe->Draw();

    double chi2  = xframe->chiSquare("sum","data");
    double meanf  = mean.getVal();
    double meanfe  = mean.getError();
    double sigmaf1  = sigma1.getVal();
    double sigmaf2  = sigma2.getVal();
    double bkgf  = polysig.getVal();
    double sigf1  = sig1.getVal();
    double sigf2  = sig2.getVal();
    double sigwf1  = sigf1 /(sigf1 +sigf2 );
    double sigwf2  = sigf2 /(sigf1 +sigf2 );
    double c1 = a.getVal();
    double c2 = b.getVal();

    double sigmaf  = sqrt(sigmaf1 **2*sigwf1  + sigmaf2 **2*sigwf2 );
    double massmin  = meanf  - 2*sigmaf ;
    double massmax  = meanf  + 2*sigmaf ;

    int nmin  =  inputHist  ->GetXaxis()->FindBin(massmin );
    int nmax  =  inputHist  ->GetXaxis()->FindBin(massmax );
    int anmin  =  inputHist  ->GetXaxis()->FindBin(0.44);
    int anmax  =  inputHist  ->GetXaxis()->FindBin(0.56);

    double awyh1  =  inputHist  ->Integral(anmin ,nmin );
    double awyh2  =  inputHist  ->Integral(nmax ,anmax );
    double awyh  = awyh1  + awyh2 ;
    double totyh  =  inputHist  ->Integral(nmin ,nmax );

    x.setRange("cut",massmin ,massmax );
    RooAbsReal* ibkg  = poly.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig1  = gaus1.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig2  = gaus2.createIntegral(x,NormSet(x),Range("cut"));
    double ibkgf  = ibkg ->getVal();
    double bkgfe  = polysig.getError();
    double isig1f  = isig1 ->getVal();
    double isig2f  = isig2 ->getVal();

    double bkgy  = ibkgf *bkgf ;
    double bkgye  = ibkgf *bkgfe ;
    double sigy1  = isig1f *sigf1 ;
    double sigy2  = isig2f *sigf2 ;
    double sigy  = sigy1  + sigy2 ;
    double toty  = bkgy  + sigy ;

    double abkgy  = (1-ibkgf )*bkgf ;
    double asigy1  = (1-isig1f )*sigf1 ;
    double asigy2  = (1-isig2f )*sigf2 ;
    double asigy  = asigy1  + asigy2 ;
    double awy  = abkgy  + asigy ;

    double sigfrac  = sigy /toty ;
    double bkgfrac  = bkgy /toty ;

    double sigyh  = totyh  - bkgy ;
    double sigfrach  = sigyh /totyh ;
    double bkgfrach  = bkgy /totyh ;

    double signif  = sigyh / sqrt( totyh );

    return sigyh;
    
}

//Similarly for Lambda;
double la_YieldCal( TH1F* inputHist ){

//starting to do RooFit:

    RooRealVar x("x","mass",1,1.2);
    RooDataHist data("data","dataset",x, inputHist );
    RooPlot* xframe = x.frame(160);

    data.plotOn(xframe, Name("data"));

    RooRealVar mean("mean","mean",1.115,1.11,1.12);
    RooRealVar sigma1("sigma1","sigma1",0.003,0.001,0.01);
    RooRealVar sigma2("sigma2","sigma2",0.003,0.001,0.01);
    RooRealVar sig1("sig1","signal1",10,0,10000000);
    RooRealVar sig2("sig2","signal2",10,0,10000000);
    RooRealVar a("a","a",0,-100000,100000);
    RooRealVar b("b","b",0,-100000,100000);
    RooRealVar cp("cp","cp",0,-100000,100000);
    RooRealVar d("d","d",0,-100000,100000);

    RooRealVar f("f","f",0,-100000,100000);
    RooRealVar g("g","g",0,-100000,100000);
    RooRealVar h("h","h",0,-100000,100000);
    RooRealVar k("k","k",0,-100000,100000);

    RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
    RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
    RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
    RooRealVar polysig("polysig","polysig",10,0,10000000);
    RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

    x.setRange("cut",1.10,1.14);

    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));
    sum.fitTo(data,Range("cut"));

    sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
    sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));

    xframe->Draw();

    double chi2  = xframe->chiSquare("sum","data");
    double meanf  = mean.getVal();
    double meanfe  = mean.getError();
    double sigmaf1  = sigma1.getVal();
    double sigmaf2  = sigma2.getVal();
    double bkgf  = polysig.getVal();
    double sigf1  = sig1.getVal();
    double sigf2  = sig2.getVal();
    double sigwf1  = sigf1 /(sigf1 +sigf2 );
    double sigwf2  = sigf2 /(sigf1 +sigf2 );
    double c1 = a.getVal();
    double c2 = b.getVal();

    double sigmaf  = sqrt(sigmaf1 **2*sigwf1  + sigmaf2 **2*sigwf2 );
    double massmin  = meanf  - 2*sigmaf ;
    double massmax  = meanf  + 2*sigmaf ;

    int nmin  = inputHist->GetXaxis()->FindBin(massmin );
    int nmax  = inputHist->GetXaxis()->FindBin(massmax );
    int anmin  = inputHist->GetXaxis()->FindBin(1.0);
    int anmax  = inputHist->GetXaxis()->FindBin(1.2);

    double awyh1  = inputHist->Integral(anmin ,nmin );
    double awyh2  = inputHist->Integral(nmax ,anmax );
    double awyh  = awyh1  + awyh2 ;
    double totyh  = inputHist->Integral(nmin ,nmax );

    x.setRange("cut",massmin ,massmax );
    RooAbsReal* ibkg  = poly.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig1  = gaus1.createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig2  = gaus2.createIntegral(x,NormSet(x),Range("cut"));
    double ibkgf  = ibkg ->getVal();
    double bkgfe  = polysig.getError();
    double isig1f  = isig1 ->getVal();
    double isig2f  = isig2 ->getVal();

    double bkgy  = ibkgf *bkgf ;
    double bkgye  = ibkgf *bkgfe ;
    double sigy1  = isig1f *sigf1 ;
    double sigy2  = isig2f *sigf2 ;
    double sigy  = sigy1  + sigy2 ;
    double toty  = bkgy  + sigy ;

    double abkgy  = (1-ibkgf )*bkgf ;
    double asigy1  = (1-isig1f )*sigf1 ;
    double asigy2  = (1-isig2f )*sigf2 ;
    double asigy  = asigy1  + asigy2 ;
    double awy  = abkgy  + asigy ;

    double sigfrac  = sigy /toty ;
    double bkgfrac  = bkgy /toty ;

    double sigyh  = totyh  - bkgy ;
    double sigfrach  = sigyh /totyh ;
    double bkgfrach  = bkgy /totyh ;

    //double signif  = sigyh / sqrt( totyh );

    return sigyh;

}

void jetConeEfficiency_Regit_check1(){

    //define the path to data:
	TFile* file = new TFile("~/2014Research/ROOT_file/HiRegit/RegitPYTHIA80_JUNE29_ppCuts_updates3_2014.root");
    TTree* theTree = ( TTree* )file->Get("v0analyzerNew/PFJet"); 
    TTree* KshortTree = ( TTree* )file->Get("v0analyzerNew/v0_Kshort");
    TTree* LambdaTree = ( TTree* )file->Get("v0analyzerNew/v0_Lambda"); 
    TTree* genParticleTree = ( TTree* )file->Get("v0analyzerNew/GenParticle");

    //v0analyzerNew is for the V0 inputTag from the pp standard tracking;
    //v0analyzerHI is for Regit tracking;

    //RECO variables:

        float jet_pt[3000];
        float jet_eta[3000];
        float jet_phi[3000];
        int nJets;

        int nK0short;
        int nLambda;

        float ks_eta[3000];
        float ks_phi[3000];
        float ks_pt[3000];
        float ks_px[3000];
        float ks_py[3000];
        float ks_pz[3000];
        float ks_dau1_dzos[3000];
        float ks_dau2_dzos[3000];
        float ks_dau1_dxyos[3000];
        float ks_dau2_dxyos[3000];
        float ks_dau1_nhit[3000];
        float ks_dau2_nhit[3000];
        float ks_dau1_eta[3000];
        float ks_dau2_eta[3000];
        float ks_dau1_phi[3000];
        float ks_dau2_phi[3000];
        float ks_agl[3000];
        float ks_dlos[3000];
        float ks_mass[3000];
       
        float la_px[3000];
        float la_py[3000];
        float la_pz[3000];
        float la_eta[3000];
        float la_phi[3000];
        float la_pt[3000];
        float la_dau1_dzos[3000];
        float la_dau2_dzos[3000];
        float la_dau1_dxyos[3000];
        float la_dau2_dxyos[3000];
        float la_dau1_nhit[3000];
        float la_dau2_nhit[3000];
        float la_dau1_eta[3000];
        float la_dau2_eta[3000];
        float la_dau1_phi[3000];
        float la_dau2_phi[3000];
        float la_agl[3000];
        float la_dlos[3000];
        float la_mass[3000];

        theTree->SetBranchAddress( "nJets", &nJets );
        theTree->SetBranchAddress( "jet_pt", &jet_pt );
        theTree->SetBranchAddress( "jet_eta", &jet_eta );
        theTree->SetBranchAddress( "jet_phi", &jet_phi );
        

        KshortTree->SetBranchAddress( "N", &nK0short );
        KshortTree->SetBranchAddress( "ks_phi", &ks_phi );
        KshortTree->SetBranchAddress( "ks_eta", &ks_eta );
        KshortTree->SetBranchAddress( "ks_pt", &ks_pt );
        KshortTree->SetBranchAddress( "ks_dau1_dzos", &ks_dau1_dzos );
        KshortTree->SetBranchAddress( "ks_dau2_dzos", &ks_dau2_dzos );
        KshortTree->SetBranchAddress( "ks_dau1_dxyos", &ks_dau1_dxyos );
        KshortTree->SetBranchAddress( "ks_dau2_dxyos", &ks_dau2_dxyos );
        KshortTree->SetBranchAddress( "ks_dau1_nhit", &ks_dau1_nhit );
        KshortTree->SetBranchAddress( "ks_dau2_nhit", &ks_dau2_nhit );
        KshortTree->SetBranchAddress( "ks_dau1_eta", &ks_dau1_eta );
        KshortTree->SetBranchAddress( "ks_dau2_eta", &ks_dau2_eta );
        KshortTree->SetBranchAddress( "ks_dau1_phi", &ks_dau1_phi );
        KshortTree->SetBranchAddress( "ks_dau2_phi", &ks_dau2_phi );
        KshortTree->SetBranchAddress( "ks_agl", &ks_agl );
        KshortTree->SetBranchAddress( "ks_dlos", &ks_dlos );    
        KshortTree->SetBranchAddress( "ks_px", &ks_px );
        KshortTree->SetBranchAddress( "ks_py", &ks_py );
        KshortTree->SetBranchAddress( "ks_pz", &ks_pz );
        KshortTree->SetBranchAddress( "ks_mass", &ks_mass );
        
        LambdaTree->SetBranchAddress( "N1", &nLambda );
        LambdaTree->SetBranchAddress( "la_px", &la_px );
        LambdaTree->SetBranchAddress( "la_py", &la_py );
        LambdaTree->SetBranchAddress( "la_pz", &la_pz );
        LambdaTree->SetBranchAddress( "la_phi", &la_phi );
        LambdaTree->SetBranchAddress( "la_eta", &la_eta );
        LambdaTree->SetBranchAddress( "la_pt", &la_pt );
        LambdaTree->SetBranchAddress( "la_dau1_dzos", &la_dau1_dzos );
        LambdaTree->SetBranchAddress( "la_dau2_dzos", &la_dau2_dzos );
        LambdaTree->SetBranchAddress( "la_dau1_dxyos", &la_dau1_dxyos );
        LambdaTree->SetBranchAddress( "la_dau2_dxyos", &la_dau2_dxyos );
        LambdaTree->SetBranchAddress( "la_dau1_nhit", &la_dau1_nhit );
        LambdaTree->SetBranchAddress( "la_dau2_nhit", &la_dau2_nhit );
        LambdaTree->SetBranchAddress( "la_dau1_eta", &la_dau1_eta );
        LambdaTree->SetBranchAddress( "la_dau2_eta", &la_dau2_eta );
        LambdaTree->SetBranchAddress( "la_dau1_phi", &la_dau1_phi );
        LambdaTree->SetBranchAddress( "la_dau2_phi", &la_dau2_phi );
        LambdaTree->SetBranchAddress( "la_agl", &la_agl );
        LambdaTree->SetBranchAddress( "la_dlos", &la_dlos );
        LambdaTree->SetBranchAddress( "la_mass", &la_mass );
    
    //GEN variables:

        int nGen;
        int pdg[3000];
        float genP_pt[3000];
        float genP_eta[3000];
        float genP_phi[3000];
        int mid[3000];

        genParticleTree->SetBranchAddress("nGen", &nGen );
        genParticleTree->SetBranchAddress("pdg", &pdg );
        genParticleTree->SetBranchAddress("genP_mid", &mid);
        genParticleTree->SetBranchAddress("genP_pt", &genP_pt );
        genParticleTree->SetBranchAddress("genP_eta", &genP_eta );
        genParticleTree->SetBranchAddress("genP_phi", &genP_phi );

        theTree->AddFriend( KshortTree );
        theTree->AddFriend( LambdaTree );
        theTree->AddFriend( genParticleTree );

    int maxEvent = theTree->GetEntries();
    int nEvent = 1000;
 
//==================
//define Histograms:
//==================
   
        TH1F* k0Hist[10];
        TH1F* laHist[10];
    
        stringstream KSname;
        stringstream LAname;

        for (int nhist = 0; nhist < 10; nhist++){

            KSname << "K^{0}_{s} dauTrack #DeltaR w.r.t jet axis_";
            KSname << nhist + 1;

            LAname << "#Lambda/#bar{#Lambda} dauTrack #DeltaR w.r.t jet axis__";
            LAname << nhist + 1;

            k0Hist[nhist] = new TH1F( KSname.str().c_str(),KSname.str().c_str(),200,0,1);
            k0Hist[nhist]->SetXTitle("#DeltaR ");
            k0Hist[nhist]->SetYTitle("#counts");

            laHist[nhist] = new TH1F( LAname.str().c_str(),LAname.str().c_str(),200,0,1);
            laHist[nhist]->SetXTitle("#DeltaR");
            laHist[nhist]->SetYTitle("#counts");
            
            KSname.str("");
            LAname.str("");
        
        }

        for (int it = 0; it < maxEvent; it++ ){

            theTree->GetEntry(it);

            cout << "let's see the nEvent: " << it << endl;

                for (int i = 0; i < nJets; i++){

                    if ( jet_pt[i] < 60 || TMath::Abs( jet_eta[i]) > 2.0 ) continue;

                    for (int x = 0; x < nK0short; x++){

                        float delta_ks_eta = (jet_eta[i]) - (ks_eta[x]);
                        float delta_ks_phi = (jet_phi[i]) - (ks_phi[x]);

                        if ( delta_ks_phi > 3.14 ){

                            float KSconeSize = sqrt((6.28 - delta_ks_phi)*(6.28 - delta_ks_phi)+(delta_ks_eta)*(delta_ks_eta));
                        }
                        else if ( delta_ks_phi < -3.14 ){

                            float KSconeSize = sqrt((6.28 + delta_ks_phi)*(6.28 + delta_ks_phi)+(delta_ks_eta)*(delta_ks_eta));

                        }
                        else{

                            float KSconeSize = sqrt((delta_ks_phi)*(delta_ks_phi)+(delta_ks_eta)*(delta_ks_eta));

                        }



                    if ( KSconeSize < 0.3 ){

                        float delta_ks1_eta = (jet_eta[i]) - (ks_dau1_eta[x]);
                        float delta_ks1_phi = (jet_phi[i]) - (ks_dau1_phi[x]);

                        if ( delta_ks1_phi > 3.14 ){

                            float KS1coneSize1 = sqrt((6.28 - delta_ks1_phi)*(6.28 - delta_ks1_phi)+(delta_ks1_eta)*(delta_ks1_eta));
                        }
                        else if ( delta_ks1_phi < -3.14 ){

                            float KS1coneSize1 = sqrt((6.28 + delta_ks1_phi)*(6.28 + delta_ks1_phi)+(delta_ks1_eta)*(delta_ks1_eta));

                        }
                        else{

                            float KS1coneSize = sqrt((delta_ks1_phi)*(delta_ks1_phi)+(delta_ks1_eta)*(delta_ks1_eta));
                        }

                        float delta_ks2_eta = (jet_eta[i]) - (ks_dau2_eta[x]);
                        float delta_ks2_phi = (jet_phi[i]) - (ks_dau2_phi[x]);

                        if ( delta_ks2_phi > 3.14 ){

                            float KS2coneSize2 = sqrt((6.28 - delta_ks2_phi)*(6.28 - delta_ks2_phi)+(delta_ks2_eta)*(delta_ks2_eta));
                        }
                        else if ( delta_ks2_phi < -3.14 ){

                            float KS2coneSize2 = sqrt((6.28 + delta_ks2_phi)*(6.28 + delta_ks2_phi)+(delta_ks2_eta)*(delta_ks2_eta));

                        }
                        else{

                            float KS2coneSize = sqrt((delta_ks2_phi)*(delta_ks2_phi)+(delta_ks2_eta)*(delta_ks2_eta));
                        }


                                if ( ks_pt[x] > 0.7 && ks_pt[x] < 1.0 ){

                                    k0Hist[0]->Fill( KS1coneSize );
                                    k0Hist[0]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 1.0 && ks_pt[x] < 1.4 ){

                                    k0Hist[1]->Fill( KS1coneSize );
                                    k0Hist[1]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 1.4 && ks_pt[x] < 1.8 ){

                                    k0Hist[2]->Fill( KS1coneSize );
                                    k0Hist[2]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 1.8 && ks_pt[x] < 2.2 ){

                                    k0Hist[3]->Fill( KS1coneSize );
                                    k0Hist[3]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 2.2 && ks_pt[x] < 2.8 ){

                                    k0Hist[4]->Fill( KS1coneSize );
                                    k0Hist[4]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 2.8 && ks_pt[x] < 3.6 ){

                                    k0Hist[5]->Fill( KS1coneSize );
                                    k0Hist[5]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 3.6 && ks_pt[x] < 4.6 ){

                                    k0Hist[6]->Fill( KS1coneSize );
                                    k0Hist[6]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 4.6 && ks_pt[x] < 6.0 ){

                                    k0Hist[7]->Fill( KS1coneSize );
                                    k0Hist[7]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 6.0 && ks_pt[x] < 9.0 ){

                                    k0Hist[8]->Fill( KS1coneSize );
                                    k0Hist[8]->Fill( KS2coneSize );
                                }

                                if ( ks_pt[x] > 9.0 && ks_pt[x] < 12.0 ){

                                    k0Hist[9]->Fill( KS1coneSize );
                                    k0Hist[9]->Fill( KS2coneSize );
                                } 
                        }   

                }
                    
                    //Finding Lambda in Jet cone:

                    for(int x = 0; x < nLambda; x++){

                        float delta_la_eta = (jet_eta[i]) - (la_eta[x]);
                        float delta_la_phi = (jet_phi[i]) - (la_phi[x]);

                        if ( delta_la_phi > 3.14 ){

                            float LAconeSize = sqrt((6.28 - delta_la_phi)*(6.28 - delta_la_phi)+(delta_la_eta)*(delta_la_eta));
                        }
                        else if ( delta_la_phi < -3.14 ){

                            float LAconeSize = sqrt((6.28 + delta_la_phi)*(6.28 + delta_la_phi)+(delta_la_eta)*(delta_la_eta));

                        }
                        else{

                            float LAconeSize = sqrt((delta_la_phi)*(delta_la_phi)+(delta_la_eta)*(delta_la_eta));

                        }

                    if ( LAconeSize < 0.3 ){

                        float delta_la1_eta = (jet_eta[i]) - (la_dau1_eta[x]);
                        float delta_la1_phi = (jet_phi[i]) - (la_dau1_phi[x]);

                        if ( delta_la1_phi > 3.14 ){

                            float LA1coneSize = sqrt((6.28 - delta_la1_phi)*(6.28 - delta_la1_phi)+(delta_la1_eta)*(delta_la1_eta));
                        }
                        else if ( delta_la1_phi < -3.14 ){

                            float LA1coneSize = sqrt((6.28 + delta_la1_phi)*(6.28 + delta_la1_phi)+(delta_la1_eta)*(delta_la1_eta));

                        }
                        else{

                            float LA1coneSize = sqrt((delta_la1_phi)*(delta_la1_phi)+(delta_la1_eta)*(delta_la1_eta));
                        }

                        float delta_la2_eta = (jet_eta[i]) - (la_dau2_eta[x]);
                        float delta_la2_phi = (jet_phi[i]) - (la_dau2_phi[x]);

                        if ( delta_la2_phi > 3.14 ){

                            float LA2coneSize = sqrt((6.28 - delta_la2_phi)*(6.28 - delta_la2_phi)+(delta_la2_eta)*(delta_la2_eta));
                        }
                        else if ( delta_la2_phi < -3.14 ){

                            float LA2coneSize = sqrt((6.28 + delta_la2_phi)*(6.28 + delta_la2_phi)+(delta_la2_eta)*(delta_la2_eta));

                        }
                        else{

                            float LA2coneSize = sqrt((delta_la2_phi)*(delta_la2_phi)+(delta_la2_eta)*(delta_la2_eta));
                        }


                                if ( la_pt[x] > 0.7 && la_pt[x] < 1.0 ){

                                    laHist[0]->Fill( LA1coneSize );
                                    laHist[0]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 1.0 && la_pt[x] < 1.4 ){

                                    laHist[1]->Fill( LA1coneSize );
                                    laHist[1]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 1.4 && la_pt[x] < 1.8 ){

                                    laHist[2]->Fill( LA1coneSize );
                                    laHist[2]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 1.8 && la_pt[x] < 2.2 ){

                                    laHist[3]->Fill( LA1coneSize );
                                    laHist[3]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 2.2 && la_pt[x] < 2.8 ){

                                    laHist[4]->Fill( LA1coneSize );
                                    laHist[4]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 2.8 && la_pt[x] < 3.6 ){

                                    laHist[5]->Fill( LA1coneSize );
                                    laHist[5]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 3.6 && la_pt[x] < 4.6 ){

                                    laHist[6]->Fill( LA1coneSize );
                                    laHist[6]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 4.6 && la_pt[x] < 6.0 ){

                                    laHist[7]->Fill( LA1coneSize );
                                    laHist[7]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 6.0 && la_pt[x] < 9.0 ){

                                    laHist[8]->Fill( LA1coneSize );
                                    laHist[8]->Fill( LA2coneSize );
                                }

                                if ( la_pt[x] > 9.0 && la_pt[x] < 12.0 ){

                                    laHist[9]->Fill( LA1coneSize );
                                    laHist[9]->Fill( LA2coneSize );
                                } 

                    }
       
                } 
                    
            }  
                
        }

//pT range;

    char* pTrange[10];
    pTrange[0] = "0.7 < pT < 1.0";
    pTrange[1] = "1.0 < pT < 1.4";
    pTrange[2] = "1.4 < pT < 1.8";
    pTrange[3] = "1.8 < pT < 2.2";
    pTrange[4] = "2.2 < pT < 2.8";
    pTrange[5] = "2.8 < pT < 3.6";
    pTrange[6] = "3.6 < pT < 4.6";
    pTrange[7] = "4.6 < pT < 6.0";
    pTrange[8] = "6.0 < pT < 9.0";
    pTrange[9] = "9.0 < pT < 12.0";

    TCanvas* c1 = new TCanvas();
    c1->Print("k0_Daughter_DeltaR_0.3.pdf[");
        for (int r = 0; r < 10; r++){

            k0Hist[r]->SetStats(kFALSE);
            k0Hist[r]->Draw();
            TLatex* tlx=new TLatex(0.6, 0.8, pTrange[r] );
            tlx->SetNDC(kTRUE); // <- use NDC coordinate
            tlx->SetTextSize(0.06);
            tlx->Draw("same");
            c1->Print("k0_Daughter_DeltaR_0.3.pdf");

        }
    c1->Print("k0_Daughter_DeltaR_0.3.pdf]");

    TCanvas* c2 = new TCanvas();
    c2->Print("la_Daughter_DeltaR_0.3.pdf[");
    for (int w = 0; w < 10; w++){

        laHist[w]->SetStats(kFALSE);
        laHist[w]->Draw();
        TLatex* tlx1=new TLatex(0.6, 0.8, pTrange[w] );
        tlx1->SetNDC(kTRUE); // <- use NDC coordinate
        tlx1->SetTextSize(0.06);
        tlx1->Draw("same");
        c2->Print("la_Daughter_DeltaR_0.3.pdf");

    }
    c2->Print("la_Daughter_DeltaR_0.3.pdf]");
    
}
