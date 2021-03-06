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

void HijingEfficiency(){

    TFile* file = new TFile("~/Desktop/HIJING_June21_2014.root");
    TTree* theTree = ( TTree* )file->Get("ana/v0_Kshort");
    theTree->AddFriend("ana/v0_Lambda");
    theTree->AddFriend("ana/GenParticle");

//RECO vairiables:

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
        float la_agl[3000];
        float la_dlos[3000];
        float la_mass[3000];

        theTree->SetBranchAddress( "N", &nK0short );
        theTree->SetBranchAddress( "ks_phi", &ks_phi );
        theTree->SetBranchAddress( "ks_eta", &ks_eta );
        theTree->SetBranchAddress( "ks_pt", &ks_pt );
        theTree->SetBranchAddress( "ks_dau1_dzos", &ks_dau1_dzos );
        theTree->SetBranchAddress( "ks_dau2_dzos", &ks_dau2_dzos );
        theTree->SetBranchAddress( "ks_dau1_dxyos", &ks_dau1_dxyos );
        theTree->SetBranchAddress( "ks_dau2_dxyos", &ks_dau2_dxyos );
        theTree->SetBranchAddress( "ks_dau1_nhit", &ks_dau1_nhit );
        theTree->SetBranchAddress( "ks_dau2_nhit", &ks_dau2_nhit );
        theTree->SetBranchAddress( "ks_agl", &ks_agl );
        theTree->SetBranchAddress( "ks_dlos", &ks_dlos );    
        theTree->SetBranchAddress( "ks_px", &ks_px );
        theTree->SetBranchAddress( "ks_py", &ks_py );
        theTree->SetBranchAddress( "ks_pz", &ks_pz );
        theTree->SetBranchAddress( "ks_mass", &ks_mass );
        
        theTree->SetBranchAddress( "N1", &nLambda );
        theTree->SetBranchAddress( "la_px", &la_px );
        theTree->SetBranchAddress( "la_py", &la_py );
        theTree->SetBranchAddress( "la_pz", &la_pz );
        theTree->SetBranchAddress( "la_phi", &la_phi );
        theTree->SetBranchAddress( "la_eta", &la_eta );
        theTree->SetBranchAddress( "la_pt", &la_pt );
        theTree->SetBranchAddress( "la_dau1_dzos", &la_dau1_dzos );
        theTree->SetBranchAddress( "la_dau2_dzos", &la_dau2_dzos );
        theTree->SetBranchAddress( "la_dau1_dxyos", &la_dau1_dxyos );
        theTree->SetBranchAddress( "la_dau2_dxyos", &la_dau2_dxyos );
        theTree->SetBranchAddress( "la_dau1_nhit", &la_dau1_nhit );
        theTree->SetBranchAddress( "la_dau2_nhit", &la_dau2_nhit );
        theTree->SetBranchAddress( "la_agl", &la_agl );
        theTree->SetBranchAddress( "la_dlos", &la_dlos );
        theTree->SetBranchAddress( "la_mass", &la_mass );

//GEN variables:

        int nGen;
        int pdg[50000];
        int mid[50000];
        float genP_pt[50000];
        float genP_eta[50000];
        float genP_phi[50000];

        theTree->SetBranchAddress("nGen", &nGen );
        theTree->SetBranchAddress("pdg", &pdg );
        theTree->SetBranchAddress("genP_mid", &mid );
        theTree->SetBranchAddress("genP_pt", &genP_pt );
        theTree->SetBranchAddress("genP_eta", &genP_eta );
        theTree->SetBranchAddress("genP_phi", &genP_phi );
      
    int maxEvent = theTree->GetEntries();
    int nEvent = 1000000;

//define Histograms:

//---------------
//define pt bins:
//---------------

        double ptbins[] = {0.7,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0};
        double ptbinwidth[10] = {0.3,0.4,0.4,0.6,0.6,0.8,1.0,1.4,3.0,3.0};

        TH1F* h3 = new TH1F("h3","gen_K0short pT",10,ptbins);
        TH1F* h4 = new TH1F("h4","gen_Lambda pT",10,ptbins);

        TH1F* k0Hist[10];
        TH1F* laHist[10];

        TH1F* genk0Hist[10];
        TH1F* genlaHist[10];
       
        stringstream KSname;
        stringstream LAname;

        stringstream genKSname;
        stringstream genLAname;

        for (int nhist = 0; nhist < 10; nhist++){

            KSname << "k0Hist_";
            KSname << nhist + 1;

            LAname << "laHist_";
            LAname << nhist + 1;

            genKSname << "genk0PtHist_";
            genKSname << nhist + 1;

            genLAname << "genlaPtHist_";
            genLAname << nhist + 1;

            k0Hist[nhist] = new TH1F( KSname.str().c_str(),KSname.str().c_str(),240,0.44,0.56);
            k0Hist[nhist]->SetXTitle("mass (GeV/c^{2})");
            k0Hist[nhist]->SetYTitle("#counts");

            laHist[nhist] = new TH1F( LAname.str().c_str(),LAname.str().c_str(),240,1.08,1.16);
            laHist[nhist]->SetXTitle("mass (GeV/c^{2})");
            laHist[nhist]->SetYTitle("#counts");

            genk0Hist[nhist] = new TH1F( genKSname.str().c_str(),genKSname.str().c_str(),240,0.44,0.56);
            genk0Hist[nhist]->SetXTitle("mass (GeV/c^{2})");
            genk0Hist[nhist]->SetYTitle("#counts");

            genlaHist[nhist] = new TH1F( genLAname.str().c_str(),genLAname.str().c_str(),240,1.08,1.16);
            genlaHist[nhist]->SetXTitle("mass (GeV/c^{2})");
            genlaHist[nhist]->SetYTitle("#counts");
            

            KSname.str("");
            LAname.str("");

            genKSname.str("");
            genLAname.str("");
        
        }


    for (int it = 0; it < nEvent; it++){

        theTree->GetEntry(it);

        cout << "#events: " << it << endl;

        for (int x = 0; x < nK0short; x++){

            if ( ks_pt[x] > 0.7 && ks_pt[x] < 1.0 ){

                k0Hist[0]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 1.0 && ks_pt[x] < 1.4 ){

                k0Hist[1]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 1.4 && ks_pt[x] < 1.8 ){

                k0Hist[2]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 1.8 && ks_pt[x] < 2.2 ){

                k0Hist[3]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 2.2 && ks_pt[x] < 2.8 ){

                k0Hist[4]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 2.8 && ks_pt[x] < 3.6 ){

                k0Hist[5]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 3.6 && ks_pt[x] < 4.6 ){

                k0Hist[6]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 4.6 && ks_pt[x] < 6.0 ){

                k0Hist[7]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 6.0 && ks_pt[x] < 9.0 ){

                k0Hist[8]->Fill( ks_mass[x] );
            }

            if ( ks_pt[x] > 9.0 && ks_pt[x] < 12.0 ){

                k0Hist[9]->Fill( ks_mass[x] );
            }
        }

        for (int x = 0; x < nLambda; x++){

            if ( la_pt[x] > 0.7 && la_pt[x] < 1.0 ){

                laHist[0]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 1.0 && la_pt[x] < 1.4 ){

                laHist[1]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 1.4 && la_pt[x] < 1.8 ){

                laHist[2]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 1.8 && la_pt[x] < 2.2 ){

                laHist[3]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 2.2 && la_pt[x] < 2.8 ){

                laHist[4]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 2.8 && la_pt[x] < 3.6 ){

                laHist[5]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 3.6 && la_pt[x] < 4.6 ){

                laHist[6]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 4.6 && la_pt[x] < 6.0 ){

                laHist[7]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 6.0 && la_pt[x] < 9.0 ){

                laHist[8]->Fill( la_mass[x] );
            }

            if ( la_pt[x] > 9.0 && la_pt[x] < 12.0 ){

                laHist[9]->Fill( la_mass[x] );
            }
        }

        for (int y = 0; y < nGen; y++){

            if ( pdg[y] == 310 ){

                h3->Fill( genP_pt[y] );

            }

            if ( TMath::Abs(pdg[y]) == 3122 && TMath::Abs(mid[x]) != 3322 && TMath::Abs(mid[x]) != 3312 && TMath::Abs(mid[x]) != 3324 && TMath::Abs(mid[x]) != 3314 && TMath::Abs(mid[x]) != 3334 ){

                h4->Fill( genP_pt[y] );
            }

        }

    }

    double ks_yield[10];
    double la_yield[10];
    double genKS_yield[10];
    double genLA_yield[10];

//Obtain the yield from the pt distribution for genKS AND genLA:

    for (int i = 0; i < 10; i++){

        genKS_yield[i] = h3->GetBinContent(i+1);
        genLA_yield[i] = h4->GetBinContent(i+1);

    } 
 
//store K0short and Lambda pt RAW mass distribution:

    TCanvas* r1 = new TCanvas();

    r1->Print("K0short_June26_HIJING.pdf[");
    for (int p = 0; p < 10; p++){

        ks_yield[p] = ks_YieldCal( k0Hist[p] );
        r1->Print("K0short_June26_HIJING.pdf");

    }
    r1->Print("K0short_June26_HIJING.pdf]");


    TCanvas* r2 = new TCanvas();

    r2->Print("Lambda_June26_HIJING.pdf[");
    for(int p1 = 0; p1 < 10; p1++){

        la_yield[p1] = la_YieldCal( laHist[p1] );
        r2->Print("Lambda_June26_HIJING.pdf");

    }
    r2->Print("Lambda_June26_HIJING.pdf]");


//-----------------
//drawing histogram:
//-----------------

//K0short:

    TCanvas* s1 = new TCanvas("s1","s1",600,800);
    TH1F* h1 = new TH1F("h1","Kshort efficiency",10,ptbins);

    double temp[10] = {0,0,0,0,0,0,0,0,0,0};

    for (int it = 0 ; it < 10; it++){

    //replace ptbinwith[it] with genKS_yield[it]:

        temp[it] = ks_yield[it]/genKS_yield[it];
        h1->SetBinContent(it+1, temp[it] );
    }

    h1->SetMarkerStyle(21);
    h1->SetMarkerColor(kBlue);
    h1->SetAxisRange(0,0.3,"Y");
    h1->SetXTitle("P^{}_{T,V0} (GeV/c)");
    h1->SetYTitle("Efficiency");
    h1->Draw("P");

//Lambda:

    TH1F* h2 = new TH1F("h2","Lambda efficiency",10,ptbins);

    double temp1[10] = {0,0,0,0,0,0,0,0,0,0};

    for (int is = 0; is < 10; is++){

    //replace ptbinwidth[is] with genLA_yield[is]:

        temp1[is] = la_yield[is]/genLA_yield[is];
        h2->SetBinContent(is+1, temp1[is] );
    }

    h2->SetMarkerStyle(21);
    h2->SetMarkerColor(kRed);
    h2->SetXTitle("P^{}_{T,V0} (GeV/c)");
    h2->SetYTitle("Efficiency");
    h2->Draw("Psame");


    TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
    w1->AddEntry(h1,"K^{0}_{s}");
    w1->AddEntry(h2,"#Lambda/#bar{#Lambda}");
    w1->Draw("same");

    
}