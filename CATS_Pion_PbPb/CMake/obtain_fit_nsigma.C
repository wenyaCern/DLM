void readFile(TString filename, TFile*& file);
void nSigma(TH1F* h, TF1* fit, TH1F*& hnSig);
void drawDaranSigma(TCanvas* c, TPad*& pad1, TPad*& pad2, int iCanva, TH1F* h, TF1* fit, TF1* baseline, TH1F* nSigma, int centLow, int centUp, double mTLow, double mTUp);
void SetStyle(Bool_t graypalette=kFALSE);

TFile* fileIn;
int nQnBins = 10;

double mTBin[] = {0.3, 0.4};
const int nMtBins = sizeof(mTBin)/sizeof(mTBin[0]) - 1;
int centBin[] = {30, 50};
const int nCentBins = sizeof(centBin)/sizeof(centBin[0]) - 1;

void obtain_fit_nsigma(int icent = 0, int imt = 0, bool qnDivided=false){
	// qn-integral
	

		TString filename = "fit_results_nonGauss_fixfrac_resoultion_test.root";
		readFile(filename, fileIn);
		TH1F* hData = (TH1F*)fileIn->Get("copyhisto");
		TH1F* hDataClone = (TH1F*)hData->Clone();
		TF1* fFit = (TF1*)fileIn->Get("CATS_fit");
		// TF1* fBaseline = (TF1*)fileIn->Get("baseline");
		TH1F* hnSig;

		nSigma(hDataClone, fFit, hnSig);

		SetStyle();
		TPad* pad1;
		TPad* pad2;
		TCanvas* c = new TCanvas("c","Fit with Pull",900, 600);
		drawDaranSigma(c, pad1, pad2, 1, hDataClone, fFit, nullptr, hnSig, centBin[0], centBin[1], mTBin[1], mTBin[2]);
	
}

void readFile(TString filename, TFile*& file){
	file = TFile::Open(filename);
}

void nSigma(TH1F* h, TF1* fit, TH1F*& hnSig){
	int nbins = h->GetNbinsX();
	TH1F* hPull = new TH1F("hpull",";k^{*} (GeV/c); n#sigma", h->GetNbinsX(), h->GetXaxis()->GetBinLowEdge(1), h->GetXaxis()->GetBinUpEdge(h->GetNbinsX()));

	for (int i=1; i<=nbins; i++) {
	    double O  = h->GetBinContent(i);
	    double err= h->GetBinError(i);
	    double x  = h->GetBinCenter(i);
	    double E  = fit->Eval(x);
	
	    double pull = (err>0) ? (O - E)/err : 0;
	    hPull->SetBinContent(i, pull);
	}
	hnSig = hPull;
}

void drawDaranSigma(TCanvas* c, TPad*& pad1, TPad*& pad2, int iCanva, TH1F* h, TF1* fit, TF1* baseline, TH1F* nSigma, int centLow, int centUp, double mTLow, double mTUp){
	c->cd(iCanva);
    pad1 = new TPad(Form("pad1_%d",iCanva),"",0,0.25,1,1.0);
    pad1->SetPad(0,0.25,1,1.0);
    pad1->SetBottomMargin(0.01);

    pad2 = new TPad(Form("pad2_%d",iCanva),"",0,0.0,1,0.25);
    pad2->SetPad(0,0.0,1,0.25);
    pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.4);

    pad1->Draw();
    pad2->Draw();

	pad1->cd();
	h->SetMarkerStyle(20);
	h->SetMarkerSize(1.0);
	h->SetMarkerColor(kBlack);
	h->SetLineColor(kBlack);
	h->GetYaxis()->SetRangeUser(0.45, 1.35);
	h->GetYaxis()->SetTitle("C(k^{*})");
	h->GetYaxis()->SetTitleOffset(0.5);
	h->Draw("p"); 
	fit->SetLineColor(kRed);
	fit->SetLineWidth(3); 
	fit->Draw("Same");
	// baseline->SetLineColor(kGreen-2); 
	// baseline->SetLineWidth(3);
	// baseline->Draw("Same");

	TLegend* lg = new TLegend(0.54, 0.33, 0.9, 0.6,"Pb#minusPb, #sqrt{#it{s}_{NN}} = 5.36 TeV" );
 	lg->SetFillStyle(0);
 	lg->SetBorderSize(0);
 	lg->SetTextSize(0.04);
 	lg->AddEntry(h,Form("%d-%d %% p#minusp", centLow, centUp),"apl");
 	lg->AddEntry(fit,"CATS fit","l");
 	// lg->AddEntry(baseline,"Baseline","l");
 	lg->AddEntry((TObject*)0,Form("%.2f<m_{T}<%.2f GeV/c^{2}", mTLow, mTUp),"");
 	if (iCanva==1) lg->Draw();

    pad2->cd();
    nSigma->SetStats(0);
    nSigma->SetMarkerSize(0.01);
    nSigma->SetMarkerColor(kWhite);
    nSigma->SetLineColor(kWhite);
    nSigma->GetYaxis()->SetTitle("n#sigma");
    nSigma->GetYaxis()->SetTitleOffset(0.23); 
    nSigma->GetYaxis()->SetRangeUser(-12,12);
    nSigma->GetXaxis()->SetRangeUser(0,300);
    nSigma->GetYaxis()->SetTitleSize(0.18);
    nSigma->GetYaxis()->SetLabelSize(0.1);
    nSigma->GetXaxis()->SetTitleSize(0.18);
    nSigma->GetXaxis()->SetLabelSize(0.13);
    nSigma->SetFillStyle(1001);
    nSigma->SetFillColorAlpha(kRed, 0.35);
    nSigma->Draw("e3");
}

void SetStyle(Bool_t graypalette) {
  cout << "Setting style!" << endl;
  
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  if(graypalette) gStyle->SetPalette(8,0);
  else gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameLineWidth(1);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadColor(10);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(kRed);
  gStyle->SetFuncWidth(2);
  gStyle->SetFuncColor(kGreen);
  gStyle->SetLineWidth(2);
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetLabelOffset(0.01,"y");
  gStyle->SetLabelOffset(0.01,"x");
  gStyle->SetLabelColor(kBlack,"xyz");
  gStyle->SetTitleSize(0.05,"xyz");
  gStyle->SetTitleOffset(1.15,"y");
  gStyle->SetTitleOffset(1.15,"x");
  gStyle->SetTitleFillColor(kWhite);
  gStyle->SetTextSizePixels(26);
  gStyle->SetTextFont(42);
  //  gStyle->SetTickLength(0.04,"X");  gStyle->SetTickLength(0.04,"Y"); 

  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);
  //  gStyle->SetFillColor(kWhite);
  gStyle->SetLegendFont(42);
}