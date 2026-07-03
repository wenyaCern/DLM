#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> SplitTab(const std::string& line)
{
  std::vector<std::string> columns;
  std::stringstream ss(line);
  std::string item;
  while(std::getline(ss, item, '\t')) columns.push_back(item);
  return columns;
}

double ToDouble(const std::string& text)
{
  if(text=="NA" || text.empty()) return 0.0;
  return std::stod(text);
}

double MtCenter(const std::string& mtBin)
{
  const std::size_t pos = mtBin.find("to");
  if(pos==std::string::npos) return 0.0;
  return 0.5*(std::stod(mtBin.substr(0,pos)) + std::stod(mtBin.substr(pos+2)));
}

bool KeepPoint(const double value, const double error)
{
  if(value<=0.0) return false;
  if(error<0.0) return false;
  if(error>5.0) return false;
  if(error>2.0*value) return false;
  return true;
}

struct FitPoint {
  double x = 0.0;
  double y = 0.0;
  double ey = 0.0;
};

TGraphErrors* MakeGraph(const std::vector<FitPoint>& points,
                        const int color,
                        const int marker,
                        const TString& name)
{
  std::vector<FitPoint> sortedPoints = points;
  std::sort(sortedPoints.begin(), sortedPoints.end(),
            [](const FitPoint& lhs, const FitPoint& rhs) {
              return lhs.x < rhs.x;
            });

  TGraphErrors* graph = new TGraphErrors(sortedPoints.size());
  graph->SetName(name);
  graph->SetLineColor(color);
  graph->SetMarkerColor(color);
  graph->SetMarkerStyle(marker);
  graph->SetMarkerSize(1.1);
  graph->SetLineWidth(2);

  for(int i=0; i<int(sortedPoints.size()); ++i){
    graph->SetPoint(i, sortedPoints[i].x, sortedPoints[i].y);
    graph->SetPointError(i, 0.0, sortedPoints[i].ey);
  }
  return graph;
}

void DrawPanel(const std::map<std::string, std::vector<FitPoint>>& values,
               const TString& yTitle,
               const TString& frameTitle,
               const std::map<std::string, int>& colors,
               const std::map<std::string, int>& markers,
               TLegend* legend,
               std::vector<TGraphErrors*>& ownedGraphs)
{
  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle(frameTitle + ";m_{T} (GeV/#it{c}^{2});" + yTitle);

  for(const auto& item : values){
    const std::string& cent = item.first;
    TGraphErrors* graph = MakeGraph(item.second,
                                    colors.at(cent),
                                    markers.at(cent),
                                    frameTitle + "_" + cent.c_str());
    mg->Add(graph, "LP");
    ownedGraphs.push_back(graph);
    if(legend) legend->AddEntry(graph, ("cent " + cent).c_str(), "lp");
  }

  mg->Draw("A");
  mg->GetXaxis()->SetLimits(0.25, 0.85);
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleSize(0.045);
  mg->GetYaxis()->SetTitleSize(0.045);
  mg->GetXaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetLabelSize(0.04);
  gPad->SetGrid(1,1);
}

} // namespace

void plot_rside_qgauss_lognormal_vs_mT(
    TString tableFileName = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS/source_shape_fit_parameters3D_Rside.tsv",
    TString outputPrefix = "/Users/wuwenya/Documents/code/femto_flow/CATS/CATS_install/CATS_3D/EPOS_data/scan_EPOS/rside_qgauss_lognormal_vs_mT_by_cent_clean")
{
  std::ifstream input(tableFileName.Data());
  if(!input.is_open()){
    std::cerr << "ERROR: cannot open " << tableFileName << std::endl;
    return;
  }

  std::map<std::string, std::vector<FitPoint>> qRadius;
  std::map<std::string, std::vector<FitPoint>> qNu;
  std::map<std::string, std::vector<FitPoint>> logMedian;
  std::map<std::string, std::vector<FitPoint>> logSigma;

  std::string line;
  bool isHeader = true;
  while(std::getline(input,line)){
    if(line.empty()) continue;
    if(isHeader){
      isHeader = false;
      continue;
    }

    const std::vector<std::string> col = SplitTab(line);
    if(col.size()<25) continue;

    const std::string& cent = col[0];
    const std::string& mtBin = col[1];
    const std::string& model = col[2];
    const double mt = MtCenter(mtBin);

    if(model=="normGammaMixtureQGaussian"){
      const double r = ToDouble(col[9]);
      const double er = ToDouble(col[10]);
      const double nu = ToDouble(col[23]);
      const double enu = ToDouble(col[24]);
      if(KeepPoint(r,er)) qRadius[cent].push_back({mt,r,er});
      if(KeepPoint(nu,enu)) qNu[cent].push_back({mt,nu,enu});
    }
    else if(model=="normLogNormalRadiusMixture"){
      const double median = ToDouble(col[19]);
      const double emedian = ToDouble(col[20]);
      const double sigma = ToDouble(col[21]);
      const double esigma = ToDouble(col[22]);
      if(KeepPoint(median,emedian)) logMedian[cent].push_back({mt,median,emedian});
      if(KeepPoint(sigma,esigma)) logSigma[cent].push_back({mt,sigma,esigma});
    }
  }

  gStyle->SetOptStat(0);
  TCanvas* canvas = new TCanvas("cRsideQGaussLogNormalMt",
                                "Rside source-shape parameters vs mT",
                                1300, 1000);
  canvas->Divide(2,2);

  const std::map<std::string, int> colors = {
      {"0to10", kBlack},
      {"10to20", kRed + 1},
      {"20to30", kBlue + 1},
      {"30to40", kGreen + 2},
      {"40to50", kMagenta + 1}
  };
  const std::map<std::string, int> markers = {
      {"0to10", 20},
      {"10to20", 21},
      {"20to30", 22},
      {"30to40", 23},
      {"40to50", 33}
  };

  std::vector<TGraphErrors*> ownedGraphs;
  TLegend* legend = new TLegend(0.58, 0.62, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

  canvas->cd(1);
  DrawPanel(qRadius,
            "R_{side} q-Gaussian R (fm)",
            "R_{side}: Gamma-mixture / q-Gaussian",
            colors, markers, legend, ownedGraphs);

  canvas->cd(2);
  DrawPanel(qNu,
            "R_{side} q-Gaussian #nu",
            "R_{side}: Gamma-mixture / q-Gaussian",
            colors, markers, nullptr, ownedGraphs);

  canvas->cd(3);
  DrawPanel(logMedian,
            "R_{side} log-normal median R (fm)",
            "R_{side}: log-normal radius mixture",
            colors, markers, nullptr, ownedGraphs);

  canvas->cd(4);
  DrawPanel(logSigma,
            "R_{side} log-normal #sigma_{ln R}",
            "R_{side}: log-normal radius mixture",
            colors, markers, nullptr, ownedGraphs);
  legend->Draw();

  canvas->SaveAs(outputPrefix + ".pdf");
  canvas->SaveAs(outputPrefix + ".png");

  TFile outputFile(outputPrefix + ".root", "RECREATE");
  canvas->Write();
  for(TGraphErrors* graph : ownedGraphs) graph->Write();
  outputFile.Close();

  std::cout << "Saved " << outputPrefix << ".pdf/.png/.root" << std::endl;
}
