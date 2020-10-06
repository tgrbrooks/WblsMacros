#ifndef PLOTTING_H
#define PLOTTING_H

#include "Target.h"

//std::vector<TString> files = {"gd_h2o_10k", "gd_wbls_1pct", "wbls1pct", "wbls3pct", "wbls5pct"};
std::vector<TString> files = {"gd_water", "wbls_1pc_gd", "wbls_1pc", "wbls_3pc", "wbls_5pc"};
//std::vector<TString> files = {"gd_h2o_pos", "wbls_1pc_gd_pos", "wbls_1pc_pos", "wbls_3pc_pos", "wbls_5pc_pos"};
//std::vector<TString> files = {"gd_h2o_tl208", "gd_wbls_1pct_tl208", "wbls1pct_tl208", "wbls3pct_tl208", "wbls5pct_tl208"};
std::vector<TString> names = {"Gd-H_{2}O", "Gd-WbLS(1%)", "WbLS(1%)", "WbLS(3%)", "WbLS(5%)"};
//std::vector<TString> names = {"n9", "n100", "seln9", "n400"};
TString dir = "FitPlots/";
std::vector<int> colours = {46, 33, 38, 42, 30, 49};


// Plot all target materials on one plot of 1D histograms
void PlotTargets1D(std::vector<Target*> targets, TString variable, int bins, double min, double max){

  TString name = variable;
  TCanvas *canvas = new TCanvas(name, name, 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

  THStack *stack = new THStack("hs"+name, "");
  TString xtitle;
  TString ytitle;

  for(size_t i = 0; i < targets.size(); i++){

    TH1F* hist = targets[i]->Hist1D(variable, bins, min, max);
    hist->SetLineColor(colours[i]);
    hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
    hist->SetYTitle("Percentage of events");
    hist->SetLineWidth(3);

    stack->Add(hist);
    xtitle = hist->GetXaxis()->GetTitle();
    ytitle = hist->GetYaxis()->GetTitle();

    legend->AddEntry(hist, names[i], "l");
  }

  stack->Draw("nostack hist");
  stack->GetXaxis()->SetTitle(xtitle);
  stack->GetYaxis()->SetTitle(ytitle);
  canvas->Modified();
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("SAME");
  canvas->SaveAs(dir+name+".png");

  // clean up
  stack->GetHists()->Delete();
  delete canvas, legend, stack;

} // PlotTargets1D


// Plot 2D histogram for each target material
void PlotTargets2D(std::vector<Target*> targets, TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

  for(size_t i = 0; i < targets.size(); i++){
    TString name = variable1+"_"+variable2+"_"+files[i];
    TCanvas *canvas = new TCanvas(name, name, 900, 600);
    canvas->SetMargin(0.12, 0.12, 0.14, 0.1);

    TH2F* hist = targets[i]->Hist2D(variable1, bins1, min1, max1, variable2, bins2, min2, max2);

    hist->Draw("COLZ");
    canvas->SaveAs(dir+name+".png");

    // clean up
    delete canvas;
    delete hist;
  }
} // PlotTargets2D


void PlotCuts1D(std::vector<Target*> targets, TString variable, int bins, double min, double max, std::vector<int> scale){

  TString name = variable + "_cut";
  TCanvas *canvas = new TCanvas(name, name, 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

  THStack *stack = new THStack("hs"+name, "");
  TString xtitle;
  TString ytitle;

  for(size_t i = 0; i < targets.size(); i++){

    TH1F* hist = targets[i]->Cut1D(variable, bins, min, max);
    hist->SetLineColor(colours[i]);
    hist->SetMarkerColor(colours[i]);
    hist->Scale(1./scale[i]);
    hist->SetYTitle("Percentage of triggers cut");
    hist->SetLineWidth(3);

    stack->Add(hist);
    xtitle = hist->GetXaxis()->GetTitle();
    ytitle = hist->GetYaxis()->GetTitle();

    legend->AddEntry(hist, names[i], "l");
  }

  stack->Draw("nostack p E1 x0");
  stack->GetHistogram()->SetMarkerColor(kWhite);
  canvas->Update();
  stack->Draw("nostack Lhist same");
  stack->GetXaxis()->SetTitle(xtitle);
  stack->GetYaxis()->SetTitle(ytitle);
  stack->GetYaxis()->SetRangeUser(0.6, 1);
  canvas->Modified();
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("SAME");
  canvas->SaveAs(dir+name+".png");

  // clean up
  stack->GetHists()->Delete();
  delete canvas, legend, stack;

} // PlotTargets1D


// Plot the arithmetic mean of a variable as a function of another for all target materials
void PlotTargetsMeans(std::vector<Target*> targets, TString variable1, int bins, double min, double max, TString variable2, double trunc=-1){

  TString name = "mean_"+variable2+"_vs_" + variable1;
  TCanvas *canvas = new TCanvas(name, name, 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

  THStack *stack = new THStack("hs"+name, "");
  TString xtitle;
  TString ytitle;

  for(size_t i = 0; i < targets.size(); i++){
    TH1F* hist = targets[i]->MeanStdHist(variable1, bins, min, max, variable2, trunc);
    hist->SetLineColor(colours[i]);
    hist->SetMarkerColor(colours[i]);
    hist->SetLineWidth(3);

    stack->Add(hist);
    xtitle = hist->GetXaxis()->GetTitle();
    ytitle = hist->GetYaxis()->GetTitle();

    legend->AddEntry(hist, names[i], "p");
  }

  stack->Draw("nostack p E1 x0");
  stack->GetHistogram()->SetMarkerColor(kWhite);
  canvas->Update();
  stack->Draw("nostack Lhist same");
  stack->GetXaxis()->SetTitle(xtitle);
  stack->GetYaxis()->SetTitle(ytitle);
  canvas->Modified();
  
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("SAME");
  canvas->SaveAs(dir+name+".png");

  // clean up
  stack->GetHists()->Delete();
  delete canvas, legend, stack;
} // PlotTargetsMeans


// Plot the cumulative mean of a variable as a function of another for all target materials
void PlotTargetsCumulative(std::vector<Target*> targets, TString variable1, int bins, double min, double max, TString variable2, double percentage, double error){

  TString name = "cdf_"+variable2+"_vs_" + variable1;
  TCanvas *canvas = new TCanvas(name, name, 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

  bool first = true;
  for(size_t i = 0; i < targets.size(); i++){
    TGraphAsymmErrors* graph = targets[i]->Cumulative(variable1, bins, min, max, variable2, percentage, error);
    graph->SetLineColor(colours[i]);
    graph->SetMarkerColor(colours[i]);
    graph->SetLineWidth(3);
    if (first){
      graph->Draw("apl");
      graph->GetYaxis()->SetRangeUser(0, 180);
      first = false;
    }
    else graph->Draw("pl same");

    legend->AddEntry(graph, names[i], "p");
  }
  
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("SAME");
  canvas->SaveAs(dir+name+".png");

  // clean up
  delete canvas, legend;
} // PlotTargetsMeans


// Plot the bias and resolution (from Gaussian fits to binned distributions) of one variable in terms of another for all target materials
void PlotTargetsBiasRes(std::vector<Target*> targets, TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

  TString name = variable2+"_vs_" + variable1;
  THStack *bias_stack = new THStack("bhs"+name, "");
  THStack *res_stack = new THStack("rhs"+name, "");
  TLegend *bias_legend = new TLegend(0.18,0.91,0.92,0.98);
  TLegend *res_legend = new TLegend(0.18,0.91,0.92,0.98);
  TString bxtitle, bytitle, rxtitle, rytitle;

  for(size_t i = 0; i < targets.size(); i++){
    std::pair<TH1F*, TH1F*> bias_res = targets[i]->BiasRes(variable1, bins1, min1, max1, variable2, bins2, min2, max2);

    bias_res.first->SetLineColor(colours[i]);
    bias_res.first->SetMarkerColor(colours[i]);
    bias_res.first->SetLineWidth(3);
    bias_stack->Add(bias_res.first);
    bias_legend->AddEntry(bias_res.first, names[i], "p");
    bxtitle = bias_res.first->GetXaxis()->GetTitle();
    bytitle = bias_res.first->GetYaxis()->GetTitle();

    bias_res.second->SetLineColor(colours[i]);
    bias_res.second->SetMarkerColor(colours[i]);
    bias_res.second->SetLineWidth(3);
    res_stack->Add(bias_res.second);
    res_legend->AddEntry(bias_res.second, names[i], "p");
    rxtitle = bias_res.second->GetXaxis()->GetTitle();
    rytitle = bias_res.second->GetYaxis()->GetTitle();

  }

  TCanvas *bias_canvas = new TCanvas(name+"_bias", "", 900, 600);
  bias_canvas->SetMargin(0.12, 0.07, 0.14, 0.1);
  bias_stack->Draw("nostack p E1 x0");
  bias_stack->GetHistogram()->SetMarkerColor(kWhite);
  bias_canvas->Update();
  bias_stack->Draw("nostack Lhist same");
  bias_stack->GetXaxis()->SetTitle(bxtitle);
  bias_stack->GetYaxis()->SetTitle(bytitle);
  bias_canvas->Modified();
  TLine *l_bias = new TLine(min1, 0, max1, 0);
  l_bias->SetLineStyle(9);
  l_bias->Draw("same");
  bias_legend->SetNColumns(bias_legend->GetNRows() * bias_legend->GetNColumns());
  bias_legend->Draw("same");
  bias_canvas->SaveAs(dir+name+"_bias.png");

  TCanvas *res_canvas = new TCanvas(name+"_res", "", 900, 600);
  res_canvas->SetMargin(0.12, 0.07, 0.14, 0.1);
  res_stack->Draw("nostack p E1 x0");
  res_stack->GetHistogram()->SetMarkerColor(kWhite);
  res_canvas->Update();
  res_stack->Draw("nostack Lhist same");
  res_stack->GetXaxis()->SetTitle(rxtitle);
  res_stack->GetYaxis()->SetTitle(rytitle);
  res_canvas->Modified();

  res_legend->SetNColumns(res_legend->GetNRows() * res_legend->GetNColumns());
  res_legend->Draw("same");
  res_canvas->SaveAs(dir+name+"_resolution.png");

  // clean up
  bias_stack->GetHists()->Delete();
  res_stack->GetHists()->Delete();
  delete bias_canvas, bias_legend, res_canvas, res_legend, bias_stack, res_stack, l_bias;
} // PlotTargetsBiasRes


// Plot the bias and resolution (from Gaussian fits to binned distributions) of one variable in terms of another for all target materials
void PlotTargetsFraction(std::vector<Target*> targets, TString variable1, double min1, TString variable2, int bins, double min2, double max2){

  TString name = variable2+"_vs_" + variable1;
  THStack *stack = new THStack("fhs"+name, "");
  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);
  TString xtitle, ytitle;

  for(size_t i = 0; i < targets.size(); i++){
    TH1F* pass = targets[i]->Fraction(variable1, min1, variable2, bins, min2, max2);

    pass->SetLineColor(colours[i]);
    pass->SetMarkerColor(colours[i]);
    pass->SetMarkerSize(2);
    pass->SetLineWidth(3);
    stack->Add(pass);
    legend->AddEntry(pass, names[i], "p");
    xtitle = pass->GetXaxis()->GetTitle();
    ytitle = pass->GetYaxis()->GetTitle();

  }

  TCanvas *canvas = new TCanvas(name+"_frac", "", 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);
  stack->Draw("nostack p");
  stack->GetHistogram()->SetMarkerColor(kWhite);
  canvas->Update();
  stack->Draw("nostack Lhist same");
  stack->GetXaxis()->SetTitle(xtitle);
  stack->GetYaxis()->SetTitle(ytitle);
  canvas->Modified();
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("same");
  canvas->SaveAs(dir+name+"_frac.png");

  // clean up
  stack->GetHists()->Delete();
  delete canvas, legend, stack;
} // PlotTargetsBiasRes


// Plot the mean and standard deviation of a gaussian fit of one variable binned in another for all target materials
void PlotTargetsMeanFit(std::vector<Target*> targets, TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

  TString name = "mean_fit_"+variable2+"_vs_" + variable1;
  THStack *stack = new THStack("hs"+name, "");
  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);
  TString xtitle, ytitle;

  for(size_t i = 0; i < targets.size(); i++){
    std::pair<TH1F*, TH1F*> bias_res = targets[i]->BiasRes(variable1, bins1, min1, max1, variable2, bins2, min2, max2);

    bias_res.first->SetLineColor(colours[i]);
    bias_res.first->SetMarkerColor(colours[i]);
    bias_res.first->SetLineWidth(3);
    for(size_t j = 1; j <= bias_res.first->GetNbinsX(); j++){
      bias_res.first->SetBinError(j, bias_res.second->GetBinContent(j));
    }
    stack->Add(bias_res.first);
    legend->AddEntry(bias_res.first, names[i], "p");
    xtitle = bias_res.first->GetXaxis()->GetTitle();
    ytitle = bias_res.first->GetYaxis()->GetTitle();

    // clean up
    delete bias_res.second;
  }

  TCanvas *canvas = new TCanvas(name, "", 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);
  stack->Draw("nostack p E1 x0");
  stack->GetHistogram()->SetMarkerColor(kWhite);
  canvas->Update();
  stack->Draw("nostack L hist same");
  stack->GetXaxis()->SetTitle(xtitle);
  stack->GetYaxis()->SetTitle(ytitle);
  canvas->Modified();

  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("same");
  canvas->SaveAs(dir+name+".png");

  // clean up
  stack->GetHists()->Delete();
  delete canvas, legend, stack;
  //for(size_t i = 0; i < v_bias.size(); i++) delete v_bias[i];
} // PlotTargetsMeanFit
  

// Plot 1D histograms for both fitters on one plot
void PlotFitters1D(std::vector<Target*> targets, TString variable, int bins, double min, double max){

  for(size_t i = 0; i < targets.size(); i++){
    TString name = variable+"_"+files[i];
    TCanvas *canvas = new TCanvas(name, name, 900, 600);
    canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

    TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

    THStack *stack = new THStack("hs"+name, "");

    TH1F* bhist = targets[i]->Hist1D("bonsai_"+variable, bins, min, max);
    bhist->SetLineColor(colours[0]);
    bhist->Scale(1./bhist->Integral(0, bhist->GetNbinsX()+1));
    bhist->SetYTitle("Percentage of events");
    bhist->SetLineWidth(3);
    stack->Add(bhist);
    legend->AddEntry(bhist, "bonsai", "l");

    TH1F* chist = targets[i]->Hist1D("centroid_"+variable, bins, min, max);
    chist->SetLineColor(colours[1]);
    chist->Scale(1./chist->Integral(0, chist->GetNbinsX()+1));
    chist->SetYTitle("Percentage of events");
    chist->SetLineWidth(3);
    stack->Add(chist);
    legend->AddEntry(chist, "centroid", "l");

    stack->Draw("nostack HIST");
    stack->GetXaxis()->SetTitle(bhist->GetXaxis()->GetTitle());
    stack->GetYaxis()->SetTitle(bhist->GetYaxis()->GetTitle());
    canvas->Modified();
    legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
    legend->Draw("SAME");
    canvas->SaveAs(dir+name+".png");

    // clean up
    delete canvas, legend, bhist, chist, stack;
  }
}


// Plot the arithmetic mean of a variable as a function of another for both fitters
void PlotFittersMeans(Target* target, int i, TString variable1, int bins, double min, double max, TString variable2, double trunc=-1){

  TString name = "mean_"+variable2+"_vs_" + variable1 + "_" + files[i];
  TCanvas *canvas = new TCanvas(name, name, 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

  THStack *stack = new THStack("hs"+name, "");

  TH1F* bhist = target->MeanStdHist(variable1, bins, min, max, "bonsai_"+variable2, trunc);
  bhist->SetLineColor(colours[0]);
  bhist->SetMarkerColor(colours[0]);
  bhist->SetLineWidth(3);
  stack->Add(bhist);
  legend->AddEntry(bhist, "bonsai", "p");

  TH1F* chist = target->MeanStdHist(variable1, bins, min, max, "centroid_"+variable2, trunc);
  chist->SetLineColor(colours[1]);
  chist->SetMarkerColor(colours[1]);
  chist->SetLineWidth(3);
  stack->Add(chist);
  legend->AddEntry(chist, "centroid", "p");

  stack->Draw("nostack p E1 x0");
  stack->GetHistogram()->SetMarkerColor(kWhite);
  canvas->Update();
  stack->Draw("nostack L hist same");
  stack->GetXaxis()->SetTitle(bhist->GetXaxis()->GetTitle());
  stack->GetYaxis()->SetTitle(bhist->GetYaxis()->GetTitle());
  canvas->Modified();
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("SAME");
  canvas->SaveAs(dir+name+".png");

  // clean up
  delete canvas, legend, bhist, chist, stack;

} // PlotFittersMean


void PlotFittersMeans(std::vector<Target*> targets, TString variable1, int bins, double min, double max, TString variable2, double trunc=-1){

  for(size_t i = 0; i < targets.size(); i++){
    PlotFittersMeans(targets[i], i, variable1, bins, min, max, variable2, trunc);
  }
}


// Plot the arithmetic mean of a variable as a function of another for both fitters
void PlotFittersCumulative(Target* target, int i, TString variable1, int bins, double min, double max, TString variable2, double percentage, double error){

  TString name = "cdf_"+variable2+"_vs_" + variable1 + "_" + files[i];
  TCanvas *canvas = new TCanvas(name, name, 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

  TGraphAsymmErrors* bgraph = target->Cumulative(variable1, bins, min, max, "bonsai_"+variable2, percentage, error);
  bgraph->SetLineColor(colours[0]);
  bgraph->SetMarkerColor(colours[0]);
  bgraph->SetLineWidth(3);
  legend->AddEntry(bgraph, "bonsai", "p");

  TGraphAsymmErrors* cgraph = target->Cumulative(variable1, bins, min, max, "centroid_"+variable2, percentage, error);
  cgraph->SetLineColor(colours[1]);
  cgraph->SetMarkerColor(colours[1]);
  cgraph->SetLineWidth(3);
  legend->AddEntry(cgraph, "centroid", "p");

  bgraph->Draw("apl");
  //bgraph->GetYaxis()->SetRangeUser(0, 3);
  cgraph->Draw("pl same");
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("SAME");
  canvas->SaveAs(dir+name+".png");

  // clean up
  delete canvas, legend, bgraph, cgraph;

} // PlotFittersCumulative


void PlotFittersCumulative(std::vector<Target*> targets, TString variable1, int bins, double min, double max, TString variable2, double percentage, double error){

  for(size_t i = 0; i < targets.size(); i++){
    PlotFittersCumulative(targets[i], i, variable1, bins, min, max, variable2, percentage, error);
  }
}


// Plot the bias and resolution (from Gaussian fits to binned distributions) of one variable in terms of another for both fitters
void PlotFittersBiasRes(Target* target, int i, TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

  TString name = variable2+"_vs_" + variable1+"_"+files[i];
  THStack *bias_stack = new THStack("bhs"+name, "");
  THStack *res_stack = new THStack("rhs"+name, "");
  TLegend *bias_legend = new TLegend(0.18,0.91,0.92,0.98);
  TLegend *res_legend = new TLegend(0.18,0.91,0.92,0.98);

  std::pair<TH1F*, TH1F*> bbias_res = target->BiasRes(variable1, bins1, min1, max1, "bonsai_"+variable2, bins2, min2, max2);

  bbias_res.first->SetLineColor(colours[0]);
  bbias_res.first->SetMarkerColor(colours[0]);
  bbias_res.first->SetLineWidth(3);
  bias_stack->Add(bbias_res.first);
  bias_legend->AddEntry(bbias_res.first, "bonsai", "p");

  bbias_res.second->SetLineColor(colours[0]);
  bbias_res.second->SetMarkerColor(colours[0]);
  bbias_res.second->SetLineWidth(3);
  res_stack->Add(bbias_res.second);
  res_legend->AddEntry(bbias_res.second, "bonsai", "p");

  std::pair<TH1F*, TH1F*> cbias_res = target->BiasRes(variable1, bins1, min1, max1, "centroid_"+variable2, bins2, min2, max2);

  cbias_res.first->SetLineColor(colours[1]);
  cbias_res.first->SetMarkerColor(colours[1]);
  cbias_res.first->SetLineWidth(3);
  bias_stack->Add(cbias_res.first);
  bias_legend->AddEntry(cbias_res.first, "centroid", "p");

  cbias_res.second->SetLineColor(colours[1]);
  cbias_res.second->SetMarkerColor(colours[1]);
  cbias_res.second->SetLineWidth(3);
  res_stack->Add(cbias_res.second);
  res_legend->AddEntry(cbias_res.second, "centroid", "p");

  TCanvas *bias_canvas = new TCanvas(name+"_bias", "", 900, 600);
  bias_canvas->SetMargin(0.12, 0.07, 0.14, 0.1);
  bias_stack->Draw("nostack p E1 x0");
  bias_stack->GetHistogram()->SetMarkerColor(kWhite);
  bias_canvas->Update();
  bias_stack->Draw("nostack Lhist same");
  bias_stack->GetXaxis()->SetTitle(bbias_res.first->GetXaxis()->GetTitle());
  bias_stack->GetYaxis()->SetTitle(bbias_res.first->GetYaxis()->GetTitle());
  bias_canvas->Modified();

  TLine *l_bias = new TLine(min1, 0, max1, 0);
  l_bias->SetLineStyle(9);
  l_bias->Draw("same");
  bias_legend->SetNColumns(bias_legend->GetNRows() * bias_legend->GetNColumns());
  bias_legend->Draw("same");
  bias_canvas->SaveAs(dir+name+"_bias.png");

  TCanvas *res_canvas = new TCanvas(name+"_res", "", 900, 600);
  res_canvas->SetMargin(0.12, 0.07, 0.14, 0.1);
  res_stack->Draw("nostack p E1 x0");
  res_stack->GetHistogram()->SetMarkerColor(kWhite);
  res_canvas->Update();
  res_stack->Draw("nostack Lhist same");
  res_stack->GetXaxis()->SetTitle(bbias_res.second->GetXaxis()->GetTitle());
  res_stack->GetYaxis()->SetTitle(bbias_res.second->GetYaxis()->GetTitle());
  res_canvas->Modified();

  res_legend->SetNColumns(res_legend->GetNRows() * res_legend->GetNColumns());
  res_legend->Draw("same");
  res_canvas->SaveAs(dir+name+"_resolution.png");

  // clean up
  delete bias_canvas, bias_legend, res_canvas, res_legend, bbias_res.first, bbias_res.second, cbias_res.first, cbias_res.second, bias_stack, res_stack;//, l_bias;

} // PlotFittersBiasRes


// Plot the bias and resolution (from Gaussian fits to binned distributions) of one variable in terms of another for both fitters
void PlotFittersBiasRes(std::vector<Target*> targets, TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

  for(size_t i = 0; i < targets.size(); i++){
    PlotFittersBiasRes(targets[i], i, variable1, bins1, min1, max1, variable2, bins2, min2, max2);
  }

}


// Plot the mean and standard deviation of a gaussian fit of one variable binned in another for both fitters
void PlotFittersMeanFit(Target* target, int i, TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

  TString name = variable2+"_vs_" + variable1 + "_" + files[i];
  THStack *stack = new THStack("hs"+name, "");
  TLegend *legend = new TLegend(0.18,0.91,0.92,0.98);

  std::pair<TH1F*, TH1F*> bonsai = target->BiasRes(variable1, bins1, min1, max1, "bonsai_"+variable2, bins2, min2, max2);

  bonsai.first->SetLineColor(colours[0]);
  bonsai.first->SetMarkerColor(colours[0]);
  bonsai.first->SetLineWidth(3);
  for(size_t j = 1; j <= bonsai.first->GetNbinsX(); j++){
    bonsai.first->SetBinError(j, bonsai.second->GetBinContent(j));
  }
  stack->Add(bonsai.first);
  legend->AddEntry(bonsai.first, "bonsai", "p");

  std::pair<TH1F*, TH1F*> centroid = target->BiasRes(variable1, bins1, min1, max1, "centroid_"+variable2, bins2, min2, max2);

  centroid.first->SetLineColor(colours[1]);
  centroid.first->SetMarkerColor(colours[1]);
  centroid.first->SetLineWidth(3);
  for(size_t j = 1; j <= centroid.first->GetNbinsX(); j++){
    centroid.first->SetBinError(j, centroid.second->GetBinContent(j));
  }
  stack->Add(centroid.first);
  legend->AddEntry(centroid.first, "centroid", "p");
  
  TCanvas *canvas = new TCanvas(name, "", 900, 600);
  canvas->SetMargin(0.12, 0.07, 0.14, 0.1);

  stack->Draw("nostack p E1 x0");
  stack->GetHistogram()->SetMarkerColor(kWhite);
  canvas->Update();
  stack->Draw("nostack Lhist same");
  stack->GetXaxis()->SetTitle(bonsai.first->GetXaxis()->GetTitle());
  stack->GetYaxis()->SetTitle(bonsai.first->GetYaxis()->GetTitle());
  canvas->Modified();

  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("same");
  canvas->SaveAs(dir+name+".png");

  // clean up
  delete canvas, legend, bonsai.first, bonsai.second, centroid.first, centroid.second, stack;

} // PlotFittersMeanFit


void PlotFittersMeanFit(std::vector<Target*> targets, TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

  for(size_t i = 0; i < targets.size(); i++){

    PlotFittersMeanFit(targets[i], i, variable1, bins1, min1, max1, variable2, bins2, min2, max2);

  }

}

#endif
