//std::vector<TString> files = {"gd_h2o_10k", "gd_wbls_1pct", "wbls1pct", "wbls3pct", "wbls5pct"};
std::vector<TString> files = {"wbls5pct"};
std::vector<TString> names = {"Gd-H_{2}O", "Gd-WbLS(1%)", "WbLS(1%)", "WbLS(3%)", "WbLS(5%)"};
std::vector<int> colours = {46, 33, 38, 42, 30, 49};


void SaveFig(TString name, std::vector<TH1F*> hists){
  TCanvas *canvas = new TCanvas(name, name, 900, 600);
  canvas->SetLeftMargin(0.12);
  canvas->SetBottomMargin(0.14);
  canvas->SetTopMargin(0.1);
  TLegend *legend = new TLegend(0.12,0.91,0.92,0.98);
  int ind = 0;
  double hmax = 0;
  for(auto const& hist : hists){
    hist->SetLineColor(colours[ind]);
    hist->Scale(1./hist->Integral(0, hist->GetNbinsX()+1));
    hist->SetYTitle("Percentage of events");
    hist->SetLineWidth(4);
    if(ind == 0) hist->Draw("HIST");
    else hist->Draw("HIST SAME");
    int max_bin = hist->GetMaximumBin();
    if((hist->GetBinContent(max_bin) + hist->GetBinError(max_bin)) > hmax){
      hmax = hist->GetBinContent(max_bin) + hist->GetBinError(max_bin);
    }
    legend->AddEntry(hist, names[ind], "l");
    ind++;
  }
  hists[0]->GetYaxis()->SetRangeUser(0, 1.1*hmax);
  canvas->Modified();
  legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
  legend->Draw("SAME");
  canvas->SaveAs("Plots/"+name+".png");
}


void Save2D(TString name, std::vector<TH2F*> hists){
  int ind = 0;
  for(auto const& hist : hists){
    TCanvas *canvas = new TCanvas(name+files[ind], name, 900, 600);
    canvas->SetLeftMargin(0.12);
    canvas->SetBottomMargin(0.14);
    canvas->SetRightMargin(0.12);
    hist->Draw("colz");
    canvas->SaveAs("Plots/"+name+"_"+files[ind]+".png");
    ind++;
  }
}


void PlotBiasRes(TString name, std::vector<TH2F*> hists, int nbins){
  
  std::vector<TH1F*> v_bias; 
  std::vector<TH1F*> v_res; 

  int ind = 0;
  double xmin = 0;
  double xmax = 0;

  double bmax = 0;
  double rmax = 0;

  for(auto const& hist : hists){

    int nbinsx = hist->GetNbinsX();
    xmin = hist->GetXaxis()->GetXmin();
    xmax = hist->GetXaxis()->GetXmax();

    TH1F *h_bias = new TH1F("h_bias_"+name+files[ind], "", nbins, xmin, xmax);
    h_bias->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
    h_bias->GetYaxis()->SetTitle("Bias [m]");

    TH1F *h_res = new TH1F("h_res_"+name+files[ind], "", nbins, xmin, xmax);
    h_res->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
    h_res->GetYaxis()->SetTitle("Resolution [m]");

    for(int i = 1; i <= nbins; i++){

      TString hBin_name = Form("hBin_%i", i);
      int firstbin = nbinsx*((double)(i-1)/(double)nbins) + 1;
      int lastbin = nbinsx*((double)i/(double)nbins);
      TH1D* hBin = hist->ProjectionY(hBin_name+name+files[ind], firstbin, lastbin);

      //if(hBin->GetEntries() < 50) continue;

      TF1 *gaus = new TF1("mg", "gaus", -2, 2);
      hBin->Fit(gaus, "RQE");

      double bias = gaus->GetParameter(1);
      double bias_e = gaus->GetParError(1);
      double res = gaus->GetParameter(2);
      double res_e = gaus->GetParError(2);

      // If error too large ignore
      if(bias_e > 10 || res_e > 10) continue;

      h_bias->SetBinContent(i, bias);
      h_bias->SetBinError(i, bias_e);
      if(std::abs(bias) > bmax) bmax = std::abs(bias);

      h_res->SetBinContent(i, res);
      h_res->SetBinError(i, res_e);
      if(std::abs(res) > rmax) rmax = std::abs(res);

    }
    
    v_bias.push_back(h_bias);
    v_res.push_back(h_res);

    ind++;
  }

  TCanvas *bias_canvas = new TCanvas(name+"_bias", "", 900, 600);
  bias_canvas->SetLeftMargin(0.12);
  bias_canvas->SetBottomMargin(0.14);
  bias_canvas->SetTopMargin(0.1);
  TLegend *bias_legend = new TLegend(0.12,0.91,0.92,0.98);
  for(size_t i = 0; i < v_bias.size(); i++){
    v_bias[i]->SetMarkerColor(colours[i]);
    v_bias[i]->SetLineColor(colours[i]);
    if(i==0){ 
      v_bias[i]->Draw("p X0");
      v_bias[i]->GetYaxis()->SetRangeUser(-1.1*bmax, 1.1*bmax);
      bias_canvas->Modified();
    }
    else v_bias[i]->Draw("p X0 same");
    bias_legend->AddEntry(v_bias[i], names[i], "p");
  }
  TLine *l_bias = new TLine(xmin, 0, xmax, 0);
  l_bias->SetLineStyle(9);
  l_bias->Draw("same");
  bias_legend->SetNColumns(bias_legend->GetNRows() * bias_legend->GetNColumns());
  bias_legend->Draw("same");
  bias_canvas->SaveAs("Plots/"+name+"_bias.png");

  TCanvas *res_canvas = new TCanvas(name+"_res", "", 900, 600);
  res_canvas->SetLeftMargin(0.12);
  res_canvas->SetBottomMargin(0.14);
  res_canvas->SetTopMargin(0.1);
  TLegend *res_legend = new TLegend(0.12,0.91,0.92,0.98);
  for(size_t i = 0; i < v_res.size(); i++){
    v_res[i]->SetMarkerColor(colours[i]);
    v_res[i]->SetLineColor(colours[i]);
    if(i==0){ 
      v_res[i]->Draw("p X0");
      v_res[i]->GetYaxis()->SetRangeUser(0, 1.1*rmax);
      res_canvas->Modified();
    }
    else v_res[i]->Draw("p X0 same");
    res_legend->AddEntry(v_res[i], names[i], "p");
  }
  res_legend->SetNColumns(res_legend->GetNRows() * res_legend->GetNColumns());
  res_legend->Draw("same");
  res_canvas->SaveAs("Plots/"+name+"_resolution.png");

}

void PlotSwitch(TString name, std::vector<TH2F*> hists_centroid, std::vector<TH2F*> hists_bonsai, int nbins){
  
  for(size_t ind = 0; ind < hists_centroid.size(); ind++){

    double max = 0;

    int nbinsx = hists_centroid[ind]->GetNbinsX();
    double xmin = hists_centroid[ind]->GetXaxis()->GetXmin();
    double xmax = hists_centroid[ind]->GetXaxis()->GetXmax();

    TH1F *h_meanb = new TH1F("h_meanb_"+name+files[ind], "", nbins, xmin, xmax);
    h_meanb->GetXaxis()->SetTitle(hists_centroid[ind]->GetXaxis()->GetTitle());
    h_meanb->GetYaxis()->SetTitle("Mean vertex distance [m]");
    TH1F *h_meanc = new TH1F("h_meanc_"+name+files[ind], "", nbins, xmin, xmax);
    h_meanc->GetXaxis()->SetTitle(hists_centroid[ind]->GetXaxis()->GetTitle());
    h_meanc->GetYaxis()->SetTitle("Mean vertex distance [m]");

    for(int i = 1; i <= nbins; i++){

      TString hBin_name = Form("hBin_%i", i);
      int firstbin = nbinsx*((double)(i-1)/(double)nbins) + 1;
      int lastbin = nbinsx*((double)i/(double)nbins);

      TH1D* hBinc = hists_centroid[ind]->ProjectionY(hBin_name+name+files[ind]+"c", firstbin, lastbin);
      double totc = 0;
      double numc = 0;
      for(int cbin = 1; cbin < hBinc->GetNbinsX(); cbin++){
        totc += hBinc->GetBinCenter(cbin) * hBinc->GetBinContent(cbin);
        numc += hBinc->GetBinContent(cbin);
      }
      double meanc = totc/numc;

      TH1D* hBinb = hists_bonsai[ind]->ProjectionY(hBin_name+name+files[ind]+"b", firstbin, lastbin);
      double totb = 0;
      double numb = 0;
      for(int bbin = 1; bbin < hBinb->GetNbinsX(); bbin++){
        totb += hBinb->GetBinCenter(bbin) * hBinb->GetBinContent(bbin);
        numb += hBinb->GetBinContent(bbin);
      }
      double meanb = totb/numb;

      h_meanb->SetBinContent(i, meanb);
      h_meanc->SetBinContent(i, meanc);

    }

    TCanvas *canvas = new TCanvas(name+files[ind]+"_switch", "", 900, 600);
    canvas->SetLeftMargin(0.12);
    canvas->SetBottomMargin(0.14);
    canvas->SetTopMargin(0.1);
    TLegend *legend = new TLegend(0.12,0.91,0.92,0.98);
    h_meanb->SetMarkerColor(colours[0]);
    h_meanb->SetLineColor(colours[0]);
    h_meanb->Draw("p");
    h_meanb->GetYaxis()->SetRangeUser(0, 15);
    legend->AddEntry(h_meanb, "bonsai", "p");
    h_meanc->SetMarkerColor(colours[1]);
    h_meanc->SetLineColor(colours[1]);
    h_meanc->Draw("p same");
    legend->AddEntry(h_meanc, "centroid", "p");
    legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
    legend->Draw("same");
    canvas->SaveAs("Plots/"+name+"_"+files[ind]+"_switch.png");

  }

}

TH1F* ToMeanHist(std::map<double, std::vector<double>> energy_dist, TString name, int bins, double min, double max){
  int bin = 0;
  TH1F *h_mean = new TH1F("h_mean"+name, "", bins, min, max);
  h_mean->GetXaxis()->SetTitle("Energy [MeV]");
  h_mean->GetYaxis()->SetTitle("Mean vertex distance [m]");
  for(auto const& energy_dists : energy_dist){
    double mean = 0;
    for(auto const& dist : energy_dists.second){
      mean += dist;
    }
    mean /= energy_dists.second.size();
    double std_dev = 0;
    for(auto const& dist : energy_dists.second){
      std_dev += std::pow(dist - mean, 2);
    }
    std_dev = std::sqrt(std_dev/(energy_dists.second.size()-1));
    h_mean->SetBinContent(bin, mean);
    h_mean->SetBinError(bin, std_dev);
    bin ++;
  }
  return h_mean;
}

void PlotMeans(TString name, std::vector<TH1F*> h_bonsai, std::vector<TH1F*> h_centroid){
  
  for(size_t ind = 0; ind < h_bonsai.size(); ind++){

    TCanvas *canvas = new TCanvas(name+files[ind]+"_switch", "", 900, 600);
    canvas->SetLeftMargin(0.12);
    canvas->SetBottomMargin(0.14);
    canvas->SetTopMargin(0.1);

    TLegend *legend = new TLegend(0.12,0.91,0.92,0.98);

    h_bonsai[ind]->SetMarkerColor(colours[0]);
    h_bonsai[ind]->SetLineColor(colours[0]);
    h_bonsai[ind]->Draw("p");
    h_bonsai[ind]->GetYaxis()->SetRangeUser(0, 15);
    legend->AddEntry(h_bonsai[ind], "bonsai", "p");

    h_centroid[ind]->SetMarkerColor(colours[1]);
    h_centroid[ind]->SetLineColor(colours[1]);
    h_centroid[ind]->Draw("p same");
    legend->AddEntry(h_centroid[ind], "centroid", "p");

    legend->SetNColumns(legend->GetNRows() * legend->GetNColumns());
    legend->Draw("same");
    canvas->SaveAs("Plots/"+name+"_"+files[ind]+"_switch.png");
  }
}
