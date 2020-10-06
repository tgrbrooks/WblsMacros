#ifndef TARGET_H
#define TARGET_H

class Target
{
  public:

  TString name;
  double pmt_bound_R;
  double pmt_bound_Z;

  // One entry for each event (first MC particle and first trigger)
  std::vector<double> mc_energy;
  std::vector<double> mc_time;
  std::vector<double> mc_cherenkov;
  std::vector<double> mc_scintillation;
  std::vector<double> mc_x;
  std::vector<double> mc_y;
  std::vector<double> mc_z;
  std::vector<double> mc_u;
  std::vector<double> mc_v;
  std::vector<double> mc_w;

  std::vector<double> pmt_hits;
  std::vector<double> inner_hits;
  std::vector<double> veto_hits;
  std::vector<double> reco_pe;
  std::vector<double> reco_time;
  std::vector<double> n_triggers;

  std::vector<double> bonsai_x;
  std::vector<double> bonsai_y;
  std::vector<double> bonsai_z;
  std::vector<double> bonsai_u;
  std::vector<double> bonsai_v;
  std::vector<double> bonsai_w;
  std::vector<double> n9;
  std::vector<double> good_dir;
  std::vector<double> centroid_x;
  std::vector<double> centroid_y;
  std::vector<double> centroid_z;

  std::map<TString, TString> labels;

  Target(TString n, double r, double z)
  {
    name = n;
    pmt_bound_R = r;
    pmt_bound_Z = z;
    labels["mc_x"] = "X^{mc} [m]";
    labels["mc_y"] = "Y^{mc} [m]";
    labels["mc_z"] = "Z^{mc} [m]";
    labels["mc_u"] = "Dir X^{mc}";
    labels["mc_v"] = "Dir Y^{mc}";
    labels["mc_w"] = "Dir Z^{mc}";
    labels["mc_theta"] = "Dir #theta^{mc} [deg]";
    labels["mc_phi"] = "Dir #phi^{mc} [deg]";
    labels["mc_wall"] = "Distance to Wall^{mc} [m]";
    labels["mc_r"] = "R^{mc} [m]";
    labels["mc_angle"] = "#theta^{mc} [deg]";
    labels["mc_energy"] = "E^{mc} [MeV]";
    labels["mc_time"] = "t^{mc} [ns]";
    labels["mc_cherenkov"] = "Cherenkov #gamma";
    labels["mc_scintillation"] = "Scintillation #gamma";
    labels["pmt_hits"] = "PMT Hits";
    labels["inner_hits"] = "Inner PMT Hits";
    labels["veto_hits"] = "Veto PMT Hits";
    labels["reco_pe"] = "PE^{reco}";
    labels["reco_time"] = "t^{reco} [ns]";
    labels["dtime"] = "#Delta t^{reco-mc} [ns]";
    labels["n_triggers"] = "Triggers";
    labels["bonsai_x"] = "X^{bonsai} [m]";
    labels["bonsai_y"] = "Y^{bonsai} [m]";
    labels["bonsai_z"] = "Z^{bonsai} [m]";
    labels["bonsai_u"] = "Dir X^{bonsai}";
    labels["bonsai_v"] = "Dir Y^{bonsai}";
    labels["bonsai_w"] = "Dir Z^{bonsai}";
    labels["bonsai_theta"] = "Dir #theta^{bonsai} [deg]";
    labels["bonsai_dtheta"] = "#Delta Dir #theta^{bonsai-mc} [deg]";
    labels["bonsai_phi"] = "Dir #phi^{bonsai} [deg]";
    labels["bonsai_dphi"] = "#Delta Dir #phi^{bonsai-mc} [deg]";
    labels["bonsai_wall"] = "Distance to Wall^{bonsai} [m]";
    labels["bonsai_dwall"] = "#Delta Distance to Wall^{bonsai-mc} [m]";
    labels["bonsai_r"] = "R^{bonsai} [m]";
    labels["bonsai_dx"] = "#Delta X^{bonsai-mc} [m]";
    labels["bonsai_dy"] = "#Delta Y^{bonsai-mc} [m]";
    labels["bonsai_dz"] = "#Delta Z^{bonsai-mc} [m]";
    labels["bonsai_ddir"] = "Dir Angle^{bonsai-mc} [deg]";
    labels["bonsai_dv"] = "#Delta Dir Y^{bonsai-mc} [m]";
    labels["bonsai_dw"] = "#Delta Dir Z^{bonsai-mc} [m]";
    labels["bonsai_dr"] = "#Delta R^{bonsai-mc} [m]";
    labels["bonsai_dangle"] = "#Delta #theta^{bonsai-mc} [deg]";
    labels["bonsai_dist"] = "Vertex Distance^{bonsai-mc} [m]";
    labels["bonsai_angle"] = "#theta^{bonsai} [deg]";
    labels["n9"] = "N9";
    labels["good_dir"] = "Direction Goodness";
    labels["centroid_x"] = "X^{centroid} [m]";
    labels["centroid_y"] = "Y^{centroid} [m]";
    labels["centroid_z"] = "Z^{centroid} [m]";
    labels["centroid_r"] = "R^{centroid} [m]";
    labels["centroid_angle"] = "#theta^{centroid} [m]";
    labels["centroid_dx"] = "#Delta X^{centroid-mc} [m]";
    labels["centroid_dy"] = "#Delta Y^{centroid-mc} [m]";
    labels["centroid_dz"] = "#Delta Z^{centroid-mc} [m]";
    labels["centroid_wall"] = "Distance to Wall^{centroid} [m]";
    labels["centroid_dwall"] = "#Delta Distance to Wall^{centroid-mc} [m]";
    labels["centroid_dr"] = "#Delta R^{centroid-mc} [m]";
    labels["centroid_dangle"] = "#Delta #theta^{centroid-mc} [m]";
    labels["centroid_dist"] = "Vertex Distance^{centroid-mc} [m]";
  }

  void SetBoundR(double r){
    pmt_bound_R = r;
  }

  void SetBoundZ(double z){
    pmt_bound_Z = z;
  }

  double Distance(int i, TString fitter){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(fitter == "bonsai"){
      if(bonsai_x[i] <= -99.99) return -99999;
      double dx = bonsai_x[i] - mc_x[i];
      double dy = bonsai_y[i] - mc_y[i];
      double dz = bonsai_z[i] - mc_z[i];
      double dist = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2));
      return dist;
    }
    if(bonsai_x[i] <= -99.99) return -99999;
    double dx = centroid_x[i] - mc_x[i];
    double dy = centroid_y[i] - mc_y[i];
    double dz = centroid_z[i] - mc_z[i];
    return std::sqrt(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2));
  }

  double WallDist(int i, TString fitter="none"){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(bonsai_x[i] <= -99.99) return -99999;
    double r = Radius(i, fitter);
    if(fitter == "bonsai"){
      return std::min(pmt_bound_R - r, pmt_bound_Z - std::abs(bonsai_z[i]));
    }
    else if(fitter == "centroid"){
      return std::min(pmt_bound_R - r, pmt_bound_Z - std::abs(centroid_z[i]));
    }
    return std::min(pmt_bound_R - r, pmt_bound_Z - std::abs(mc_z[i]));
  }

  double Radius(int i, TString fitter="none"){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(fitter == "bonsai"){
      return std::sqrt(std::pow(bonsai_x[i], 2) + std::pow(bonsai_y[i], 2));
    }
    else if(fitter == "centroid"){
      return std::sqrt(std::pow(centroid_x[i], 2) + std::pow(centroid_y[i], 2));
    }
    return std::sqrt(std::pow(mc_x[i], 2) + std::pow(mc_y[i], 2));
  }

  double Angle(int i, TString fitter="none"){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(fitter == "bonsai"){
      return std::atan2(bonsai_y[i], bonsai_x[i])*180./TMath::Pi();
    }
    else if(fitter == "centroid"){
      return std::atan2(centroid_y[i], centroid_x[i])*180./TMath::Pi();
    }
    return std::atan2(mc_y[i], mc_x[i])*180./TMath::Pi();
  } 

  double DirAngle(int i){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(bonsai_u[i] == -99999) return -99999;
    TVector3 mc_unit(mc_u[i], mc_v[i], mc_w[i]);
    TVector3 unit(bonsai_u[i], bonsai_v[i], bonsai_w[i]);
    return mc_unit.Angle(unit)*180./TMath::Pi();
  }

  double DirTheta(int i, TString fitter="none"){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(fitter == "bonsai"){
      double r = sqrt(std::pow(bonsai_u[i],2)+std::pow(bonsai_v[i],2)+std::pow(bonsai_w[i],2));
      return std::acos(bonsai_w[i]/r)*180./TMath::Pi();
    }
    double r = sqrt(std::pow(mc_u[i],2)+std::pow(mc_v[i],2)+std::pow(mc_w[i],2));
    return std::acos(mc_w[i]/r)*180./TMath::Pi();
  }

  double DirPhi(int i, TString fitter="none"){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(fitter == "bonsai"){
      return std::atan2(bonsai_v[i], bonsai_u[i])*180./TMath::Pi();
    }
    return std::atan2(mc_v[i], mc_u[i])*180./TMath::Pi();
  }
  
  double Diff(int i, TString variable1, TString variable2){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(Get(i, variable1) <= -99999 || Get(i, variable2) <= -99999) return -99999;
    if((variable1.Contains("angle") && variable2.Contains("angle")) ||
       (variable1.Contains("theta") && variable2.Contains("theta")) ||
       (variable1.Contains("phi") && variable2.Contains("phi"))){
      double a1 = Get(i, variable1)*TMath::Pi()/180.;
      double a2 = Get(i, variable2)*TMath::Pi()/180.;
      return std::atan2(std::sin(a1-a2), std::cos(a1-a2))*180./TMath::Pi();
    }
    return Get(i, variable1) - Get(i, variable2);
  }

  double Truncation(TString variable, double trunc){
    std::vector<double> var;
    for(size_t i = 0; i < mc_x.size(); i++){
      double v = Get(i, variable);
      if(v <= -99999) continue;
      var.push_back(Get(i, variable));
    }
    std::sort(var.begin(), var.end());
    int bin = std::floor((1.-trunc/100.)*var.size());
    return var[bin];
  }

  double Get(int i, TString variable){
    if(i >= mc_x.size() || i < 0) return -99999;
    if(variable == "mc_x") return mc_x[i];
    if(variable == "mc_y") return mc_y[i];
    if(variable == "mc_z") return mc_z[i];
    if(variable == "mc_u") return mc_u[i];
    if(variable == "mc_v") return mc_v[i];
    if(variable == "mc_w") return mc_w[i];
    if(variable == "mc_theta") return DirTheta(i);
    if(variable == "mc_phi") return DirPhi(i);
    if(variable == "mc_wall") return WallDist(i);
    if(variable == "mc_r") return Radius(i);
    if(variable == "mc_angle") return Angle(i);
    if(variable == "mc_energy") return mc_energy[i];
    if(variable == "mc_time") return mc_z[i];
    if(variable == "mc_cherenkov") return mc_cherenkov[i];
    if(variable == "mc_scintillation") return mc_scintillation[i];
    if(variable == "pmt_hits") return pmt_hits[i];
    if(variable == "inner_hits") return inner_hits[i];
    if(variable == "veto_hits") return veto_hits[i];
    if(variable == "reco_pe") return reco_pe[i];
    if(variable == "reco_time") return reco_time[i];
    if(variable == "dtime") return Diff(i, "reco_time", "mc_time");
    if(variable == "n_triggers") return n_triggers[i];
    if(variable == "bonsai_x") return bonsai_x[i];
    if(variable == "bonsai_y") return bonsai_y[i];
    if(variable == "bonsai_z") return bonsai_z[i];
    if(variable == "bonsai_u") return bonsai_u[i];
    if(variable == "bonsai_v") return bonsai_v[i];
    if(variable == "bonsai_w") return bonsai_w[i];
    if(variable == "bonsai_theta") return DirTheta(i, "bonsai");
    if(variable == "bonsai_dtheta") return Diff(i, "bonsai_theta", "mc_theta");
    if(variable == "bonsai_phi") return DirPhi(i, "bonsai");
    if(variable == "bonsai_dphi") return Diff(i, "bonsai_phi", "mc_phi");
    if(variable == "bonsai_wall") return WallDist(i, "bonsai");
    if(variable == "bonsai_r") return Radius(i, "bonsai");
    if(variable == "bonsai_angle") return Angle(i, "bonsai");
    if(variable == "bonsai_dx") return Diff(i, "bonsai_x", "mc_x");
    if(variable == "bonsai_dy") return Diff(i, "bonsai_y", "mc_y");
    if(variable == "bonsai_dz") return Diff(i, "bonsai_z", "mc_z");
    if(variable == "bonsai_ddir") return DirAngle(i);
    if(variable == "bonsai_dr") return Diff(i, "bonsai_r", "mc_r");
    if(variable == "bonsai_dwall") return Diff(i, "bonsai_wall", "mc_wall");
    if(variable == "bonsai_dangle") return Diff(i, "bonsai_angle", "mc_angle");
    if(variable == "bonsai_dist") return Distance(i, "bonsai");
    if(variable == "n9") return n9[i];
    if(variable == "good_dir") return good_dir[i];
    if(variable == "centroid_x") return centroid_x[i];
    if(variable == "centroid_y") return centroid_y[i];
    if(variable == "centroid_z") return centroid_z[i];
    if(variable == "centroid_wall") return WallDist(i, "centroid");
    if(variable == "centroid_r") return Radius(i, "centroid");
    if(variable == "centroid_angle") return Angle(i, "centroid");
    if(variable == "centroid_dx") return Diff(i, "centroid_x", "mc_x");
    if(variable == "centroid_dy") return Diff(i, "centroid_y", "mc_y");
    if(variable == "centroid_dz") return Diff(i, "centroid_z", "mc_z");
    if(variable == "centroid_dwall") return Diff(i, "centroid_wall", "mc_wall");
    if(variable == "centroid_dr") return Diff(i, "centroid_r", "mc_r");
    if(variable == "centroid_dangle") return Diff(i, "centroid_angle", "mc_angle");
    if(variable == "centroid_dist") return Distance(i, "centroid");
    return -99999;
  }


  TH1F* Hist1D(TString variable, int bins, double min, double max){
    TH1F *hist = new TH1F("hist1d"+name+variable, "", bins, min, max);
    for(size_t i = 0; i < mc_x.size(); i++){
      hist->Fill(Get(i, variable));
    }
    hist->SetXTitle(labels[variable]);
    return hist;
  }

  std::pair<double, double> FitHist1D(TString variable, int bins, double min, double max){
    TH1F *hist = new TH1F("fithist1d"+name+variable, "", bins, min, max);
    for(size_t i = 0; i < mc_x.size(); i++){
      hist->Fill(Get(i, variable));
    }
    TF1 *gaus = new TF1("mg", "gaus", min, max);
    hist->Fit(gaus, "RQE");
    double bias = gaus->GetParameter(1);
    double res = gaus->GetParameter(2);
    delete hist;
    return std::make_pair(bias, res);
  }

  double Count(TString variable, double max){
    int count = 0;
    for(size_t i = 0; i < mc_x.size(); i++){
      double var = Get(i, variable);
      if(var < max){
        count++;
      }
    }
    return count;
  }

  TH1F* Cut1D(TString variable, int bins, double min, double max){
    TH1F *hist = new TH1F("cut1d"+name+variable, "", bins, min, max);
    for(size_t i = 0; i < bins; i++){
      double cut = min + (i+0.5)*(max - min)/bins;
      hist->SetBinContent(i+1, Count(variable, cut));
    }
    hist->SetXTitle("Cut on "+labels[variable]);
    return hist;
  }

  TH2F* Hist2D(TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){
    TH2F *hist = new TH2F("hist2d"+name+variable1+variable2, "", bins1, min1, max1, bins2, min2, max2);
    for(size_t i = 0; i < mc_x.size(); i++){
      hist->Fill(Get(i, variable1), Get(i, variable2));
    }
    hist->SetXTitle(labels[variable1]);
    hist->SetYTitle(labels[variable2]);
    return hist;
  }


  std::pair<TH1F*, TH1F*> BiasRes(TString variable1, int bins1, double min1, double max1, TString variable2, int bins2, double min2, double max2){

    TH2F* hist = Hist2D(variable1, bins1, min1, max1, variable2, bins2, min2, max2);

    TH1F* bias = new TH1F("bias"+name+variable1+variable2, "", bins1, min1, max1);
    bias->GetXaxis()->SetTitle(labels[variable1]);
    bias->GetYaxis()->SetTitle("Bias "+labels[variable2]);

    TH1F* res = new TH1F("res"+name+variable1+variable2, "", bins1, min1, max1);
    res->GetXaxis()->SetTitle(labels[variable1]);
    res->GetYaxis()->SetTitle("Resolution "+labels[variable2]);

    for(int i = 1; i <= bins1; i++){
      TH1D* hbin = hist->ProjectionY("bin"+name+variable1+variable2, i, i+1);
      if(hbin->GetEntries() < 20) continue;
      
      TF1 *gaus = new TF1("mg", "gaus", min2,max2);
      gaus->SetParLimits(1, min2, max2);
      hbin->Fit(gaus, "RQE");
      if(gaus->GetParError(1) > 3) continue;
      bias->SetBinContent(i, gaus->GetParameter(1));
      bias->SetBinError(i, gaus->GetParError(1));
      res->SetBinContent(i, gaus->GetParameter(2));
      res->SetBinError(i, gaus->GetParError(2));
    }

    delete hist;
    return std::make_pair(bias, res);
  }


  TH1F* Fraction(TString variable1, double min1, TString variable2, int bins, double min2, double max2){

    TH1F* pass = new TH1F("pass"+name+variable1+variable2, "", bins, min2, max2);
    pass->GetXaxis()->SetTitle(labels[variable2]);
    pass->GetYaxis()->SetTitle("Fraction with "+labels[variable1]+" > 9");

    TH1F* tot = new TH1F("fail"+name+variable1+variable2, "", bins, min2, max2);

    for(size_t i = 0; i < mc_x.size(); i++){
      double var1 = Get(i, variable1);
      double var2 = Get(i, variable2);
      if(var1 <= -99999) continue;
      if(var1 > min1) pass->Fill(var2);
      tot->Fill(var2);
    }
    pass->Divide(tot);
    return pass;
  }


  std::pair<double, double> MeanStd(TString variable, double min=0, double max=0){

    double mean = 0;
    int count = 0;
    for(size_t i = 0; i < mc_x.size(); i++){
      double var = Get(i, variable);
      if(var <= -99999) continue;
      if(min != 0 && max != 0 && (var < min || var > max)) continue;
      mean += var;
      count++;
    }
    mean /= count;
    double std_dev = 0;
    for(size_t i = 0; i < mc_x.size(); i++){
      double var = Get(i, variable);
      if(var <= -99999) continue;
      if(min != 0 && max != 0 && (var < min || var > max)) continue;
      std_dev += std::pow(var - mean, 2);
    }
    std_dev = std::sqrt(std_dev/(count-1));
    return std::make_pair(mean, std_dev);
  }


  TH1F* MeanStdHist(TString variable1, int bins, double min, double max, TString variable2, double trunc=-1){

    TH1F *hist = new TH1F("mean"+name+variable1+variable2, "", bins, min, max);
    hist->GetXaxis()->SetTitle(labels[variable1]);
    hist->GetYaxis()->SetTitle("Mean "+labels[variable2]);

    double trunc_val = -1;
    if(trunc != -1) trunc_val = Truncation(variable2, trunc);

    double width = (max - min)/bins;
    for(int b = 1; b <= bins; b++){
      double mean = 0;
      int count = 0;
      double std_dev = 0;
      for(size_t i = 0; i < mc_x.size(); i++){
        double v1 = Get(i, variable1);
        double v2 = Get(i, variable2);
        if(v1 >= min+(b-1)*width && v1 < min+b*width){
          if(v2 <= -99999) continue;
          if(trunc_val != -1 && v2 > trunc_val) continue;
          mean += v2;
          count++;
        }
      }
      mean /= count;
      for(size_t i = 0; i < mc_x.size(); i++){
        double v1 = Get(i, variable1);
        double v2 = Get(i, variable2);
        if(v1 >= min + (b-1)*width && v1 < min + b*width){
          if(v2 <= -99999) continue;
          if(trunc_val != -1 && v2 > trunc_val) continue;
          std_dev += std::pow(v2 - mean, 2);
        }
      }
      std_dev = std::sqrt(std_dev/(count-1));
      if(count == 0){
        mean = 0;
        std_dev = 0;
      }
      if(count == 1) std_dev = 0;
      hist->SetBinContent(b, mean);
      hist->SetBinContent(b, std_dev);
    }

    return hist;
  }

  TGraphAsymmErrors* Cumulative(TString variable1, int bins, double min, double max, TString variable2, double percent, double err){

    double width = (max - min)/bins;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ey_min;
    std::vector<double> ey_max;
    for(int b = 1; b <= bins; b++){
      std::vector<double> vals;
      for(size_t i = 0; i < mc_x.size(); i++){
        double v1 = Get(i, variable1);
        double v2 = Get(i, variable2);
        if(v1 >= min+(b-1)*width && v1 < min+b*width){
          if(v2 <= -99999) continue;
          vals.push_back(v2);
        }
      }
      if(vals.size()==0) continue;
      std::sort(vals.begin(), vals.end());
      int entry = std::floor(percent*vals.size());
      int entry_min = std::floor((percent-err)*vals.size());
      int entry_max = std::floor((percent+err)*vals.size());
      x.push_back(min+(b-0.5)*width);
      y.push_back(vals[entry]);
      ey_min.push_back(vals[entry]-vals[entry_min]);
      ey_max.push_back(vals[entry_max]-vals[entry]);
    }
    std::vector<double> empty(0, x.size());
    TGraphAsymmErrors *graph = new TGraphAsymmErrors(x.size(), &x[0], &y[0], &empty[0], &empty[0], &ey_min[0], &ey_max[0]);
    graph->GetXaxis()->SetTitle(labels[variable1]);
    graph->GetYaxis()->SetTitle(Form("%.1f CDF "+labels[variable2], percent));

    return graph;
  }

  // FIXME Bad, doesn't work
  void CorrectBias(){

    //std::pair<TH1F*, TH1F*> centroid_z_bias = BiasRes("centroid_z", 10, -7, 7, "centroid_dz", 15, -5, 5);
    std::pair<TH1F*, TH1F*> centroid_r_bias = BiasRes("mc_r", 7, 0, 7, "centroid_r", 7, 0, 7);
    TF1 *pol1 = new TF1("p", "pol1");
    centroid_r_bias.first->Fit(pol1, "RQE");
    double m = pol1->GetParameter(0);
    double c = pol1->GetParameter(1);
    //std::pair<TH1F*, TH1F*> bonsai_z_bias = BiasRes("bonsai_z", 10, -7, 7, "bonsai_dz", 15, -5, 5);
    //std::pair<TH1F*, TH1F*> bonsai_r_bias = BiasRes("bonsai_r", 7, 0, 7, "bonsai_dr", 15, -5, 5);

    std::cout<<"number of events = "<<mc_x.size()<<"\n";
    for(size_t i = 0; i < mc_x.size(); i++){
      if(bonsai_z[i] == -99.99) continue;
      //centroid_z[i] = centroid_z[i] - centroid_z_bias.first->Interpolate(centroid_z[i]);
      std::cout<<"Initial (x, y, z, r, a) = ("<<centroid_x[i]<<", "<<centroid_y[i]<<", "<<centroid_z[i]<<", "<<Radius(i, "centroid")<<", "<<Angle(i, "centroid")<<")\n";
      std::cout<<"Correction = "<<centroid_r_bias.first->Interpolate(Radius(i, "centroid"))<<"\n";
      //double centroid_r = Radius(i, "centroid") - centroid_r_bias.first->Interpolate(Radius(i, "centroid"));
      double centroid_r = (Radius(i, "centroid") - c)/m;
      double centroid_angle = Angle(i, "centroid");
      centroid_x[i] = centroid_r * std::cos(centroid_angle);
      centroid_y[i] = centroid_r * std::sin(centroid_angle);
      std::cout<<"Final (x, y, z, r, a) = ("<<centroid_x[i]<<", "<<centroid_y[i]<<", "<<centroid_z[i]<<", "<<Radius(i, "centroid")<<", "<<Angle(i, "centroid")<<")\n";
      //bonsai_z[i] = bonsai_z[i] - bonsai_z_bias.first->Interpolate(bonsai_z[i]);
      //double bonsai_r = Radius(i, "bonsai") - bonsai_r_bias.first->Interpolate(Radius(i, "bonsai"));
      //double bonsai_angle = Angle(i, "bonsai");
      //bonsai_x[i] = bonsai_r * std::cos(bonsai_angle);
      //bonsai_y[i] = bonsai_r * std::sin(bonsai_angle);
    }

    delete centroid_r_bias.first, centroid_r_bias.second;
  }


};

#endif
