#include "Target.h"
#include "Plotting.h"

// Plots needed
// Truth:
// * True kinetic energy spectrum
// * True time spectrum
// * True number of cherenkov photons
// * True number of scintillation photons
// * True cherenkov/scintillation ratio
// Trigger:
// * Number of PMT hits (vs R, Z, KE)
// * Total reconstructed PE (vs R, Z, KE)
// Reconstruction:
// * Difference between true and reco vertex X, Y, Z, R, theta (bias and resolution vs self, KE, PE)
// * Distance between true and reco vertex (mean distance vs KE, PE)

void fitterAna()
{

  gROOT->ProcessLine(".x ~/UniversalFiles/tomStyle.C");

  std::vector<Target*> targets;
  std::vector<int> vec_no_trigger;
  std::vector<int> vec_trigger;

  for(int f = 0; f < files.size(); f++){

    //FIXME how to enter bounds before 
    Target *target = new Target(files[f], 6.7, 6.7);
    int no_trigger = 0;
    int trigger = 0;

    // Get files in the directory and loop over them
    TString dirname = "/data/kneale/production/"+files[f]+"/positrons/";
    TString bdirname = "/data/brooks/wbls/"+files[f]+"_noang/";
    /*if(f != 0){
      bdirname = "/data/brooks/wbls/"+files[f]+"_noang/";
      std::cout<<"Gets here\n";
    }*/
    std::cout<<f<<": "<<dirname<<" "<<bdirname<<"\n";
    TSystemDirectory dir(dirname, dirname); 
    TList *files = dir.GetListOfFiles(); 
    if (!files) continue;
    TSystemFile *file; 
    TString fname; 
    TString bfname; 
    TIter next(files); 
    while ((file=(TSystemFile*)next())) { 
      fname = file->GetName(); 
      if(fname == "positrons.root") continue;
      if(fname == "run1589498224813277892.root") continue;
      bfname = fname.ReplaceAll(".root", "_bonsai.root");
      fname.ReplaceAll("_bonsai.root", ".root");
      if (file->IsDirectory() || !fname.EndsWith(".root")) continue;
      std::cout << (dirname+fname).Data() << " " << (bdirname+bfname).Data() << std::endl; 

      // Get PMT info
      TFile* rf = new TFile(dirname+fname);
      //TFile* rf = new TFile("/data/kneale/production/gd_water/run1589492333848056055.root");
      TTree* run_tree=(TTree*) rf->Get("runT");
      RAT::DS::Run  *run=new RAT::DS::Run();
      run_tree->SetBranchAddress("run", &run);
      run_tree->GetEntry(0);
      RAT::DS::PMTInfo *pmtinfo = run->GetPMTInfo();
      double pmt_bound_R = 0;
      double pmt_bound_Z = 0;
      for(int i = 0; i < pmtinfo->GetPMTCount(); i++){
        if(pmtinfo->GetType(i) != 1) continue;
        TVector3 pos = pmtinfo->GetPosition(i);
        if(pos.X()/1000. > pmt_bound_R) pmt_bound_R = pos.X()/1000.;
        if(pos.Z()/1000. > pmt_bound_Z) pmt_bound_Z = pos.Z()/1000.;
      }
      target->SetBoundR(pmt_bound_R);
      target->SetBoundZ(pmt_bound_Z);

      // Create an event reader from the input file
      RAT::DSReader reader(dirname+fname);
      //RAT::DSReader reader("/data/kneale/production/gd_water/run1589492333848056055.root");

      // Loop through the events
      int num_events = reader.GetT()->GetEntries();
      std::vector<int> saved_events;
      int num_cent = 0;
      for(int i = 0; i < num_events; i++){

        RAT::DS::Root *event = reader.NextEvent();

        // Get the monte-carlo information
        RAT::DS::MC *mc = event->GetMC();
        int num_particles = mc->GetMCParticleCount();
        // Get the number of triggered events
        int num_triggers = event->GetEVCount();

        // Ignore if no particles or triggers
        if(num_particles < 1) continue;
        // Get the initial position of the first particle
        double x_mc = mc->GetMCParticle(0)->GetPosition().X()/1000.;
        double y_mc = mc->GetMCParticle(0)->GetPosition().Y()/1000.;
        double z_mc = mc->GetMCParticle(0)->GetPosition().Z()/1000.;
        double r_mc = std::sqrt(std::pow(x_mc, 2)+std::pow(y_mc, 2));
        // Ignore if true vertex not in the fiducial volume
        if(std::abs(r_mc) > 6.7 || std::abs(z_mc) > 6.7) continue;
        //if(mc->GetMCParticle(0)->GetKE() < 2) continue;
        if(num_triggers < 1){ 
          no_trigger++;
          continue;
        }
        num_cent++;

        

        // FIXME try to fix the mess made by simulating both positrons and neutrons
        // Skip events where the trigger time does not match the true time
        //if(std::abs(event->GetEV(0)->GetCalibratedTriggerTime()-mc->GetMCParticle(0)->GetTime())>100) continue;
        // Skip events with low energy positron and only one trigger
        //if(mc->GetMCParticle(0)->GetKE() < 0.5 && event->GetEVCount() == 1) continue;
        RAT::DS::EV *trig = event->GetEV(0);

        // Get the PMTs involved in the trigger and add up the PE
        int num_pmts = trig->GetPMTCount();
        int inner_hits = 0;
        int veto_hits = 0;
        double totPE = 0;
        for(int k = 0; k < num_pmts; k++){
          RAT::DS::PMT *pmt = trig->GetPMT(k);
          totPE += pmt->GetCharge()/1.6;
          int id = pmt->GetID();
          if(pmtinfo->GetType(id) == 1) inner_hits++;
          else veto_hits++;
        }
        //if(inner_hits < 9){ 
        //  no_trigger++;
        //  continue;
        //}
        trigger++;
        saved_events.push_back(i);

        // Save truth info
        target->mc_x.push_back(x_mc);
        target->mc_y.push_back(y_mc);
        target->mc_z.push_back(z_mc);
        TVector3 mom = mc->GetMCParticle(0)->GetMomentum().Unit();
        target->mc_u.push_back(mom.X());
        target->mc_v.push_back(mom.Y());
        target->mc_w.push_back(mom.Z());
        target->mc_time.push_back(mc->GetMCParticle(0)->GetTime());
        target->mc_energy.push_back(mc->GetMCParticle(0)->GetKE());
        target->mc_cherenkov.push_back(mc->GetMCSummary()->GetNumCerenkovPhoton());
        target->mc_scintillation.push_back(mc->GetMCSummary()->GetNumScintPhoton());

        
        // Save trigger info
        target->pmt_hits.push_back(num_pmts);
        target->inner_hits.push_back(inner_hits);
        target->veto_hits.push_back(veto_hits);
        target->reco_pe.push_back(totPE);
        target->reco_time.push_back(trig->GetCalibratedTriggerTime());
        target->n_triggers.push_back(num_triggers);

        // Save centroid fitter info
        RAT::DS::Centroid *centroid = trig->GetCentroid();
        target->centroid_x.push_back(centroid->GetPosition().X()/1000.);
        target->centroid_y.push_back(centroid->GetPosition().Y()/1000.);
        target->centroid_z.push_back(centroid->GetPosition().Z()/1000.);

      }

      // Create an event reader from the input file
      TFile b_file(bdirname+bfname);
      //TFile b_file("/data/brooks/wbls/gd_water/run1589492333848056055_bonsai.root");
      TTreeReader b_reader("data", &b_file);

      TTreeReaderValue<double> x(b_reader, "x");
      TTreeReaderValue<double> y(b_reader, "y");
      TTreeReaderValue<double> z(b_reader, "z");
      TTreeReaderValue<double> u(b_reader, "u");
      TTreeReaderValue<double> v(b_reader, "v");
      TTreeReaderValue<double> w(b_reader, "w");
      TTreeReaderValue<double> n9(b_reader, "n9");
      TTreeReaderValue<int> mcid(b_reader, "mcid");
      TTreeReaderValue<int> subid(b_reader, "subid");

      int num_entries = 0;
      // Loop through the events
      while(b_reader.Next()){

        if(std::find(saved_events.begin(), saved_events.end(), *mcid) == saved_events.end()) continue;
        if(*subid != 0) continue;
        num_entries++;

        target->bonsai_x.push_back(*x/1000.);
        target->bonsai_y.push_back(*y/1000.);
        target->bonsai_z.push_back(*z/1000.);
        target->bonsai_u.push_back(*u);
        target->bonsai_v.push_back(*v);
        target->bonsai_w.push_back(*w);
        target->n9.push_back(*n9);

      }
    }

    targets.push_back(target);
    vec_no_trigger.push_back(no_trigger);
    vec_trigger.push_back(trigger);

  }

  // Apply bias corrections
  //for(size_t i = 0; i < targets.size(); i++){
    //targets[1]->CorrectBias();
  //}

  double min_e = 2;//6
  double max_e = 6;//6
  // Make plots
  // Truth:
  // * True kinetic energy spectrum
  PlotTargets1D(targets, "mc_energy", 20, min_e, max_e);
  // * True time spectrum
  PlotTargets1D(targets, "mc_time", 20, -8, 8);
  // * True number of cherenkov photons
  PlotTargets1D(targets, "mc_cherenkov", 30, 0, 4000);//6000
  // * True number of scintillation photons
  PlotTargets1D(targets, "mc_scintillation", 40, 0.0001, 2000);//4500
  // * Vertex position
  PlotTargets1D(targets, "mc_x", 20, -10, 10);
  PlotTargets1D(targets, "mc_y", 20, -10, 10);
  PlotTargets1D(targets, "mc_z", 20, -10, 10);
  PlotTargets1D(targets, "mc_r", 20, 0, 7);
  PlotTargets1D(targets, "mc_angle", 20, -180, 180);
  PlotTargets1D(targets, "mc_wall", 20, 0, 6);
  // * Vertex direction
  PlotTargets1D(targets, "mc_u", 20, -1, 1);
  PlotTargets1D(targets, "mc_v", 20, -1, 1);
  PlotTargets1D(targets, "mc_w", 20, -1, 1);
  // * True cherenkov/scintillation ratio
  // Trigger:
  // * Number of PMT hits (vs R, Z, KE)
  PlotTargets1D(targets, "inner_hits", 20, 0, 170);
  PlotTargetsMeans(targets, "mc_r", 6, 0, 7, "inner_hits");
  PlotTargetsMeans(targets, "mc_z", 6, -7, 7, "inner_hits");
  PlotTargetsMeans(targets, "mc_energy", 6, 0, max_e, "inner_hits");
  PlotTargetsFraction(targets, "inner_hits", 9, "mc_energy", 80, 0, 4);
  // * Total reconstructed PE (vs R, Z, KE)
  PlotTargets1D(targets, "reco_pe", 30, 0, 120);
  PlotTargetsMeans(targets, "mc_r", 6, 0, 7, "reco_pe");
  PlotTargetsMeans(targets, "mc_z", 6, -7, 7, "reco_pe");
  PlotTargetsMeans(targets, "mc_energy", 6, 0, max_e, "reco_pe");
  // * Time difference between trigger and particle
  PlotTargets1D(targets, "dtime", 30, -20, 800);
  // * Number of triggers
  PlotTargets1D(targets, "n_triggers", 7, 0, 7);
  // Reconstruction:
  // * Energy
  PlotTargets2D(targets, "mc_energy", 12, 0, max_e, "reco_pe", 20, 0, 140);
  PlotTargets2D(targets, "mc_energy", 12, 0, max_e, "inner_hits", 20, 0, 200);
  PlotTargets2D(targets, "mc_energy", 12, 0, max_e, "n9", 20, 0, 40);
  // * Difference between true and reco vertex X, Y, Z, R, theta (bias and resolution vs self, KE, PE)
  PlotTargets1D(targets, "bonsai_dx", 30, -2.5, 2.5);
  PlotTargets1D(targets, "bonsai_dy", 30, -2.5, 2.5);
  PlotTargets1D(targets, "bonsai_dz", 30, -2.5, 2.5);
  PlotTargetsBiasRes(targets, "mc_energy", 10, min_e, max_e, "bonsai_dx", 15, -2.5, 2.5);
  PlotTargetsBiasRes(targets, "mc_energy", 10, min_e, max_e, "bonsai_dy", 15, -2.5, 2.5);
  PlotTargetsBiasRes(targets, "mc_energy", 10, min_e, max_e, "bonsai_dz", 15, -2.5, 2.5);
  PlotTargetsBiasRes(targets, "mc_energy", 10, min_e, max_e, "bonsai_dr", 15, -2.5, 2.5);
  PlotTargetsBiasRes(targets, "bonsai_z", 10, -7, 7, "bonsai_dz", 15, -5, 5);
  PlotTargetsBiasRes(targets, "mc_z", 10, -7, 7, "bonsai_dz", 15, -5, 5);
  PlotTargets1D(targets, "bonsai_dr", 30, -2.5, 2.5);
  PlotTargetsBiasRes(targets, "bonsai_r", 7, 0, 7, "bonsai_dr", 15, -5, 5);
  PlotTargetsBiasRes(targets, "mc_r", 7, 0, 7, "bonsai_dr", 15, -5, 5);
  PlotTargets1D(targets, "bonsai_dangle", 30, -80, 80);
  PlotTargetsBiasRes(targets, "bonsai_angle", 10, -180, 180, "bonsai_dangle", 15, -40, 40);
  PlotTargetsBiasRes(targets, "mc_angle", 10, -180, 180, "bonsai_dangle", 15, -40, 40);
  PlotTargets1D(targets, "bonsai_dwall", 30, -5, 5);
  PlotTargetsBiasRes(targets, "bonsai_wall", 6, 0, 5, "bonsai_dwall", 15, -5, 5);
  PlotTargetsBiasRes(targets, "mc_wall", 6, 0, 5, "bonsai_dwall", 15, -5, 5);
  PlotTargets2D(targets, "mc_energy", 20, 0, max_e, "bonsai_dwall", 30, -5, 5);
  PlotTargets1D(targets, "centroid_dx", 30, -5, 5);
  PlotTargets1D(targets, "centroid_dy", 30, -5, 5);
  PlotTargets1D(targets, "centroid_dz", 30, -5, 5);
  PlotTargetsBiasRes(targets, "centroid_z", 10, -7, 7, "centroid_dz", 15, -5, 5);
  PlotTargetsBiasRes(targets, "mc_z", 10, -7, 7, "centroid_dz", 15, -5, 5);
  PlotTargets1D(targets, "centroid_dr", 30, -5, 5);
  PlotTargetsBiasRes(targets, "centroid_r", 7, 0, 7, "centroid_dr", 15, -5, 5);
  PlotTargetsBiasRes(targets, "mc_r", 7, 0, 7, "centroid_dr", 15, -5, 5);
  PlotTargets1D(targets, "centroid_dangle", 30, -80, 80);
  PlotTargetsBiasRes(targets, "centroid_angle", 10, -180, 180, "centroid_dangle", 15, -40, 40);
  PlotTargetsBiasRes(targets, "mc_angle", 10, -180, 180, "centroid_dangle", 15, -40, 40);
  PlotTargets1D(targets, "centroid_dwall", 30, -5, 5);
  PlotTargetsBiasRes(targets, "centroid_wall", 5, 0, 5, "centroid_dwall", 15, -5, 5);
  PlotTargetsBiasRes(targets, "mc_wall", 5, 0, 5, "centroid_dwall", 15, -5, 5);
  PlotTargets2D(targets, "mc_energy", 20, 0, max_e, "centroid_dwall", 30, -5, 5);
  //Temp
  PlotTargets2D(targets, "mc_energy", 6, 0, max_e, "centroid_dx", 30, -10, 10);
  PlotTargets2D(targets, "mc_energy", 6, 0, max_e, "centroid_dy", 30, -10, 10);
  PlotTargets2D(targets, "mc_energy", 6, 0, max_e, "centroid_dz", 30, -10, 10);
  PlotTargets2D(targets, "mc_energy", 6, 0, max_e, "bonsai_dx", 30, -10, 10);
  PlotTargets2D(targets, "mc_energy", 6, 0, max_e, "bonsai_dy", 30, -10, 10);
  PlotTargets2D(targets, "mc_energy", 6, 0, max_e, "bonsai_dz", 30, -10, 10);
  // * Vertex direction
  PlotTargets1D(targets, "bonsai_ddir", 20, 0, 180);
  PlotTargetsMeans(targets, "mc_energy", 6, 0, max_e, "bonsai_ddir");
  // * Distance between true and reco vertex (mean distance vs KE, PE)
  PlotTargets1D(targets, "bonsai_dist", 30, 0, 8);
  PlotTargets2D(targets, "mc_energy", 20, 0, max_e, "bonsai_dist", 30, 0, 8);
  //PlotTargetsMeans(targets, "mc_energy", 6, 0, max_e, "bonsai_dist", 10);
  PlotTargetsMeanFit(targets, "mc_energy", 6, 0, max_e, "bonsai_dist", 30, 0, 10);
  //PlotTargetsMeans(targets, "reco_pe", 6, 0, 100, "bonsai_dist");
  PlotTargetsMeanFit(targets, "reco_pe", 6, 0, 100, "bonsai_dist", 30, 0, 10);
  PlotTargets1D(targets, "centroid_dist", 30, 0, 8);
  PlotTargets2D(targets, "mc_energy", 20, 0, max_e, "centroid_dist", 30, 0, 8);
  PlotTargetsMeanFit(targets, "mc_energy", 6, 0, max_e, "centroid_dist", 30, 0, 10);
  PlotTargetsMeanFit(targets, "reco_pe", 6, 0, 100, "centroid_dist", 30, 0, 10);
  // Comparison:
  // * Mean and standard deviation from Gaussian fit
  PlotFittersMeanFit(targets, "mc_energy", 5, 0, max_e, "dist", 30, 0, 10);
  PlotFittersMeanFit(targets[0], 0, "reco_pe", 5, 0, 40, "dist", 30, 0, 10);
  PlotFittersMeanFit(targets[1], 1, "reco_pe", 5, 0, 60, "dist", 30, 0, 10);
  PlotFittersMeanFit(targets[2], 2, "reco_pe", 5, 0, 60, "dist", 30, 0, 10);
  PlotFittersMeanFit(targets[3], 3, "reco_pe", 5, 0, 100, "dist", 30, 0, 10);
  PlotFittersMeanFit(targets[4], 4, "reco_pe", 5, 0, 140, "dist", 30, 0, 10);
  // * Arithmetic mean and standard deviation
  PlotFittersMeans(targets, "mc_energy", 6, 0, max_e, "dist", 5);
  PlotFittersMeans(targets[0], 0, "reco_pe", 5, 0, 40, "dist", 5);
  PlotFittersMeans(targets[1], 1, "reco_pe", 5, 0, 60, "dist", 5);
  PlotFittersMeans(targets[2], 2, "reco_pe", 5, 0, 60, "dist", 5);
  PlotFittersMeans(targets[3], 3, "reco_pe", 5, 0, 100, "dist", 5);
  PlotFittersMeans(targets[4], 4, "reco_pe", 5, 0, 140, "dist", 5);
  // * 50% cumulative distance with +/- 10% errors
  PlotFittersCumulative(targets, "mc_energy", 6, 0, max_e, "dist", 0.5, 0.1);
  PlotTargetsCumulative(targets, "mc_energy", 6, 0, max_e, "bonsai_dist", 0.5, 0.1);
  PlotTargetsCumulative(targets, "reco_pe", 12, 0, 120, "bonsai_dist", 0.5, 0.1);
  PlotFittersCumulative(targets[0], 0, "reco_pe", 5, 0, 40, "dist", 0.5, 0.1); //40
  PlotFittersCumulative(targets[1], 1, "reco_pe", 5, 0, 120, "dist", 0.5, 0.1); //60
  PlotFittersCumulative(targets[2], 2, "reco_pe", 5, 0, 60, "dist", 0.5, 0.1); //60
  PlotFittersCumulative(targets[3], 3, "reco_pe", 5, 10, 60, "dist", 0.5, 0.1); //100
  PlotFittersCumulative(targets[4], 4, "reco_pe", 5, 20, 80, "dist", 0.5, 0.1); //140
  // * Resolution of R, Z, angle vertex reconstruction
  PlotFittersBiasRes(targets, "mc_r", 5, 0, 7, "dr", 15, -5, 5);
  PlotFittersBiasRes(targets, "mc_wall", 6, 0, 4, "dwall", 15, -5, 5);
  PlotFittersBiasRes(targets, "mc_energy", 6, 0, max_e, "dwall", 15, -5, 5);
  PlotFittersBiasRes(targets, "mc_z", 10, -7, 7, "dz", 15, -5, 5);
  PlotFittersBiasRes(targets, "mc_angle", 10, -180, 180, "dangle", 15, -40, 40);
  PlotFittersBiasRes(targets[0], 0, "reco_pe", 6, 0, 40, "dwall", 15, -5, 5);
  PlotFittersBiasRes(targets[1], 1, "reco_pe", 6, 0, 60, "dwall", 15, -5, 5);
  PlotFittersBiasRes(targets[2], 2, "reco_pe", 6, 0, 60, "dwall", 15, -5, 5);
  PlotFittersBiasRes(targets[3], 3, "reco_pe", 6, 0, 90, "dwall", 15, -5, 5);
  PlotFittersBiasRes(targets[4], 4, "reco_pe", 6, 0, 140, "dwall", 15, -5, 5);

  for(int i = 0; i < files.size(); i++){
    std::cout<<files[i]<<" no trigger = "<<vec_no_trigger[i]<<" trigger = "<<vec_trigger[i]<<" percent triggered = "<<(double)vec_trigger[i]/(vec_trigger[i]+vec_no_trigger[i])<<"\n";
  }

}
