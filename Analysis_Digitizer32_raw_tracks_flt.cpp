// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include "TRandom.h"
#include "TMath.h"
#include "TVirtualFFT.h"
#include <unordered_map>
#include <memory>
#include <limits> // For numeric limits

// Standard includes
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>

// To compile:
// clang++ -std=c++17 -O3 -Wall -Wextra    `root-config --cflags --libs`   Analysis_Digitizer32_raw_tracks_flt.cpp -o Analysis_Digitizer32_raw_tracks_flt

// Lorenzo includes
#include "Analysis_Digitizer16.h"

using namespace std;

// Fable5 edition of Amplitudes(): same algorithm, but the TF1 and TGraph are
// created once and reused. Building a formula-based TF1 at every call costs
// far more than the fit itself on ~20 points.
void AmplitudesF5(Int_t camp, Float_t amp[], Float_t timeS[], Float_t RMSbck,
                  Int_t t_start, Int_t t_stop,
                  Float_t *pamp, Float_t *tamp, Float_t *chi)
{
  int j_max = 0;
  int j_start = 0;
  int n_j = 0;
  // Limit of gaus fitting around j_max
  const float AmpFrac = 0.6;
  double max = -100000;
  double tmax = -100000;

  double npbck = 0;
  double bck = 0;
  for (int j = 0; j < t_start; j++)
    {
      bck += amp[j];
      npbck++;
    }
  if (npbck > 0) bck /= npbck;

  // Maximum calculation
  for (int j = t_start; j < t_stop; j++)
    {
      if (amp[j] - bck > max)
	{
	  max = amp[j] - bck;
	  tmax = (timeS[j+1]*amp[j+1] + timeS[j]*amp[j] + timeS[j-1]*amp[j-1]) /
	         (amp[j+1] + amp[j] + amp[j-1]);
	  j_max = j;
	}
    }

  if (max < 5*RMSbck)
    {
      *pamp = max;
      *tamp = tmax;
      *chi = -1;
      return;
    }

  // Find limits at AmpFrac: from the maximum backwards
  for (int j = j_max; j > 0; j--)
    {
      if (amp[j] - bck < AmpFrac*max)
	{
	  j_start = j;
	  break;
	}
    }

  // from the lower limit, up the maximum, down to the limit on the falling edge
  for (int j = j_start+1; j < camp-1; j++)
    {
      n_j++;
      if (amp[j] - bck < AmpFrac*max)
	break;
    }

  static TF1    *f1 = new TF1("gausfitF5","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",-1.e6,1.e6);
  static TGraph *g  = new TGraph();

  g->Set(n_j);
  for (int j = 0; j < n_j; j++)
    g->SetPoint(j, timeS[j+j_start], amp[j+j_start] - bck);

  f1->SetParameters(max, tmax, (timeS[j_start]-timeS[j_start+n_j])/2.);
  g->Fit(f1, "QN", "", timeS[j_start]-1, timeS[j_start+n_j-1]+1);

  *pamp = f1->GetParameter(0);
  *tamp = f1->GetParameter(1);
  *chi  = f1->GetChisquare()/f1->GetNDF();

  if (*chi > 50 || fabs((*pamp - max)/max) > 0.2)
    {
      *pamp = max;
      *tamp = tmax;
      *chi = -2;
    }
}

// Fable5 edition of mobileAVG(): running sum, O(N) independent of navg.
// navg<=1 (including the negative values used to flag the FFT filter) copies
// the input unchanged instead of indexing out of bounds.
void mobileAVGF5(Int_t camp, Float_t amp[], Int_t navg, Float_t m_amp[])
{
  const int half = (navg-1)/2;

  if (navg <= 1 || camp < 2*half+2)
    {
      for (int j = 0; j < camp; j++) m_amp[j] = amp[j];
      return;
    }

  double sum = 0;
  for (int k = 0; k < 2*half+1; k++) sum += amp[k];
  m_amp[half] = sum/float(navg);

  for (int j = half+1; j < camp-half-1; j++)
    {
      sum += amp[j+half] - amp[j-half-1];
      m_amp[j] = sum/float(navg);
    }

  for (int j = 0; j < half; j++)             m_amp[j] = m_amp[half];
  for (int j = camp-half-1; j < camp; j++)   m_amp[j] = m_amp[camp-half-2];
}

int main()
{

  bool showFFT=false;
  bool doFit = false;
  bool FWaveform = true;

  TFile *OutputFile = new TFile("RunDigitizer32_raw_tracks_flt.root","recreate");
  OutputFile->SetCompressionAlgorithm(ROOT::kLZ4); // faster writes, slightly larger file

  ifstream InputCARD("Input_Folder_Digitizer32_raw_tracks.txt");

  if(!InputCARD)
    {
      cout << "Error: could not find the InputCARD" << endl;
      return (0);
    }

  // read the active number of channels
  int ntmp;
  int np_acq; // defined reading the file
  int np_offset = 0; // 500; //100;
  int np_Max = 1000;

  nMCP = 7;

  bool Interpolation = true;

  event = 0;

  // Levels array
  const int NArrayL = 11;
  float ThresholdL[NArrayL] = {10, 15,20,30,40,50,60,80,100,200,300};
  float TValuesL[NArrayL];
  int NTValuesL[NArrayL];
  int NArr = 0;

  // Levels trail array
  const int TrailNArrayL = 12;
  float TrailThresholdL[TrailNArrayL] = {10,20,30,40,50,60,80,100,120,140,160,180};
  float TrailTValuesL[TrailNArrayL];
  int TrailNTValuesL[TrailNArrayL];
  int TrailNArr = 0;

  // fraction array
  const int NArrayF = 12;
  float ThresholdF[NArrayF] = {0.05, 0.1,0.15,0.20,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9};
  float TValuesF[NArrayF];
  int NTValuesF[NArrayF];

  // fraction trail array
  const int TrailNArrayF = 5;
  float TrailThresholdF[TrailNArrayF] = {0.1,0.3,0.5,0.7,0.9};
  float TrailTValuesF[TrailNArrayF];
  int TrailNTValuesF[TrailNArrayF];

  // Active-channel mask: skips the per-sample analysis loops for channels
  // with no real data, to avoid wasting time on empty waveforms.
  int AcCh[36] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,   // ch 0-16: waveforms used
		  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,        // ch 17-31: unused
		  1,1,1,1};                             // ch 32-35: trigger

  //samples in t
  InputCARD >> pip >> ntmp;
  nchro = ntmp;

  InputCARD >> pip >> ntmp;
  FWaveform = ntmp;

  if (FWaveform) cout << "Acquisition stores waveforms" << endl;
  else
    cout << "Acquisition does not store waveforms" << endl;

  TTree *OutTree = new TTree("Analysis","Analysis");
  OutTree->SetDirectory(OutputFile);

  // Analysis branches
  OutTree->Branch("ntrig",&ntrig,"ntrig/I");
  OutTree->Branch("event",&event,"event/I");
  OutTree->Branch("nrun",&nrun,"nrun/I");
  OutTree->Branch("nchro",&nchro,"nchro/I");

  OutTree->Branch("xtrk1",&xtrk1,"xtrk1/F");
  OutTree->Branch("ytrk1",&ytrk1,"ytrk1/F");
  OutTree->Branch("xtrk2",&xtrk2,"xtrk2/F");
  OutTree->Branch("ytrk2",&ytrk2,"ytrk2/F");
  OutTree->Branch("trkchi2",&trkchi2,"trkchi2/F");


  OutTree->Branch("t_I",t_I,"t_I[nchro]/F");
  OutTree->Branch("t_II",t_II,"t_II[nchro]/F");

  OutTree->Branch("area_pos",area_pos,"area_pos[nchro]/F");
  OutTree->Branch("area",area,"area[nchro]/F");
  OutTree->Branch("ampl",ampl,"ampl[nchro]/F");
  OutTree->Branch("nampl",nampl,"nampl[nchro]/F");
  OutTree->Branch("pol",pol,"pol[nchro]/I");
  OutTree->Branch("t_max_sample",t_max_sample,"t_max_sample[nchro]/F");
  OutTree->Branch("ampl_sample",ampl_sample,"ampl_sample[nchro]/F");
  OutTree->Branch("ampl_chi2",ampl_chi2,"ampl_chi2[nchro]/F");
  OutTree->Branch("dVdt3070",dVdt3070,"dVdt3070[nchro]/F");
  OutTree->Branch("dVdt1030",dVdt1030,"dVdt1030[nchro]/F");
  OutTree->Branch("dVdt2080",dVdt2080,"dVdt2080[nchro]/F");

  OutTree->Branch("ToT20",ToT20,"ToT20[nchro]/F");
  OutTree->Branch("ToT30",ToT30,"ToT30[nchro]/F");
  OutTree->Branch("ToT40",ToT40,"ToT40[nchro]/F");

  OutTree->Branch("bck",bck,"bck[nchro]/F");
  OutTree->Branch("max_bck_before",max_bck_before,"max_bck_before[nchro]/F");
  OutTree->Branch("rms_bck_before",rms_bck_before,"rms_bck_before[nchro]/F");
  OutTree->Branch("max_bck_after",max_bck_after,"max_bck_after[nchro]/F");
  OutTree->Branch("rms_bck_after",rms_bck_after,"rms_bck_after[nchro]/F");


  OutTree->Branch("t_cent",t_cent,"t_cent[nchro]/F");

  OutTree->Branch("t_level10",t_level10,"t_level10[nchro]/F");
  OutTree->Branch("t_level15",t_level15,"t_level15[nchro]/F");
  OutTree->Branch("t_level20",t_level20,"t_level20[nchro]/F");
  OutTree->Branch("t_level30",t_level30,"t_level30[nchro]/F");
  OutTree->Branch("t_level40",t_level40,"t_level40[nchro]/F");
  OutTree->Branch("t_level50",t_level50,"t_level50[nchro]/F");
  OutTree->Branch("t_level60",t_level60,"t_level60[nchro]/F");
  OutTree->Branch("t_level80",t_level80,"t_level80[nchro]/F");
  OutTree->Branch("t_level100",t_level100,"t_level100[nchro]/F");
  OutTree->Branch("t_level200",t_level200,"t_level200[nchro]/F");
  OutTree->Branch("t_level300",t_level300,"t_level300[nchro]/F");

  OutTree->Branch("trail_t_level10",trail_t_level10,"trail_t_level10[nchro]/F");
  OutTree->Branch("trail_t_level20",trail_t_level20,"trail_t_level20[nchro]/F");
  OutTree->Branch("trail_t_level30",trail_t_level30,"trail_t_level30[nchro]/F");
  OutTree->Branch("trail_t_level40",trail_t_level40,"trail_t_level40[nchro]/F");
  OutTree->Branch("trail_t_level50",trail_t_level50,"trail_t_level50[nchro]/F");
  OutTree->Branch("trail_t_level60",trail_t_level60,"trail_t_level60[nchro]/F");
  OutTree->Branch("trail_t_level80",trail_t_level80,"trail_t_level80[nchro]/F");
  OutTree->Branch("trail_t_level100",trail_t_level100,"trail_t_level100[nchro]/F");
  OutTree->Branch("trail_t_level120",trail_t_level120,"trail_t_level120[nchro]/F");
  OutTree->Branch("trail_t_level140",trail_t_level140,"trail_t_level140[nchro]/F");
  OutTree->Branch("trail_t_level160",trail_t_level160,"trail_t_level160[nchro]/F");
  OutTree->Branch("trail_t_level180",trail_t_level180,"trail_t_level180[nchro]/F");


  OutTree->Branch("cfd05",cfd05,"cfd05[nchro]/F");
  OutTree->Branch("cfd10",cfd10,"cfd10[nchro]/F");
  OutTree->Branch("cfd15",cfd15,"cfd15[nchro]/F");
  OutTree->Branch("cfd20",cfd20,"cfd20[nchro]/F");
  OutTree->Branch("cfd25",cfd25,"cfd25[nchro]/F");
  OutTree->Branch("cfd30",cfd30,"cfd30[nchro]/F");
  OutTree->Branch("cfd40",cfd40,"cfd40[nchro]/F");
  OutTree->Branch("cfd50",cfd50,"cfd50[nchro]/F");
  OutTree->Branch("cfd60",cfd60,"cfd60[nchro]/F");
  OutTree->Branch("cfd70",cfd70,"cfd70[nchro]/F");
  OutTree->Branch("cfd80",cfd80,"cfd80[nchro]/F");
  OutTree->Branch("cfd90",cfd90,"cfd90[nchro]/F");

  OutTree->Branch("Ncfd05",Ncfd05,"Ncfd05[nchro]/I");
  OutTree->Branch("Ncfd10",Ncfd10,"Ncfd10[nchro]/I");
  OutTree->Branch("Ncfd15",Ncfd15,"Ncfd15[nchro]/I");
  OutTree->Branch("Ncfd20",Ncfd20,"Ncfd20[nchro]/I");
  OutTree->Branch("Ncfd25",Ncfd25,"Ncfd25[nchro]/I");
  OutTree->Branch("Ncfd30",Ncfd30,"Ncfd30[nchro]/I");
  OutTree->Branch("Ncfd40",Ncfd40,"Ncfd40[nchro]/I");
  OutTree->Branch("Ncfd50",Ncfd50,"Ncfd50[nchro]/I");
  OutTree->Branch("Ncfd60",Ncfd60,"Ncfd60[nchro]/I");
  OutTree->Branch("Ncfd70",Ncfd70,"Ncfd70[nchro]/I");
  OutTree->Branch("Ncfd80",Ncfd80,"Ncfd80[nchro]/I");
  OutTree->Branch("Ncfd90",Ncfd90,"Ncfd90[nchro]/I");

  OutTree->Branch("trail_Ncfd10",trail_Ncfd10,"trail_Ncfd10[nchro]/I");
  OutTree->Branch("trail_Ncfd30",trail_Ncfd30,"trail_Ncfd30[nchro]/I");
  OutTree->Branch("trail_Ncfd50",trail_Ncfd50,"trail_Ncfd50[nchro]/I");
  OutTree->Branch("trail_Ncfd70",trail_Ncfd70,"trail_Ncfd70[nchro]/I");
  OutTree->Branch("trail_Ncfd90",trail_Ncfd90,"trail_Ncfd90[nchro]/I");


  OutTree->Branch("t_rms3",t_rms3,"t_rms3[nchro]/F");
  OutTree->Branch("t_rms5",t_rms5,"t_rms5[nchro]/F");
  OutTree->Branch("trail_cfd10",trail_cfd10,"trail_cfd10[nchro]/F");
  OutTree->Branch("trail_cfd30",trail_cfd30,"trail_cfd30[nchro]/F");
  OutTree->Branch("trail_cfd50",trail_cfd50,"trail_cfd50[nchro]/F");
  OutTree->Branch("trail_cfd70",trail_cfd70,"trail_cfd70[nchro]/F");
  OutTree->Branch("trail_cfd90",trail_cfd90,"trail_cfd90[nchro]/F");

  OutTree->Branch("samples",samples,"samples[nchro]/I");
  OutTree->Branch("t_max",t_max,"t_max[nchro]/F");


  if(doFit)
    {
      OutTree->Branch("t_zero",t_zero,"t_zero[nchro]/F");
      OutTree->Branch("rise_lin0",rise_lin0,"rise_lin0[nchro]/F");
      OutTree->Branch("rise_lin1",rise_lin1,"rise_lin1[nchro]/F");
      OutTree->Branch("rise_lin_chi2",&rise_lin_chi2,"rise_lin_chi2[nchro]/F");
      OutTree->Branch("rise_exp0",rise_exp0,"rise_exp0[nchro]/F");
      OutTree->Branch("rise_exp1",rise_exp1,"rise_exp1[nchro]/F");
      OutTree->Branch("rise_exp_chi2",rise_exp_chi2,"rise_exp_chi2[nchro]/F");
    }

  cout << "nchro = " << nchro << endl;
  if (FWaveform)
    {

      for(i=0;i<nchro;i++)
	{

	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "samples_" << i;
	  leafl << "samples_" << i <<"/I";
	  OutTree->Branch(leaf.str().c_str(),&samples[i],leafl.str().c_str());
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "samplesrec_" << i;
	  leafl << "samplesrec_" << i <<"/I";
	  OutTree->Branch(leaf.str().c_str(),&samplesrec[i],leafl.str().c_str());
	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();

	  leaf << "amp" << i;
	  leafl << "amp" << i <<"[samplesrec_" << i << "]/F";
	  OutTree->Branch(leaf.str().c_str(),&amprec[i][0],leafl.str().c_str());

	  leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	  leaf << "m_amp" << i;
	  leafl << "m_amp" << i <<"[samplesrec_" << i << "]/F";
	  OutTree->Branch(leaf.str().c_str(),&m_amprec[i][0],leafl.str().c_str());
	  cout << " =  " << leafl.str().c_str() << endl;

	  if(showFFT)
	    {
	      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	      leaf << "FFT_abs" << i;
	      leafl << "FFT_abs" << i <<"[samplesrec_" << i << "]/F";
	      OutTree->Branch(leaf.str().c_str(),&FFT_abs[i][0],leafl.str().c_str());

	      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	      leaf << "FFT_real" << i;
	      leafl << "FFT_real" << i <<"[samplesrec_" << i << "]/F";
	      OutTree->Branch(leaf.str().c_str(),&FFT_real[i][0],leafl.str().c_str());

	      leaf.str("");leaf.clear();leafl.str("");leafl.clear();
	      leaf << "FFT_comp" << i;
	      leafl << "FFT_comp" << i <<"[samplesrec_" << i << "]/F";
	      OutTree->Branch(leaf.str().c_str(),&FFT_comp[i][0],leafl.str().c_str());
	    }

	}

      OutTree->Branch("time",timerec,"time[samples_0]/F");

    }

 cout << "Branch settings done" << endl;


  //Loop on different runs
  while(1)
    {
      pip="";
      toffee="";
      if(InputCARD.eof())
	break;

      //Read run variables; a failed read (missing -1 sentinel) also ends the loop
      if(!(InputCARD >> pip >> nrun >> pip >> ntrig >> pip >> nchro >> pip >> MaxEvt >> pip >> nMCP))
	break;
      if (nrun == -1 ) break;
      cout << endl << endl << "Run " << nrun << "\t " << nchro << " channels" << " Events to be analized = " << MaxEvt - ntrig << endl;
      cout << " MCP on ch = " << nMCP << endl;

      // Read tracker text file
      InputCARD >> toffee >> TrkFile;

      std::ifstream Filein;
      Filein.open(TrkFile, std::ios::in);

      if(!Filein.is_open())
	{
	  std::cout << "Tracker file failed: " << TrkFile << std::endl;
	  exit(1);
	}
      else { std::cout << "Opened file " << TrkFile << std::endl;
             std::string _hdr; std::getline(Filein, _hdr); } // skip CSV header

      //Note: ntrig is the number of the FIRST trigger of the run

      //Inizialization of data run variables
      ntrig=ntrig-1;

      ostringstream convert;
      convert << nrun ;
      String_nrun = convert.str();

      //Inizialization of program run variables
      running=1;

      //Inizialization of channel variables for this run
      InputCARD >> toffee >> InputNA;
      snprintf(title,sizeof(title),"%s",InputNA.c_str());
      InputCARD >> pip >> nme;

      cout << endl << " File name:" << endl << title << " averaging " << nme << " points " << endl;

      for (int c=0;c<NCHRO;c++) FreqCut[c] = 0;

      if (nme<0)
	{
	  FreqCut[0] = -nme*1.e6;
	  FreqCut[2] = -nme*1.e6;
	  FreqCut[3] = -nme*1.e6;
	  FreqCut[1] = 0;
	  nmedia[1] = 1;
	  cout << endl << "Channel 0  with Freq. Cut = " << FreqCut[0]*1e-6 << " MHz" <<endl;
	  cout << "Channel 1  with  Freq. Cut  = " <<  FreqCut[1]*1e-6 << " MHz"  <<endl;
	}

      // Single digitizer — TTreeReader with float arrays
      TFile *f = new TFile(title);
      TTree *tree = (TTree*)f->Get("wfm");

      TTreeReader reader(tree);
      TTreeReaderValue<Float_t> wfmSize(reader, "size");
      TTreeReaderArray<float> w0(reader,  "w0");
      TTreeReaderArray<float> w1(reader,  "w1");
      TTreeReaderArray<float> w2(reader,  "w2");
      TTreeReaderArray<float> w3(reader,  "w3");
      TTreeReaderArray<float> w4(reader,  "w4");
      TTreeReaderArray<float> w5(reader,  "w5");
      TTreeReaderArray<float> w6(reader,  "w6");
      TTreeReaderArray<float> w7(reader,  "w7");
      TTreeReaderArray<float> w8(reader,  "w8");
      TTreeReaderArray<float> w9(reader,  "w9");
      TTreeReaderArray<float> w10(reader, "w10");
      TTreeReaderArray<float> w11(reader, "w11");
      TTreeReaderArray<float> w12(reader, "w12");
      TTreeReaderArray<float> w13(reader, "w13");
      TTreeReaderArray<float> w14(reader, "w14");
      TTreeReaderArray<float> w15(reader, "w15");
      TTreeReaderArray<float> w16(reader, "w16");
      TTreeReaderArray<float> w17(reader, "w17");
      TTreeReaderArray<float> w18(reader, "w18");
      TTreeReaderArray<float> w19(reader, "w19");
      TTreeReaderArray<float> w20(reader, "w20");
      TTreeReaderArray<float> w21(reader, "w21");
      TTreeReaderArray<float> w22(reader, "w22");
      TTreeReaderArray<float> w23(reader, "w23");
      TTreeReaderArray<float> w24(reader, "w24");
      TTreeReaderArray<float> w25(reader, "w25");
      TTreeReaderArray<float> w26(reader, "w26");
      TTreeReaderArray<float> w27(reader, "w27");
      TTreeReaderArray<float> w28(reader, "w28");
      TTreeReaderArray<float> w29(reader, "w29");
      TTreeReaderArray<float> w30(reader, "w30");
      TTreeReaderArray<float> w31(reader, "w31");
      TTreeReaderArray<float> trg0(reader, "trg0");
      TTreeReaderArray<float> trg1(reader, "trg1");
      TTreeReaderArray<float> trg2(reader, "trg2");
      TTreeReaderArray<float> trg3(reader, "trg3");

      // Channel index -> waveform branch, replaces the 36-way if-else per sample
      TTreeReaderArray<float>* wf[NCHRO] = {
	&w0,&w1,&w2,&w3,&w4,&w5,&w6,&w7,&w8,&w9,&w10,&w11,&w12,&w13,&w14,&w15,
	&w16,&w17,&w18,&w19,&w20,&w21,&w22,&w23,&w24,&w25,&w26,&w27,&w28,&w29,&w30,&w31,
	&trg0,&trg1,&trg2,&trg3};

      Long64_t nentries = tree->GetEntries();

      cout << "Number of events = " << nentries << endl;
      // Valid entries are 0..nentries-1: clamp MaxEvt so SetEntry never runs past the end
      if (MaxEvt == -1 || MaxEvt > nentries-1) MaxEvt = nentries-1;

      bool SamplesChecked = false;

      emptycount=0;

      // here we are looping on each event
      while(running)
	{

	  ntrig++;
	  event++;


	  if (MaxEvt > 0 && ntrig >MaxEvt) goto Write;
	  if(event==1)
	    cout << "Event = " << event << " at trigger = " << ntrig << endl;
	  else if(ntrig%200==0)
	    cout << "Event = " << event << " at trigger = " << ntrig << endl;

	  if (reader.SetEntry(ntrig) != TTreeReader::kEntryValid)
	    {
	      cout << "TTreeReader: entry " << ntrig << " not valid, closing the run" << endl;
	      goto Write;
	    }

	  for(i=0;i<nchro;i++)
	    {
	      nmedia[i] = nme;
	      samples[i]=0;
	      samplesrec[i]=0;
	      AMax[i]=0;
	      TMax[i]=0;
	      NMax[i]=0;
	      t_max[i] = -1000;
	      max_bck_after[i]= -1000;
	      rms_bck_after[i] = -1000;
	      max_bck_before[i]= -1000;
	      rms_bck_before[i] = -1000;
	      ampl[i] = -1000;
	      nampl[i] = 1000;
	      ampl_chi2[i] = -1000;
	      bck[i] = -1000;
	      Fakebck[i] = 0;
	      polarity[i] = 0;
	      pol[i] = 0;
	      DTCh0i[i] = 0;
	      t_max_sample[i] = 0;
	      ampl_sample[i] = 0;
	      area[i] = 0;
	      area_pos[i] = 0;
	      dVdt1030[i] = 0; dVdt3070[i] = 0; dVdt2080[i] = 0;
	      ToT20[i] = 0; ToT30[i] = 0; ToT40[i] = 0;
	      t_cent[i] = 0; t_rms3[i] = 0; t_rms5[i] = 0;
	      t_level10[i]=0; t_level15[i]=0; t_level20[i]=0; t_level30[i]=0;
	      t_level40[i]=0; t_level50[i]=0; t_level60[i]=0; t_level80[i]=0;
	      t_level100[i]=0; t_level200[i]=0; t_level300[i]=0;
	      trail_t_level10[i]=0; trail_t_level20[i]=0; trail_t_level30[i]=0;
	      trail_t_level40[i]=0; trail_t_level50[i]=0; trail_t_level60[i]=0;
	      trail_t_level80[i]=0; trail_t_level100[i]=0; trail_t_level120[i]=0;
	      trail_t_level140[i]=0; trail_t_level160[i]=0; trail_t_level180[i]=0;
	      cfd05[i]=0; cfd10[i]=0; cfd15[i]=0; cfd20[i]=0; cfd25[i]=0; cfd30[i]=0;
	      cfd40[i]=0; cfd50[i]=0; cfd60[i]=0; cfd70[i]=0; cfd80[i]=0; cfd90[i]=0;
	      Ncfd05[i]=0; Ncfd10[i]=0; Ncfd15[i]=0; Ncfd20[i]=0; Ncfd25[i]=0; Ncfd30[i]=0;
	      Ncfd40[i]=0; Ncfd50[i]=0; Ncfd60[i]=0; Ncfd70[i]=0; Ncfd80[i]=0; Ncfd90[i]=0;
	      trail_cfd10[i]=0; trail_cfd30[i]=0; trail_cfd50[i]=0; trail_cfd70[i]=0; trail_cfd90[i]=0;
	      trail_Ncfd10[i]=0; trail_Ncfd30[i]=0; trail_Ncfd50[i]=0; trail_Ncfd70[i]=0; trail_Ncfd90[i]=0;
	    }

	  chmax = -10;
	  Amax = -10;

	  if (ntrig>1) {
	    std::string _ln; std::getline(Filein, _ln);
	    std::istringstream _ss(_ln); std::string _t;
	    float _v[8]={}; int _nv=0;
	    while(std::getline(_ss,_t,'\t') && _nv<8)
	      try{_v[_nv++]=std::stof(_t);}catch(...){}
	    ntrk=(int)_v[0]; xtrk1=_v[1]; ytrk1=_v[2]; xtrk2=_v[3]; ytrk2=_v[4]; trkchi2=_v[7];
	    if (xtrk1==0 && ytrk1==0) {
	      std::cout << "Error: track coordinates x=y=0 at trigger " << ntrig
	                << " -- wrong tracker file format?" << std::endl;
	      exit(1);
	    }
	  }

	  //	  cout << ntrk << " " << xtrk1 << " " << ytrk1 << endl;
	  if (xtrk1<-800) goto NoRoot;


	  for(i=0;i<nchro;i++)
	    {
	      if (!AcCh[i]) continue;

	      np=0;
	      nprec=0;

	      if (i ==0 && !SamplesChecked)
		{
		  SamplesChecked = true;
		  treesize = (int)(*wfmSize);
		  cout << "Number of samples in the acquisition = " << treesize << endl;

		  np_acq = treesize;

		  cout << "Number of samples in the analysis = " << np_Max-np_offset << " starting at samples " << np_offset << endl;
		  if (np_Max>np_acq)
		    {
		      cout << "Number of samples in analysis larger than in the acquisition = " << np_Max  << " > " << np_acq << endl;
		      cout << "The program STOPS" << endl;
		      return 1;
		    }

		}

	      {
		TTreeReaderArray<float> &wch = *wf[i];
		for (int npp =np_offset;npp<np_Max;npp++)
		  {
		    timeS[np] =  0.2*(npp);
		    amp[i][np] = wch[npp]*ADCmV;
		    np++;
		  }
	      }

	      DT = (timeS[10]-timeS[9]); // delta T in nanosecond

	      samples[i]=np-1;

	      AMax[i]=0;

	      for (np = 0;np<samples[i]-10;np++)
		{

		  if ( fabs(amp[i][np]-amp[i][0])>AMax[i])
		    {
		      AMax[i] = fabs(amp[i][np]-amp[i][0]);
		      TMax[i] = timeS[np];
		      NMax[i] = np;
		    }

		}
	      if (AMax[i] > 0)
		polarity[i] = (amp[i][NMax[i]]-amp[i][0])/AMax[i];
	      else
		polarity[i] = 1; // flat waveform: keep the sign neutral

	      pol[i] =  polarity[i];

	      if (ntrig == 2  && i == 0) cout << "Number of points: " << np  << " with a time step of " << DT << " [ns] " << endl;

	      if (ntrig <20  ) cout << "Event: " << event << " channel " << i << " has maximum value of " << AMax[i] << " [mV] at  "
				    <<  TMax[i] << " [ns] " << " on sample " << NMax[i] << " polarity = " << polarity[i] << endl;

	      if ( TMax[i] ==0 ||  AMax[i]>10000)
		{
		  cout << "Maximum at 0 ns, or signal too large ==> something must be wrong. The program skips this channel" << endl;
		  continue;
		}


	    } // end of preparation loop

	  AMaxSignal = 0;
	  for(i=0;i<nchro;i++)
	       {
		 t_max_sample[i] = TMax[i];
		 ampl_sample[i] = AMax[i];

		 if (i!=nMCP && i<32 &&  AMax[i] > AMaxSignal)
		   {
		     AMaxSignal = AMax[i];
		     TMaxSignal = TMax[i];
		     ChMaxSignal = i;
		   }
	       }
	  if (ntrig <20  )
	    cout <<" MaxSignal on channel " << ChMaxSignal << " Amplitude = " <<  AMaxSignal << " at t = " <<   TMaxSignal << endl;

	  for(i=0;i<nchro;i++)
	    {
	      // to do: rivedere il timing dei canali 16-31
	      // rivedere il timing 32-35, l'ampiezza viene molto piccola
	      //
	      t_I[i] = 50;
	      t_II[i] =70;

	      if (i == nMCP)
		{
		  t_I[i] = 60;
		  t_II[i] = 80;
		}

	      if (i == 32 || i == 33 || i == 34 || i == 35)
		{
		  t_I[i] = 140;
		  t_II[i] = 200;
		}

	    }

	  for(i=0;i<nchro;i++)
	       {
	      if (!AcCh[i]) continue;

	      if (FreqCut[i] == 0)
		{
		  mobileAVGF5(samples[i],&amp[i][0],nmedia[i],&m_amp[i][0]);
		}
	      else
		{
		  if (i==0)
		    for(int l=0;l<samples[i];l++)
		      freq[l]=(1.e9/DT)*(1./samples[i])*l;

		  FFTrealtocomplex(samples[i],&amp[i][0],&FFT_real[i][0],&FFT_comp[i][0],&FFT_abs[i][0]);
		  for(int l=0;l<samples[i];l++)
		    {
		      if (freq[l]< FreqCut[i])
			{
			  FFT_real_cut[i][l] = FFT_real[i][l];
			  FFT_comp_cut[i][l] = FFT_comp[i][l];
			}
		      else
			{
			  FFT_real_cut[i][l] = 0;
			  FFT_comp_cut[i][l] = 0;
			}
		    }
		  // Single inverse FFT, after the full spectrum has been cut
		  FFTcomplextoreal(samples[i],&FFT_real_cut[i][0],&FFT_comp_cut[i][0],&m_amp[i][0]);
		}

	      // measured in the last 10 ns
	      BackgroundGPT2(&m_amp[i][0], samples[i]-10./DT,samples[i],&bck[i],&max_bck_after[i],&rms_bck_after[i]);

	      // measured in the first 5 ns
	      BackgroundGPT2(&m_amp[i][0], 0, 5./DT-np_offset,&bck[i],&max_bck_before[i],&rms_bck_before[i]);

	      nprec = 0;
	      for(j=0;j<samples[i];j++)
		{
		  if(j%CAMPFACT==0)
		    {
		      timerec[nprec]=timeS[j];
		      amprec[i][nprec]=amp[i][j];
		      m_amprec[i][nprec]=(m_amp[i][j]-bck[i])*polarity[i];
		      nprec++;
		    }
		}

	      samplesrec[i]=nprec;

	     //negative amplitudes: nAmplitudes() compares against timeS[] so it
	     //takes the window in nanoseconds, not in samples
	      nAmplitudes(samples[i],&m_amprec[i][0],timeS,rms_bck_before[i],t_I[i],t_II[i],&nampl[i],&t_max[i],&ampl_chi2[i]);

	     t_max[i] = -1000; // t_max is for the positive amplitudes

	     AmplitudesF5(samples[i],&m_amprec[i][0],timeS,rms_bck_before[i],t_I[i]/DT,t_II[i]/DT,&ampl[i],&t_max[i],&ampl_chi2[i]);

	    }

	  for(i=0;i<nchro;i++)
	    {
	      if (!AcCh[i]) continue;

	      NArr = 0;
	      for (int kl = 0; kl<NArrayL; kl++)
		{
		  NTValuesL[kl] = 0;
		  TValuesL[kl] = 0;
		  if ( ThresholdL[kl]<ampl[i] ) NArr++;
		}

	      TrailNArr = 0;
	      for (int kl = 0; kl<TrailNArrayL; kl++)
		{
		  TrailTValuesL[kl] = 0;
		  TrailNTValuesL[kl] = 0;
		  if ( TrailThresholdL[kl]<ampl[i] )  TrailNArr++;
		}

	      // CFD arrays must be cleared as well: the CFD functions run only for
	      // ampl>5 and would otherwise leak the previous channel's values
	      for (int kl = 0; kl<NArrayF; kl++)
		{
		  TValuesF[kl] = 0;
		  NTValuesF[kl] = 0;
		}

	      for (int kl = 0; kl<TrailNArrayF; kl++)
		{
		  TrailTValuesF[kl] = 0;
		  TrailNTValuesF[kl] = 0;
		}

	       if (NArr>0)      Levtime_N_arrayGPT2( &NTValuesL[0], &TValuesL[0],  NArr, &ThresholdL[0], &m_amprec[i][0],DT,NMax[i],Fakebck[i],Interpolation);
	       if (TrailNArr>0) TrailLevtime_N_arrayGPT2(&TrailTValuesL[0],  TrailNArr, &TrailThresholdL[0], samples[i], &m_amprec[i][0],DT,NMax[i],Fakebck[i],Interpolation);

	       if (ampl[i]>5) ConstFractime_N_arrayGPT2( &NTValuesF[0], &TValuesF[0],  NArrayF, &ThresholdF[0],  &m_amprec[i][0],DT,Fakebck[i],ampl[i],NMax[i],Interpolation);
	       if (ampl[i]>5) TrailConstFractime_N_arrayGPT2(&TrailNTValuesF[0], &TrailTValuesF[0],  TrailNArrayF, &TrailThresholdF[0], samples[i],  &m_amprec[i][0],DT,Fakebck[i],ampl[i],NMax[i],Interpolation);

	       // Fill Level Values
	       t_level10[i] = TValuesL[0]+DTCh0i[i];
	       t_level15[i] = TValuesL[1]+DTCh0i[i];
	       t_level20[i] = TValuesL[2]+DTCh0i[i];
	       t_level30[i] = TValuesL[3]+DTCh0i[i];
	       t_level40[i] = TValuesL[4]+DTCh0i[i];
	       t_level50[i] = TValuesL[5]+DTCh0i[i];
	       t_level60[i] = TValuesL[6]+DTCh0i[i];
	       t_level80[i] = TValuesL[7]+DTCh0i[i];
	       t_level100[i] = TValuesL[8]+DTCh0i[i];
	       t_level200[i] = TValuesL[9]+DTCh0i[i];
	       t_level300[i] = TValuesL[10]+DTCh0i[i];

	       trail_t_level10[i] = TrailTValuesL[0]+DTCh0i[i];
	       trail_t_level20[i] = TrailTValuesL[1]+DTCh0i[i];
	       trail_t_level30[i] = TrailTValuesL[2]+DTCh0i[i];
	       trail_t_level40[i] = TrailTValuesL[3]+DTCh0i[i];
	       trail_t_level50[i] = TrailTValuesL[4]+DTCh0i[i];
	       trail_t_level60[i] = TrailTValuesL[5]+DTCh0i[i];
	       trail_t_level80[i] = TrailTValuesL[6]+DTCh0i[i];
	       trail_t_level100[i] = TrailTValuesL[7]+DTCh0i[i];
	       trail_t_level120[i] = TrailTValuesL[8]+DTCh0i[i];
	       trail_t_level140[i] = TrailTValuesL[9]+DTCh0i[i];
	       trail_t_level160[i] = TrailTValuesL[10]+DTCh0i[i];
	       trail_t_level180[i] = TrailTValuesL[11]+DTCh0i[i];

	       // Fill Constant Fraction Values
	       ToT20[i] = trail_t_level20[i] -  t_level20[i];
	       ToT30[i] = trail_t_level30[i] -  t_level30[i];
	       ToT40[i] = trail_t_level40[i] -  t_level40[i];

	       Ncfd05[i] = NTValuesF[0];
	       Ncfd10[i] = NTValuesF[1];
	       Ncfd15[i] = NTValuesF[2];
	       Ncfd20[i] = NTValuesF[3];
	       Ncfd25[i] = NTValuesF[4];
	       Ncfd30[i] = NTValuesF[5];
	       Ncfd40[i] = NTValuesF[6];
	       Ncfd50[i] = NTValuesF[7];
	       Ncfd60[i] = NTValuesF[8];
	       Ncfd70[i] = NTValuesF[9];
	       Ncfd80[i] = NTValuesF[10];
	       Ncfd90[i] = NTValuesF[11];

	       cfd05[i] = TValuesF[0]+DTCh0i[i];
	       cfd10[i] = TValuesF[1]+DTCh0i[i];
	       cfd15[i] = TValuesF[2]+DTCh0i[i];
	       cfd20[i] = TValuesF[3]+DTCh0i[i];
	       cfd25[i] = TValuesF[4]+DTCh0i[i];
	       cfd30[i] = TValuesF[5]+DTCh0i[i];
	       cfd40[i] = TValuesF[6]+DTCh0i[i];
	       cfd50[i] = TValuesF[7]+DTCh0i[i];
	       cfd60[i] = TValuesF[8]+DTCh0i[i];
	       cfd70[i] = TValuesF[9]+DTCh0i[i];
	       cfd80[i] = TValuesF[10]+DTCh0i[i];
	       cfd90[i] = TValuesF[11]+DTCh0i[i];

	       trail_cfd10[i] = TrailTValuesF[0]+DTCh0i[i];
	       trail_cfd30[i] = TrailTValuesF[1]+DTCh0i[i];
	       trail_cfd50[i] = TrailTValuesF[2]+DTCh0i[i];
	       trail_cfd70[i] = TrailTValuesF[3]+DTCh0i[i];
	       trail_cfd90[i] = TrailTValuesF[4]+DTCh0i[i];

	       trail_Ncfd10[i] = TrailNTValuesF[0];
	       trail_Ncfd30[i] = TrailNTValuesF[1];
	       trail_Ncfd50[i] = TrailNTValuesF[2];
	       trail_Ncfd70[i] = TrailNTValuesF[3];
	       trail_Ncfd90[i] = TrailNTValuesF[4];

	       // Compute area between bck-pulse
	         Charges(&m_amprec[i][0],DT,t_I[i]/DT,t_II[i]/DT,&area[i],&area_pos[i],rms_bck_before[i]);

	       if (( (cfd10[i]/DT) <samples[i] &&   (cfd20[i]/DT) < samples[i] && (cfd30[i]/DT)<samples[i]
		     &&  (cfd70[i]/DT) <samples[i]&&   (cfd80[i]/DT) <samples[i]) &&
		   ( (cfd10[i]/DT) >0 &&   (cfd20[i]/DT) >0 && (cfd30[i]/DT)>0
		     &&  (cfd70[i]/DT) >0 &&   (cfd80[i]/DT) >0))
		 {

		   if ((cfd30[i]- cfd10[i]) != 0)
		     dVdt1030[i] =( m_amprec[i][ int (cfd30[i]/DT)]- m_amprec[i][int (cfd10[i]/DT)])/(cfd30[i]- cfd10[i]);

		   if  ((cfd70[i]- cfd30[i]) != 0)
		     dVdt3070[i] =( m_amprec[i][ int (cfd70[i]/DT)]- m_amprec[i][int (cfd30[i]/DT)])/(cfd70[i]- cfd30[i]);

		   if  ((cfd80[i]- cfd20[i]) != 0)
		     dVdt2080[i] =( m_amprec[i][ int (cfd80[i]/DT)]- m_amprec[i][int (cfd20[i]/DT)])/(cfd80[i]- cfd20[i]);
		 }

	       if(doFit)
		 {
		   riselinfit(samples[i],&m_amprec[i][0],timeS, Fakebck[i], t_level20[i],t_level40[i], &rise_lin0[i], &rise_lin1[i],&rise_lin_chi2[i]) ;
		   if (rise_lin1[i] != 0) t_zero[i]= (-rise_lin0[i])/rise_lin1[i];
		   riseexpfit(samples[i],&m_amprec[i][0],timeS, Fakebck[i], cfd20[i]-5,cfd20[i], &rise_exp0[i], &rise_exp1[i],&rise_exp_chi2[i]) ;
		 }

	    } // end loop on channels


	  OutTree->Fill();

	NoRoot:
	  continue;


	}  // while{running}

 Write:

      f->Close();
      delete f;
      continue;

    }

   OutputFile->cd();
   OutTree->Write();
   OutputFile->Close();

   return(EXIT_SUCCESS);
}
