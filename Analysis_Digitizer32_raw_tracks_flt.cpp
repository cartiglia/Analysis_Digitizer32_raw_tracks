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

//#include <TThread.h>
//#include <ROOT/TTreeProcessorMT.hxx>


// To compile:
// clang++ -std=c++17 -O3 -Wall -Wextra    `root-config --cflags --libs`   Analysis_Digitizer32_raw_tracks_flt.cpp -o Analysis_Digitizer32_raw_tracks_flt

// Lorenzo includes
#include "Analysis_Digitizer16.h"

using namespace std;

int main()
{

  //bool showFFT=false;
  bool showFFT=false;
  bool doFit = false;
  bool FWaveform = true;

  TFile *OutputFile = new TFile("RunDigitizer32_raw_tracks_flt.root","recreate");

  ifstream InputCARD("Input_Folder_Digitizer32_raw_tracks.txt");

  if(!InputCARD)
    {
      cout << "Error: could not find the InputCARD" << endl;
      return (0);
    }

  // read the active number of channels
  int ntmp;
  int np_acq; // defined reading the file
  int  Nth = 0;
  int np_offset = 0; // 500; //100;
  int np_Max = 1000;

  nMCP = 7;

  float Totie, Totih, Totieg, Totihg;
  bool FToffeeCh[4] = {0,0,0,0};

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

  int AcCh[36] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
		  1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0};

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

	  OutTree->Branch("time",timerec,"time[samples_0]/F");

	}

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

 cout << "Branch settings done" << endl;


  //Loop on different runs
  while(1)
    {
      pip="";
      toffee="";
      if(InputCARD.eof())
	break;

      //Read run variables
      InputCARD >> pip >> nrun >> pip >> ntrig >> pip >> nchro >> pip >> MaxEvt >> pip >> nMCP;
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
      sprintf(title,"%s",InputNA.c_str());
      InputCARD >> pip >> nme;

      cout << endl << " File name:" << endl << title << " averaging " << nme << " points " << endl;

      FreqCut[0] = 0;
      FreqCut[1] = 0;
      FreqCut[2] = 0;
      FreqCut[3] = 0;

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

      Long64_t nentries = tree->GetEntries();

      cout << "Number of events = " << nentries << endl;
      if (MaxEvt == -1) MaxEvt = nentries;
      string mystr;

      string RunBit;
      string RunRoot;
      RunRoot = String_nrun[0];


      emptycount=0;

      float pippo;


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

	  reader.SetEntry(ntrig);

	  for(i=0;i<nchro;i++)
	    {
	      nmedia[i] = nme;
	      samples[i]=0;
	      AMax[i]=0;
	      NMax[i]=0;
	      RunBit = String_nrun[i+1];
	      t_max[i] = -1000;
	      max_bck_after[i]= -1000;
	      rms_bck_after[i] = -1000;
	      max_bck_before[i]= -1000;
	      rms_bck_before[i] = -1000;
	      ampl[i] = -1000;
	      nampl[i] = 1000;
	      ampl_chi2[i] = -1000;
	      samples[i]=0;
	      bck[i] = -1000;
	      Fakebck[i] = 0;
	      polarity[i] = 0;

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
	       DTCh0i[i] = 0;
	      if (ntmp>0) sprintf(number,"%d", ntrig);
	      if (ntmp<0) sprintf(number,"%d", ntrig);

	      tempsum=0;
	      np=0;
	      nprec=0;
	      Totie = 0;
	      Totih = 0;
	      Totieg = 0;
	      Totihg = 0;
	      WF2gain = 0;
	      treesize = (int)(*wfmSize);

	      if (i ==0 && ntrig ==0)
		{
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

	      for (int npp =np_offset;npp<np_Max;npp++)
		{
		  timeS[np] =  0.2*(npp);

		  if      (i ==  0) amp[i][np] = w0[npp]  *ADCmV;
		  else if (i ==  1) amp[i][np] = w1[npp]  *ADCmV;
		  else if (i ==  2) amp[i][np] = w2[npp]  *ADCmV;
		  else if (i ==  3) amp[i][np] = w3[npp]  *ADCmV;
		  else if (i ==  4) amp[i][np] = w4[npp]  *ADCmV;
		  else if (i ==  5) amp[i][np] = w5[npp]  *ADCmV;
		  else if (i ==  6) amp[i][np] = w6[npp]  *ADCmV;
		  else if (i ==  7) amp[i][np] = w7[npp]  *ADCmV;
		  else if (i ==  8) amp[i][np] = w8[npp]  *ADCmV;
		  else if (i ==  9) amp[i][np] = w9[npp]  *ADCmV;
		  else if (i == 10) amp[i][np] = w10[npp] *ADCmV;
		  else if (i == 11) amp[i][np] = w11[npp] *ADCmV;
		  else if (i == 12) amp[i][np] = w12[npp] *ADCmV;
		  else if (i == 13) amp[i][np] = w13[npp] *ADCmV;
		  else if (i == 14) amp[i][np] = w14[npp] *ADCmV;
		  else if (i == 15) amp[i][np] = w15[npp] *ADCmV;
		  else if (i == 16) amp[i][np] = w16[npp] *ADCmV;
		  else if (i == 17) amp[i][np] = w17[npp] *ADCmV;
		  else if (i == 18) amp[i][np] = w18[npp] *ADCmV;
		  else if (i == 19) amp[i][np] = w19[npp] *ADCmV;
		  else if (i == 20) amp[i][np] = w20[npp] *ADCmV;
		  else if (i == 21) amp[i][np] = w21[npp] *ADCmV;
		  else if (i == 22) amp[i][np] = w22[npp] *ADCmV;
		  else if (i == 23) amp[i][np] = w23[npp] *ADCmV;
		  else if (i == 24) amp[i][np] = w24[npp] *ADCmV;
		  else if (i == 25) amp[i][np] = w25[npp] *ADCmV;
		  else if (i == 26) amp[i][np] = w26[npp] *ADCmV;
		  else if (i == 27) amp[i][np] = w27[npp] *ADCmV;
		  else if (i == 28) amp[i][np] = w28[npp] *ADCmV;
		  else if (i == 29) amp[i][np] = w29[npp] *ADCmV;
		  else if (i == 30) amp[i][np] = w30[npp] *ADCmV;
		  else if (i == 31) amp[i][np] = w31[npp] *ADCmV;
		  else if (i == 32) amp[i][np] = trg0[npp]*ADCmV;
		  else if (i == 33) amp[i][np] = trg1[npp]*ADCmV;
		  else if (i == 34) amp[i][np] = trg2[npp]*ADCmV;
		  else if (i == 35) amp[i][np] = trg3[npp]*ADCmV;

		  np++;
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
	      polarity[i] = (amp[i][NMax[i]]-amp[i][0])/AMax[i];

	      pol[i] =  polarity[i];

	      if (ntrig == 2  && i == 0) cout << "Number of points: " << np  << " with a time step of " << DT << " [ns] " << endl;

	      if (ntrig <20  ) cout << "Event: " << event << " channel " << i << " has maximum value of " << AMax[i] << " [mV] at  "
				    <<  TMax[i] << " [ns] " << " on sample " << NMax[i] << " polarity = " << polarity[i] << endl;

	      if ( TMax[i] ==0 ||  AMax[i]>10000)
		{
		  cout << "Maximum at 0 ns, or signal too large ==> something must be wrong. The program skips this event" << endl;
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
	      t_I[i] =  t_max_sample[ChMaxSignal]-2; 
	      //	        t_II[i] =  t_max_sample[nMCP]-1;
	      t_II[i] = t_max_sample[ChMaxSignal]+3; //  80;


	      t_I[i] = 50;
	      //	        t_II[i] =  t_max_sample[nMCP]-1;
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

	      if (FreqCut[i] == 0)
		{
		  mobileAVG(samples[i],&amp[i][0],nmedia[i],&m_amp[i][0]);
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
		      FFTcomplextoreal(samples[i],&FFT_real_cut[i][0],&FFT_comp_cut[i][0],&m_amp[i][0]);
		    }
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

	     //negative amplitudes
	      nAmplitudes(samples[i],&m_amprec[i][0],timeS,rms_bck_before[i],t_I[i]/DT,t_II[i]/DT,&nampl[i],&t_max[i],&ampl_chi2[i]);

	     t_max[i] = -1000; // t_max is for the positive amplitudes

	     Amplitudes(samples[i],&m_amprec[i][0],timeS,rms_bck_before[i],t_I[i]/DT,t_II[i]/DT,&ampl[i],&t_max[i],&ampl_chi2[i]);

	     //	     cout << t_I[i]/DT << " " << t_II[i]/DT  << endl;
	    }

	  //WARNING, to be changed for new configuration
	  jump = 0;
	  for(i=1000;i<35;i++)
	    {
	      if (AcCh[i] && ampl[i] > 10)
		{
		  jump++;
		  break;
		}

	    }

	  //	  if (!jump) goto NoRoot;

	  for(i=0;i<nchro;i++)
	    {

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
		  if ( TrailThresholdL[kl]<ampl[i] )  TrailNArr++;
		}

	      for (int kl = 0;  kl<NArrayL; kl++)
		{
		  TValuesL[kl]=0;
		  NTValuesL[kl]=0;
		}

	      for (int kl = 0; kl<TrailNArrayL; kl++)
		{
		  TrailTValuesL[kl]=0;
		  TrailNTValuesL[kl]=0;
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


	  if(running==0 && emptycount<100)
	    {
	      running=1;
	      continue;
	    }
	  if(running==0 && emptycount>=100)
	    {
	      break;
	    }


	  OutTree->Fill();

	NoRoot:
	  continue;


	}  // while{running}

 Write:

      tree->Delete();
      continue;

    }

   OutputFile->cd();
   OutTree->Write();
   OutputFile->Close();

   return(EXIT_SUCCESS);
}
