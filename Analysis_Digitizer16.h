#ifndef WAVE_ANALYZER_H
#define WAVE_ANALYZER_H

#define NCHRO 36
#define MAX_CAMP 20000
#define CAMPFACT 1
#define MAX_CAMPREC MAX_CAMP/CAMPFACT+100

char title[200];
char title1[200];

Int_t ntrk;
Float_t trkx1, trky1, trkx2,trky2;

Int_t osc[NCHRO];
Int_t chan[NCHRO];
Int_t chmax,TRKev;
Float_t Amax;
Float_t bck[NCHRO];
Float_t Fakebck[NCHRO];
Float_t  TRKx[120000], TRKy[120000], Nothing;
Float_t xtrk1;
Float_t ytrk1;
Float_t xtrk2;
Float_t ytrk2;
Float_t trkchi2;
Float_t x_pos;
Float_t y_pos;
Float_t ADCmV = 0.25;

int nMCP = 0;

TString  TrkFile;

int jump = 0;
int Flag =0;

Float_t t_I[NCHRO];
Float_t t_II[NCHRO];
Int_t camp[NCHRO];
Int_t pol[NCHRO];
Int_t camprec[NCHRO];
Int_t samples[NCHRO];
Int_t samplesrec[NCHRO];
Float_t FreqCut[NCHRO];
Float_t t_max_sample[NCHRO];
Float_t ampl_sample[NCHRO];
Float_t DTCh0i[NCHRO];

Float_t ToT20[NCHRO];
Float_t ToT30[NCHRO];
Float_t ToT40[NCHRO];

Float_t DT;

Float_t timeS[MAX_CAMP]; 
Float_t Itot[NCHRO][MAX_CAMP];
Float_t amp[NCHRO][MAX_CAMP];
Float_t m_amp[NCHRO][MAX_CAMP];
Float_t der_amp[NCHRO][MAX_CAMP];
Float_t amprec[NCHRO][MAX_CAMPREC];
Float_t m_amprec[NCHRO][MAX_CAMPREC];
Float_t d_amprec[NCHRO][MAX_CAMP];
Float_t FFT_real[NCHRO][MAX_CAMP];
Float_t FFT_comp[NCHRO][MAX_CAMP];
Float_t FFT_abs[NCHRO][MAX_CAMP];

Float_t FFT_real_cut[NCHRO][MAX_CAMP];
Float_t FFT_comp_cut[NCHRO][MAX_CAMP];
Float_t FFT_abs_cut[NCHRO][MAX_CAMP];

Float_t freq[MAX_CAMPREC];
Float_t timerec[MAX_CAMPREC];

Float_t XCross;
Float_t YCross;

Float_t AMaxSignal;
Float_t TMaxSignal;

Int_t ChMaxSignal;

  //Analysis variables

Float_t day;
Float_t month;
Float_t year;
Float_t hour;
Float_t minute;

Float_t MaxDeriv;
Float_t sec[NCHRO];
Float_t ampl[NCHRO];
Float_t nampl[NCHRO];
Float_t area[NCHRO];
Float_t area_pos[NCHRO];
Float_t dVdt1030[NCHRO];
Float_t dVdt3070[NCHRO];
Float_t dVdt2080[NCHRO];

  Float_t totcha[NCHRO];
  Float_t maxS[NCHRO];
  Float_t maxD[NCHRO];
  Float_t max_bck_before[NCHRO];
  Float_t rms_bck_before[NCHRO];
  Float_t max_bck_after[NCHRO];
  Float_t rms_bck_after[NCHRO];
  Float_t charge[NCHRO];
  Int_t strip[NCHRO];
  Float_t FEcharge[NCHRO];
  Float_t FEchargeerr[NCHRO];
  Float_t totcharge[NCHRO];
  Float_t p_charge;
  Float_t p_chargeerr;
  Int_t sat[NCHRO];
  Float_t maxamp;
  Int_t maxstrip;
  Float_t eff_dis08[NCHRO];
  Float_t eff_dis10[NCHRO];
  Float_t eff_dis15[NCHRO];
  Float_t eff_dis20[NCHRO];
  Float_t eff_dis30[NCHRO];
  Float_t eff_dis40[NCHRO];
  Float_t eff_dis50[NCHRO];
  Float_t eff_dis60[NCHRO];
  Float_t eff_dis70[NCHRO];
Float_t eff_dis100[NCHRO];
  Float_t eff_dis150[NCHRO];
  Float_t eff_dis200[NCHRO];
  Float_t eff_dis250[NCHRO];
  Float_t eff_dis300[NCHRO];
  Float_t eff_dis400[NCHRO];
  Float_t eff_dis500[NCHRO];
  Float_t eff_dis800[NCHRO];
  Float_t eff_dis1000[NCHRO];
  Float_t eff_dis2000[NCHRO];
  Float_t eff_dis5000[NCHRO];
  Float_t eff_dis10000[NCHRO];
  Float_t t_level10[NCHRO];
  Float_t t_level15[NCHRO];
  Float_t t_level18[NCHRO];
  Float_t t_level20[NCHRO];
  Float_t t_level30[NCHRO];
  Float_t t_level40[NCHRO];
  Float_t t_level50[NCHRO];
  Float_t t_level60[NCHRO];
  Float_t t_level80[NCHRO];
  Float_t t_level100[NCHRO];
  Float_t t_level200[NCHRO];
  Float_t t_level300[NCHRO];

 
  Float_t trail_t_level10[NCHRO];
  Float_t trail_t_level15[NCHRO];
  Float_t trail_t_level18[NCHRO];
  Float_t trail_t_level20[NCHRO];
  Float_t trail_t_level30[NCHRO];
  Float_t trail_t_level40[NCHRO];
  Float_t trail_t_level50[NCHRO];
  Float_t trail_t_level60[NCHRO];
  Float_t trail_t_level80[NCHRO];
  Float_t trail_t_level100[NCHRO];
  Float_t trail_t_level120[NCHRO];
  Float_t trail_t_level140[NCHRO];
 Float_t trail_t_level160[NCHRO];
  Float_t trail_t_level180[NCHRO];
  Float_t trail_t_level200[NCHRO];
  Float_t trail_t_level300[NCHRO];

  Float_t cfd05[NCHRO];
  Float_t cfd10[NCHRO];
  Float_t cfd15[NCHRO];
  Float_t cfd20[NCHRO];
  Float_t cfd25[NCHRO];
  Float_t cfd30[NCHRO];
  Float_t cfd40[NCHRO];
  Float_t cfd50[NCHRO];
  Float_t cfd60[NCHRO];
  Float_t cfd70[NCHRO];
  Float_t cfd80[NCHRO];
  Float_t cfd90[NCHRO];

  Int_t Ncfd05[NCHRO];
  Int_t Ncfd10[NCHRO];
  Int_t Ncfd15[NCHRO];
  Int_t Ncfd20[NCHRO];
  Int_t Ncfd25[NCHRO];
  Int_t Ncfd30[NCHRO];
  Int_t Ncfd40[NCHRO];
  Int_t Ncfd50[NCHRO];
  Int_t Ncfd60[NCHRO];
  Int_t Ncfd70[NCHRO];
  Int_t Ncfd80[NCHRO];
  Int_t Ncfd90[NCHRO];

  Int_t trail_Ncfd10[NCHRO];
  Int_t trail_Ncfd30[NCHRO];
  Int_t trail_Ncfd50[NCHRO];
  Int_t trail_Ncfd70[NCHRO];
  Int_t trail_Ncfd90[NCHRO];

  Float_t dist[NCHRO];
  Int_t ind[NCHRO];
  Int_t XPad[NCHRO];
  Int_t YPad[NCHRO];
  Float_t XDist[NCHRO];
  Float_t YDist[NCHRO];

  Float_t trail_cfd10[NCHRO];
  Float_t trail_cfd30[NCHRO];
  Float_t trail_cfd50[NCHRO];
  Float_t trail_cfd70[NCHRO];
  Float_t trail_cfd90[NCHRO];
 


  Int_t Nt_level30[NCHRO];
  Int_t Ntrail_t_level30[NCHRO];
  Float_t t_cent30[NCHRO];
  Float_t t_centMax[NCHRO];
  Float_t t_cent[NCHRO];

  Float_t t_rms3[NCHRO];
  Float_t t_rms5[NCHRO];
  Float_t t_max[NCHRO];
Float_t t_zero[NCHRO];
Float_t fitchi[NCHRO];
Float_t ampl_chi2[NCHRO];
Float_t rise_lin0[NCHRO];
Float_t rise_lin1[NCHRO];
Float_t rise_lin_chi2[NCHRO];
  Float_t rise_exp0[NCHRO];
  Float_t rise_exp1[NCHRO];
Float_t rise_exp_chi2[NCHRO];
Float_t AMax[NCHRO];
Float_t TMax[NCHRO];
int NMax[NCHRO];
int polarity[NCHRO];


  // Analysis branches
Float_t WF2gain;
Int_t Wafer;
Int_t LGain;
Int_t FType;
Int_t DMetal;
Int_t DPitch;

Float_t WF2xpos;
Float_t WF2angle;
Float_t TotCharge;
Int_t MaxEvt;
Float_t XCord;
Float_t YCord;
Float_t XAmpl;
Float_t YAmpl;
Float_t XTime;
Float_t YTime;
Float_t XCirc;
Float_t YCirc;
Float_t XArea;
Float_t YArea;


Int_t CS;
Int_t eff;
Float_t noise=3.;
Float_t bar0;
Float_t bar1;
Float_t A0;
Float_t A1;
Float_t delta0;
Float_t delta1;
Float_t chinorm0;
Float_t chinorm1;


  // Float_t gum_amp[MAX_CAMP];


Int_t nchro; // = NCHRO;  //Maximum number of Read out channels to be present in the rootple
//Int_t nchActive = NCHRO;  //Number of active Read out channels, to be filled written in the InputFolder
Int_t nchro_active;
bool doFFT = false;  
bool showFFT = false;

std::string pip;
std::string toffee;
std::string InputNAME[1000];
std::string InputNA;
int i;
int j;
int k;
int np;
int nprec;


Float_t hv;
Float_t hvcorr;
Float_t temp;
Float_t pres;
Float_t TTCfd30;
Float_t TTCfd50;
Float_t TCfd20;
Float_t TCfd30;
Float_t TCfd40;
Float_t TCfd50;
Float_t TTmax;

Int_t nrun=0;
std::string String_nrun;
Int_t ntrig=0;
Int_t ds=0;
Int_t event=0;
Int_t CAmax=0;


TFile *OutputFile;
TTree *OutTree;

std::ostringstream leaf(std::ostringstream::out);
std::ostringstream leafl(std::ostringstream::out);

  // Program variables
std::ostringstream datafile;

  float Tpeak = 0;
  float BallDef = 0;

int nmedia[NCHRO];
int nme;

  float tempsum;
  int running=1;
  int emptycount=0;
  float time0temp;
  int segsize;
  char number[10];
  char date[12];
  char timeofday[8];
  float FFT_ratio;



   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t        i_timestamp;
   Float_t        i_current;
   std::vector<float>  *pos;
   std::vector<double>  *t0;
   std::vector<double>  *t1;
   std::vector<double>  *t2;
   std::vector<double>  *t3;
   std::vector<double>  *t4;
   std::vector<double>  *t5;
   std::vector<double>  *t6;
   std::vector<double>  *t7;
   std::vector<double>  *t8;



TBranch        *b_i_timestamp;   //!
TBranch        *b_i_current;   //!


TBranch        *b_pos;   //!
TBranch        *b_t0;   //!
TBranch        *b_t1;   //!
TBranch        *b_t2;   //!
TBranch        *b_t3;   //!
TBranch        *b_t4;   //!
TBranch        *b_t5;   //!
TBranch        *b_t6;   //!
TBranch        *b_t7;   //!
TBranch        *b_t8;   //!


 // Declaration of leaf types for digitizer
   Float_t        bias;
   Float_t        Dig_freq;
   Float_t        treesize;

   std::vector<double>  *w0;
   std::vector<double>  *w1;
   std::vector<double>  *w2;
   std::vector<double>  *w3;
   std::vector<double>  *w4;
   std::vector<double>  *w5;
   std::vector<double>  *w6;
   std::vector<double>  *w7;
   std::vector<double>  *w8;
   std::vector<double>  *w9;
   std::vector<double>  *w10;
   std::vector<double>  *w11;
   std::vector<double>  *w12;
   std::vector<double>  *w13;
   std::vector<double>  *w14;
   std::vector<double>  *w15;
   std::vector<double>  *trg0;
   std::vector<double>  *trg1;







Float_t   fxtrk1;
Float_t   fytrk1;
Float_t   fxtrk2;
Float_t   fytrk2;
Float_t  fchi2trk;

   // List of branches
   TBranch        *b_bias;   //!
   TBranch        *b_freq;   //!
   TBranch        *b_size;   //!
   TBranch        *b_w0;   //!
   TBranch        *b_w1;   //!
   TBranch        *b_w2;   //!
   TBranch        *b_w3;   //!
   TBranch        *b_w4;   //!
   TBranch        *b_w5;   //!
   TBranch        *b_w6;   //!
   TBranch        *b_w7;   //!
   TBranch        *b_w8;   //!
   TBranch        *b_w9;   //!
   TBranch        *b_w10;   //!
   TBranch        *b_w11;   //!
   TBranch        *b_w12;   //!
   TBranch        *b_w13;   //!
   TBranch        *b_w14;   //!
   TBranch        *b_w15;   //!
   TBranch        *b_trg0;   //!
   TBranch        *b_trg1;   //!
   TBranch        *b_fxtrk1;   //!
   TBranch        *b_fytrk1;   //!
   TBranch        *b_fxtrk2;   //!
   TBranch        *b_fytrk2;   //!
   TBranch        *b_fchi2trk;   //!




   std::vector<double>  *w0_1;
   std::vector<double>  *w1_1;
   std::vector<double>  *w2_1;
   std::vector<double>  *w3_1;
   std::vector<double>  *w4_1;
   std::vector<double>  *w5_1;
   std::vector<double>  *w6_1;
   std::vector<double>  *w7_1;
   std::vector<double>  *w8_1;
   std::vector<double>  *w9_1;
   std::vector<double>  *w10_1;
   std::vector<double>  *w11_1;
   std::vector<double>  *w12_1;
   std::vector<double>  *w13_1;
   std::vector<double>  *w14_1;
   std::vector<double>  *w15_1;
   std::vector<double>  *trg0_1;
   std::vector<double>  *trg1_1;

   TBranch        *b_w0_1;   //!
   TBranch        *b_w1_1;   //!
   TBranch        *b_w2_1;   //!
   TBranch        *b_w3_1;   //!
   TBranch        *b_w4_1;   //!
   TBranch        *b_w5_1;   //!
   TBranch        *b_w6_1;   //!
   TBranch        *b_w7_1;   //!
   TBranch        *b_w8_1;   //!
   TBranch        *b_w9_1;   //!
   TBranch        *b_w10_1;   //!
   TBranch        *b_w11_1;   //!
   TBranch        *b_w12_1;   //!
   TBranch        *b_w13_1;   //!
   TBranch        *b_w14_1;   //!
   TBranch        *b_w15_1;   //!
   TBranch        *b_trg0_1;   //!
   TBranch        *b_trg1_1;   //!



//Convert little endian int to big endian int
unsigned int intToBigEndianess(const unsigned int x)
{
  return  ( x >> 24 ) |  // Move first byte to the end,
    ( ( x << 8 ) & 0x00FF0000 ) | // move 2nd byte to 3rd,
    ( ( x >> 8 ) & 0x0000FF00 ) | // move 3rd byte to 2nd,
    ( x << 24 ); // move last byte to start.
}


//returns the float value corresponding to a given bit represention of a scalar int value or vector of int values
float intBitsToFloat(const int x)
{
  union 
  {
    float f;  // assuming 32-bit IEEE 754 single-precision
    int i;    // assuming 32-bit 2's complement int
  } 
  u;
  
  u.i = x;
  return u.f;
}


//Always run background function before doing further analyses

//Compute average amplitude for background
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, max background amplitude)
void Background(Float_t amp[],Int_t Nt0, Int_t Nt1, Float_t *bck, Float_t *maxbck, Float_t *rmsbck)
{

  //  double min=1000000;
  double max=-1000000;
  float sumamp=0;
  float sumampsq=0;
  int pointsavg=0;
  float meanbck;
  *maxbck = 0;

  // TH1D hlite("hlite","hlite",camp,min,max);
  *bck=0;
  if(Nt1>Nt0)
    {
      //
      // std::cout << Nt0 << " " << Nt1 << std::endl;
      for(int j=Nt0;j<Nt1;j++)
	{
	  //	  std::cout << ampl[j] << std::endl;
	  pointsavg++;
	  sumamp += amp[j];
	  sumampsq += amp[j]*amp[j];
	}
      meanbck=sumamp/float(pointsavg);
      *rmsbck=sqrt((1./float(pointsavg))*(sumampsq-float(pointsavg)*meanbck*meanbck));
      
      //*bck=hlite.GetBinCenter(hlite.GetMaximumBin());
      *bck=meanbck;
      
      for(int j=Nt0;j<Nt1;j++)
	{

	  if(std::fabs(amp[j]-meanbck)>max)
	    {
	      //	 std::cout << amp[j] << std::endl;
	      max=fabs(amp[j]-meanbck);
	    }
	      
	}
      
      *maxbck=max;
      //       std::cout << *bck << std::endl;
    }
  else
    {
      std::cout << "Cannot calculate bck: initial sample = " << Nt0 << " final sample =  " << Nt1 << std::endl;
    }
  return;
}



void BackgroundGPT2(Float_t amp[], Int_t Nt0, Int_t Nt1, Float_t *bck, Float_t *maxbck, Float_t *rmsbck) {
    // Initialize variables
    double maxDeviation = -std::numeric_limits<double>::infinity();
    double sumAmp = 0.0;
    double sumAmpSq = 0.0;
    int pointsCount = 0;

    // Check the validity of input range
    if (Nt1 <= Nt0) {
        std::cerr << "Cannot calculate background: initial sample = " << Nt0 << " final sample = " << Nt1 << std::endl;
        return;
    }

    // Accumulate sums for mean and RMS calculations
    for (int j = Nt0; j < Nt1; j++) {
        sumAmp += amp[j];
        sumAmpSq += amp[j] * amp[j];
        pointsCount++;
    }

    // Calculate mean background
    double meanBck = sumAmp / pointsCount;
    *bck = meanBck;

    // Calculate RMS background
    *rmsbck = std::sqrt((sumAmpSq / pointsCount) - meanBck * meanBck);

    // Find the maximum deviation from the background
    for (int j = Nt0; j < Nt1; j++) {
        double deviation = std::fabs(amp[j] - meanBck);
        if (deviation > maxDeviation) {
            maxDeviation = deviation;
        }
    }

    *maxbck = maxDeviation;
}

//Inversione degli impulsi negativi
void Inversion(Int_t camp, Float_t neg_amp[])
{
  for(int j=0;j<camp;j++)
	{
	  neg_amp[j]=-neg_amp[j];
        }
}


//DEVELOPEMENT
// Compute maximum amplitude and peak time for waveform using a gaussian fit to the peak
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse maximum amplitude (averaged), pulse maximum amplitude)
void Amplitudes(Int_t camp, Float_t amp[],  Float_t timeS[],   Float_t RMSbck, Int_t t_start,Int_t t_stop,  Float_t *pamp,  Float_t *tamp,  Float_t *chi)
{


  Int_t np=0;
  //  int cont=0;
  int j_max= 0;
  int j_start = 0;
  int n_j = 0;
  // Limit of gaus fitting around j_max
  float AmpFrac = 0.6;
  // double tempmean;
  // double tempsig;
  double_t   max=-100000;
  double_t   tmax=-100000;

  double tempminch=0;
  double bck=0;

  for(int j=0;j<t_start;j++)
    {
      bck += amp[j];
      tempminch++;    
    }
  
  bck *=1/tempminch;
 

  //Maximum calculation
  for(int j=t_start;j<t_stop;j++)
    {

      np++;
      if(amp[j]-bck>max)
	{
	  max=amp[j]-bck;
	  tmax=(timeS[j+1]*amp[j+1]+timeS[j]*amp[j]+timeS[j-1]*amp[j-1])/(amp[j+1]+amp[j]+amp[j-1]);
	  j_max=j;
	}
    }


  if (max<5*RMSbck)
    {
      //     std::cout << max  << " " << RMSbck << std::endl;
      *pamp= max;
      *tamp= tmax;
      *chi = -1;
    } 
  else
    {
      //Fit amplitude calculation
      // Find limits at AmpFrac
      
      //from the maximum backwards 
      for(int j=j_max;j>0;j--)
	{
	  if(amp[j]-bck<AmpFrac*(max))
	    {
	      j_start=j;
	      break;
	    }
	}
      
      //from the minimum limit, up the maximum, and then down till the limit on the falling edge  
      for(int j=j_start+1;j<camp-1;j++)
	{
	  n_j++;
	  if(amp[j]-bck<AmpFrac*(max))
	    {
	      break;
	    }
	}
      
      
      //  std::cout << j_start << " " <<j_max << " Number of points in the fit= " << n_j << std::endl; 
      
      double temp0;
      double temp1;
      double temp2;
      TF1 *f1 = new TF1("gausfit","[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))",timeS[j_start-1],timeS[j_start+n_j+1]);
      // TGraph g(n_j);
      TGraph *g = new TGraph();
      g->Set(n_j);
      for(int j=0;j<n_j;j++)
	{
	  g->SetPoint(j,timeS[j+j_start],amp[j+j_start]-bck);
	}
      
      temp0=max;
      temp1=tmax;
      temp2=(timeS[j_start]-timeS[j_start+n_j])/2.;
      f1->SetParameters(temp0,temp1,temp2);
      g->Fit("gausfit","QN","",timeS[j_start]-1,timeS[j_start+n_j-1]+1);
      //  g.Fit("gausfit","QN","",timeS[j_start]-1,timeS[j_start+n_j-1]+1);
      
      // std::cout <<  max - f1->GetParameter(0) << " chi2 " << f1->GetChisquare()/f1->GetNDF() <<  std::endl;
      *pamp=f1->GetParameter(0);
      *tamp=f1->GetParameter(1);
      *chi=f1->GetChisquare()/f1->GetNDF();
      
      if ( *chi > 50 || fabs( (*pamp-max)/max)>0.2)
	{

	  //	  std::cout << *tamp<< " " <<  *chi << " " << fabs( (*pamp-max)/max) <<std::endl;
	  *pamp= max;
	  *tamp= tmax;
	  *chi = -2;
	}
      

      g->Delete();
      f1->Delete();

    }

  return;
}


void AmplitudesGPT(Int_t camp, Float_t amp[], Float_t timeS[], Float_t bck, Float_t t_I, Float_t t_II, Float_t *pamp, Float_t *tamp, Float_t *chi) {
    Int_t np = 0;
    int j_max = 0;
    int j_start = 0;
    int n_j = 0;
    float AmpFrac = 0.6;
    double max = -100000;
    double tmax = -100000;

    // Maximum calculation
    for (int j = 0; j < camp - 1; j++) {
        if (timeS[j] > t_I && timeS[j] < t_II) {
            np++;
            if (amp[j] - bck > max) {
                max = amp[j] - bck;
                tmax = (timeS[j + 1] + timeS[j]) / 2.;
                j_max = j;
            }
        }
    }

    // Find limits at AmpFrac
    for (int j = j_max; j > 0; j--) {
        if (amp[j] - bck < AmpFrac * max) {
            j_start = j;
            break;
        }
    }

    for (int j = j_start + 1; j < camp - 1; j++) {
        n_j++;
        if (amp[j] - bck < AmpFrac * max) {
            break;
        }
    }

    TF1 *f1 = new TF1("gausfit", "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))", timeS[j_start - 1], timeS[j_start + n_j + 1]);
    TGraph *g = new TGraph();
    g->Set(n_j);

    for (int j = 0; j < n_j; j++) {
        g->SetPoint(j, timeS[j + j_start], amp[j + j_start] - bck);
    }

    double temp0 = max;
    double temp1 = tmax;
    double temp2 = (timeS[j_start] - timeS[j_start + n_j]) / 2.;
    f1->SetParameters(temp0, temp1, temp2);
    g->Fit("gausfit", "QN", "", timeS[j_start] - 1, timeS[j_start + n_j - 1] + 1);

    *pamp = f1->GetParameter(0);
    *tamp = f1->GetParameter(1);
    *chi = f1->GetChisquare() / f1->GetNDF();

    if (*tamp < t_I || *tamp > t_II || *chi > 50 || fabs((*pamp - max) / max) > 0.2) {
        *pamp = max;
        *tamp = tmax;
        *chi = -1;
    }

    // Use ROOT's memory management
    delete g;
    delete f1;
}



void AmplitudesGPT2(Int_t camp, Float_t amp[], Float_t timeS[], Float_t bck, Float_t t_I, Float_t t_II, Float_t *pamp, Float_t *tamp, Float_t *chi) {
    Int_t np = 0;
    int j_max = 0;
    int j_start = 0;
    int n_j = 0;
    float AmpFrac = 0.6;
    double max = -100000;
    double tmax = -100000;

    // Maximum calculation
    for (int j = 0; j < camp - 1; j++) {
        if (timeS[j] > t_I && timeS[j] < t_II) {
            np++;
            if (amp[j] - bck > max) {
                max = amp[j] - bck;
                tmax = (timeS[j + 1] + timeS[j]) / 2.;
                j_max = j;
            }
        }
    }

    // Find limits at AmpFrac
    for (int j = j_max; j > 0; j--) {
        if (amp[j] - bck < AmpFrac * max) {
            j_start = j;
            break;
        }
    }

    for (int j = j_start + 1; j < camp - 1; j++) {
        n_j++;
        if (amp[j] - bck < AmpFrac * max) {
            break;
        }
    }

    double temp0 = max;
    double temp1 = tmax;
    double temp2 = (timeS[j_start] - timeS[j_start + n_j]) / 2.;

    // Using unique_ptr for automatic memory management
    std::unique_ptr<TF1> f1(new TF1("gausfit", "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))", timeS[j_start - 1], timeS[j_start + n_j + 1]));
    std::unique_ptr<TGraph> g(new TGraph());
    g->Set(n_j);

    for (int j = 0; j < n_j; j++) {
        g->SetPoint(j, timeS[j + j_start], amp[j + j_start] - bck);
    }

    f1->SetParameters(temp0, temp1, temp2);
    g->Fit("gausfit", "QN", "", timeS[j_start] - 1, timeS[j_start + n_j - 1] + 1);

    *pamp = f1->GetParameter(0);
    *tamp = f1->GetParameter(1);
    *chi = f1->GetChisquare() / f1->GetNDF();

    if (*tamp < t_I || *tamp > t_II || *chi > 50 || fabs((*pamp - max) / max) > 0.2) {
        *pamp = max;
        *tamp = tmax;
        *chi = -1;
    }
    // No need for explicit delete, smart pointers handle it
}



void nAmplitudesGPT2(Int_t camp, Float_t amp[], Float_t timeS[], Float_t bck, Float_t t_I, Float_t t_II, Float_t *pamp, Float_t *tamp, Float_t *chi) {
    Int_t j_max = -1;
    double max = -1e5;
    double tmax = -1e5;
    float AmpFrac = 0.6;

    // Find the maximum within the specified time range
    for (int j = 0; j < camp - 1; j++) {
        if (timeS[j] > t_I && timeS[j] < t_II) {
            double adjusted_amp = -(amp[j] - bck);
            if (adjusted_amp > max) {
                max = adjusted_amp;
                tmax = (timeS[j + 1] + timeS[j]) / 2.0;
                j_max = j;
            }
        }
    }

    // If no maximum found, exit early
    if (j_max == -1) {
        *pamp = 0;
        *tamp = 0;
        *chi = -1;
        return;
    }

    // Find the start and end for fitting
    int j_start = 0;
    for (int j = j_max; j > 0; j--) {
        if (-(amp[j] - bck) < AmpFrac * max) {
            j_start = j;
            break;
        }
    }

    int n_j = 0;
    for (int j = j_start + 1; j < camp - 1; j++) {
        n_j++;
        if (-(amp[j] - bck) < AmpFrac * max) {
            break;
        }
    }

    if (n_j == 0) {
        *pamp = -max;
        *tamp = tmax;
        *chi = -1;
        return;
    }

    // Using smart pointers for TF1 and TGraph to avoid memory leaks
    auto f1 = std::make_unique<TF1>("gausfit", "[0]*exp(-(x-[1])*(x-[1])/(2*[2]*[2]))", timeS[j_start - 1], timeS[j_start + n_j + 1]);
    auto g = std::make_unique<TGraph>();
    g->Set(n_j);

    for (int j = 0; j < n_j; j++) {
        g->SetPoint(j, timeS[j + j_start], -(amp[j + j_start] - bck));
    }

    double temp0 = max;
    double temp1 = tmax;
    double temp2 = (timeS[j_start] - timeS[j_start + n_j]) / 2.0;
    f1->SetParameters(temp0, temp1, temp2);
    g->Fit("gausfit", "QN", "", timeS[j_start - 1], timeS[j_start + n_j - 1]);

    *pamp = -f1->GetParameter(0);
    *tamp = f1->GetParameter(1);
    *chi = f1->GetChisquare() / f1->GetNDF();

    // Check if the fit results are valid
    if (*tamp < t_I || *tamp > t_II || *chi > 50 || fabs((*pamp - max) / max) > 0.2) {
        *pamp = -max;
        *tamp = tmax;
        *chi = -1;
    }
}

void nAmplitudes(Int_t camp, Float_t amp[], Float_t timeS[],  Float_t RMSbck, Float_t t_I,Float_t t_II, Float_t *pamp, Float_t *tamp, Float_t *chi)
{
  
  Int_t np=0;
  //  int cont=0;
  int j_max= 0;
  int j_start = 0;
  int n_j = 0;
  // Limit of parabola fitting around j_max
  float AmpFrac = 0.6;
  // double tempmean;
  // double tempsig;
  double_t   max=0;
  double_t   tmax=-100000;

  //Maximum calculation
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_I && timeS[j]<t_II)
	{
	  np++;
	  if(amp[j]<max)
	    {
	      max= amp[j];
	      tmax=(timeS[j+1]*amp[j+1]+timeS[j]*amp[j]+timeS[j-1]*amp[j-1])/(amp[j+1]+amp[j]+amp[j-1]);
	      j_max=j;
	    }
	}
    }
    *pamp = max;
    *tamp = tmax;

  return;
}


//Amplitude without the Fit
void Amplitudes_NF( Float_t amp[], Float_t timeS[],  Float_t bck, Int_t NMax, Float_t *max,Float_t *tmax)
{
  

  *max=-100000;
  *tmax=-100000;
  *tmax=(timeS[NMax-1]*amp[NMax-1]+timeS[NMax]*amp[NMax]+timeS[NMax+1]*amp[NMax+1])/(amp[NMax-1]+amp[NMax]+amp[NMax+1]);
  *max=(amp[NMax-1]+amp[NMax]+amp[NMax+1]-3*bck)/3;
 
  return;
}



//DEVELOPEMENT
// Compute maximum, amplitude and peak time for waveform
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse maximum amplitude (averaged), pulse maximum amplitude)
void Neg_Amplitudes(Int_t camp, Float_t amp[], Float_t timeS[], Float_t t_I, Float_t t_II, Float_t bck, Float_t *pamp, Float_t *pmax, Float_t *tmax)
{
  *pamp = 0;
  Int_t np=0;
  // int cont=0;
  double maxtemp=-100000;
  //  double tempmean;
  //  double tempsig;
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_I && timeS[j]<t_II)
	{
	  np++;
	  if(-(amp[j]-bck)>maxtemp)
	    {
	      maxtemp=-(amp[j]-bck);
	      *tmax=(timeS[j+1]+timeS[j])/2.;
	    }
	}
    }
  *pmax=maxtemp;

  /*
  TGraph *g = new TGraph(np);
  TF1 *f= new TF1("gauss","[0]*exp(-((x-[1])*(x-[1]))/(2.*[2]*[2]))",t_I,t_II);
  f->SetParameters(*pmax,(t_II+t_I)/2.,(t_II-t_I)/2.);
  for(int j=0;j<camp-1;j++)
    {
      if(timeS[j]>t_I && timeS[j]<t_II)
	{
	  cont++;
	  g->SetPoint(cont,timeS[j],amp[j]-bck);
	}
    }
  g->Fit("gauss","QN","",t_I,t_II);
  
  tempmean=f->GetParameter(1);
  tempsig=f->GetParameter(2);
  g->Fit("gauss","QN","",tempmean-tempsig,tempmean+tempsig);
  
  *pamp=f->GetParameter(0);
  */
  return;
}

// Calibration using curve [0](1-e^(x/[1]))
// Set the function with the proper front end each time
// (Parameter zero, parameter 1, amplitude of amplified pulse, corresponding charge/amplitude of input pulse, error on input charge/amplitude, error on output amplitude - must be known)
void Calibration(Float_t cal_low_par0,Float_t cal_low_par1,Float_t cal_high_par0,Float_t cal_high_par1, Float_t max,Float_t *charge, Float_t *chargeerr, Float_t amperr)
{
  TF1 f_low("cal_low","[0]*log(1.-x/[1])",0.,cal_low_par1-0.1);
  TF1 f_high("cal_high","[0]*log(1.-pow((x/[1]),2))",0.,cal_high_par1-0.1);
  f_low.SetParameters(cal_low_par0,cal_low_par1);
  f_high.SetParameters(cal_high_par0,cal_high_par1);
  if(max<444)
    {
      *charge=f_low(max);
      *chargeerr=(-cal_low_par0/(cal_low_par1-(max)))*amperr;
    }
  
  if(max>=444 && max<cal_high_par1)
    {
      *charge=f_high(max);
      *chargeerr=(-cal_high_par0)*(1./(1.-pow(1.-max/cal_high_par1,2)))*2.*(1.-max/cal_high_par1)*amperr;
    }
  if(max>=cal_high_par1)
    {
      *charge=1200.;
      *chargeerr=1000.;
    }
  
  return;
}

// Calibration using curve [0](1-e^(x/[1]))
// Set the function with the proper front end each time
// (Parameter zero, parameter 1, amplitudesof amplified pulse, corresponding charge/amplitude of input pulse, error on input charge/amplitude, error on output amplitude - must be known)
void Calibration_ip(Float_t cal_A0,Float_t cal_x0,Float_t cal_B0,Float_t max,Float_t *charge, Float_t *chargeerr, Float_t amperr)
{
  TF1 f("cal_ip","[1]+pow([0]/(x-[2]),0.25)",-10.,900);
  f.SetParameters(cal_A0,cal_x0,cal_B0);
  
  *charge=f(max);
  *chargeerr=-0.25*(pow(cal_A0/(max-cal_B0),1.25))*amperr/cal_A0;
    
   
  if(max>=cal_B0)
    {
      *charge=1500.;
      *chargeerr=1000.;
    }
  
  return;
}


// Integrates the pulse to give charge - readout impedance must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, pulse charge, total charge, readout resistance)
void Charges(Float_t amp[], Float_t Dt, Int_t t_start, Int_t t_stop, Float_t *area, Float_t *area_pos, double RMSbck)
{
  double tempminch=0;
  double tempminch2=0;
  double bck=0;
  double TMax = 0;
  double AMax = 0;
  int NFlag = 0;
  for(int j=0;j<t_start;j++)
    {
      bck += amp[j];
      tempminch++;
      
    }
  
  bck *=1/tempminch;
  tempminch=0;
  
  //  std::cout << "start - stop " << t_start << " - " << t_stop << std::endl;
  for(int j=t_start;j<t_stop;j++)
    {
      if( amp[j]> AMax)
	{
	  TMax = j;
	  AMax = amp[j];
	}
    }
  for(int j=t_start;j<t_stop;j++)
    {
      //   if there is a signal, do not sum the undershoot
      if (AMax> 3*RMSbck)
	{
	  if (j<=TMax)
	    tempminch2 += amp[j]-bck;
	  else  if (j>TMax && amp[j]-bck > 0 && NFlag ==0)
	    {
	      tempminch2 += amp[j]-bck;
	    }
	  else if (j>TMax && amp[j]-bck < 0)
	    {
	      NFlag = 1;
	    }
	}
      else
	tempminch2 += amp[j]-bck;
      
      tempminch += amp[j]-bck;
    }
      *area_pos=tempminch2*Dt;
      *area=tempminch*Dt;
      
  return;
}


//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, threshold to comupute efficiency, efficiency)
void Efficiencies(Int_t camp, Float_t amp[], Float_t timeS[], Float_t t_I, Float_t t_II, Float_t bck,Float_t threshold, Float_t *eff)
{
  *eff=0;
  for(int j=0;j<camp;j++)
    {
      if(timeS[j]<t_I)
	continue;
      if(timeS[j]>t_I && timeS[j]<t_II)
	{
	  if((amp[j]-bck-threshold)>0 && threshold>0)
	    {
	      *eff=1;
	      break;
	    }
	  if((amp[j]-bck-threshold)<0 && threshold<0)
	    {
	      *eff=1;
	      break;
	    }
	}
    }
  return;
}

// Time at which a fixed threshold is passed. Makes a linear fit on 6 points and inverts the equation
// time  = (threshold-q)/q
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, threshold to measure time, time at which the threshold is passed)
void Thrtime_old(Int_t camp, Float_t amp[], Float_t timeS[], Float_t t_I, Float_t t_II, Float_t bck,Float_t threshold, Float_t *thrtime)
{
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=0;j<camp;j++)
    {
      if(timeS[j]<t_I)
	continue;
      if(timeS[j]>t_II)
	continue;
      if(timeS[j]>t_I && timeS[j]<t_II)
	{
	  //	  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
	  if((amp[j]-bck-threshold)>0 && threshold>0)
	    {
	      //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	      //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	      //std::cout << "Der = " << Der << std::endl;
	      //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //   std::cout << "th1 = " << *thrtime << std::endl;
	      for ( int kk = 0;kk<npoints;kk++)
		{
		  g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
		  //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
		  //  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
		}

	      break;

	    }
	  //	  if((amp[j]-bck-threshold)<0 && threshold<0)
	  // {
	  //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
	}
    }
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
    //    Der= lin->GetParameter(1);
  if (lin->GetParameter(1)>0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *thrtime = 0;
  // std::cout << "th1bis = " << *thrtime << std::endl;
  g->Delete();
  lin->Delete();
  return;
}

// Time at which a threshold relative to the amplitude of the pulse is passed. Amplitude of the pulse must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, relative threshold to measure time, time at which the threshold is passed)


void ConstFractime(Float_t amp[], Float_t timeS[], Float_t bck,Float_t threshold, Float_t max, Int_t NMax,Float_t *thrtime)
{
  
  
  *thrtime=0;

  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  //  *thrtime=0;
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  //  for(int j=0;j<camp;j++)
    for(int j=NMax;j>1;j--)
    {
      //      if(timeS[j]>t_max || timeS[j]<t_I )
      //	continue;
      // std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  //	  for ( int kk = 0;kk<npoints;kk++)
	  // {
	  //   g->SetPoint(kk,(amp[j-npoints/2+kk]-bck)/max,timeS[j-npoints/2+kk]);
	      //   std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << (amp[j-npoints/2+kk]-bck)/max << " " << threshold << std::endl;
	  // }
	  C2 = (amp[j+1]-bck)/max;
	  C1 = (amp[j]-bck)/max;
	  //  std::cout << "time = " << timeS[j]  <<  " " << timeS[j+1]<<  " C1 = " << C1 << " C2 = " << C2 <<std::endl;
	  if ( (C2-C1) != 0) TCrossing =  timeS[j]+(timeS[j+1]-timeS[j])*(threshold-C1)/(C2-C1);
	  break;
	  
	}
	
    }
  //  std::cout << "th1 = " << *thrtime << std::endl;
    //  TF1 *lin = new TF1("lin","[0]+[1]*x");
    // g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  
    //  if (lin->GetParameter(1) !=0)
    //  *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
    //  *thrtime= lin->GetParameter(0) +threshold*lin->GetParameter(1);
    *thrtime= TCrossing;
    //else
    // *thrtime = 0;
    
  //  std::cout << "th fit = " << *thrtime <<  " Th frac = " << TCrossing << " "  << " C1 = " << C1 << " C2 = " << C2 << std::endl;
    // g->Delete();
    // lin->Delete();
  return;
}


// Time at which a threshold relative to the amplitude of the pulse is passed. Amplitude of the pulse must be known
//(record lenght, amplitude vector, time vector, upper limit for background, upper limit for pulse, average background amplitude, relative threshold to measure time, time at which the threshold is passed)
void Zerofittime(Int_t camp, Float_t amp[], Float_t timeS[], Float_t t_I, Float_t t_II, Float_t bck,Float_t threshold1, Float_t threshold2, Float_t t_max,Float_t *thrtime)
{
  t_II = 0;
  *thrtime=0;
  double N=0;
  double sumV=0;
  double sumT=0;
  double sumVT=0;
  double sumTT=0;
  for(int j=camp-1;j>=0;j--)
    {
      if(timeS[j]>t_max)
	continue;
      if(timeS[j]<t_I)
	continue;
      if(amp[j]-bck>threshold1 && amp[j]-bck<threshold2 && timeS[j]<t_max && N<3)
	{
	  N=N+1;
	  sumV=sumV+amp[j]-bck;
	  sumT=sumT+timeS[j];
	  sumVT=sumVT+(amp[j]-bck)*timeS[j];
	  sumTT=sumTT+timeS[j]*timeS[j];
	}
    }
  if(N>0)
    {
      double Vm=(N*sumVT-sumV*sumT)/(sumTT-sumT*sumT);
      double V0=(sumV-Vm*sumT)/N;
      *thrtime=-V0/Vm;
    }
  if(N<2)
    *thrtime=-99;
  return;
}

// Time at which a fixed threshold is passed (descending)
//
void Trailtime(Int_t camp, Float_t amp[], Float_t timeS[], Int_t NMax, Float_t bck,Float_t threshold, Float_t *thrtime)
{
  //  std::cout << "NMAX = " << NMax << std::endl;
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j<camp;j++)
    {
     
      //  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if((amp[j]-bck-threshold)<0 && (amp[j-1]-bck-threshold)>0 &&  threshold>0)
	{
	  //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	  //std::cout << "Der = " << Der << std::endl;
	  //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //	       std::cout << "th1 = " << timeS[j] << std::endl;
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
	      //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
	      //	  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
	    }
	  
	  break;
	  
	}
      //	  if((amp[j]-bck-threshold)<0 && threshold<0)
      // {
      //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
    }
    
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  //    Der= lin->GetParameter(1);
  if (lin->GetParameter(1)<0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *thrtime = 0;
  //  std::cout << "th1bis = " << *thrtime << std::endl;
  g->Delete();
  lin->Delete();
  return;
}


// Time at which a fixed threshold is passed (going to the maximum)
// The search starts from the maximum and goes backwards
void Thrtime( Float_t amp[], Float_t timeS[], Int_t NMax, Float_t bck,Float_t threshold, Float_t *thrtime)
{
  //  std::cout << "NMAX = " << NMax << std::endl;
  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j>0;j--)
    {
     
      //std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if((amp[j]-bck-threshold)<0 && (amp[j+1]-bck-threshold)>0 &&  threshold>0)
	{
	  //	      *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //	      Der = fabs(  (threshold - (amp[j-1]-bck))/((amp[j]-amp[j-1])));
	  //std::cout << "Der = " << Der << std::endl;
	  //*thrtime=timeS[j-1]+Der*fabs(timeS[j]-timeS[j-1]);
	      //	       std::cout << "th1 = " << timeS[j] << std::endl;
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],amp[j-npoints/2+kk]-bck);
	      //	  g->SetPoint(kk,amp[j-npoints/2+kk]-bck,timeS[j-npoints/2+kk]);
	      //	  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << amp[j-npoints/2+kk]-bck << std::endl;
	    }
	  
	  break;
	  
	}
      //	  if((amp[j]-bck-threshold)<0 && threshold<0)
      // {
      //   *thrtime=(timeS[j]+timeS[j-1])/2.;
	  //   break;
	  // }
    }
    
  TF1 *lin = new TF1("lin","[0]+[1]*x");
  g->Fit("lin","QN","goff");
  //  g->Fit("lin");
  //    Der= lin->GetParameter(1);
  // if (lin->GetParameter(1)<0)
    *thrtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
    //else
    // *thrtime = 0;
    //    std::cout << "th1bis = " << *thrtime << std::endl;
    g->Delete();
    lin->Delete();
  return;
}


// Time at which a constant fraction threshold is passed (descending)
//(record lenght, amplitude vector, time vector, average background amplitude, relative threshold to measure time, starting time for the search (example: time of the matimum), time of the trailing edge, isteresys with respect to the leading edge)
void TrailConstFractime(Int_t camp, Float_t amp[], Float_t timeS[],Float_t bck,Float_t threshold, Float_t max, Int_t NMax,Float_t *trailtime)
{
  *trailtime=0;

  TGraph *g = new TGraph();
  int npoints=20;
  g->Set(npoints);	      
  //  *thrtime=0;
  //  float Der = 0;
  for(int j=NMax;j<camp;j++)
    {
      //	  std::cout << "sample = " << j << " " << (amp[j]-bck) << " " << threshold << " time = " << timeS[j] << std::endl;
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  for ( int kk = 0;kk<npoints;kk++)
	    {
	      g->SetPoint(kk,timeS[j-npoints/2+kk],(amp[j-npoints/2+kk]-bck)/max);
	      //		  std::cout << kk << " " <<  timeS[j-npoints/2+kk] << " " << (amp[j-npoints/2+kk]-bck)/max << " " << threshold << std::endl;
	    }
	  
	  break;
	  
	}
      
    }
    // std::cout << "th1 = " << *trailtime << std::endl;

    // TF1 *quad = new TF1("quad","[0]+[1]*x+[2]*x*x");
    TF1 *lin = new TF1("lin","[0]+[1]*x");
    g->Fit("lin","QN","goff");
    //g->Fit("lin");
    // g->Fit("quad");
  if (lin->GetParameter(1) !=0)
    *trailtime=(threshold - lin->GetParameter(0))/lin->GetParameter(1);
  else
    *trailtime = 0;
    
  //  std::cout << "th2 = " << *trailtime << std::endl;
  g->Delete();
  lin->Delete();

  
  return;
}


// Says if the event suffered from saturation of the waveform readout or not. Must be calibrated each time to evaluate the sensitivity of the routine
//(record lenght, amplitude vector, time vector,saturation flag)
void Saturation(Int_t camp, Float_t amp[], Float_t timeS[], Float_t bck, Int_t *satur)
{
  *satur=0;
  timeS = 0;
  int cont=0;
  int tempcont=0;
  double value=100000;
  for(int j=1;j<camp;j++)
    {
      if(fabs(amp[j]-amp[j-1])>0.0001)
	{
	  tempcont=0;
	}
      if(fabs(amp[j]-amp[j-1])<0.0001)
	{
	  tempcont++;
	  if(tempcont>cont)
	    {
	      cont=tempcont;
	      value=amp[j];
	    }
	}
    }
  if(fabs(value-bck)>1. && cont>int(camp/500)+5)
    *satur=1;
  return;
}

// Computed the charge centroid using the A0/cosh((x-x0)/delta0) function.
//(Charge, Error of charge, number of RO channels, distribution width, distribution width 2nd method, Distribution const, Distribution const 2nd method, Centroid , Centroid 2nd method , ChiSquare, ChiSquare 2nd method,pitch, Cluster size, threshold)
void Centroid(Float_t ch[], Float_t cherr[],Int_t nchro, Float_t *delta0, Float_t *delta1, Float_t *A0, Float_t *A1, Float_t *x0, Float_t *x1, Float_t *chinorm0, Float_t *chinorm1, Float_t pitch, Int_t CS, Float_t threshold)
{
  TF1 cint("cint",Form("[0]*atan(exp((x+%f-[1])/[2]))-[0]*atan(exp((x-%f-[1])/[2]))",pitch/2.,pitch/2.),-100.,100.);
  TF1 cdis("cdis","[0]*(1./cosh((x-[1])/[2]))",-100.,100.);
  TF1 cintf("cintf",Form("[0]*atan(exp((x+%f-[1])/[2]))-[0]*atan(exp((x-%f-[1])/[2]))",pitch/2.,pitch/2.),-100.,100.);
  TGraphErrors gbar(CS);
  double tmpmax=-100;
  double tmpbar=0;
  int k=0;
  *delta0=-1;
  *delta1=-1;
  *A0=-1;
  *A1=-1;
  *x0=-10;
  *x1=-10;
  *chinorm0=-1;
  *chinorm1=-1;
  cintf.FixParameter(2,3.9);
  if(CS>2)
    {
      for(int i=0;i<nchro;i++)
	{
	  if(ch[i]>tmpmax)
	    {
	      tmpmax=ch[i];
	      tmpbar=double(i)*pitch+pitch/2.;
	    }
	}
      
      cint.SetParameters(tmpmax,tmpbar,3.5);
      cint.SetParLimits(1,-pitch*5.,double(nchro)*pitch+pitch*5.);
      cintf.SetParameters(tmpmax,tmpbar,3.9);
      cintf.FixParameter(2,3.9);
      cintf.SetParLimits(1,-pitch*5.,double(nchro)*pitch+pitch*5.);

      for(int i=0;i<nchro;i++)
	{
	  if(ch[i]>threshold)
	    {
	      k++;
	      gbar.SetPoint(k,double(i)*pitch+pitch/2.,ch[i]);
	      gbar.SetPointError(k,0.,cherr[i]);
	    }
	}
      
      gbar.Fit("cint","QN","",0.,nchro*pitch+pitch/2.);
      *delta0=cint.GetParameter(2);
      *x0=cint.GetParameter(1);
      *A0=cint.GetParameter(0);
      *chinorm0=cint.GetChisquare()/cint.GetNDF();
      gbar.Fit("cintf","QN","",0.,nchro*pitch+pitch/2.);
      *delta1=cintf.GetParameter(2);
      *x1=cintf.GetParameter(1);
      *A1=cintf.GetParameter(0);
      *chinorm1=cintf.GetChisquare()/cintf.GetNDF();
    }
  if(CS==2)
    {
      *x0=0;
      float tempnum=0;
      float tempden=0;
      float bar0temp=0;
      float tempor=0;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      bar0temp=bar0temp+(float(j)*pitch+pitch/2.);
	    }
	}
      bar0temp=bar0temp/2.;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      tempnum=tempnum+ch[j]*(float(j)*pitch+pitch/2.-bar0temp);
	      tempden=tempden+ch[j];
	      //cout << FEcharge[j] << "\t";
	    }
	}
      tempor=tempnum/tempden;
      *x0=bar0temp+7.5*tempor;  //7.5 conversion factor from charge centroid to mm
    }
  if(CS==1)
    {
      *x0=0;
      for(int j=0;j<nchro;j++)
	{
	  if(ch[j]>threshold)
	    {
	      *x0=*x0+(float(j)*pitch+pitch/2.);
	    }
	}
    }
  return;
}


//  Computes the slope of the rise of a signal as an exponential
void riseexpfit(Int_t camp, Float_t amp[], Float_t timeS[],Float_t bck,Float_t t_low,  Float_t t_high, Float_t *rise_exp0, Float_t *rise_exp1, Float_t *chi)
{
  TGraph *g = new TGraph();
  int npoints=0;
  int nnpoints=0;
  int Nlow=0;

  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)  continue;
      if(timeS[i]>t_high) break;
      npoints++;
      if (npoints == 1) 
	{
	  Nlow = i;
	}
    }

  g->Set(npoints);
  nnpoints=0;
  for(int i=Nlow;i<Nlow+npoints-1;i++)
    {
      g->SetPoint(nnpoints,timeS[i]-t_low,amp[i]-bck);
	    //      g->SetPoint(nnpoints,timeS[i],amp[i]-bck);
	    //  std::cout << timeS[i]<< " " << amp[i]-bck << std::endl;
      nnpoints++;
    }
  
  TF1 *exp = new TF1("exp","exp([0])*exp([1]*x) +[2]");
  

  g->Draw();


  exp->SetParLimits(1, 0.,5.);
  //  exp->SetParLimits(2, -5.,5.);
  exp->SetParameters(-1,1.,bck);
  
  g->Fit("exp","QN","goff",0.,t_high-t_low);
  // g->Fit("exp"," "," ",0.,t_high-t_low);
  *rise_exp0=exp->GetParameter(0);
  *rise_exp1=exp->GetParameter(1);
 
  *chi=exp->GetChisquare()/exp->GetNDF();

  //  std::cout << " rise = " << *rise_exp1 << " chi2=  " << *chi << std::endl;
  exp->Delete();
  g->Delete();

  return;
}


void riselinfit(Int_t camp, Float_t amp[], Float_t timeS[],Float_t bck,Float_t t_low,  Float_t t_high, Float_t *rise_lin0, Float_t *rise_lin1, Float_t *chi)
{
  TGraph *g = new TGraph();
  int npoints=0;
  int nnpoints=0;
  int Nlow=0;
  double_t amp_low =0;
  double_t amp_high =0;
  //  t_low = t_high-0.2;
  // std::cout << " t low " << t_low << " t_high " << t_high << std::endl;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)  continue;
      if(timeS[i]>t_high) break;
      npoints++;
      if (npoints == 1) 
	{
	  Nlow = i;
	}
    }
  amp_low = amp[Nlow]-bck;
  amp_high = amp[npoints]-bck;
  g->Set(npoints);
  nnpoints=0;
 
  for(int i=Nlow;i<Nlow+npoints-1;i++)
    {
      g->SetPoint(nnpoints,timeS[i],amp[i]-bck);
      //  std::cout << nnpoints << " \t" << timeS[i]<< " \t" <<amp[i]-bck << " \n";
      nnpoints++;
    }
  
  TF1 *lin = new TF1("lin","[0]+[1]*x",t_low,t_high);
  // TF1 *exp = new TF1("exp","[0]*exp([1]*x)",t_low,t_high);

  lin->SetParameters(amp_low-t_low*(amp_high-amp_low)/(t_high-t_low),(amp_high-amp_low)/(t_high-t_low));
  g->Draw();
  g->Fit("lin","QN","goff",t_low,t_high);
  // g->Fit("lin"," ","g",t_low,t_high);
  *rise_lin0=lin->GetParameter(0);
  *rise_lin1=lin->GetParameter(1);
  *chi=lin->GetChisquare()/lin->GetNDF();
  lin->Delete();
  g->Delete();

  return;
}
// Computes the mobile average of the signal
// (record lenght,       amplitude vector, time vector, number of points for the mobile average (odd number), averaged amplitude vector)
void mobileAVG(Int_t camp, Float_t amp[], Int_t navg, Float_t m_amp[])
{
  float tempsum=0;
  for(int j=int((navg-1)/2);j<camp-int((navg-1)/2)-1;j++)
    {
      tempsum=0;
      for(int k=-int((navg-1)/2);k<int((navg-1)/2)+1;k++)
	{
	  tempsum=tempsum+amp[j+k];
	}
      m_amp[j]=tempsum/float(navg);
    }
  for(int j=0;j<int((navg-1)/2);j++)
    {
      m_amp[j]=m_amp[int((navg-1)/2)];
    }
  for(int j=camp-int((navg-1)/2)-1;j<camp;j++)
    {
      m_amp[j]=m_amp[camp-int((navg-1)/2)-2];
    }

  return;
}

void Centroid_N( Float_t TIMEUNIT, Int_t Nmin, Int_t Nmax, Float_t Itot[], double *t_centr)
{
  float Num =  0;// BB RC in seconds
  float Den = 0; // BB Tau in seconds

  //tmp = new double[2*camp];//output voltage
  for (int i=Nmin;i<Nmax;i++)
    {
      Num += Itot[i]*TIMEUNIT*i;
      Den += Itot[i];
      //      tmp[i]=0;
    }
  * t_centr =  Num/Den;
  // std::cout  <<  t_centr  << " " <<  Num <<" " <<Den << "\n";
  return;
}


void CSA_NA62( Float_t TIMEUNIT, Int_t camp, Float_t Itot[],Float_t ShaperOut_Q[], Int_t *IMaxSh)
{

  *IMaxSh = 7./TIMEUNIT; //TimeUnit in [ns]
  double tr = 5.6E-9 ; 
  double til = 1.8E-9 ;
  double tac = 47.E-9 ;
  double tau = 0;
  double tt = 0;
  double qtot = 0; 
  double  Qdif_Shaper = 0;
  double sh_max = 0;
  TIMEUNIT *=1e-9; // in seconds

  for (int i=0;i< *IMaxSh;i++)  ShaperOut_Q[i] = 0;
  for (int i=0;i<camp;i++)
    {
      qtot +=  Itot[i]* TIMEUNIT*1e-6;

      if (fabs(Itot[i]) == 0) continue;
      //   std::cout <<  "i = " << i <<  " Itot[i] = " << Itot[i] << " \n";
      Qdif_Shaper  = Itot[i]*1e-6*(TIMEUNIT);// This is the charge on the Cap      

      for (int ll = 0; ll<*IMaxSh-i;ll++)
	{

  // NA62
  
	  tt = ll*TIMEUNIT;
	  tau = 0;
	  ShaperOut_Q[i+ll] += Qdif_Shaper*(-1.*tac * tr)
	    * ( ((-1. * exp(-(tt-tau)/tac) * tac) / (pow((tac-til),2.) * (tac -tr)))+((exp(-(tt-tau)/tr) * tr) / (pow((til-tr),2.) * (tac -tr))) +
		((exp(-(tt-tau)/til) * (pow(til,3.) + (tt-tau)*(tac-til)*(til-tr) -tac*til*tr)) / (pow((tac - til),2.) * til*pow((til - tr),2.))) ) ;

	  if (fabs(ShaperOut_Q[i+ll]) > fabs(sh_max))  sh_max  = ShaperOut_Q[i+ll]; 
	}

    }
  for (int ll = 0; ll<*IMaxSh;ll++)
    {
      ShaperOut_Q[ll] *=  75*1e15*qtot/sh_max; // [mV/fQ *Q/Q]
      // cout << ShaperOut_V[i] << endl;
    }
      

  // std::cout << qtot << " " << sh_max << " " << 75.*qtot/sh_max << "\n";
  // for (int ll = 0; ll<5*camp-i;ll++) std::cout  << " ll = " << ll << " time =  " << ll*TIMEUNIT*1e9 << " " << ShaperOut_Q[ll] << "\n";
  
  return;
}
void BBA( Float_t TIMEUNIT, Float_t Cdet, Float_t BBImp, Float_t BBGain, double_t BBABW, Int_t camp, Float_t Itot[],Float_t BBGraph[],Int_t *IMaxSh)
{

  *IMaxSh = 2./TIMEUNIT; //Number of time steps to be computed 
  double tau_BB_RC = 1.0e-12*BBImp*Cdet; // BB RC in seconds
  double tau_BB_BW = 0.35/(1.0e9*BBABW)/2.2; // BB Tau in seconds
  TIMEUNIT *=1e-9; // in seconds
  double Qi = 0;
  double   *VonR_BB_RC = 0;

  
  VonR_BB_RC=new double[*IMaxSh];//output voltage
  
  for (int i=0;i<*IMaxSh;i++)
    {
      BBGraph[i]=0;
      VonR_BB_RC[i] = 0;

    }
  
  for (int i=0;i<camp;i++)
    {
            
      if (i == 0) Qi = Itot[0]*1e-6*(TIMEUNIT);
      else      Qi  = (Itot[i]-Itot[i-1])*1e-6*(TIMEUNIT);// This is the charge on the Cap added in step i
      
      //     pippo += Idif/(1e-6*(TIMEUNIT));
      //	  std::cout << " i-1 = " << i-1 << " " << Itot[i-1] <<  " " << pippo << std::endl;
      //   std::cout << " i = " << i << " " <<  Qi << " " << Itot[i] << "\n";
      for (int ll = 0; ll<*IMaxSh-i;ll++)
	{
	  
 //Voltage on Rinput, this is the discharge of Qi on Rinput
	  
	  VonR_BB_RC[i+ll] += Qi/(1e-12*Cdet)*TMath::Exp(-ll*TIMEUNIT/tau_BB_RC);	  
	  BBGraph[i+ll] += 1e+3*BBGain*VonR_BB_RC[i]*(1.-TMath::Exp( (float) -ll*TIMEUNIT/tau_BB_BW));
	  
	}
    }



  delete[] VonR_BB_RC;
  return;
}
void CSA( Float_t TimeUnit, Float_t TransImp, Float_t TauFall, double_t TauRise, Int_t camp, Float_t Itot[],Float_t CSAout[], Int_t *IMaxSh)
{

  *IMaxSh = 12./TimeUnit; //TimeUnit in [ns]
  int IMax;
  float  PreAmp_Q[*IMaxSh];
  float CSAMax = 0;
  float Qtot = 0;

  for(int k=0;k<*IMaxSh;k++)
    {
      PreAmp_Q[k]=0.0;
      CSAout[k] = 0.0;
    }


  IMax = camp; // In case Itot never goes to zero  
  for(int j=0;j<camp;j++)
    {
      if ( j>2 && Itot[j] == 0)
	{
	  IMax = j;
	  break;
	}
    }
  
  for(int i=1;i<IMax;i++)
    {
      
      PreAmp_Q[i]= Itot[i]*TimeUnit;
      for (int ll = 0; ll<*IMaxSh-i;ll++)  // valid only up to IMaxSh 
	{
	  CSAout[i+ll] += TauFall/(TauFall+TauRise)*(PreAmp_Q[i])*
	    (TMath::Exp(-ll*TimeUnit/TauFall)-TMath::Exp(-ll*TimeUnit/TauRise)); // [Q] HS eq 4.3 This is the shaper	  
	  if (fabs(CSAout[i+ll]) > fabs(CSAMax))  CSAMax = CSAout[i+ll]; 
      	}
	  

      Qtot += Itot[i]*TimeUnit;
      
     }

  for (int ll = 0; ll<*IMaxSh;ll++)
    {
      CSAout[ll] *= fabs(TransImp*Qtot/CSAMax);
      //std::cout  <<  "CSA = " << CSAout[ll] << std::endl;
    }
  //   std::cout  << " Qtot = " << Qtot <<  " CSAMax  = " << CSAMax << std::endl;
  // for (int ll = 0; ll<IMaxSh;ll++)  // valid only up to IMaxSh 
  // std::cout <<  " Total Charge = " << Qtot << " Tpeak = " << Tpeak << " Ball. Deficit = " << BallDef << " CSAMax = "<< CSAMax << std::endl;       

  
  return;
}


//
//Calcola la derivata a quattro punti
// (record lenght, amplitude vector,reshaping vector)
void Derivative(Float_t TimeUnit,Int_t camp, Float_t amp[],Float_t resh[])
{

  TGraph *g = new TGraph();
  int npoints=3;
  int nnpoints=0;
  // int Nlow=0;
  g->Set(2*npoints);
 
  for(int j=npoints;j<camp-npoints;j++)
	{
	   nnpoints=0;
	  for (int kk = j-npoints; kk< j+npoints; kk++)
	    {
	      g->SetPoint(nnpoints,nnpoints*TimeUnit,amp[kk]);
	      //	      std::cout << "Fit = " <<  nnpoints << " \t" << nnpoints*TimeUnit << " \t" <<amp[kk] << " \n";
	      nnpoints++;
	    }	    
	  // resh[j]=(amp[j+10]-amp[j-10])/(2*npoints*TimeUnit);
	  TF1 *lin = new TF1("lin","[0]+[1]*x",0,2*nnpoints*TimeUnit);
	  g->Fit("lin","QN","goff",0,2*nnpoints*TimeUnit);
	  //*rise_lin1=lin->GetParameter(1);
	  resh[j]=lin->GetParameter(1);
	  //	  std::cout << "j = " << j << " "  << resh[j] << " " << *rise_lin1 << std::endl;
	  lin->Delete();
	}
   g->Delete();

return;
}



// Generates FFT of a real vector (inreal) of lenght camp. Real part of output FFT is outre, immaginary part is outim. Module output is outmod.
void FFTrealtocomplex(Int_t camp, Float_t inreal[], Float_t outre[], Float_t outim[], Float_t outmod[])
{
  Double_t re;
  Double_t im;
  Double_t *in = new Double_t[camp];
  for(int i=0;i<camp;i++)
    {
      *(&in[i])=inreal[i];
    }
  TVirtualFFT *fftr2c = TVirtualFFT::FFT(1, &camp, "R2C");
  fftr2c->SetPoints(in);
  fftr2c->Transform();
  for (Int_t i=0; i<camp; i++)
    {
      fftr2c->GetPointComplex(i, re, im);
      outre[i]=re;
      outim[i]=im;
      outmod[i]=sqrt(re*re+im*im);
    }
  fftr2c->Delete();
  delete[] in;

  return;
}


// Generates the inverse FFT from a complex vector of lenght camp (two vectors: inre real part, inim immaginary part). The output vector is outreal.
void FFTcomplextoreal(Int_t camp, Float_t inre[], Float_t inim[], Float_t outreal[])
{
  Double_t *imma = new Double_t[camp];
  Double_t *real = new Double_t[camp];
  
  for(int i=0;i<camp;i++)
    {
      *(&real[i])=inre[i];
      *(&imma[i])=inim[i];
    }
  TVirtualFFT *fftc2r = TVirtualFFT::FFT(1, &camp, "C2R");
  fftc2r->SetPointsComplex(real,imma);
  fftc2r->Transform();
  for (Int_t i=0; i<camp; i++)
    {
      outreal[i]=fftc2r->GetPointReal(i)/double(camp);
    }
  fftc2r->Delete();
  delete[] imma;
  delete[] real;

  return;
}


// Redines coordinates correction distorsion of centroid. polinomial correction of coordinates (expansion at center of strips and compressio at borders)
//(centroid variable address, curve parameter)
void positionbias(Float_t *bar,Float_t sca)
{
  sca = 0;
  float barfz=*bar;
  float xloop=(1/100000.)*(int((100000.)*(barfz-4))%800000);
  
  *bar=((1/(70.+16.))*((70.)*(xloop-4)+pow((xloop-4),3))+(barfz-xloop-4)+8);

  return;
}


//Function to fit the waveform with Gumbel function
void Gumbel_fit(Int_t camp, Float_t amp[], Float_t timeS[],Float_t bck,Float_t t_low,Float_t t_up,Float_t t_amp, Float_t t_width, Float_t max, Float_t *bck_fit, Float_t *t_amp_fit, Float_t *t_width_fit, Float_t *amp_fit)
{
  TF1 *fgum = new TF1("gumbel","[0]+[3]*(exp(-(x-[1])/[2]+exp(-(x-[1])/[2])))",t_low,t_up);
  fgum->SetParameters(bck,t_amp,t_width,max);

  TGraph *g = new TGraph();
  int npoints=0;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)
	continue;
      if(timeS[i]>t_up)
	break;
      npoints++;
    }
  g->Set(npoints);
  npoints=0;
  for(int i=0;i<camp;i++)
    {
      if(timeS[i]<t_low)
	continue;
      if(timeS[i]>t_up)
	break;
      g->SetPoint(npoints,timeS[i],amp[i]);
      npoints++;
    }

  g->Draw();
  g->Fit("gumbel","QN","goff",t_low,t_up);
  *bck_fit=fgum->GetParameter(0);
  *t_amp_fit=fgum->GetParameter(1);
  *t_width_fit=fgum->GetParameter(2);
  *amp_fit=fgum->GetParameter(3);

  g->Delete();
  fgum->Delete();

  return;
}

double string_to_double( const std::string& s )
 {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 } 





void ConstFractime_N(Float_t amp[], Float_t DT, Float_t bck,Float_t threshold, Float_t max, Int_t NMax,Float_t *thrtime, bool Interpolation, Int_t  *Nth)
{    
  *thrtime=0;

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  //  for(int j=0;j<camp;j++)
    for(int j=NMax;j>1;j--)
    {
      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j+1]-bck)/max; // Above Vth; too late
	      C1 = (amp[j]-bck)/max;  // Below Vth; too early
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1); //time[j] is too late, so I need to subtract a fraction of time bin	  
	    }
	  else   TCrossing =  timeS[j];

	  break;	  
	}
	
    }
    *thrtime= TCrossing;

  return;
}




void Trailtime_N(Int_t camp, Float_t amp[], Float_t DT, Int_t NMax, Float_t bck,Float_t threshold, Float_t *thrtime, bool Interpolation, Int_t  *Nth)
{
  
  *thrtime=0;
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  for(int j=NMax;j<camp;j++)
    {
      if((amp[j]-bck-threshold)<0 && (amp[j-1]-bck-threshold)>0 &&  threshold>0)
	{
	    if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j-1]-bck); //above Vth, too early
	      C1 = (amp[j]-bck);  //below Vth, too late	      
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1);	  
	    }
	  else   TCrossing =  timeS[j];

	  break;
	  
	}
    }
    
  *thrtime= TCrossing;
  return;
}

void TrailConstFractime_N(Int_t camp, Float_t amp[], Float_t DT,Float_t bck,Float_t threshold, Float_t max, Int_t NMax,Float_t *trailtime, bool Interpolation, Int_t  *Nth)
{
  *trailtime=0;
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  for(int j=NMax;j<camp;j++)
    {

      if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	{
	  if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j-1]-bck)/max; //above Vth, too early
	      C1 = (amp[j]-bck)/max;  // below vth, too late
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1); //timeS[j] is too late, I need to subtract  a bit of DT
	      //	      std:: cout << " C1 = " << C1 << " C2 = " << C2 << " threshold = " << threshold << std:: endl;
	      //
	      // std:: cout <<  timeS[j-1] << "  threshold - C1 = " << threshold - C1 << " C2 - C1 = " << C2-C1 << " TCrossing = " << TCrossing << std:: endl; 
	    }
	  else   TCrossing =  timeS[j];
	  
	  break;	  
	}
      
    }
  *trailtime= TCrossing;  
  return;
}


#endif


void Levtime_N( Float_t amp[], Float_t DT, Int_t NMax, Float_t bck,Float_t threshold, Float_t *thrtime, bool Interpolation, Int_t  *Nth)
{

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  for(int j=NMax;j>0;j--)
    {
 
      if((amp[j]-bck-threshold)<0 && (amp[j+1]-bck-threshold)>0 &&  threshold>0)
	{
	  if (Interpolation)
	    {
	      *Nth = j;
	      C2 = (amp[j+1]-bck); // above Vth
	      C1 = (amp[j]-bck); //Below Vth
	      if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1);  //time[j] is too early, so I need to add a fraction of time bin  
	    }
	  else   TCrossing =  timeS[j];

	  break;
	  
	}
    }
    *thrtime= TCrossing;    

  return;
}


void Levtime_N_array(int NTValues[], Float_t TValues[], Int_t NArray, Float_t Threshold[], Float_t amp[], Float_t DT, Int_t NMax, Float_t bck, bool Interpolation)
{

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;

  NArray--;
  NTValues[NArray] = -100;
  TValues[NArray] = -100;
  int Nth;
  for(int j=NMax;j>0;j--)
    {
      if (NArray < 0) break;
      for(int k = NArray;k > -1;k--)
	{
	  threshold = Threshold[NArray]; 
	  if((amp[j]-bck-threshold)<0 && (amp[j+1]-bck-threshold)>0 &&  threshold>0)
	    {
	      if (Interpolation)
		{
		  Nth = j;
		  C2 = (amp[j+1]-bck); // above Vth
		  C1 = (amp[j]-bck); //Below Vth
		  if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1);  //time[j] is too early, so I need to add a fraction of time bin
		  
		  // std::cout << "Crossing 2-points: = "<< TCrossing << std::endl;
		  // std::cout << "Level = "<<  threshold << std::endl;
		  
		  //	      TF1 *f1 = new TF1("f1","[0]+[1]*x");
		  if (1>2) // dofit
		    {
		      TF1 *f1 = new TF1("f1","pol2");
		      TGraph *g = new TGraph();
		      g->Set(7);
		      for(int b=0;b<7;b++)
			{
			  g->SetPoint(b,amp[j-2+b]-bck,timeS[j-2+b]);
			  //  std::cout <<amp[j-2+b]-bck << " " << timeS[j-2+b]<<std::endl;
			}
		      g->Fit("f1","QN");
		      TCrossing = f1->GetParameter(0)+ threshold*f1->GetParameter(1)+ threshold*threshold*f1->GetParameter(2);
		      //  std::cout << "Parameter Crossing fits: = "<<f1->GetParameter(0) << " " << f1->GetParameter(1) <<  " " << f1->GetParameter(2)<< std::endl;
		      // std::cout << "Crossing fits: = "<< TCrossing << std::endl;
		      g->Delete();
		      f1->Delete();
		    }
		}
	  
	      else   TCrossing =  timeS[j];
	  
	      NTValues[NArray] = Nth;
	      TValues[NArray] = TCrossing ;
	  
	      //std::cout << " Level th = " <<   Threshold[NArray]  << " at time " <<  TValues[NArray] <<  " at index = " << NTValues[NArray] << std::endl;
	      NArray--;

	      //while(amp[j]-bck<Threshold[NArray]) j++;
	    

	  
	    }
	  else
	    break;
	}
      
    }
  
  return;
}


void Levtime_N_arrayGPT2(int NTValues[], Float_t TValues[], Int_t NArray, Float_t Threshold[], Float_t amp[], Float_t DT, Int_t NMax, Float_t bck, bool Interpolation) {
    float TCrossing = 0;
    NArray--;  // Adjust NArray to be zero-based
    NTValues[NArray] = -100;
    TValues[NArray] = -100;

    for (int j = NMax; j > 0; j--) {
        if (NArray < 0) break;  // Break if no more thresholds to process

        for (int k = NArray; k >= 0; k--) {
            Float_t threshold = Threshold[k];  // Correct threshold selection

            // Check for crossing the threshold
            if ((amp[j] - bck - threshold) < 0 && (amp[j + 1] - bck - threshold) > 0 && threshold > 0) {
                if (Interpolation) {
                    double C2 = (amp[j + 1] - bck);  // Above threshold
                    double C1 = (amp[j] - bck);      // Below threshold

                    // Calculate the exact crossing time using linear interpolation
                    if ((C2 - C1) != 0) {
                        TCrossing = j * DT + DT * (threshold - C1) / (C2 - C1);  
                    }

                    // Optional higher-order fit (typically not necessary for crossing points)
 
                } else {
                    TCrossing = j * DT;  // Use direct crossing time if no interpolation
                }

                NTValues[k] = j;
                TValues[k] = TCrossing;
                NArray--;  // Move to next threshold
            }
        }
    }
}



void TrailLevtime_N_array(Float_t TValues[], Int_t NArray, Float_t Threshold[], Int_t camp, Float_t amp[], Float_t DT, Int_t NMax, Float_t bck, bool Interpolation)
{

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;

  NArray--;
  TValues[NArray] = -100;
  int Nth;
  for(int j=NMax;j<camp;j++)
    {
       if (NArray < 0) break;
       for(int k = NArray;k > -1;k--)
	 {
	   
	   threshold = Threshold[NArray]; 
	   if((amp[j]-bck-threshold)<0 && (amp[j-1]-bck-threshold)>0 &&  threshold>0)
	     {
	       if (Interpolation)
		 {
		   Nth = j;
		   C2 = (amp[j-1]-bck); // above Vth
		   C1 = (amp[j]-bck); //Below Vth
		   if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1);  //time[j] is too early, so I need to add a fraction of time bin  
		 }
	       else   TCrossing =  timeS[j];
	       
	       TValues[NArray] = TCrossing ;
	       
	       NArray--;
	       
	       //      while(amp[j]-bck<Threshold[NArray]) j--;

	  
	     }
	   else
	     break;
	 }
       
    }
  
  return;
}

void TrailLevtime_N_arrayGPT2(Float_t TValues[], Int_t NArray, Float_t Threshold[], Int_t camp, Float_t amp[], Float_t DT, Int_t NMax, Float_t bck, bool Interpolation) {
    float TCrossing = 0.0;
    NArray--;  // Adjust to zero-based index
    TValues[NArray] = -100;  // Initialize the last element

    for (int j = NMax; j < camp; j++) {
        if (NArray < 0) break;  // Exit loop if no more thresholds to process

        for (int k = NArray; k >= 0; k--) {
            double threshold = Threshold[k];

            // Check for crossing the threshold from above
            if ((amp[j] - bck - threshold) < 0 && (amp[j - 1] - bck - threshold) > 0 && threshold > 0) {
                if (Interpolation) {
                    double C2 = (amp[j - 1] - bck);  // Above threshold
                    double C1 = (amp[j] - bck);      // Below threshold

                    // Calculate the exact crossing time using linear interpolation
                    if ((C2 - C1) != 0) {
                        TCrossing = j * DT - DT * (threshold - C1) / (C2 - C1);
                    } else {
                        TCrossing = j * DT;
                    }
                } else {
                    TCrossing = j * DT;  // Use direct crossing time if no interpolation
                }

                TValues[k] = TCrossing;
                NArray--;  // Move to the next threshold

                // Move to the next threshold without further checks for this crossing
                break;
            }
        }
    }
}

void ConstFractime_N_array(int NTValues[], Float_t TValues[], Int_t NArray, Float_t Threshold[], Float_t amp[], Float_t DT, Float_t bck, Float_t max, Int_t NMax, bool Interpolation)
{    
 
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;
  int Nth;
  // double thrtime;

  //  for(int j=0;j<camp;j++)
  
  NArray--;
  NTValues[NArray] = -100;
  TValues[NArray] = -100;
    for(int j=NMax;j>1;j--)
    {
      if (NArray < 0) break;
      for(int k = NArray;k > -1;k--)
	{
	  threshold = Threshold[NArray]; 
	  if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	    {
	      if (Interpolation)
		{
		  Nth = j;
		  C2 = (amp[j+1]-bck)/max; // Above Vth; too late
		  C1 = (amp[j]-bck)/max;  // Below Vth; too early
		  if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1); //time[j] is too late, so I need to subtract a fraction of time bin	  
		}
	      else   TCrossing =  timeS[j];
	      
	      NTValues[NArray] = Nth;
	      TValues[NArray] = TCrossing ;
	  
	      //	      std::cout << " Frac th = " <<   Threshold[NArray]  << " at time " <<  TValues[NArray] <<  " at index = " << NTValues[NArray] << std::endl;
	      NArray--;
	    
	    }
	  else
	    break;
	}

    }
  
    //    *ANTValues = NNTValues[0];
    // std::cout << " &NNTValues[0]  " <<  &NNTValues[0] << std::endl;
    // std::cout << " * ANTValues " <<   *ANTValues  << std::endl;
    // std::cout << " ANTValues " <<   ANTValues  << std::endl;
    // std::cout << " &ANTValues " <<   &ANTValues  << std::endl;
    return;
}

void ConstFractime_N_arrayGPT2(int NTValues[], Float_t TValues[], Int_t NArray, Float_t Threshold[], Float_t amp[], Float_t DT, Float_t bck, Float_t max, Int_t NMax, bool Interpolation)
{    
 
// pippo 

  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;
  int Nth;
  // double thrtime;

  //  for(int j=0;j<camp;j++)
  
 
//  NTValues[NArray] = -100;
//  TValues[NArray] = -100;
  for(int k = 0;k < NArray;k++)
    {
      threshold = Threshold[k];
//	std::cout << Threshold[NArray] << std::endl ;
      for(int j=NMax;j>1;j--)
	{

	  if(((amp[j]-bck)/max-threshold)<0 && threshold>0) // from above, below Vth
	    {
	      if (Interpolation)
		{
		  Nth = j;
		  C2 = (amp[j+1]-bck)/max; // Above Vth; too late
		  C1 = (amp[j]-bck)/max;  // Below Vth; too early
		  if ( (C2-C1) != 0) TCrossing =  timeS[j]+DT*(threshold-C1)/(C2-C1); //time[j] is too late, so I need to subtract a fraction of time bin	  
		}
	      else   TCrossing =  timeS[j];
	      
	      NTValues[k] = Nth;
	      TValues[k] = TCrossing ;
//		std::cout <<  Threshold[k] << " " << NTValues[k]  << " " << TValues[k] << std::endl;
		break;
	    }
	}
    }
}


void TrailConstFractime_N_array(int NTValues[], Float_t TValues[], Int_t NArray, Float_t Threshold[], Int_t camp, Float_t amp[], Float_t DT, Float_t bck, Float_t max, Int_t NMax, bool Interpolation)
{    
 
  float TCrossing = 0;
  float C1 = 0;
  float C2 = 0;
  float threshold;
  int Nth;
  // double thrtime;

  //  for(int j=0;j<camp;j++)
  
//  NArray--;
//  TValues[NArray] = -100;
//  NTValues[NArray] = -100;
  
    for(int j=NMax;j<camp;j++)
    {
//      if (NArray < 0) break;
      for(int k = 0;k < NArray;k++)
	{
	  threshold = Threshold[k]; 
	  if(((amp[j]-bck)/max-threshold)<0 && threshold>0)
	    {
	      if (Interpolation)
		{
		  Nth = j;
		  C2 = (amp[j-1]-bck)/max; // Above Vth; too late
		  C1 = (amp[j]-bck)/max;  // Below Vth; too early
		  if ( (C2-C1) != 0) TCrossing =  timeS[j]-DT*(threshold-C1)/(C2-C1); //time[j] is too late, so I need to subtract a fraction of time bin	  
		}
	      else   TCrossing =  timeS[j];
	      TValues[k] = TCrossing ;
	      NTValues[k] = Nth;
	      
	  //	  std::cout << " Frac th = " <<   Threshold[NArray]  << " at time " <<  TValues[NArray] <<  " at index = " << NTValues[NArray] << std::endl;
	      NArray--;	  
	    }
	  else
	    break;
	}
	
    }
  
    return;
}




void TrailConstFractime_N_arrayGPT2(int NTValues[], Float_t TValues[], Int_t NArray, Float_t Threshold[], Int_t camp, Float_t amp[], Float_t DT, Float_t bck, Float_t max, Int_t NMax, bool Interpolation) {    
    double TCrossing = 0.0;
    int Nth = -100;

    // Ensure NArray is valid and zero-based
    NArray--;
    TValues[NArray] = -100;
    NTValues[NArray] = -100;

    // Loop from NMax (peak) to the end of the array (camp)

 for(int k = NArray;k > -1;k--)
    {
      double  threshold = Threshold[k];
      
      for (int j = NMax; j < camp; j++) {
        

	//        for (int k = NArray; k >= 0; k--)
	//  {
	//   double threshold = Threshold[k];

            // Check if the amplitude crosses the threshold
            if (((amp[j] - bck) / max) < threshold && threshold > 0)
	      {
                if (Interpolation) {
                    Nth = j;
                    double C2 = (amp[j - 1] - bck) / max; // Above threshold; too late
                    double C1 = (amp[j] - bck) / max;    // Below threshold; too early

                    if ((C2 - C1) != 0) {
                        TCrossing = j * DT - DT * (threshold - C1) / (C2 - C1);  // Calculate crossing time using interpolation
                    } else {
                        TCrossing = j * DT;
                    }
                } else {
                    TCrossing = j * DT;  // Use direct time if no interpolation
                }

                TValues[k] = TCrossing;
                NTValues[k] = Nth;

                // Move to the next threshold
		//   NArray--;
                break;
            }
        }
    }
}

void calculateThreeCircleIntersectionN(Int_t xc[], Int_t yc[], Float_t rc[], Float_t am[])
			  
{
  double a, dx, dy, d, h, rx, ry;
  double point2_x, point2_y;
  double EPSILON = 20;
  int i,j;
  XCross =-100;
  YCross = -100;
  for (i = 0; i < 3 ; i++)
    {
       for (j = i; j < 4; j++)
         if (am[j] > am[i])
	   {
	     //  cout << i << " " << j << endl;
	     temp = am[i];
	     am[i] = am[j];
	     am[j] = temp;

	     
	     temp = xc[i];
	     xc[i] = xc[j];
	     xc[j] = temp;


	     
	     temp = rc[i];
	     rc[i] = rc[j];
	     rc[j] = temp;
	     
	     	     
	     temp = yc[i];
	     yc[i] = yc[j];
	     yc[j] = temp;
	     /**/
	     
	     
	   }
    }
  //  std::cout << "pippo " << am[0] << " " << am[1] << " " << am[2] << " " << am[3] <<std::endl;
 
  

    /* dx and dy are the vertical and horizontal distances between
    * the circle centers.
    */
    dx = xc[1] - xc[0];
    dy = yc[1] - yc[0];

    /* Determine the straight-line distance between the centers. */
    d = sqrt((dy*dy) + (dx*dx));

    //    std::cout << "distance between centers: " << d << std::endl;
    /* Check for solvability. */
    if (d-EPSILON > (rc[0] + rc[1]))
    {
        /* no solution. circles do not intersect. */
        return;
    }
    if (d+EPSILON < fabs(rc[0] - rc[1]))
    {
        /* no solution. one circle is contained in the other */
        return;
    }

    /* 'point 2' is the point where the line through the circle
    * intersection points crosses the line between the circle
    * centers.
    */

    /* Determine the distance from point 0 to point 2. */
    a = ((rc[0]*rc[0]) - (rc[1]*rc[1]) + (d*d)) / (2.0 * d) ;

    //  std::cout << "point 2  " <<  a << "," << r0 <<  "," << r1 << "," << d << std::endl;
    /* Determine the coordinates of point 2. */
    point2_x = xc[0] + (dx * a/d);
    point2_y = yc[0] + (dy * a/d);
    
    //    std::cout << "point 2  " <<  point2_x << x0 <<  "," << dx*a/d <<  std::endl;
    /* Determine the distance from point 2 to either of the
    * intersection points.
    */
    h = sqrt((rc[0]*rc[0]) - (a*a));

    /* Now determine the offsets of the intersection points from
    * point 2.
    */
    rx = -dy * (h/d);
    ry = dx * (h/d);

    /* Determine the absolute intersection points. */
    double intersectionPoint1_x = point2_x + rx;
    double intersectionPoint2_x = point2_x - rx;
    double intersectionPoint1_y = point2_y + ry;
    double intersectionPoint2_y = point2_y - ry;

    //    std::cout << "INTERSECTION Circle1 AND Circle2: (" <<  intersectionPoint1_x << "," << intersectionPoint1_y<< ") AND ("<< intersectionPoint2_x << "," << intersectionPoint2_y <<")" << std::endl;
    //    XCross =(intersectionPoint2_x+intersectionPoint1_x)/2;
    // YCross = (intersectionPoint2_y+intersectionPoint1_y)/2;      
    
    /* Lets determine if circle 3 intersects at either of the above intersection points. */
    dx = intersectionPoint1_x - xc[2];
    dy = intersectionPoint1_y - yc[2];
    double d1 = sqrt((dy*dy) + (dx*dx));

    dx = intersectionPoint2_x - xc[2];
    dy = intersectionPoint2_y - yc[2];
    double d2 = sqrt((dy*dy) + (dx*dx));

    if(fabs(d1 - rc[2]) < EPSILON) {
      XCross =intersectionPoint1_x;
      YCross = intersectionPoint1_y;
     // std::cout << "INTERSECTION Circle1 AND Circle2 AND Circle3: (" <<  XCross << "," << YCross << ")"<< std::endl;


    }
    else if(fabs(d2 - rc[2]) < EPSILON) {
      XCross =intersectionPoint2_x;
      YCross = intersectionPoint2_y;      
      //  std::cout << "INTERSECTION Circle1 AND Circle2 AND Circle3: (" <<  XCross << "," << YCross << ")"<< std::endl;
    }
    else {
      //        Log.d("INTERSECTION Circle1 AND Circle2 AND Circle3:", "NONE");
    }
    return;
}
