#include <iostream>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include <TROOT.h>
#include <TPad.h>
#include <TMapFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>

typedef unsigned int u_int32_t;
u_int32_t *data=(u_int32_t*)malloc(128);
u_int32_t *buf= (u_int32_t*)malloc(8);
Int_t bufsize = 0;
Int_t data0[9][32];
Int_t data1[9][34];
Int_t data2[9][1];
Int_t cldata[1000][9][34];

Int_t scaler[33],scaler_temp[33],
      scint[33]={1,2,3,4,5,6,7,8,9,10,
		             11,12,13,14,15,16,17,18,19,20,
                 21,22,23,24,25,26,27,28,29,30,
                 31,32,33},
      scal_ch[33] ={ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
		                 11,12,13,14,15,16,17,18,19,20,
		                 21,22,23,24,25,26,27,28,29,30,
		                 31,32,33};
  Int_t si13[48],si15[64],ge5[32],ge11[32],qscint[14],tscint[14],rear[4],scalertt[32];
  const int nrch = 33;
  Int_t DADC1[nrch],DADC2[nrch],DADC3[nrch],DADC4[nrch],DADC5[nrch],DADC6[nrch],DQDC1[nrch],DTDC1[nrch];
  Int_t clevent=0;
  Int_t count = 0;

	Float_t madc_offset[7][34];
	Float_t madc_gain[7][34];
	Float_t madc_factor[7][34];

	TTree *day1data=new TTree("day1data","detector strips");
	TTree *rawdata=new TTree("rawdata","channels");

  TMapFile *mfile;

//   TH2S *hits = new TH2S("HitsOnStrip","Spectrum on each strip",70,0.,70., 8192,0.,8192.);
	double mean[67];
	Int_t hit[67]={0,0,0,0,0};
	Int_t strip[67]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                     20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
                     37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
                     54,55,56,57,58,59,60,61,62,63,64,65,66,67};
	Int_t channel[142]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                     20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
                     37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
                     54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,
		     71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,
		     91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,
		     108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,
		     124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,
		     140,141,142};
//   double Si_13_strip[49],Si_15_strip[65],Ge_5mm_strip[33],Ge_11mm_strip[33];

	Int_t Si_13_strip[49]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                     20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
                     37,38,39,40,41,42,43,44,45,46,47,48,49};
	Int_t Si_15_strip[65]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                     20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,
                     37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,
                     54,55,56,57,58,59,60,61,62,63,64,65};
	Int_t Ge_5mm_strip[33]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                     20,21,22,23,24,25,26,27,28,29,30,31,32,33};
	Int_t Ge_11mm_strip[33]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
                     20,21,22,23,24,25,26,27,28,29,30,31,32,33};

	Int_t nbinsx_madc = 8192,
	      nbinsx_mqdc = 4096,
	      nbinsx_mtdc = 4096,
	      nbinsx_v785 = 4096;
	Int_t xlow_madc = 0, xup_madc = 8192,
	      xlow_mqdc = 0, xup_mqdc = 4096,
	      xlow_mtdc = 0, xup_mtdc = 75536,
	      xlow_v785 = 0, xup_v785 = 4096;

    	Int_t nbinsx_si13 = 5000;
	Int_t nbinsx_si15 = 12000;
	Int_t nbinsx_ge5 = 35000;
	Int_t nbinsx_ge11 = 70000;
	Double_t xlow_si13 = 0., xup_si13 = 5.;
	Double_t xlow_si15 = 0, xup_si15 = 12;
	Double_t xlow_ge5 = 0, xup_ge5 = 35;
	Double_t xlow_ge11 = 0, xup_ge11 = 70;

 	const int m = 64,l = 48, n=32;
	TH1F *h_ADC1[n], *h_ADC2[n], *h_ADC3[n],*h_ADC4[n], *h_ADC5[n],*h_ADC6[n], *h_QDC1[n],*h_TDC1[n];
/*	TH1F *rawh_si1[l], *rawh_si2[m], *rawh_ge1[n],*rawh_ge2[n];
	for(Int_t j=0;j<n;j++)
	{
        h_ADC1[j] = new TH1F(Form("ADC1#%d",j+1),Form("Si1_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);
        h_ADC2[j] = new TH1F(Form("ADC2#%d",j+1),Form("Si1_Strip%02d",j+32),nbinsx_madc, xlow_madc, xup_madc);
        h_ADC3[j] = new TH1F(Form("ADC3#%d",j+1),Form("Si2_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);
        h_ADC4[j] = new TH1F(Form("ADC4#%d",j+1),Form("Si2_Strip%02d",j+32),nbinsx_madc, xlow_madc, xup_madc);
        h_ADC5[j] = new TH1F(Form("ADC5#%d",j+1),Form("Ge1_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);
        h_ADC6[j] = new TH1F(Form("ADC6#%d",j+1),Form("Ge2_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);

	}
*/


   TH1F *timestampADC1 = new TH1F("TimeStampADC1","TS_ADC1",1000,0,8000);
   TH1F *timestampADC2 = new TH1F("TimeStampADC2","TS_ADC2",1000,0,8000);
   TH1F *timestampADC3 = new TH1F("TimeStampADC3","TS_ADC3",1000,0,8000);
   TH1F *timestampADC4 = new TH1F("TimeStampADC4","TS_ADC4",1000,0,8000);
   TH1F *timestampADC5 = new TH1F("TimeStampADC5","TS_ADC5",1000,0,8000);
   TH1F *timestampADC6 = new TH1F("TimeStampADC6","TS_ADC5",1000,0,8000);
   TH1F *timestampQDC1 = new TH1F("TimeStampQDC1","TS_QDC1",1000,0,8000);
   TH1F *timestampTDC1 = new TH1F("TimeStampTDC1","TS_TDC1",1000,0,8000);

   TH1F *h101 = new TH1F("ADC1_ch1","Si_13_Strip1",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h102 = new TH1F("ADC1_ch2","Si_13_Strip2",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h103 = new TH1F("ADC1_ch3","Si_13_Strip3",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h104 = new TH1F("ADC1_ch4","Si_13_Strip4",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h105 = new TH1F("ADC1_ch5","Si_13_Strip5",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h106 = new TH1F("ADC1_ch6","Si_13_Strip6",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h107 = new TH1F("ADC1_ch7","Si_13_Strip7",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h108 = new TH1F("ADC1_ch8","Si_13_Strip8",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h109 = new TH1F("ADC1_ch9","Si_13_Strip9",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h110 = new TH1F("ADC1_ch10","Si_13_Strip10",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h111 = new TH1F("ADC1_ch11","Si_13_Strip11",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h112 = new TH1F("ADC1_ch12","Si_13_Strip12",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h113 = new TH1F("ADC1_ch13","Si_13_Strip13",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h114 = new TH1F("ADC1_ch14","Si_13_Strip14",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h115 = new TH1F("ADC1_ch15","Si_13_Strip15",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h116 = new TH1F("ADC1_ch16","Si_13_Strip16",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h117 = new TH1F("ADC1_ch17","Si_13_Strip17",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h118 = new TH1F("ADC1_ch18","Si_13_Strip18",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h119 = new TH1F("ADC1_ch19","Si_13_Strip19",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h120 = new TH1F("ADC1_ch20","Si_13_Strip20",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h121 = new TH1F("ADC1_ch21","Si_13_Strip21",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h122 = new TH1F("ADC1_ch22","Si_13_Strip22",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h123 = new TH1F("ADC1_ch23","Si_13_Strip23",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h124 = new TH1F("ADC1_ch24","Si_13_Strip24",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h125 = new TH1F("ADC1_ch25","Si_13_Strip25",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h126 = new TH1F("ADC1_ch26","Si_13_Strip26",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h127 = new TH1F("ADC1_ch27","Si_13_Strip27",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h128 = new TH1F("ADC1_ch28","Si_13_Strip28",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h129 = new TH1F("ADC1_ch29","Si_13_Strip29",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h130 = new TH1F("ADC1_ch30","Si_13_Strip30",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h131 = new TH1F("ADC1_ch31","Si_13_Strip31",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h132 = new TH1F("ADC1_ch32","Si_13_Strip32",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h201 = new TH1F("ADC2_ch1","Si_13_Strip33-34",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h202 = new TH1F("ADC2_ch2","Si_13_Strip35-36",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h203 = new TH1F("ADC2_ch3","Si_13_Strip37-38",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h204 = new TH1F("ADC2_ch4","Si_13_Strip39-40",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h205 = new TH1F("ADC2_ch5","Si_13_Strip41-42",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h206 = new TH1F("ADC2_ch6","Si_13_Strip43-44",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h207 = new TH1F("ADC2_ch7","Si_13_Strip45-46",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h208 = new TH1F("ADC2_ch8","Si_13_Strip47-48",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h209 = new TH1F("ADC2_ch9","Si_13_Strip49-50",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h210 = new TH1F("ADC2_ch10","Si_13_Strip51-52",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h211 = new TH1F("ADC2_ch11","Si_13_Strip53-54",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h212 = new TH1F("ADC2_ch12","Si_13_Strip55-56",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h213 = new TH1F("ADC2_ch13","Si_13_Strip57-58",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h214 = new TH1F("ADC2_ch14","Si_13_Strip59-60",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h215 = new TH1F("ADC2_ch15","Si_13_Strip61-62",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h216 = new TH1F("ADC2_ch16","Si_13_Strip63-64",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h217 = new TH1F("ADC2_ch17","Si_13_rear",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h218 = new TH1F("ADC2_ch18","Si_15_rear",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h219 = new TH1F("ADC2_ch19","Ge_5mm_rear",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h220 = new TH1F("ADC2_ch20","Ge_11mm_rear",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h221 = new TH1F("ADC2_ch21","Si_13_bipolar",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h222 = new TH1F("ADC2_ch22","Si_15_bipolar",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h223 = new TH1F("ADC2_ch23","Ge_5mm_bipolar",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h224 = new TH1F("ADC2_ch24","Ge_11mm_bipolar",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h225 = new TH1F("ADC2_ch25","",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h226 = new TH1F("ADC2_ch26","",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h227 = new TH1F("ADC2_ch27","",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h228 = new TH1F("ADC2_ch28","",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h229 = new TH1F("ADC2_ch29","",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h230 = new TH1F("ADC2_ch30","",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h231 = new TH1F("ADC2_ch31","",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h232 = new TH1F("ADC2_ch32","",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h301 = new TH1F("ADC3_ch1","Si_15_Strip1",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h302 = new TH1F("ADC3_ch2","Si_15_Strip2",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h303 = new TH1F("ADC3_ch3","Si_15_Strip3",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h304 = new TH1F("ADC3_ch4","Si_15_Strip4",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h305 = new TH1F("ADC3_ch5","Si_15_Strip5",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h306 = new TH1F("ADC3_ch6","Si_15_Strip6",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h307 = new TH1F("ADC3_ch7","Si_15_Strip7",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h308 = new TH1F("ADC3_ch8","Si_15_Strip8",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h309 = new TH1F("ADC3_ch9","Si_15_Strip9",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h310 = new TH1F("ADC3_ch10","Si_15_Strip10",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h311 = new TH1F("ADC3_ch11","Si_15_Strip11",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h312 = new TH1F("ADC3_ch12","Si_15_Strip12",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h313 = new TH1F("ADC3_ch13","Si_15_Strip13",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h314 = new TH1F("ADC3_ch14","Si_15_Strip14",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h315 = new TH1F("ADC3_ch15","Si_15_Strip15",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h316 = new TH1F("ADC3_ch16","Si_15_Strip16",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h317 = new TH1F("ADC3_ch17","Si_15_Strip17",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h318 = new TH1F("ADC3_ch18","Si_15_Strip18",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h319 = new TH1F("ADC3_ch19","Si_15_Strip19",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h320 = new TH1F("ADC3_ch20","Si_15_Strip20",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h321 = new TH1F("ADC3_ch21","Si_15_Strip21",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h322 = new TH1F("ADC3_ch22","Si_15_Strip22",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h323 = new TH1F("ADC3_ch23","Si_15_Strip23",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h324 = new TH1F("ADC3_ch24","Si_15_Strip24",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h325 = new TH1F("ADC3_ch25","Si_15_Strip25",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h326 = new TH1F("ADC3_ch26","Si_15_Strip26",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h327 = new TH1F("ADC3_ch27","Si_15_Strip27",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h328 = new TH1F("ADC3_ch28","Si_15_Strip28",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h329 = new TH1F("ADC3_ch29","Si_15_Strip29",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h330 = new TH1F("ADC3_ch30","Si_15_Strip30",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h331 = new TH1F("ADC3_ch31","Si_15_Strip31",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h332 = new TH1F("ADC3_ch32","Si_15_Strip32",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h401 = new TH1F("ADC4_ch1","Si_15_Strip33",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h402 = new TH1F("ADC4_ch2","Si_15_Strip34",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h403 = new TH1F("ADC4_ch3","Si_15_Strip35",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h404 = new TH1F("ADC4_ch4","Si_15_Strip36",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h405 = new TH1F("ADC4_ch5","Si_15_Strip37",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h406 = new TH1F("ADC4_ch6","Si_15_Strip38",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h407 = new TH1F("ADC4_ch7","Si_15_Strip39",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h408 = new TH1F("ADC4_ch8","Si_15_Strip40",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h409 = new TH1F("ADC4_ch9","Si_15_Strip41",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h410 = new TH1F("ADC4_ch10","Si_15_Strip42",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h411 = new TH1F("ADC4_ch11","Si_15_Strip43",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h412 = new TH1F("ADC4_ch12","Si_15_Strip44",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h413 = new TH1F("ADC4_ch13","Si_15_Strip45",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h414 = new TH1F("ADC4_ch14","Si_15_Strip46",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h415 = new TH1F("ADC4_ch15","Si_15_Strip47",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h416 = new TH1F("ADC4_ch16","Si_15_Strip48",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h417 = new TH1F("ADC4_ch17","Si_15_Strip49",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h418 = new TH1F("ADC4_ch18","Si_15_Strip50",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h419 = new TH1F("ADC4_ch19","Si_15_Strip51",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h420 = new TH1F("ADC4_ch20","Si_15_Strip52",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h421 = new TH1F("ADC4_ch21","Si_15_Strip53",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h422 = new TH1F("ADC4_ch22","Si_15_Strip54",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h423 = new TH1F("ADC4_ch23","Si_15_Strip55",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h424 = new TH1F("ADC4_ch24","Si_15_Strip56",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h425 = new TH1F("ADC4_ch25","Si_15_Strip57",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h426 = new TH1F("ADC4_ch26","Si_15_Strip58",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h427 = new TH1F("ADC4_ch27","Si_15_Strip59",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h428 = new TH1F("ADC4_ch28","Si_15_Strip60",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h429 = new TH1F("ADC4_ch29","Si_15_Strip51",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h430 = new TH1F("ADC4_ch30","Si_15_Strip62",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h431 = new TH1F("ADC4_ch31","Si_15_Strip63",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h432 = new TH1F("ADC4_ch32","Si_15_Strip64",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h501 = new TH1F("ADC5_ch1","Ge_5mm_Strip1-2",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h502 = new TH1F("ADC5_ch2","Ge_5mm_Strip3-4",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h503 = new TH1F("ADC5_ch3","Ge_5mm_Strip5-6",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h504 = new TH1F("ADC5_ch4","Ge_5mm_Strip7-8",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h505 = new TH1F("ADC5_ch5","Ge_5mm_Strip9-10",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h506 = new TH1F("ADC5_ch6","Ge_5mm_Strip11-13",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h507 = new TH1F("ADC5_ch7","Ge_5mm_Strip14-16",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h508 = new TH1F("ADC5_ch8","Ge_5mm_Strip17-19",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h509 = new TH1F("ADC5_ch9","Ge_5mm_Strip20-22",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h510 = new TH1F("ADC5_ch10","Ge_5mm_Strip23-25",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h511 = new TH1F("ADC5_ch11","Ge_5mm_Strip26-28",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h512 = new TH1F("ADC5_ch12","Ge_5mm_Strip29-31",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h513 = new TH1F("ADC5_ch13","Ge_5mm_Strip31-34",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h514 = new TH1F("ADC5_ch14","Ge_5mm_Strip35-37",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h515 = new TH1F("ADC5_ch15","Ge_5mm_Strip38-40",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h516 = new TH1F("ADC5_ch16","Ge_5mm_Strip41-43",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h517 = new TH1F("ADC5_ch17","Ge_5mm_Strip44-46",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h518 = new TH1F("ADC5_ch18","Ge_5mm_Strip47-48",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h519 = new TH1F("ADC5_ch19","Ge_5mm_Strip49-50",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h520 = new TH1F("ADC5_ch20","Ge_5mm_Strip51-52",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h521 = new TH1F("ADC5_ch21","Ge_5mm_Strip53-54",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h522 = new TH1F("ADC5_ch22","Ge_5mm_Strip55-56",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h523 = new TH1F("ADC5_ch23","Ge_5mm_Strip57-58",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h524 = new TH1F("ADC5_ch24","Ge_5mm_Strip59",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h525 = new TH1F("ADC5_ch25","Ge_5mm_Strip60",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h526 = new TH1F("ADC5_ch26","Ge_5mm_Strip61",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h527 = new TH1F("ADC5_ch27","Ge_5mm_Strip62",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h528 = new TH1F("ADC5_ch28","Ge_5mm_Strip63",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h529 = new TH1F("ADC5_ch29","Ge_5mm_Strip64",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h530 = new TH1F("ADC5_ch30","Ge_5mm_Strip65",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h531 = new TH1F("ADC5_ch31","Ge_5mm_Strip66",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h532 = new TH1F("ADC5_ch32","Ge_5mm_Strip67",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h601 = new TH1F("ADC6_ch1","Ge_11mm_Strip1-2",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h602 = new TH1F("ADC6_ch2","Ge_11mm_Strip3-4",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h603 = new TH1F("ADC6_ch3","Ge_11mm_Strip5-6",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h604 = new TH1F("ADC6_ch4","Ge_11mm_Strip7-8",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h605 = new TH1F("ADC6_ch5","Ge_11mm_Strip9-10",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h606 = new TH1F("ADC6_ch6","Ge_11mm_Strip11-12",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h607 = new TH1F("ADC6_ch7","Ge_11mm_Strip13-14",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h608 = new TH1F("ADC6_ch8","Ge_11mm_Strip15-16",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h609 = new TH1F("ADC6_ch9","Ge_11mm_Strip17-18",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h610 = new TH1F("ADC6_ch10","Ge_11mm_Strip19-20",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h611 = new TH1F("ADC6_ch11","Ge_11mm_Strip21-22",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h612 = new TH1F("ADC6_ch12","Ge_11mm_Strip23-24",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h613 = new TH1F("ADC6_ch13","Ge_11mm_Strip25-26",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h614 = new TH1F("ADC6_ch14","Ge_11mm_Strip27-28",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h615 = new TH1F("ADC6_ch15","Ge_11mm_Strip29-30",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h616 = new TH1F("ADC6_ch16","Ge_11mm_Strip31-32",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h617 = new TH1F("ADC6_ch17","Ge_11mm_Strip33-34",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h618 = new TH1F("ADC6_ch18","Ge_11mm_Strip35-36",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h619 = new TH1F("ADC6_ch19","Ge_11mm_Strip37-38",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h620 = new TH1F("ADC6_ch20","Ge_11mm_Strip39-40",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h621 = new TH1F("ADC6_ch21","Ge_11mm_Strip41-42",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h622 = new TH1F("ADC6_ch22","Ge_11mm_Strip43-44",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h623 = new TH1F("ADC6_ch23","Ge_11mm_Strip45-46",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h624 = new TH1F("ADC6_ch24","Ge_11mm_Strip47-48",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h625 = new TH1F("ADC6_ch25","Ge_11mm_Strip49-50",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h626 = new TH1F("ADC6_ch26","Ge_11mm_Strip51-52",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h627 = new TH1F("ADC6_ch27","Ge_11mm_Strip53-54",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h628 = new TH1F("ADC6_ch28","Ge_11mm_Strip55-56",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h629 = new TH1F("ADC6_ch29","Ge_11mm_Strip57-58",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h630 = new TH1F("ADC6_ch30","Ge_11mm_Strip59-60",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h631 = new TH1F("ADC6_ch31","Ge_11mm_Strip61-63",nbinsx_madc, xlow_madc, xup_madc);
   TH1F *h632 = new TH1F("ADC6_ch32","Ge_11mm_Strip64-67",nbinsx_madc, xlow_madc, xup_madc);

   TH1F *h701 = new TH1F("QDC1_ch1","Q_Big_Up",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h702 = new TH1F("QDC1_ch2","Q_Big_Down",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h703 = new TH1F("QDC1_ch3","Q_Small_Up1",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h704 = new TH1F("QDC1_ch4","Q_Small_Up2",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h705 = new TH1F("QDC1_ch5","Q_Small_Down1",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h706 = new TH1F("QDC1_ch6","Q_Small_Down2",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h707 = new TH1F("QDC1_ch7","Q_Small_Left1",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h708 = new TH1F("QDC1_ch8","Q_Small_Left2",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h709 = new TH1F("QDC1_ch9","Q_Small_Right1",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h710 = new TH1F("QDC1_ch10","Q_Small_Right2",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h711 = new TH1F("QDC1_ch11","Q_Big_Left1",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h712 = new TH1F("QDC1_ch12","Q_Big_Left2",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h713 = new TH1F("QDC1_ch13","Q_Big_Right1",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h714 = new TH1F("QDC1_ch14","Q_Big_Right2",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h715 = new TH1F("QDC1_ch15","h715",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h716 = new TH1F("QDC1_ch16","h716",nbinsx_mqdc, xlow_mqdc, xup_mqdc);

   TH1F *h717 = new TH1F("QDC1_ch17","h717",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h718 = new TH1F("QDC1_ch18","h718",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h719 = new TH1F("QDC1_ch19","h719",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h720 = new TH1F("QDC1_ch20","h720",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h721 = new TH1F("QDC1_ch21","h721",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h722 = new TH1F("QDC1_ch22","h722",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h723 = new TH1F("QDC1_ch23","h723",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h724 = new TH1F("QDC1_ch24","h724",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h725 = new TH1F("QDC1_ch25","h725",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h726 = new TH1F("QDC1_ch26","h726",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h727 = new TH1F("QDC1_ch27","h727",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h728 = new TH1F("QDC1_ch28","h728",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h729 = new TH1F("QDC1_ch29","h729",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h730 = new TH1F("QDC1_ch30","h730",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h731 = new TH1F("QDC1_ch31","h731",nbinsx_mqdc, xlow_mqdc, xup_mqdc);
   TH1F *h732 = new TH1F("QDC1_ch32","h732",nbinsx_mqdc, xlow_mqdc, xup_mqdc);

   TH1F *h801 = new TH1F("TDC1_ch1","T_Big_Up",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h802 = new TH1F("TDC1_ch2","T_Big_Down",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h803 = new TH1F("TDC1_ch3","T_Small_Up1",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h804 = new TH1F("TDC1_ch4","T_Small_Up2",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h805 = new TH1F("TDC1_ch5","T_Small_Down1",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h806 = new TH1F("TDC1_ch6","T_Small_Down2",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h807 = new TH1F("TDC1_ch7","T_Small_Left1",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h808 = new TH1F("TDC1_ch8","T_Small_Left2",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h809 = new TH1F("TDC1_ch9","T_Small_Right1",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h810 = new TH1F("TDC1_ch10","T_Small_Right2",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h811 = new TH1F("TDC1_ch11","T_Big_Left1",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h812 = new TH1F("TDC1_ch12","T_Big_Left2",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h813 = new TH1F("TDC1_ch13","T_Big_Right1",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h814 = new TH1F("TDC1_ch14","T_Big_Right2",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h815 = new TH1F("TDC1_ch15","h815",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h816 = new TH1F("TDC1_ch16","h816",nbinsx_mtdc, xlow_mtdc, xup_mtdc);

   TH1F *h817 = new TH1F("TDC1_ch17","h817",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h818 = new TH1F("TDC1_ch18","h818",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h819 = new TH1F("TDC1_ch19","h819",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h820 = new TH1F("TDC1_ch20","h820",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h821 = new TH1F("TDC1_ch21","h821",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h822 = new TH1F("TDC1_ch22","h822",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h823 = new TH1F("TDC1_ch23","h823",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h824 = new TH1F("TDC1_ch24","h824",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h825 = new TH1F("TDC1_ch25","h825",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h826 = new TH1F("TDC1_ch26","h826",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h827 = new TH1F("TDC1_ch27","h827",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h828 = new TH1F("TDC1_ch28","h828",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h829 = new TH1F("TDC1_ch29","h829",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h830 = new TH1F("TDC1_ch30","h830",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h831 = new TH1F("TDC1_ch31","h831",nbinsx_mtdc, xlow_mtdc, xup_mtdc);
   TH1F *h832 = new TH1F("TDC1_ch32","h832",nbinsx_mtdc, xlow_mtdc, xup_mtdc);

	TH2S *scaler_hist = new TH2S("Scaler","Scaler_hist",33,0,33,1000,0,10000000);
	TH2S *Si_13_hits = new TH2S("Si_13_strip","Spectrum on Si_13",49,0.,49., 1024,0.,8192.);
	TH2S *Si_15_hits = new TH2S("Si_15_Strip","Spectrum on Si_15",65,0.,65., 1024,0.,8192.);
	TH2S *Ge_5mm_hits = new TH2S("Ge_5mm_Strip","Spectrum on Ge_5mm",33,0.,33., 1024,0.,8192.);
	TH2S *Ge_11mm_hits = new TH2S("Ge_11mm_Strip","Spectrum on Ge_11mm",33,0.,33., 1024,0.,8192.);

	TH2S *Scint_qhits = new TH2S("ScintillatorsQ","Scint_QDC",33,0,33,1024,0,4096);
	TH2S *Scint_thits = new TH2S("ScintillatorsT","Scint_TDC",33,0,33,4096,0,75536);
