/*
* energy_spectrum.C
* created on Jul 3, 2013
* Author: Huagen
*
* modified on Dec.3, 2013 to involve scaler data
*
* implementating tree data on Jul.20, 2013 by Huagen
*
*/

#include <iostream>
#include <fstream>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <TPad.h>
#include <TMapFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TStopwatch.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TRandom.h"
#include "TLegend.h"

#include "madc_decoding.h"
//#include "histogram_init.h"

#include "Riostream.h"
#include <iomanip>   // format manipulation


using namespace std;

char fname[80];
//TFile *file = new TFile("test.root","RECREATE");

void hist_filling(){

    for(Int_t j=0;j<32;j++)
    {
        h_ADC1[j] = new TH1F(Form("MADC1_%d",j+1),Form("Si1_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);}
 for(Int_t j=0;j<32;j++)
    {
        h_ADC2[j] = new TH1F(Form("MADC2_%d",j+1),Form("Si1_Strip%02d",j+32),nbinsx_madc, xlow_madc, xup_madc);}
 for(Int_t j=0;j<32;j++)
    {
        h_ADC3[j] = new TH1F(Form("MADC3_%d",j+1),Form("Si2_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);}
 for(Int_t j=0;j<32;j++)
    {
        h_ADC4[j] = new TH1F(Form("MADC4_%d",j+1),Form("Si2_Strip%02d",j+32),nbinsx_madc, xlow_madc, xup_madc);}
 for(Int_t j=0;j<32;j++)
    {
        h_ADC5[j] = new TH1F(Form("MADC5_%d",j+1),Form("Ge1_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);}
 for(Int_t j=0;j<32;j++)
    {
        h_ADC6[j] = new TH1F(Form("MADC6_%d",j+1),Form("Ge2_Strip%02d",j+1),nbinsx_madc, xlow_madc, xup_madc);}
 for(Int_t j=0;j<32;j++)
    {
        h_QDC1[j] = new TH1F(Form("MQDC1_%d",j+1),Form("Forward%02d",j+1),nbinsx_mqdc, xlow_mqdc, xup_mqdc);}
 for(Int_t j=0;j<32;j++)
    {
        h_TDC1[j] = new TH1F(Form("MTDC1_%d",j+1),Form("Time%02d",j+1),nbinsx_mtdc, xlow_mtdc, xup_mtdc);
    }

}

void hist_writing(){
    for(Int_t j=0;j<32;j++){
        h_ADC1[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC2[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC3[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC4[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC5[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC6[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_QDC1[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_TDC1[j]->Write();}
}


/*****************************************************************************************************/
//function of decoding()
/*****************************************************************************************************/

int decoding(int f)
{
    //   printf("In decoding(): Start decoding!\n");

	Int_t res=0;
	Int_t size=0; //number of all following words

	bufsize = 2;
	const int endiantest = 0x12345678;
	//const int endiantest Source_72K_IRQ_simple_TrigHand_2015_Sep_22_16_50_08.cl= 0x78563412;

	if((res=read(f,buf,4*2))!=4*2)
	{	printf("buffer[0]=%d, buffer[1]=%d\n",buf[0],buf[1]);
		printf("In decoding(): read size: res=%d, error=%s\n", res, strerror(errno));
		return 1;
	}

	if (buf[1]==endiantest)
	{
	   size=buf[0]; cout<<"The size of cluster is "<<size<<endl;
	   if(size<8) return 1;

     clustersize->Fill(size);
	    // if(size>=8) printf("the size of cluster is %d \n",size);
	} else 	{
	   printf("In decoding(): cluster endian (buffer[1]) is not 0x12345678!\n");
	}


		if (size+2>bufsize)
		{
//		printf("size+2>bufsize\n");
          buf=(u_int32_t*)realloc(buf,(size+2)*4);
		      bufsize=size+2;
		}
  //   printf("after reallocate, bufsize is %d \n", bufsize);

	 read(f,buf+2,(size-1)*4); //buf+2 means to skip the cluster size and endian number
//	 printf("read(f,buf+2,(size-1)*4)\n");

/*
	if (buf[2]!=0)
	{
		printf("buf[2]!=0 \n");
		return 0;
	} else {
		//	printf("The buf[2] ==0\n");
	}
*/

//loop one cluster data with datasize of "size"

	    u_int32_t time0_ADC1=0,
		      time0_ADC2=0,
		      time0_ADC3=0,
		      time0_ADC4=0,
		      time0_ADC5=0,
		      time0_ADC6=0,
		      time0_QDC1=0,
		      time0_TDC1=0;

	   u_int32_t time_evt_ADC1=0,
		      time_evt_ADC2=0,
		      time_evt_ADC3=0,
		      time_evt_ADC4=0,
		      time_evt_ADC5=0,
		      time_evt_ADC6=0,
		      time_evt_QDC1=0,
	        time_evt_TDC1=0;

	int temp_ID=0, scaler_ID=0;

	int evtADC1=0, evtADC2=0,evtADC3=0, evtADC4=0, evtADC5=0, evtADC6=0, evtQDC1=0, evtTDC1=0;
  int evtendADC1=0, evtendADC2=0,evtendADC3=0, evtendADC4=0, evtendADC5=0, evtendADC6=0, evtendQDC1=0, evtendTDC1=0;

	u_int32_t end0=0,end1,end2,end3,end4,end5,end6;
  int dt1=0,dt2=0,dt3=0,dt4=0,dt5=0,dt6=0,dt7=0,dt8=0;

  //loop the cluster data with size
	for (int n=0;n<size+1;n++)
	{
		for(int k=0;k<30;k++)
               {
                //	printf(" Buf [%d] is 0x%08x \n", k, buf[k]);
	       }


/**********************************************************************
* MXDC data
***********************************************************************/

      if((buf[n]&0xFFF00000) == 0x40100000 ||(buf[n]&0xFFF00000) == 0x40200000 || (buf[n]&0xFFF00000) == 0x40300000)
//0x40000021 => 0x40000022 on Aug.14.2015
	   {
			//printf("The buf[%d]= 0x%08x is event header\n ",n,buf[n] );
                          int adcres = buf[n]>>12 & 0x7;		//printf("the ADC resolution is %d \n",adcres);
			  int nrwords = buf[n]&0xfff;		//printf("the following words are %d \n",nrwords);
                //	int id = (buf[n]>>16)&0xff;		printf("the ADC id is %d \n",id);
			  int id = ((buf[n]>>16)&0xff);	//	printf("module_id=%d \n",id);
			  temp_ID = id;
			//int evtId=0;

      if(id>0){
        if(id==16){evtADC1++; eventheader1->Fill(10); }
        if(id==17){evtADC2++; eventheader2->Fill(20);}
        if(id==18){evtADC3++; eventheader3->Fill(30);}
        if(id==19){evtADC4++; eventheader4->Fill(40);}
        if(id==20){evtADC5++; eventheader5->Fill(50);}
        if(id==21){evtADC6++; eventheader6->Fill(60);}
        }


        if((buf[n+nrwords]&0xC0000000)==0xC0000000){
        //  printf("The buf[%d]= 0x%08x is event end\n", n+nrwords, buf[n+nrwords]);
        if(id==16){evtendADC1++; eventend1->Fill(10);}
        if(id==17){evtendADC2++; eventend2->Fill(20);}
        if(id==18){evtendADC3++; eventend3->Fill(30);}
        if(id==19){evtendADC4++; eventend4->Fill(40);}
        if(id==20){evtendADC5++; eventend5->Fill(50);}
        if(id==21){evtendADC6++; eventend6->Fill(60);}

        } else {cout<<"no proper event end for the data"<<endl;}

		if(id<22){//Identify ADC data by module ID
		 id = id - 15;

        for(int i=1;i<=nrwords;i++) {
		        if(id==1){
                      if((buf[n+i]&0xf4E00000)==0x04000000) {
				                     int ch= (buf[n+i]>>16) & 0x1F;
            			           data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                     data1[id][ch] = data0[id][ch];

				                     cldata[evtADC1][id][ch]=data0[id][ch];
				//cout<<"ADC1["<<id<<"]["<<ch<<"]= "<<data1[1][ch]<<endl;

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                  //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
				                     end0 = buf[n+i] &0x3FFFFFFF;
				                  //   dt1 = end0 - time_evt_ADC1;
				                     time_evt_ADC1=end0;
                             dt1 = end0 - time_evt_ADC1;
				cout<<"time stamp in the EOE of ADC1=  "<<end0<<endl;
	//			 cout<<"the ADC1 dt1= "<<dt1<<endl;
			//	cout<<"id= "<<id<<" dt= "<<dt1*62.5/1000000<<" ms"<<endl;
				// cout<<"The cluster count is            "<<count<<endl;

				                     timestampADC1->Fill(abs(dt1)+100);
				                     data2[id][0] = end0;

				                     cldata[evtADC1][id][32]=end0; //timestamp of 1 event
				                  //   evtADC1++;
			//	cout<<"evtADC1="<<evtADC1<<endl;d
                      } //event end
                      } //ADC1 data

		   else if(id==2){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
            			          data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                    data1[id][ch] = data0[id][ch];

				                    cldata[evtADC2][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                  //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
				                    end0 = buf[n+i] &0x3FFFFFFF;
				                    dt2 = end0 - time_evt_ADC1;
				                    time_evt_ADC2=end0;
	//			cout<<"time stamp in the EOE of ADC2=  "<<end0<<endl;
				 cout<<"the ADC2 dt2= "<<dt2<<endl;
			//	 cout<<"id= "<<id<<" dt= "<<dt2*62.5/1000000<<" ms"<<endl;
			//	 cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC2->Fill(abs(dt2));
				                    data2[id][0] = end0;

				                    cldata[evtADC2][id][32]=end0; //timestamp of 1 event
				                  //  evtADC2++;
			//	cout<<"evtADC2="<<evtADC2<<endl;
                      } //event end
			                } // ADC2 data
		    else if(id==3){
                      if((buf[n+i]&0xf4E00000)==0x04000000) {
				                     int ch= (buf[n+i]>>16) & 0x1F;
            			           data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                     data1[id][ch] = data0[id][ch];

				                     cldata[evtADC3][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
				                     end0 = buf[n+i] &0x3FFFFFFF;
				                     dt3 = end0 - time_evt_ADC1;
				                     time_evt_ADC3=end0;
	//			    cout<<"time stamp in the EOE of ADC3=  "<<end0<<endl;
				 cout<<"the ADC3 dt3= "<<dt3<<endl;
		//		 cout<<"id= "<<id<<" dt= "<<dt3*62.5/1000000<<" ms"<<endl;
				// cout<<"The cluster count is            "<<count<<endl;

				                     timestampADC3->Fill(abs(dt3));
				                     data2[id][0] = end0;
				                     cldata[evtADC3][id][32]=end0; //timestamp of 1 event
				                   //  evtADC3++;
		//		cout<<"evtADC3="<<evtADC3<<endl;

                      } //event end
                    } //ADC3 data
		     else if(id==4){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
            			          data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 			                      data1[id][ch] = data0[id][ch];

				                    cldata[evtADC4][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
				                    end0 = buf[n+i] &0x3FFFFFFF;
				                    dt4 = end0 - time_evt_ADC1;
				                    time_evt_ADC4=end0;
	//			cout<<"time stamp in the EOE of ADC4=  "<<end0<<endl;
				cout<<"the ADC4 dt4= "<<dt4<<endl;
				///cout<<"id= "<<id<<" dt= "<<dt4*62.5/1000000<<" ms"<<endl;
				//cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC4->Fill(abs(dt4));
				                    data2[id][0] = end0;
				                    cldata[evtADC4][id][32]=end0; //timestamp of 1 event
				               //     evtADC4++;

			                } //event end
			              } //ADC4 data
		     else if(id==5){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
            			          data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                    data1[id][ch] = data0[id][ch];

				                    cldata[evtADC5][id][ch]=data0[id][ch];
               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
				                    end0 = buf[n+i] &0x3FFFFFFF;
				                    dt5 = end0 - time_evt_ADC1;
				                    time_evt_ADC5=end0;
	//			cout<<"time stamp in the EOE of ADC5=  "<<end0<<endl;
				cout<<"the ADC5 dt5= "<<dt5<<endl;
				//cout<<"id= "<<id<<" dt= "<<dt5*62.5/1000000<<" ms"<<endl;
				//cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC5->Fill(abs(dt5));
				                    data2[id][0] = end0;
				                    cldata[evtADC5][id][32]=end0; //timestamp of 1 event
				             //       evtADC5++;

  			                } //event end
			                } //ADC5 data
		      else if(id==6){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
                      			data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                    data1[id][ch] = data0[id][ch];

				                    cldata[evtADC6][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
				                    end0 = buf[n+i] &0x3FFFFFFF;
				                    dt6 = end0 - time_evt_ADC1;
				                    time_evt_ADC6=end0;
	//			cout<<"time stamp in the EOE of ADC6=  "<<end0<<endl;
				cout<<"the ADC6 dt6= "<<dt6<<endl;
				//cout<<"id= "<<id<<" dt= "<<dt6*62.5/1000000<<" ms"<<endl;
				//cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC6->Fill(abs(dt6));
				                    data2[id][0] = end0;

				                    cldata[evtADC6][id][32]=end0; //timestamp of 1 event
				              //      evtADC6++;

			               } //event end
			             } else continue;  //module ID

		        } //loop one ADC data (nrwords)

            n += nrwords;

   	    }else if(id>31&&id<33){//Identify QDC data by module ID
		        id = id - 25;

                for(int i=1;i<=nrwords;i++){
			               if(id==7){
			                         if((buf[n+i]&0xf4E00000)==0x04000000){
			                              int ch= (buf[n+i]>>16) & 0x1F;      //cout<<" ch1 = "<<ch<<endl;
			                              data1[id][ch] = (buf[n+i]) & 0xFFF;    // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;

			                              cldata[evtQDC1][id][ch]=data1[id][ch];

			                          } else if((buf[n+i]&0xC0000000)==0xC0000000){
				                            end0 = buf[n+i] &0x3FFFFFFF;
				                            dt7 = end0 - time_evt_QDC1;
				                            time_evt_QDC1=end0;
                                    timestampQDC1->Fill(dt7*62.5/1000000);
				                            data2[id][0] = end0;

				                            cout<<"id= "<<id<<" dt= "<<dt7*62.5/1000000<<" ms"<<endl;

				                            cldata[evtQDC1][id][32]=end0; //timestamp of 1 event
				                            evtQDC1++;

				                        }//event end

				                       }else continue;
		             } //loop one QDC data

	            n += nrwords;

	      }else if(id>46&&id<50){
			       id = id-40;
                for(int i=1;i<=nrwords;i++)	{
			               if(id==8){
			                    if((buf[n+i]&0xf4C00000)==0x04000000){
			                         int ch= (buf[n+i]>>16) & 0x1F;  //    cout<<" ch2 = "<<ch<<endl;
			                         data1[id][ch] = (buf[n+i]) & 0xFFFF;  //   cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
		//	printf("data[%d][ %d ] is %d \n",id, ch, data1[id][ch]);
		//	printf("data[%d][ %d ] is 0x%08x \n",id, ch, buf[n+i]);
		//	printf("TDC data: the raw data No %d and data are  0x%08x \n",n+i,buf[n+i]);

			                         cldata[evtTDC1][id][ch]=data1[id][ch];

			                     } else if((buf[n+i]&0xC0000000)==0xC0000000){
			//cout<<"The data of buf["<<n+i<<"] is event end"<<endl;

				                       end0 = buf[n+i] &0x3FFFFFFF;
				                       dt8 = end0 - time_evt_TDC1;
				                       time_evt_TDC1=end0;

				//cout<<"time stamp in the EOE of TDC is "<<end0<<endl;
			//	 cout<<"the dt is                       "<<dt<<endl;
				                      cout<<"id="<<id<<" dt="<<dt8*62.5/1000000<<" ms"<<endl;
				                       timestampTDC1->Fill(dt8*62.5/1000000);
				                       data2[id][0] = end0;

				                       cldata[evtTDC1][id][32]=end0; //timestamp of 1 event
				                       evtTDC1++;
				                   } //event end
				               } else continue;
		             } //loop one ADC data
	         n += nrwords;
		   } else continue;

   	 }  //ADC QDC TDC data identified by module ID
	else if((buf[n]&0xFFF00000) == 0x40000000)
	{
	//		printf("The buf[%d]= 0x%08x is event header\n ",n,buf[n] );
       			  int adcres = buf[n]>>12 & 0x7;		//printf("the ADC resolution is %d \n",adcres);
			  int nrwords = buf[n]&0xfff;		//printf("the following words are %d \n",nrwords);
                //	int id = (buf[n]>>16)&0xff;		printf("the ADC id is %d \n",id);
			  int id = ((buf[n]>>16)&0xff);	//	printf("module_id=%d \n",id);
			  temp_ID = id;
			//int evtId=0;

  	    if(id>0){
          	      if(id==1){evtADC1++; eventheader1->Fill(10); }
     		          if(id==2){evtADC2++; eventheader2->Fill(20);}
                  if(id==3){evtADC3++; eventheader3->Fill(30);}
                  if(id==4){evtADC4++; eventheader4->Fill(40);}
                  if(id==5){evtADC5++; eventheader5->Fill(50);}
                  if(id==6){evtADC6++; eventheader6->Fill(60);}
                 }
            if((buf[n+nrwords]&0xC0000000)==0xC0000000){
      //    printf("The buf[%d]= 0x%08x is event end\n", n+nrwords, buf[n+nrwords]);

                  if(id==1){evtendADC1++; eventend1->Fill(10);}
                  if(id==2){evtendADC2++; eventend2->Fill(20);}
                  if(id==3){evtendADC3++; eventend3->Fill(30);}
                  if(id==4){evtendADC4++; eventend4->Fill(40);}
                  if(id==5){evtendADC5++; eventend5->Fill(50);}
                  if(id==6){evtendADC6++; eventend6->Fill(60);}
              } else {cout<<"no proper event end for the data"<<endl;}

		if(id<7){//Identify ADC data by module ID
			// id = id - 15;

                   for(int i=1;i<=nrwords;i++) {
		        if(id==1){
                      if((buf[n+i]&0xf4E00000)==0x04000000) {
				                     int ch= (buf[n+i]>>16) & 0x1F;
            			           data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                     data1[id][ch] = data0[id][ch];

				                     cldata[evtADC1][id][ch]=data0[id][ch];
				//cout<<"ADC1["<<id<<"]["<<ch<<"]= "<<data1[1][ch]<<endl;

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                //    cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
                             end1 = 0, dt1=0;
				                     end1 = buf[n+i] &0x3FFFFFFF;
				                  //   dt1 = end0 - time_evt_ADC1;
                      //    cout<<"1 time_evt_ADC1 = "<<time_evt_ADC1<<endl;
				                     time_evt_ADC1=end1;
                             dt1 = end1 - time_evt_ADC1;
                    //         cout<<"2 time_evt_ADC1 = "<<time_evt_ADC1<<endl;
		//		cout<<"time stamp in the EOE of ADC1=  "<<end0<<endl;
	//			 cout<<"the ADC1 dt1= "<<dt1<<endl;
			//	cout<<"id= "<<id<<" dt= "<<dt1*62.5/1000000<<" ms"<<endl;
				// cout<<"The cluster count is            "<<count<<endl;

				                     timestampADC1->Fill(abs(dt1)+100);
				                     data2[id][0] = end1;

				                     cldata[evtADC1][id][32]=end1; //timestamp of 1 event
				                  //   evtADC1++;
			//	cout<<"evtADC1="<<evtADC1<<endl;d
                      } //event end
                      } //ADC1 data

		   else if(id==2){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
            			          data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                    data1[id][ch] = data0[id][ch];

				                    cldata[evtADC2][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                  //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
                            end2 = 0, dt2=0;
				                    end2 = buf[n+i] &0x3FFFFFFF;
				                    dt2 = end2 - time_evt_ADC1;
				                    time_evt_ADC2=end2;
		//		cout<<"time stamp in the EOE of ADC2=  "<<end0<<endl;
			//	 cout<<"the ADC2 dt2= "<<dt2<<endl;
			//	 cout<<"id= "<<id<<" dt= "<<dt2*62.5/1000000<<" ms"<<endl;
			//	 cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC2->Fill(abs(dt2));
				                    data2[id][0] = end2;

				                    cldata[evtADC2][id][32]=end2; //timestamp of 1 event
				                  //  evtADC2++;
			//	cout<<"evtADC2="<<evtADC2<<endl;
                      } //event end
			                } // ADC2 data
		    else if(id==3){
                      if((buf[n+i]&0xf4E00000)==0x04000000) {
				                     int ch= (buf[n+i]>>16) & 0x1F;
            			           data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                     data1[id][ch] = data0[id][ch];

				                     cldata[evtADC3][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
                              end3 = 0, dt3=0;
				                     end3 = buf[n+i] &0x3FFFFFFF;
				                     dt3 = end3 - time_evt_ADC1;
				                     time_evt_ADC3=end3;
			//	    cout<<"time stamp in the EOE of ADC3=  "<<end0<<endl;
			//	 cout<<"the ADC3 dt3= "<<dt3<<endl;
		//		 cout<<"id= "<<id<<" dt= "<<dt3*62.5/1000000<<" ms"<<endl;
				// cout<<"The cluster count is            "<<count<<endl;

				                     timestampADC3->Fill(abs(dt3));
				                     data2[id][0] = end3;
				                     cldata[evtADC3][id][32]=end3; //timestamp of 1 event
				                   //  evtADC3++;
		//		cout<<"evtADC3="<<evtADC3<<endl;

                      } //event end
                    } //ADC3 data
		     else if(id==4){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
            			          data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 			                      data1[id][ch] = data0[id][ch];

				                    cldata[evtADC4][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
                            end4 = 0, dt4=0;
				                    end4 = buf[n+i] &0x3FFFFFFF;
				                    dt4 = end4 - time_evt_ADC1;
				                    time_evt_ADC4=end4;
		//		cout<<"time stamp in the EOE of ADC4=  "<<end0<<endl;
			//	cout<<"the ADC4 dt4= "<<dt4<<endl;
				///cout<<"id= "<<id<<" dt= "<<dt4*62.5/1000000<<" ms"<<endl;
				//cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC4->Fill(abs(dt4));
				                    data2[id][0] = end4;
				                    cldata[evtADC4][id][32]=end4; //timestamp of 1 event
				               //     evtADC4++;

			                } //event end
			              } //ADC4 data
		     else if(id==5){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
            			          data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                    data1[id][ch] = data0[id][ch];

				                    cldata[evtADC5][id][ch]=data0[id][ch];
               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
                            end5 = 0, dt5=0;
				                    end5 = buf[n+i] &0x3FFFFFFF;
				                    dt5 = end5 - time_evt_ADC1;
				                    time_evt_ADC5=end5;
		//		cout<<"time stamp in the EOE of ADC5=  "<<end0<<endl;
			//	cout<<"the ADC5 dt5= "<<dt5<<endl;
				//cout<<"id= "<<id<<" dt= "<<dt5*62.5/1000000<<" ms"<<endl;
				//cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC5->Fill(abs(dt5));
				                    data2[id][0] = end5;
				                    cldata[evtADC5][id][32]=end5; //timestamp of 1 event
				             //       evtADC5++;

  			                } //event end
			                } //ADC5 data
		      else if(id==6){
                      if((buf[n+i]&0xf4E00000)==0x04000000){
				                    int ch= (buf[n+i]>>16) & 0x1F;
                      			data0[id][ch] = (buf[n+i]) & 0x1FFF;   // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
 				                    data1[id][ch] = data0[id][ch];

				                    cldata[evtADC6][id][ch]=data0[id][ch];

               	      }else if((buf[n+i]&0xC0000000)==0xC0000000){
                 //  cout<<"The data of buf["<<n+i<<"] is event end"<<endl;
                            end6 = 0, dt6=0;
				                    end6 = buf[n+i] &0x3FFFFFFF;
				                    dt6 = end6 - time_evt_ADC1;
				                    time_evt_ADC6=end6;
				//cout<<"time stamp in the EOE of ADC6=  "<<end0<<endl;
				//cout<<"the ADC6 dt6= "<<dt6<<endl;
				//cout<<"id= "<<id<<" dt= "<<dt6*62.5/1000000<<" ms"<<endl;
				//cout<<"The cluster count is            "<<count<<endl;
				                    timestampADC6->Fill(abs(dt6));
				                    data2[id][0] = end6;

				                    cldata[evtADC6][id][32]=end6; //timestamp of 1 event
				              //      evtADC6++;

			               } //event end



			             } else continue;  //module ID


		        } //loop one ADC data (nrwords)

            //check the dt is abnormal to the expected;
       /*     if(abs(abs(dt2)-691)>5||abs(abs(dt3)-1555)>5||abs(abs(dt4)-2575)>5||abs(abs(dt5)-3253)>5||abs(abs(dt6)-4180)>5){
              cout<<"The time stamp is abnormal"<<endl;
              cout<<"the dt1= "<<dt1<<" TS1= "<<data2[1][0]<<" "<<time_evt_ADC1<<endl
                  <<"the dt2= "<<dt2<<" TS2= "<<data2[2][0]<<" "<<time_evt_ADC2<<endl
                  <<"the dt3= "<<dt3<<" TS3= "<<data2[3][0]<<" "<<time_evt_ADC3<<endl
                  <<"the dt4= "<<dt4<<" TS4= "<<data2[4][0]<<" "<<time_evt_ADC4<<endl
                  <<"the dt5= "<<dt5<<" TS5= "<<data2[5][0]<<" "<<time_evt_ADC5<<endl
                  <<"the dt6= "<<dt6<<" TS6= "<<data2[6][0]<<" "<<time_evt_ADC6<<endl;
            }
            */ //display the time stamp difference

            if(abs(dt2)>5000||abs(dt3)>5000||abs(dt4)>5000||abs(dt5)>5000||abs(dt6)>5000){
                   cout<<"The time stamp is abnormal"<<endl;
                   cout<<"the dt1= "<<dt1<<" TS1= "<<data2[1][0]<<" "<<time_evt_ADC1<<endl
                       <<"the dt2= "<<dt2<<" TS2= "<<data2[2][0]<<" "<<time_evt_ADC2<<endl
                       <<"the dt3= "<<dt3<<" TS3= "<<data2[3][0]<<" "<<time_evt_ADC3<<endl
                       <<"the dt4= "<<dt4<<" TS4= "<<data2[4][0]<<" "<<time_evt_ADC4<<endl
                       <<"the dt5= "<<dt5<<" TS5= "<<data2[5][0]<<" "<<time_evt_ADC5<<endl
                       <<"the dt6= "<<dt6<<" TS6= "<<data2[6][0]<<" "<<time_evt_ADC6<<endl;
                 }

            n += nrwords;

   	    }else if(id==7){//Identify QDC data by module ID
		 //       id = id - 25;

                for(int i=1;i<=nrwords;i++){
			               if(id==7){
			                         if((buf[n+i]&0xf4E00000)==0x04000000){
			                              int ch= (buf[n+i]>>16) & 0x1F;      //cout<<" ch1 = "<<ch<<endl;
			                              data1[id][ch] = (buf[n+i]) & 0xFFF;    // cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;

			                              cldata[evtQDC1][id][ch]=data1[id][ch];

			                          } else if((buf[n+i]&0xC0000000)==0xC0000000){
				                            end0 = buf[n+i] &0x3FFFFFFF;
				                            dt7 = end0 - time_evt_QDC1;
				                            time_evt_QDC1=end0;
                                    timestampQDC1->Fill(dt7*62.5/1000000);
				                            data2[id][0] = end0;

				                            cout<<"id= "<<id<<" dt= "<<dt7*62.5/1000000<<" ms"<<endl;

				                            cldata[evtQDC1][id][32]=end0; //timestamp of 1 event
				                            evtQDC1++;

				                        }//event end

				                       }else continue;
		             } //loop one QDC data

	            n += nrwords;

	      }else if(id==8){
		//	       id = id-40;
                for(int i=1;i<=nrwords;i++)	{
			               if(id==8){
			                    if((buf[n+i]&0xf4C00000)==0x04000000){
			                         int ch= (buf[n+i]>>16) & 0x1F;  //    cout<<" ch2 = "<<ch<<endl;
			                         data1[id][ch] = (buf[n+i]) & 0xFFFF;  //   cout<<" buf["<<n+i<<"] is "<<buf[ch]<<" "<<endl;
		//	printf("data[%d][ %d ] is %d \n",id, ch, data1[id][ch]);
		//	printf("data[%d][ %d ] is 0x%08x \n",id, ch, buf[n+i]);
		//	printf("TDC data: the raw data No %d and data are  0x%08x \n",n+i,buf[n+i]);

			                         cldata[evtTDC1][id][ch]=data1[id][ch];

			                     } else if((buf[n+i]&0xC0000000)==0xC0000000){
			//cout<<"The data of buf["<<n+i<<"] is event end"<<endl;

				                       end0 = buf[n+i] &0x3FFFFFFF;
				                       dt8 = end0 - time_evt_TDC1;
				                       time_evt_TDC1=end0;

				//cout<<"time stamp in the EOE of TDC is "<<end0<<endl;
			//	 cout<<"the dt is                       "<<dt<<endl;
				                      cout<<"id="<<id<<" dt="<<dt8*62.5/1000000<<" ms"<<endl;
				                       timestampTDC1->Fill(dt8*62.5/1000000);
				                       data2[id][0] = end0;

				                       cldata[evtTDC1][id][32]=end0; //timestamp of 1 event
				                       evtTDC1++;
				                   } //event end
				               } else continue;
		             } //loop one ADC data
	         n += nrwords;
		   } else continue;

	} //MXDC data with HW module ID


/**********************************************************************
* Scaler data
***********************************************************************/
		if((buf[n]==1)&&(buf[n+1]==32))
			{
			//	cout<<"**********************************************************************"<<endl;
			//	cout<<"Scaler events!"<<endl;
			//	cout<<"**********************************************************************"<<endl;

				int scaler_length = buf[n+1]&0xFFFFFFFF;
				for(int k=0;k<scaler_length;k++)
					{scaler_temp[k]=buf[n+2+k]&0xFFFFFFFF;
					 scaler[k] = scaler_temp[k];
					//	if(k<6)cout<<"The scaler_temp["<<k<<"] = "<<scaler_temp[k]<<endl;
				//cout<<"**********************************************************************"<<endl;
					}
				n += scaler_length;
				scaler_ID++;

			}//loop scaler data which means one event complete

	//cout<<"The scaler_ID is "<<scaler_ID<<endl;
	//cout<<"The temp_ID is "<<temp_ID<<endl;
	//cout<<"The n is "<<n<<endl;

	//	}while (n++<size);

/***********************************************************************************/
//Fill tree //after one events loop
/**********************************************************************************/

		//	cout<<"fill data and plot histogram"<<endl;
		//	if(temp_ID==48){
		//
//	if(scaler_ID==1){
		//	cout<<"fill data and plot histogram"<<endl;

				for(int m = 0; m<64;m++)
				{
					if(m<4){
						rear[m]=data1[2][16+m];
				//		cout<<"the rear["<<m<<"] is "<<rear[m]<<" data0[2]["<<m<<"] is "<<data0[2][16+m]<<endl;
						si13[m]=data1[1][m];
						//cout<<"The data[1]["<<m<<"] = "<<data1[1][m]<<endl;
						si15[m]=data1[3][m];
						ge5[m] =data1[5][m];
						ge11[m]=data1[6][m];
						scalertt[m]=scaler_temp[m];
						qscint[m]=data1[7][m];
						tscint[m]=data1[8][m];

						}

					else if(m>3&&m<32){
						si13[m]=data1[1][m];

						//cout<<"The data[1]["<<m<<"] = "<<data1[1][m]<<endl;
						si15[m]=data1[3][m];
						ge5[m] =data1[5][m];
						ge11[m]=data1[6][m];
						scalertt[m]=scaler_temp[m];
						qscint[m]=data1[7][m];
						tscint[m]=data1[8][m];
						}
					else if(m>31&&m<48){
						si13[m]=data1[2][m-32];
						si15[m]=data1[4][m-32];
						}
					else
						{si15[m]=data1[4][m-32];}
				}
				day1data->Fill();
			//	rawdata->Fill();
				//for(int mm=0;mm<8;mm++) cout<<"The scaler data["<<mm<<"] is "<<scalertt[mm]<<endl;



                for(int k=0;k<32;k++){
                    h_QDC1[k]->Fill(data1[7][k]);
                    h_TDC1[k]->Fill(data1[8][k]);

                }//filling histogram for each events
					      							//	for(int k=0;k<33;k++)
							//	scaler_hist->Fill(scal_ch[k-1],scaler[k]);

						for(int i=0; i<64;i++)
						{
						if(i<32) {
							Ge_5mm_hits->Fill(Ge_5mm_strip[i],data1[5][i],1);
              Ge_11mm_hits->Fill(Ge_11mm_strip[i],data1[6][i],1);
							Si_13_hits->Fill(Si_13_strip[i],data1[1][i],1);
							Si_15_hits->Fill(Si_15_strip[i],data0[3][i],1);
				//	cout<<"2D histograms are filled"<<endl;
							}
						else if(i>31&&i<48) {
							Si_13_hits->Fill(Si_13_strip[i],data1[2][i-32],1);
							Si_15_hits->Fill(Si_15_strip[i],data0[4][i-32],1);
							}
						else Si_15_hits->Fill(Si_15_strip[i],data1[4][i-32],1);
			 				}

						for(int i=0; i<33;i++)
						{
 							Scint_qhits->Fill(scint[i],data1[7][i],1);
							Scint_thits->Fill(scint[i],data1[8][i],1);
			 			}


	//		}//fill the histogram when all modules data are decoded



	   } //loop size
/*
     eventheader1->Fill(evtADC1); cout<<"Fill the histogram, evtADC1="<<evtADC1<<endl;
     eventheader2->Fill(evtADC2+10);
     eventheader3->Fill(evtADC3+20);
     eventheader4->Fill(evtADC4+30);
     eventheader5->Fill(evtADC5+40);
     eventheader6->Fill(evtADC6+50);
     eventend1->Fill(evtendADC1); cout<<"Fill the histogram, evtendADC1="<<evtendADC1<<endl;
     eventend2->Fill(evtendADC2+10);
     eventend3->Fill(evtendADC3+20);
     eventend4->Fill(evtendADC4+30);
     eventend5->Fill(evtendADC5+40);
     eventend6->Fill(evtendADC6+50);
*/

//     		cout<<"evtADC1="<<evtADC1<<endl
//	    <<"evtADC2="<<evtADC2<<endl
//	    <<"evtADC3="<<evtADC3<<endl
//      <<"evtADC5="<<evtADC5<<endl
//
    return 0;
}


/*****************************************************************************************************/
//function of openfile()
/*****************************************************************************************************/

int openfile(char *file){
    int f, length, res, count = 0;
//	char *buffer;
       printf("In openfile(): Start to read data file and decoding\n");

    ifstream is;
    is.open(file,ios::binary);
	//get length of file:
	is.seekg (0,ios::end);
	length = is.tellg();
	is.seekg (0, ios::beg);
	printf("the length of the file is s%d\n",length);

	//allocate memeory
	 hist_filling(); //define the histograms for single channel

	f=open(file,O_RDONLY);
	//f = open(file,"r");
//	printf("the f is 0x%0x\n",f);
	if(f<0) {
		printf("In openfile(): open: %s\n", strerror(errno));
		return 1;
	}else {
		printf("In openfile(): Data file is open!\n");
	}

	do {
		res = decoding(f);
		//Filling ttree with 1 cluster data
		for(int evtnr=0;evtnr<300;evtnr++){
			for(int modid=1;modid<9;modid++){
				for(int n=0;n<33;n++){
					if(modid==1&&cldata[evtnr][modid][0]>10.){DADC1[n]=cldata[evtnr][modid][n];clevent=evtnr;
							if(n<32)h_ADC1[n]->Fill(DADC1[n]);
						}//cout<<"Evtnr="<<evtnr<<endl;clevent=evtnr;}
					if(modid==2&&cldata[evtnr][modid][0]>10.){DADC2[n]=cldata[evtnr][modid][n];
							if(n<32)h_ADC2[n]->Fill(DADC2[n]);
						}
					if(modid==3&&cldata[evtnr][modid][0]>10.){DADC3[n]=cldata[evtnr][modid][n];
							if(n<32)h_ADC3[n]->Fill(DADC3[n]);
						}
					if(modid==4&&cldata[evtnr][modid][0]>10.){DADC4[n]=cldata[evtnr][modid][n];
							if(n<32)h_ADC4[n]->Fill(DADC4[n]);
						}
					if(modid==5&&cldata[evtnr][modid][0]>10.){DADC5[n]=cldata[evtnr][modid][n];
							if(n<32)h_ADC5[n]->Fill(DADC5[n]);
						}
					if(modid==6&&cldata[evtnr][modid][0]>10.){DADC6[n]=cldata[evtnr][modid][n];
							if(n<32)h_ADC6[n]->Fill(DADC6[n]);
						}
					if(modid==7&&cldata[evtnr][modid][0]>10.){DQDC1[n]=cldata[evtnr][modid][n];
							if(n<32)h_QDC1[n]->Fill(DQDC1[n]);
						}
					if(modid==8&&cldata[evtnr][modid][0]>10.){DTDC1[n]=cldata[evtnr][modid][n];
							if(n<32)h_TDC1[n]->Fill(DTDC1[n]);
						}
					//if(modid==1)cout<<"DADC1["<<evtnr<<"]["<<modid<<"]["<<n<<"]="<<DADC1[n]<<endl;

				}//module events
			}//modules

		  if(clevent==evtnr)
			{
        rawdata->Fill();
			//	cout<<"clevtnr="<<evtnr<<endl;
		//	cout<<"Tree filled"<<endl;
      } // fill rawdata
	//	cout<<"EvtLoop Nr="<<evtnr<<endl;

		} //cluster events

      //  printf("In openfile(): (decoding status) res=decoding(f) is %d \n", res);
		if(count%100 == 0) cout<<"The No of "<<count<<" clusters have been decoded"<<endl;
	} while (count++<170000&&res ==0 );

	//close(f);
       printf("In openfile():   close the file!\n");

    return 0;
} // Open file for decoding


/*****************************************************************************************************/
//main() function
/*****************************************************************************************************/

int main()
{
	TStopwatch timer;
	timer.Start();

	TApplication app("app",0,0); //draw histograms needed

 //   FILE *pFile;
  char filename[80],filename2[100];
  cout<<"Type the file to be decoded "<<endl;
  scanf("%s",&filename);

	char froot[20] = ".ADC.root";
	strcpy(filename2,filename);
	strcat(filename2,froot);

  TFile *f = new TFile(filename2,"RECREATE");

	day1data->Branch("si13",&si13,"si13[48]/I");
	day1data->Branch("si15",&si15,"si15[64]/I");
	day1data->Branch("ge5",&ge5,"ge5[32]/I");
	day1data->Branch("ge11",&ge11,"ge11[32]/I");
	day1data->Branch("rear",&rear,"rear[4]/I");
	day1data->Branch("scaler",&scalertt,"scalertt[32]/I");
	day1data->Branch("qscint",&qscint,"qscint[14]/I");
	day1data->Branch("tscint",&tscint,"tscint[14]/I");

	rawdata->Branch("DADC1",&DADC1,"DADC1[33]/I");
	rawdata->Branch("DADC2",&DADC2,"DADC2[33]/I");
	rawdata->Branch("DADC3",&DADC3,"DADC3[33]/I");
	rawdata->Branch("DADC4",&DADC4,"DADC4[33]/I");
	rawdata->Branch("DADC5",&DADC5,"DADC5[33]/I");
	rawdata->Branch("DADC6",&DADC6,"DADC6[33]/I");
	rawdata->Branch("DQDC1",&DQDC1,"DQDC1[33]/I");
	rawdata->Branch("DTDC1",&DTDC1,"DTDC1[33]/I");

//	day1data->Branch("ev"&ev,"ev/I");

   gStyle->SetFrameFillColor(10);
   gStyle->SetPadBorderMode(1);
   gStyle->SetPadColor(10);
   gStyle->SetCanvasBorderMode(1);
   gStyle->SetCanvasColor(10);
   gROOT->ForceStyle();
   gStyle->SetOptStat();

	  for(int i=0;i<33;i++){
		scaler[i]=0;
		scaler_temp[i]=0;
	  }

    openfile(filename);
    printf("the file %s has been decoded! \n",filename);

    gStyle->SetOptStat();
    gStyle->SetOptFit();

 //   hist_filling();


	  for(int i=0;i<32;i++){
	  cout<<"scaler["<<i+1<<"]      is "<<scaler[i]<<endl;
//	  cout<<"scaler_temp["<<i<<"] is "<<scaler_temp[i]<<endl;
  	}

    for(int k=0;k<32;k++){
        scaler_hist->Fill(scal_ch[k],scaler[k]);
    }

TCanvas *c71 = new TCanvas("c71","All strips spectrum",80,80,1200,600);
  c71->Divide(2,2);

	c71->cd(1); gPad->SetLogz();
				Si_13_hits->SetXTitle("Strip No");
				Si_13_hits->SetYTitle("ADC channel");
//				Si_13_hits->GetXaxis()->SetRangeUser(first_madc, last_madc);
				Si_13_hits->Draw("colz");

	c71->cd(3); gPad->SetLogz();
				Si_15_hits->SetXTitle("Strip No");
				Si_15_hits->SetYTitle("ADC channel");
//				Si_15_hits->GetXaxis()->SetRangeUser(first_madc, last_madc);
				Si_15_hits->Draw("colz");

	c71->cd(2); gPad->SetLogz();
				Ge_5mm_hits->SetXTitle("Strip No");
				Ge_5mm_hits->SetYTitle("ADC channel");
//				Ge_5mm_hits->GetXaxis()->SetRangeUser(first_madc, last_madc);
				Ge_5mm_hits->Draw("colz");

	c71->cd(4); gPad->SetLogz();
				Ge_11mm_hits->SetXTitle("Strip No");
				Ge_11mm_hits->SetYTitle("ADC channel");
//				Ge_11mm_hits->GetXaxis()->SetRangeUser(first_madc, last_madc);
				Ge_11mm_hits->Draw("colz");
//   printf("Readout  procedure has been done!\n");

TCanvas *c72 = new TCanvas("c72","EventCounts",150,150,600,600);
   c72->Divide(1,2);
   c72->cd(1);
    gPad->SetLogz();
    gPad->SetLogy();

   eventheader1->SetLineColor(kBlack);
   eventheader1->SetLineWidth(3);
   eventheader1->Draw();
   eventheader2->SetLineColor(kRed);
   eventheader2->SetLineWidth(3);
   eventheader2->Draw("same");
   eventheader3->SetLineColor(kBlue);
   eventheader3->SetLineWidth(3);
   eventheader3->Draw("same");
   eventheader4->SetLineColor(kMagenta);
   eventheader4->SetLineWidth(3);
   eventheader4->Draw("same");
   eventheader5->SetLineColor(kGreen);
   eventheader5->SetLineWidth(3);
   eventheader5->Draw("same");
   eventheader6->SetLineColor(kOrange);
   eventheader6->SetLineWidth(3);
   eventheader6->Draw("same");

   c72->cd(2);
     gPad->SetLogy();
   eventend1->SetLineColor(kBlack);
   eventend1->SetLineWidth(3);
   eventend1->Draw();
   eventend2->SetLineColor(kRed);
   eventend2->SetLineWidth(3);
   eventend2->Draw("same");
   eventend3->SetLineColor(kBlue);
   eventend3->SetLineWidth(3);
   eventend3->Draw("same");
   eventend4->SetLineColor(kGreen);
   eventend4->SetLineWidth(3);
   eventend4->Draw("same");
   eventend5->SetLineColor(kMagenta);
   eventend5->SetLineWidth(3);
   eventend5->Draw("same");
   eventend6->SetLineColor(kOrange);
   eventend6->SetLineWidth(3);
   eventend6->Draw("same");




TCanvas *c74 = new TCanvas("c74","Time Stamp of ADC1",120,120,700,800);
   c74->Divide(1,2);
   c74->cd(1);
    gPad->SetLogz();
    gPad->SetLogy();
/*   timestampADC1->SetXTitle("Time (ms)");
   timestampADC1->SetYTitle("Counts");
*/
   timestampADC1->SetLineColor(kBlack);
   timestampADC1->SetLineWidth(3);
   timestampADC1->Draw();

    timestampADC2->SetLineColor(kRed);
    timestampADC2->SetLineWidth(3);
    timestampADC2->Draw("same");

    timestampADC3->SetLineColor(kBlue);
    timestampADC3->SetLineWidth(3);
    timestampADC3->Draw("same");

     timestampADC4->SetLineColor(kGreen);
     timestampADC4->SetLineWidth(3);
     timestampADC4->Draw("same");

     timestampADC5->SetLineColor(kMagenta);
     timestampADC5->SetLineWidth(3);
     timestampADC5->Draw("same");

      timestampADC6->SetLineColor(kOrange);
      timestampADC6->SetLineWidth(3);
      timestampADC6->Draw("same");

    timestampQDC1->SetLineColor(kBlue);
    timestampQDC1->SetLineWidth(3);
//    timestampQDC1->Draw("same");

    timestampTDC1->SetLineColor(kGreen);
    timestampTDC1->SetLineWidth(3);
//    timestampTDC1->Draw("same");


    float x0=0.6, y0=0.8, x1=0.72, y1=0.9;

    TLegend* leg3 = new TLegend(x0,y0,x1,y1);
    leg3->AddEntry(timestampADC1, "ADC1","L");
    leg3->AddEntry(timestampADC2, "ADC2","L");
//    leg3->AddEntry(timestampQDC1, "QDC1","L");
//    leg3->AddEntry(timestampTDC1, "TDC1","L");
 //   leg3->AddEntry(h_ge1[2], "Strip3","L");

    leg3->SetMargin(0.15);
    leg3->SetFillColor(0);
    leg3->Draw();

    c74->cd(2);
     gPad->SetLogz();
     gPad->SetLogy();
    clustersize->SetXTitle("Size of cluster");
    clustersize->SetYTitle("Counts");
    clustersize->Draw();



    Ge_5mm_hits->Write();
    Ge_11mm_hits->Write();
	  Si_13_hits->Write();
	  Si_15_hits->Write();
	  scaler_hist->Write();
	  Scint_qhits->Write();
	  Scint_thits->Write();

    timestampADC1->Write();
    timestampADC2->Write();
    timestampADC3->Write();
    timestampADC4->Write();
    timestampADC5->Write();
    timestampADC6->Write();
    timestampQDC1->Write();
    timestampTDC1->Write();


//    hist_writing();

 for(Int_t j=0;j<32;j++){
        h_ADC1[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC2[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC3[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC4[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC5[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_ADC6[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_QDC1[j]->Write();}
 for(Int_t j=0;j<32;j++){
        h_TDC1[j]->Write();}

    f->cd();


	day1data->Print();
	day1data->Write();

	rawdata->Print();
	rawdata->Write();

     f->Write();
     f->Close();

	app.Run();

	timer.Stop();
    Double_t rtime = timer.RealTime();
    Double_t ctime = timer.CpuTime();
    printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);

   return 0;
 } // main function
