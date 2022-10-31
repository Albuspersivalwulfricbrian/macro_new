#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <RooDouble.h>
#include "like_ivashkin_wants_it.h"
#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TString.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TGraphErrors.h>
#include <TObjString.h>
#include "ChannelEntry.h"
#include "CHSH_calculator.h"
#include "Check_entanglement_cuts.h"

using namespace std;
using namespace CHSH;
using namespace CUTS;

#define UseManyRoots 0
#define PresentIntermediate 1
#define CalculateRatio 1
#define DrawTime 1
#define UseIntegralCut 1
#define UseNotEntangledPhotons 1
#define DrawToPDF 0
#define calculate_CHSH 0

/////////////////////////////////////////////////////
void Create_mini_tree()
{
    const Int_t events_divider = 1;
	gStyle->SetOptFit(1);
	//gStyle->SetOptStat(1111);
	//TString source_path = "/home/doc/entanglement/root_files_data/with_scatterer/big_file/";
	TString source_path = "/home/doc/entanglement/with_spline/entangled/";

#if UseNotEntangledPhotons
    TString result_path  = source_path + "decoh_mini_tree";
    const Bool_t entangled = kFALSE;
#else
    TString result_path = source_path + "entangled_mini_tree";
    const Bool_t entangled = kTRUE;
#endif
    TFile *mini_tree_file = new TFile(result_path+".root", "RECREATE");
    TTree *short_tree = new TTree("Signals","Signals");
    mini_tree_nrg mini_leaves;
    mini_tree_time mini_time; 
    mini_leaves.CreateBranches(short_tree);
    mini_time.CreateBranches(short_tree);

//////////////
    const Int_t left_NaI_range = 180;
    const Int_t right_NaI_range = 280;
    Float_t interval_width = 1.35;
////////////////
    Float_t low_cut[32] = {0}; Float_t high_cut[32] = {0}; 
    Float_t total_low_cut[32] = {0}; Float_t total_high_cut[32] = {0};
    Float_t sigma_i[32] = {0}; Float_t mean_int[32] = {0};
    Float_t average_scatterer_peak_position[32] = {0};
    Float_t sigma_scatterer[32] = {0};
    Float_t low_energy_cut[32] = {0};
    Float_t high_energy_cut[32] = {0};
    Float_t sum_energy[32] = {0};
    Float_t sigma_energy[32] = {0};
    Float_t low_energy_sum[32] = {0};
    Float_t high_energy_sum[32] = {0};

/////////////////////
    const Int_t histos_bins_number = 200;
    const Int_t histos_left_boarder = 1;
    const Int_t histos_right_boarder = 450;
    TH1F *true_sum_hist_left = new TH1F ("TRUE_peak_left_sum_spectrum","TRUE_peak_left_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
    TH1F *true_sum_hist_right = new TH1F ("TRUE_peak_right_sum_spectrum","TRUE_peak_right_sum_spectrum",histos_bins_number,histos_left_boarder,histos_right_boarder);
    TH1F *true_sc_sum_hist_left = new TH1F ("TRUE_peak_scatterer_left_sum_spectrum","TRUE_peak_scatterer_left_sum_spectrum",200,1,700);
    TH1F *true_sc_sum_hist_right = new TH1F ("TRUE_peak_scatterer_right_sum_spectrum","TRUE_peak_scatterer_right_sum_spectrum",200,1,700);
/////////////////////
    TH1F *hist_charge[16][16];
    //TH1F *charge_hist[32];
    Int_t NumEvents[16][16] = {0};     
    Float_t angle_arr[16];
    Float_t angle_arr_err[16];
    for (Int_t channel_number = 0; channel_number < 16; channel_number++)
    {
        for (Int_t channel_number_2 = 16; channel_number_2 < 32; channel_number_2++)    
        {
            hist_charge[channel_number][channel_number_2-16] = 
            new TH1F(Form("charge_%i_%i",channel_number,channel_number_2),
            Form("charge_%i_%i",channel_number,channel_number_2),140,10,1000);   
        }     
    }    

    for (Int_t i = 0; i <16; i++)
    {
        angle_arr[i] = (float)i*22.5;
        angle_arr_err[i] = 0.;
    }
    Float_t total_events_for_angle_diff[16] = {0};
    Float_t total_events_for_angle_diff_err[16] = {0};

///////////////////////////////////////////////
    TChain *PMT_tree = new TChain;
	PMT_tree->AddFile( (source_path + "calibrated_time.root" + "/adc64_data").Data() );
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////
    const Int_t calculate_Events_Number = PMT_tree->GetEntries()/events_divider;
    Int_t Events_for_cuts = calculate_Events_Number/1;
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////
#if PresentIntermediate    
    const Int_t total_channels = 35;
    short_energy_ChannelEntry *short_channel_info = new short_energy_ChannelEntry[total_channels];
    for(Int_t ch = 0; ch < total_channels; ch++)
        (short_channel_info+ch)->SetBranch(PMT_tree, ch);
    diff_short_energy_ChannelEntry *diff_4_cut = new diff_short_energy_ChannelEntry;
    diff_4_cut->SetBranch(PMT_tree, 34);
#else
    const Int_t total_channels = 34;
    short_energy_ChannelEntry *short_channel_info = new short_energy_ChannelEntry[total_channels];
    for(Int_t ch = 0; ch < total_channels; ch++)
        (short_channel_info+ch)->SetBranch(PMT_tree, ch);
#endif

/////////////////////////charge_cut_range_cuts

    TH1F *peak_histo[35];
    TH1F *peak_histo_w_o_cuts[35];
    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        TString hernya = Form("charge_channel_%i",channel_number);
        peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),histos_bins_number,histos_left_boarder,histos_right_boarder);
    }
    for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
    {
        PMT_tree->GetEntry(iEvent);
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            Int_t sc_number = 33;
            if (channel_number < 16) sc_number = 32;
            if (
                Apply_time_in_peak_cuts(short_channel_info, channel_number, sc_number)
            && Apply_Amplitude_Saturation_cuts(short_channel_info, channel_number)

            #if UseIntegralCut
            && short_channel_info[32].charge > 140
            && short_channel_info[32].charge < 400             
            #endif   

            #if PresentIntermediate
            && Decoherent_or_Entangled(short_channel_info, entangled) 
            #endif
            ) 
            {
                peak_histo[channel_number]->Fill(short_channel_info[channel_number].charge);
            }
        }
    }

    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        double_gauss_fit(peak_histo[channel_number], low_cut[channel_number], high_cut[channel_number]);
        hist_like_ivashkin_wants_it(peak_histo[channel_number],"Energy [keV]", "Counts");
        peak_histo[channel_number]->Write();

        if (channel_number < 16) true_sum_hist_left->Add(peak_histo[channel_number]);
        if (channel_number >=16) true_sum_hist_right->Add(peak_histo[channel_number]);
        
        #if DrawToPDF
        TCanvas *canv = new TCanvas("canv","canv");
        canv->cd();
        peak_histo[channel_number]->Draw();
        canv->SaveAs(result_path+".pdf(",".pdf");
        #endif
        delete peak_histo[channel_number];
    }
    TH1F *int_without_cuts = new TH1F("ch_34.charge","ch_34.charge",200,1,300);
    TH1F *int_with_cuts = new TH1F("ch_34.charge_with_time_cuts","ch_34.charge_with_time_cuts",200,1,300);

    TCanvas *canv_0 = new TCanvas("canv_0","canv_0");
    canv_0->Divide(2);
    canv_0->cd(1);
    true_sum_hist_left->Draw();
    canv_0->cd(2);
    true_sum_hist_right->Draw();
    canv_0->SaveAs(result_path+".pdf(",".pdf");
    
    #if DrawToPDF
    TCanvas *temp_canv = new TCanvas("temp_canv","temp_canv");
    temp_canv->Divide(2);
    temp_canv->cd(1);
    PMT_tree->Draw("channel_34.charge >> ch_34.charge","channel_34.amp > 400 && channel_34.amp < 60000");
    int_without_cuts->GetYaxis()->SetRangeUser(0,(Int_t)1.05*int_without_cuts->GetMaximum());
    temp_canv->cd(2);
    PMT_tree->Draw("channel_34.charge >> ch_34.charge_with_time_cuts","channel_34.amp > 400 && channel_34.amp < 60000 && channel_34.peak_pos - channel_32.peak_pos > -30 && channel_34.peak_pos - channel_32.peak_pos < 60");
    int_with_cuts->GetYaxis()->SetRangeUser(0,(Int_t)1.05*int_without_cuts->GetMaximum());
    temp_canv->SaveAs(result_path + ".pdf(",".pdf");
    #endif

////////////////////////energy in diffuser
    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        TString hernya = Form("charge_diffuser_ch_cut_%i",channel_number);
        peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),100,0,200);
    }
    for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
    {
        PMT_tree->GetEntry(iEvent);
        for (Int_t channel_number = 0; channel_number < 16; channel_number++)
        {
            if (short_channel_info[channel_number].charge > low_cut[channel_number]
            && short_channel_info[channel_number].charge < high_cut[channel_number]
            && Apply_time_in_peak_cuts(short_channel_info, channel_number)
            && Apply_Amplitude_Saturation_cuts(short_channel_info,channel_number)
            
            #if PresentIntermediate
            && Decoherent_or_Entangled(short_channel_info,entangled)
            #endif    
            ) 
            peak_histo[channel_number]->Fill(short_channel_info[34].charge);
        }
    }

        TString hist_name = "charge_diffuser";
        TH1F *d_hist = new TH1F (hist_name.Data(),hist_name.Data(),100,0,200);
        for (Int_t channel_number = 0; channel_number < 16; channel_number++)
        {
            d_hist->Add(peak_histo[channel_number]);
            delete peak_histo[channel_number];

        }
            hist_like_ivashkin_wants_it(d_hist,"Energy [keV]", "Counts");
            d_hist->Write();
            TCanvas *canv_1 = new TCanvas("canv","canv");
            canv_1-> cd();
            d_hist->Draw();
            canv_1->SaveAs(result_path+".pdf(",".pdf");

///////////////////////////////
/////////energy in scatterer
/////////full_spectrum
        TCanvas *canv_00 = new TCanvas("canv_00","canv_00");
        canv_00->Divide(2);
        canv_00->cd(1);
        PMT_tree->Draw("channel_32.charge >> TRUE_peak_scatterer_left_sum_spectrum"
        ,"channel_32.peak_pos - channel_33.peak_pos > -10 && channel_32.peak_pos - channel_33.peak_pos < 10"
        "&& channel_32.amp < 60000 && channel_32.amp > 400"
        "&& channel_33.amp > 400 && channel_33.amp < 60000"

#if UseNotEntangledPhotons
        "&& channel_32.peak_pos - channel_34.peak_pos > -60 && channel_32.peak_pos - channel_34.peak_pos < 60"
        "&& channel_34.amp > 400 && channel_34.amp < 60000"            
#endif

        "",0,Events_for_cuts);
        true_sc_sum_hist_left->GetYaxis()->SetRangeUser(0,(Int_t)1.05*true_sc_sum_hist_left->GetMaximum());
        canv_00->cd(2);
        PMT_tree->Draw("channel_33.charge >> TRUE_peak_scatterer_right_sum_spectrum"
        ,"channel_32.peak_pos - channel_33.peak_pos > -10 && channel_32.peak_pos - channel_33.peak_pos < 10"
        "&& channel_32.amp < 60000 && channel_32.amp > 400"
        "&& channel_33.amp > 400 && channel_33.amp < 60000"
#if UseNotEntangledPhotons
        "&& channel_32.peak_pos - channel_34.peak_pos > -70 && channel_32.peak_pos - channel_34.peak_pos < -40"
        "&& channel_34.amp > 400 && channel_34.amp < 60000"            
#endif
        "",0,Events_for_cuts);
        true_sc_sum_hist_right->GetYaxis()->SetRangeUser(0,(Int_t)1.05*true_sc_sum_hist_right->GetMaximum());

        canv_00->SaveAs(result_path + ".pdf(",".pdf");
        true_sc_sum_hist_left->Write();
        true_sc_sum_hist_right->Write();

//////if hits NaI near photopeak
    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        TString hernya = Form("charge_scatterer_ch_cut_%i",channel_number);
        peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),100,1,500);
    }

    for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
    {
        PMT_tree->GetEntry(iEvent);
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            Int_t sc_number = 0;
            if (channel_number < 16) sc_number = 32;
            if (channel_number >= 16) sc_number = 33;
            if (
                short_channel_info[channel_number].charge > low_cut[channel_number]
            && short_channel_info[channel_number].charge < high_cut[channel_number]
            && Apply_time_in_peak_cuts(short_channel_info, channel_number, sc_number)
            && Apply_Amplitude_Saturation_cuts(short_channel_info, channel_number)

            #if PresentIntermediate
            && Decoherent_or_Entangled(short_channel_info, entangled)
            #endif      
            ) 
                peak_histo[channel_number]->Fill(short_channel_info[sc_number].charge);
                if (iEvent%10000==0) cout << low_cut[channel_number] << " for ch = "<< channel_number <<endl;
        }
    }

    TString left_hist_name = "charge_ch_32";
    TString right_hist_name = "charge_ch_33";
    TH1F *right_hist = new TH1F (right_hist_name.Data(),right_hist_name.Data(),100,1,500);
    TH1F *left_hist = new TH1F (left_hist_name.Data(),left_hist_name.Data(),100,1,500);


    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        cout << "channel_number = "<<channel_number <<endl;
        double_gauss_fit(peak_histo[channel_number], low_energy_cut[channel_number], high_energy_cut[channel_number]);
        peak_histo[channel_number]->Write();
    if (channel_number < 16) left_hist->Add(peak_histo[channel_number]);
    if (channel_number >= 16) right_hist->Add(peak_histo[channel_number]);

        delete peak_histo[channel_number];

    }
    hist_like_ivashkin_wants_it(left_hist,"Energy [keV]", "Counts");
    hist_like_ivashkin_wants_it(right_hist,"Energy [keV]", "Counts");

    left_hist->Write();
    right_hist->Write();

    #if DrawToPDF
    TCanvas *canv = new TCanvas("canv","canv");
    canv->Divide(2);
    canv-> cd(1);
    left_hist->Draw();
    canv->cd(2);
    right_hist->Draw();

    canv->SaveAs(result_path+".pdf(",".pdf");
    delete canv;

    #endif

    delete right_hist;
    delete left_hist;
        
//////////////////////////////////sum_energy_histos
///////////////////////////////////

    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        sum_energy[channel_number] = 510;
        sigma_energy[channel_number] = sqrt(0.5*(pow((float)sigma_scatterer[channel_number],2)+pow((float)sigma_i[channel_number],2)));
        cout<< channel_number<<endl;
        cout<< sigma_scatterer[channel_number] << " = sigma_scatterer" <<endl; 
        cout<< sigma_energy[channel_number] << " = sigma_energy" <<endl; 
        low_energy_sum[channel_number] = sum_energy[channel_number] - interval_width*sigma_energy[channel_number];
        cout<< low_energy_cut[channel_number] << " = low_energy_cut" <<endl; 
        cout<< low_energy_sum[channel_number] << " = low_energy_sum" <<endl<<endl; 
        high_energy_sum[channel_number] = sum_energy[channel_number] + interval_width*sigma_energy[channel_number];
        cout<< high_energy_cut[channel_number] << " = high_energy_cut" <<endl; 
        cout<< high_energy_sum[channel_number] << " = high_energy_sum" <<endl<<endl;           
    }
////////////////////////////////////
//////////////////////////////////////////////////
#if DrawTime
    for (Int_t channel_number = 0; channel_number < 35; channel_number++)
    {
        Int_t sc_number = 33;
        if (channel_number == 33) sc_number = 34;
        if (channel_number < 16 || channel_number == 34) sc_number = 32;
        TString hernya = Form("Time_difference_scatterer_%i_and_channel_%i_with_energy_cut",sc_number,channel_number);
        peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),6001,-3000,3000);
        
        if (channel_number < 32)
        {
            TString fignya = hernya + "without_energy_cuts";
            peak_histo_w_o_cuts[channel_number] = new TH1F (fignya.Data(),fignya.Data(),6001,-3000,3000);                
        }
    }

    for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
    {
        PMT_tree->GetEntry(iEvent);
        for (Int_t channel_number = 0; channel_number < 35; channel_number++)
        {
            Int_t sc_number = 33;
            if (channel_number < 16 || channel_number == 34) sc_number = 32;
            if (channel_number < 32)
            {
                if (
                Non_zero_time(short_channel_info, channel_number)
                && Apply_Amplitude_Saturation_cuts(short_channel_info, channel_number)

#if PresentIntermediate
                && Decoherent_or_Entangled(short_channel_info,entangled)
#endif
                )
                {
                    peak_histo_w_o_cuts[channel_number]->Fill(short_channel_info[sc_number].peak_pos
                    -short_channel_info[channel_number].peak_pos);
                    if (short_channel_info[channel_number].charge > low_cut[channel_number]
                    && short_channel_info[channel_number].charge < high_cut[channel_number] )
                        peak_histo[channel_number]->Fill(short_channel_info[sc_number].peak_pos - short_channel_info[channel_number].peak_pos);
                }                    
            }
            else
            {
                if (channel_number == 33) sc_number = 34;
                if (
                    #if PresentIntermediate
                    Decoherent_or_Entangled(short_channel_info, entangled)&& 
                    #endif
                    Apply_Amplitude_Saturation_cuts(short_channel_info, channel_number)
                )
                {
                    if (channel_number!=34)
                        peak_histo[channel_number]->Fill(short_channel_info[sc_number].peak_pos
                        -short_channel_info[channel_number].peak_pos);
                    if (channel_number==34)
                        peak_histo[channel_number]->Fill(short_channel_info[33].peak_pos
                        -short_channel_info[34].peak_pos);                        
                }
            }
        }
    }
    for (Int_t channel_number = 0; channel_number < 35; channel_number++)
    {
        cout << "channel_number = "<<channel_number <<endl;
        #if DrawToPDF
        TCanvas *canv = new TCanvas("canv","canv");
        #endif
        if (channel_number < 32)
        {
            hist_like_ivashkin_wants_it(peak_histo_w_o_cuts[channel_number],"Time [ns]", "Counts");
            hist_like_ivashkin_wants_it(peak_histo[channel_number],"Time [ns]", "Counts");
            

            peak_histo_w_o_cuts[channel_number]->Write();
            peak_histo[channel_number]->Write();

            #if DrawToPDF
            canv->Divide(2);
            canv->cd(1);
            peak_histo_w_o_cuts[channel_number]->Draw();
            canv->cd(2);
            peak_histo[channel_number]->Draw();
            #endif
        }
        else
        {
            hist_like_ivashkin_wants_it(peak_histo[channel_number],"Time [ns]", "Counts");

            #if DrawToPDF
            canv->cd();
            peak_histo[channel_number]->Draw();
            #endif
            peak_histo[channel_number]->Write();

        }
        //if (channel_number < 33) double_gauss_fit(peak_histo[channel_number], low_time_cut[channel_number], high_time_cut[channel_number],1.5,1.5);

        #if DrawToPDF
        canv->SaveAs(result_path+".pdf(",".pdf");
        delete canv;
        #endif

        delete peak_histo[channel_number];
    }
#endif
///////////////////////////////////
    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        TString hernya = Form("charge_sum_with_scatterer_ch_cut_%i",channel_number);
//            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),70,-350,350);
        peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),70,200,900);

    }

    for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
    {
        PMT_tree->GetEntry(iEvent);
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            Int_t sc_number = 0; Short_t time_min = 0; Short_t time_max = 0; Float_t Ew = 0;
            if (channel_number < 16) {sc_number = 32; Ew = short_channel_info[34].charge;}
            if (channel_number >= 16) sc_number = 33;
            if (
                short_channel_info[channel_number].charge > low_cut[channel_number]
            && short_channel_info[channel_number].charge < high_cut[channel_number]
            && short_channel_info[sc_number].charge > low_energy_cut[channel_number]
            && short_channel_info[sc_number].charge < high_energy_cut[channel_number]
            && Apply_time_in_peak_cuts(short_channel_info, channel_number, sc_number
            //,low_time_cut[32],high_time_cut[32],low_time_cut[channel_number],high_time_cut[channel_number]
            )
            
            #if PresentIntermediate
            //&& short_channel_info[34].charge < 150
            && Decoherent_or_Entangled(short_channel_info, entangled)
            #endif          
            && Apply_Amplitude_Saturation_cuts(short_channel_info, channel_number)
            )
            {
                peak_histo[channel_number]->Fill(short_channel_info[sc_number].charge+short_channel_info[channel_number].charge+Ew);        
            }
        }
    }

    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
    {
        cout << "channel_number = "<<channel_number <<endl;

        hist_like_ivashkin_wants_it(peak_histo[channel_number],"Energy [keV]", "Counts");
        peak_histo[channel_number]->Write();

        #if DrawToPDF
        TCanvas *canv = new TCanvas("canv","canv");
        canv->cd();
        peak_histo[channel_number]->Draw();
        canv->SaveAs(result_path+".pdf(",".pdf");
        delete canv;
        #endif

        delete peak_histo[channel_number];
    }


///////////////////////////////////
//#if CalculateRatio
    for (Int_t NumEvent = 0; NumEvent < calculate_Events_Number; NumEvent++)
    {
        PMT_tree->GetEntry(NumEvent);
        Short_t counter_1 = 0; Short_t ch1 = 0;
        Short_t counter_2 = 0; Short_t ch2 = 0;
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            Int_t sc_number = 33;
            if (channel_number < 16) sc_number = 32;
            if 
            (
                Apply_time_in_peak_cuts(short_channel_info, channel_number, sc_number)
            && Apply_Amplitude_Saturation_cuts(short_channel_info, channel_number) 

            #if PresentIntermediate
            && Decoherent_or_Entangled(short_channel_info, entangled) 
            #endif
            && short_channel_info[channel_number].charge > 10
            ) 
            {
                if (sc_number == 32) {counter_1++; ch1 = channel_number;}
                if (sc_number == 33) {counter_2++; ch2 = channel_number;}
            }
        }

        if (NumEvent%100000 == 0) cout << counter_1 << " " << counter_2 << endl;
        if((counter_1 == 1 && counter_2 == 1) && (counter_1 < 2 && counter_2 < 2))
        {        
            Int_t channel_number = ch1; Int_t channel_number_2 = ch2;
            if(
            Apply_Amplitude_Saturation_cuts(short_channel_info, channel_number, channel_number_2)
            && Apply_time_in_peak_cuts(short_channel_info, channel_number, 32)
            && Apply_time_in_peak_cuts(short_channel_info, channel_number_2, 33)
            #if PresentIntermediate
            && Decoherent_or_Entangled(short_channel_info,entangled)
            && diff_cuts(diff_4_cut->min_diff,diff_4_cut->max_diff,entangled)
            #endif
            )
            {
                Float_t EdepIntermediate = 0.;
                //#if UseNotEntangledPhotons
                EdepIntermediate = short_channel_info[34].charge;
                //#endif
                if ( channel_number == ch1 && channel_number_2 == ch2) 
                {
                    mini_leaves.EdepScat0 = short_channel_info[32].charge;
                    mini_leaves.EdepScat1 = short_channel_info[33].charge;
                    mini_leaves.EdepDet0 = short_channel_info[ch1].charge;
                    mini_leaves.EdepDet1 = short_channel_info[ch2].charge;
                    mini_leaves.DetNum0 = ch1;
                    mini_leaves.DetNum1 = ch2;  
                    mini_leaves.EdepIntermediate = EdepIntermediate;
                    mini_time.TimeScat0 = short_channel_info[32].peak_pos;
                    mini_time.TimeScat1 = short_channel_info[33].peak_pos;
                    mini_time.TimeDet0 = short_channel_info[ch1].peak_pos;
                    mini_time.TimeDet1 = short_channel_info[ch2].peak_pos;
                    mini_time.TimeIntermediate = short_channel_info[34].peak_pos;
                    short_tree-> Fill();
                }
            }
        }

    }

    for (Int_t channel_number = 0; channel_number < 16; channel_number++)
    {
        for(Int_t channel_number_2 = 16; channel_number_2 < 32; channel_number_2++)
        {
            NumEvents[channel_number][channel_number_2-16] = hist_charge[channel_number][channel_number_2-16]->GetEntries();
            /*TCanvas *canv = new TCanvas("canv","canv");
            canv->cd();
            hist_charge[channel_number][channel_number_2-16]->Draw();
            canv->SaveAs(source_path+"result_ev_b_ev.pdf(",".pdf");*/
            Int_t angle = channel_number_2-channel_number-16;
            if (angle < 0) angle += 16;
            total_events_for_angle_diff[angle]+=NumEvents[channel_number][channel_number_2-16];
        }
    }

    for (Int_t i = 0; i < 16; i++) total_events_for_angle_diff_err[i] = pow(total_events_for_angle_diff[i],0.5);

    TGraphErrors * gr = new TGraphErrors(16,angle_arr,total_events_for_angle_diff, angle_arr_err, total_events_for_angle_diff_err);
	TF1 *sin_fit = new TF1 ("sin_fit","[0]+[1]*cos(3.14159265359/180*2*x)",0,360);
	gr->Fit("sin_fit","R");
    Float_t a_par = abs(sin_fit->GetParameter(0));
    Float_t b_par = abs(sin_fit->GetParameter(1));
    Float_t a_error = abs(sin_fit->GetParError(0));
    Float_t b_error = abs(sin_fit->GetParError(1));
    Float_t diff = (a_par+b_par)/(a_par-b_par);
    //Float_t diff_err = diff*sqrt(2*pow((a_error/a_par),2)+2*pow((b_error/b_par),2));
    Float_t diff_err = 2*a_par/(pow(a_par-b_par,2))*sqrt(pow((a_error*b_par/a_par),2)+pow((b_error),2));
    TCanvas *canvas = new TCanvas ( "canvas", "canvas");
    canvas->cd();
    graph_like_ivashkin_wants_it(gr,"angle [degrees]","Counts", Form("%5.1f-%5.1f*cos(2x) ratio = %4.3f +- %4.3f",a_par,b_par,diff,diff_err),1);
    gr->Draw("AP");
    canvas->SaveAs(result_path+".pdf)",".pdf");

    short_tree->Write();
    mini_tree_file->Close();

}

