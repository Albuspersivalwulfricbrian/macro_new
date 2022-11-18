#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TF1.h>
// #include <TCanvas.h>
#include "like_ivashkin_wants_it.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include "ChannelEntry.h"
#include "CHSH_class.h"
#include "Double_scattering_cuts.h"
using namespace std;
using namespace CUTS;
#define PresentIntermediate 1


/////////////////////////////////////////////////////
void Create_mini_tree_short()
{
    const Int_t events_divider = 1;
	gStyle->SetOptFit(1);
	TString source_path = "/home/doc/entanglement/double_GAGG/";

    TString result_path = source_path + "MiniDST";

    TFile *mini_tree_file = new TFile(result_path+".root", "RECREATE");
    TTree *short_tree = new TTree("Signals","Signals");
    mini_tree_nrg mini_leaves;
    mini_tree_time mini_time; 

    short_tree->Branch("MiniTree",&mini_leaves);
    short_tree->Branch("TimeTree",&mini_time);

///////////////////////////////////////////////
    TChain *PMT_tree = new TChain;
	PMT_tree->AddFile( (source_path + "calibrated_time_0.root" + "/adc64_data").Data() );
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////
    const Int_t calculate_Events_Number = PMT_tree->GetEntries()/events_divider;
    Int_t Events_for_cuts = calculate_Events_Number/1;
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////
#if PresentIntermediate    
    const Int_t total_channels = 36;

    std::array<short_energy_ChannelEntry*, total_channels> short_channel_info;
    for(Int_t ch = 0; ch < total_channels; ch++)
    {
        short_channel_info[ch] = new short_energy_ChannelEntry();
        short_channel_info[ch]->Initialize();
        PMT_tree->SetBranchAddress((TString::Format("channel_%i", ch)).Data(), &short_channel_info[ch]);
    }

    std::array<diff_short_energy_ChannelEntry*, total_channels> diff_4_cut;
    for(Int_t ch = 34; ch < total_channels; ch++)
    {
        diff_4_cut[ch] = new diff_short_energy_ChannelEntry();
        diff_4_cut[ch]->Initialize();
        PMT_tree->SetBranchAddress((TString::Format("diff_channel_%i", ch)).Data(), &diff_4_cut[ch]);
    }
#else
    const Int_t total_channels = 34;
    short_energy_ChannelEntry *short_channel_info = new short_energy_ChannelEntry[total_channels];
    for(Int_t ch = 0; ch < total_channels; ch++)
        (short_channel_info+ch)->SetBranch(PMT_tree, ch);
#endif

///////////////////////////////////

    for (Int_t NumEvent = 0; NumEvent < calculate_Events_Number; NumEvent++)
    {
        PMT_tree->GetEntry(NumEvent);
        Short_t counter_1 = 0; Short_t ch1 = 0;
        Short_t counter_2 = 0; Short_t ch2 = 0;
        Short_t EventType = -10;
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            Int_t sc_number = 33;
            if (channel_number < 16) sc_number = 32;
            if 
            (
                Apply_time_in_peak_cuts<total_channels>(short_channel_info, channel_number, sc_number)
            && Apply_Amplitude_Saturation_cuts<total_channels>(short_channel_info, channel_number) 

            #if PresentIntermediate
            && EventTypeDeterminator<total_channels>(short_channel_info, EventType) 
            #endif
            && short_channel_info[channel_number]->charge > 10
            ) 
            {
                if (sc_number == 32) {counter_1++; ch1 = channel_number;}
                if (sc_number == 33) {counter_2++; ch2 = channel_number;}
            }
        }

        if (NumEvent%100000 == 0) cout << counter_1 << " " << counter_2 << "\t" << NumEvent/calculate_Events_Number <<"%" << endl;
        EventType = -10;
        if((counter_1 == 1 && counter_2 == 1) && (counter_1 < 2 && counter_2 < 2))
        {        
            Int_t channel_number = ch1; Int_t channel_number_2 = ch2;
            if(
            Apply_Amplitude_Saturation_cuts<total_channels>(short_channel_info, channel_number, channel_number_2)
            && Apply_time_in_peak_cuts<total_channels>(short_channel_info, channel_number, 32)
            && Apply_time_in_peak_cuts<total_channels>(short_channel_info, channel_number_2, 33)
            #if PresentIntermediate
            && EventTypeDeterminator<total_channels>(short_channel_info,EventType)
            && diff_cuts(diff_4_cut, EventType)

            #endif
            )
            {
                Float_t EdepIntermediate = 0.;
                EdepIntermediate = short_channel_info[34]->charge;
                if ( channel_number == ch1 && channel_number_2 == ch2) 
                {
                    mini_leaves.EdepScat0 = short_channel_info[32]->charge;
                    mini_leaves.EdepScat1 = short_channel_info[33]->charge;
                    mini_leaves.EdepDet0 = short_channel_info[ch1]->charge;
                    mini_leaves.EdepDet1 = short_channel_info[ch2]->charge;
                    mini_leaves.DetNum0 = ch1;
                    mini_leaves.DetNum1 = ch2;  
                    mini_leaves.EdepIntermediate0 = short_channel_info[34]->charge;
                    mini_leaves.EdepIntermediate1 = short_channel_info[35]->charge;
                    mini_leaves.EventType = EventType;
                    mini_time.TimeScat0 = short_channel_info[32]->time;
                    mini_time.TimeScat1 = short_channel_info[33]->time;
                    mini_time.TimeDet0 = short_channel_info[ch1]->time;
                    mini_time.TimeDet1 = short_channel_info[ch2]->time;
                    mini_time.TimeIntermediate0 = short_channel_info[34]->time;
                    mini_time.TimeIntermediate1 = short_channel_info[35]->time;
                    short_tree-> Fill();
                }
            }
        }
    }
    short_tree->Write();
    mini_tree_file->Close();

}

