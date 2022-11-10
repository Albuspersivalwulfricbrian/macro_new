#ifndef CHECK_CUTS
#define CHECK_CUTS
#include "ChannelEntry.h"
//#include "TreeStructures.h"
namespace CUTS
{
    const Int_t channelsnumb = 36;
//Template
    template<int n> 
    Bool_t isEntangled(std::array<short_energy_ChannelEntry*, n> &short_channel_info, TString side = "left")
    {
        Bool_t flag = kFALSE;
        Int_t scnum = 32; Int_t gnum = 34;
        if (side == "right") {scnum = 33; gnum = 35;}
        if ((short_channel_info[gnum]->amp < 600 &&
            ((short_channel_info[gnum]->time - short_channel_info[scnum]->time) < 150 ||
            (short_channel_info[gnum]->time - short_channel_info[scnum]->time) > 400)) || short_channel_info[gnum]->amp < 100
            )
            flag = kTRUE;   
        return flag;   
    }
    
    template<int n> 
    Bool_t isDecoherent(std::array<short_energy_ChannelEntry*, n> &short_channel_info, TString side = "left")
    {
        Bool_t flag = kFALSE;
        Int_t scnum = 32; Int_t gnum = 34;
        if (side == "right") {scnum = 33; gnum = 35;}
        if (short_channel_info[gnum]->amp < 60000
            && short_channel_info[gnum]->amp > 400
            && (short_channel_info[gnum]->time - short_channel_info[scnum]->time) < 350
            && (short_channel_info[gnum]->time - short_channel_info[scnum]->time) > 150)
            flag = kTRUE;   
        return flag;   
    }
    template<int n> 
    Bool_t Double_Decoherent_or_Decoherent_or_Entangled(std::array<short_energy_ChannelEntry*, n> &short_channel_info, TString state = "entangled")
    {
        Bool_t flag = kFALSE;
        // if (state == kFALSE && short_channel_info[34]->amp < 60000
        //     && short_channel_info[34]->amp > 100
        //     && (short_channel_info[34]->time - short_channel_info[32]->time) < 350
        //     && (short_channel_info[34]->time - short_channel_info[32]->time) > 150)
        //     flag = kTRUE;

        // if (state == kTRUE &&
        //     short_channel_info[34]->amp < 20 &&
        //     ((short_channel_info[34]->time - short_channel_info[32]->time) < 50 ||
        //     (short_channel_info[34]->time - short_channel_info[32]->time) > 500)
        //     )
        //     flag = kTRUE;   
        if (state == "entangled") flag = isEntangled<36>(short_channel_info,"left") && isEntangled<36>(short_channel_info,"right");
        if (state == "decoherent") flag = isDecoherent<36>(short_channel_info,"left") || isDecoherent<36>(short_channel_info,"right");
        if (state == "double") flag = isDecoherent<36>(short_channel_info,"left") && isDecoherent<36>(short_channel_info,"right");

        return flag;
    }

    // Bool_t Double_Decoherent_or_Decoherent_or_Entangled(short_energy_ChannelEntry *short_channel_info, TString state = "entangled")
    // {
    //     Bool_t flag = kFALSE;
    //     if (state == kFALSE && short_channel_info[34]->amp < 60000
    //         && short_channel_info[34]->amp > 100
    //         && (short_channel_info[34]->time - short_channel_info[32]->time) < 350
    //         && (short_channel_info[34]->time - short_channel_info[32]->time) > 150)
    //         flag = kTRUE;

    //     if (state == kTRUE &&
    //         short_channel_info[34]->amp < 20 &&
    //         ((short_channel_info[34]->time - short_channel_info[32]->time) < 50 ||
    //         (short_channel_info[34]->time - short_channel_info[32]->time) > 500)
    //         )
    //         flag = kTRUE;

            
    //     return flag;
    // }    
    Bool_t diff_cuts(Short_t min_diff, TString state = "entangled")
    {
        Bool_t flag = kFALSE;
        if (state == "decoherent" && min_diff < -40)
            flag = kTRUE;

        if (state == "entangled" && min_diff > -10)
            flag = kTRUE;
        return flag;
    }

    template<int n> 
    Bool_t Non_zero_time(std::array<short_energy_ChannelEntry*, n> &short_channel_info, Int_t channel_number)
    {
        Bool_t flag = kFALSE;
        if (
        short_channel_info[channel_number]->time > 11
        &&  short_channel_info[33]->time > 11
        &&  short_channel_info[32]->time > 11) flag = kTRUE;
        return flag;
    }
    template<int n> 
    Bool_t Apply_time_in_peak_cuts(
         std::array<short_energy_ChannelEntry*, n> &short_channel_info, Int_t channel_number, 
        Int_t sc_number = 32, Int_t left_main_sc = -9, 
        Int_t right_main_sc = 9, Int_t left_ch_sc = -300, 
        Int_t right_ch_sc = -100, TString data = "new")
    {

        Bool_t flag = kFALSE;
        if (
        short_channel_info[channel_number]->time > 11//non-zero signal choosing
        &&  short_channel_info[33]->time > 11
        &&  short_channel_info[32]->time > 11//end

        && (short_channel_info[33]->time - short_channel_info[32]->time) < right_main_sc// choosing entangled pairs
        && (short_channel_info[33]->time - short_channel_info[32]->time) > left_main_sc// choosing entangled pairs
        && (short_channel_info[sc_number]->time - short_channel_info[channel_number]->time) > left_ch_sc // choosing coincidence with NaI counters
        && (short_channel_info[sc_number]->time- short_channel_info[channel_number]->time) < right_ch_sc) {flag = kTRUE;}
        return flag;
    }
    template<int n> 
    Bool_t Apply_Amplitude_Saturation_cuts(std::array<short_energy_ChannelEntry*, n> &short_channel_info, Int_t channel_number, Int_t channel_number_2 = -999)
    {

        Bool_t flag = kFALSE;
        Bool_t flag_ch_2 = kFALSE;
        if (channel_number_2==-999)
            flag_ch_2 = kTRUE;
        else if (short_channel_info[channel_number_2]->amp > 200 
                && short_channel_info[channel_number_2]->amp < 60000) flag_ch_2 = kTRUE;

        if (short_channel_info[33]->amp < 60000
        && short_channel_info[33]->amp > 200
        && short_channel_info[32]->amp < 60000     
        && short_channel_info[32]->amp > 200 
        && short_channel_info[channel_number]->amp > 200 
        && short_channel_info[channel_number]->amp < 60000
        && short_channel_info[channel_number]->charge > 1
        && flag_ch_2) flag = kTRUE;
        return flag;
    }


    void double_gauss_fit(TH1F *peak_histo, Float_t &low_cut, Float_t &high_cut, Float_t left_range_width = 1.35, Float_t right_range_width = 1.35, Float_t left_range = 150, Float_t right_range = 350)
    {
            Float_t interval_width = 1.35;

        left_range = peak_histo->GetBinCenter(peak_histo->GetMaximumBin()) - 100;
        right_range = peak_histo->GetBinCenter(peak_histo->GetMaximumBin()) + 100;

        TF1 *gaus_fit = new TF1 ("gaus_func","gaus",left_range,right_range);
        peak_histo->Fit("gaus_func","R");

        TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
        (gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
        (gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
        peak_histo->Fit("gaus_func_2","R");

        low_cut = gaus_2_fit->GetParameter(1)-left_range_width*gaus_2_fit->GetParameter(2);
        high_cut = gaus_2_fit->GetParameter(1)+right_range_width*gaus_2_fit->GetParameter(2);
    }
}

#endif