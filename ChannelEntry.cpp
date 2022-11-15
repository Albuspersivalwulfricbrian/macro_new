
#include <TTree.h>
#include "ChannelEntry.h"
#include <iostream>

using namespace std;

    void IntegralInfo::Initialize()
    {
        signal_length = 0;
        npeaks = 0;
        end_amplitude = 0;    
    }
    void short_energy_ChannelEntry::Initialize()
    {
        charge = 0.;
        time = 0.;
        amp = 0; 
        zl_rms = 0.;
        zl = 0.;
        II.Initialize();
    }


    TString short_energy_ChannelEntry::GetChName(Int_t channel_num)
    {
	    return TString::Format("channel_%i", channel_num);
    }

    Int_t short_energy_ChannelEntry::SetBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
    
    TString diff_short_energy_ChannelEntry::GetChName(Int_t channel_num)
    {
	    return TString::Format("diff_channel_%i", channel_num);
    }
    TBranch* diff_short_energy_ChannelEntry::CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "min_diff/S:min_diff_time/S:max_diff/S:max_diff_time/S");
    }

    void diff_short_energy_ChannelEntry::Initialize()
    {
        min_diff = 0;
        min_diff_time = 0;
        max_diff = 0;
        max_diff_time = 0;
    }    


    TString ChannelEntry::GetChName(Int_t channel_num)
    {
	    return TString::Format("channel_%i", channel_num);
    }

    Int_t ChannelEntry::SetBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
    
    void ChannelEntry::Initialize()
    {
        for (int i = 0; i < sizeof(wf)/sizeof(wf[0]); i++) wf[i] = 0;
        wf_size = 0;
    }

    void ChannelEntry::SplineWf()
    {
        Float_t wf1[MAX_N_SAMPLES] = {0};
        const Int_t SplineWidth = 4;
        for (Int_t i = 0; i < wf_size; i++)
        {
            Int_t il=i-SplineWidth; Int_t ir=i+SplineWidth;
            if (il<0) il=0;
            if (ir>wf_size-1) ir=wf_size-1;
            Float_t counter = 0;
            for (Int_t in = il; in <=ir; in++) {wf1[i]+=wf[in];counter++;}
            wf1[i]/=counter;
        }
        for (Int_t i = 0; i < wf_size; i++) wf[i] = wf1[i];
    }

    void ChannelEntry::DiffWf()
    {
        const Float_t Diff_window = 4;
        Short_t wf1[MAX_N_SAMPLES] = {0};
        for (Int_t i = 0; i < wf_size; i++)
        {
            Int_t il=i-Diff_window; Int_t ir=i+Diff_window;
            if (il<0) il=0;
            if (ir>wf_size-1) ir=wf_size-1;
            wf1[i]=(Short_t)((Float_t)(wf[ir]-wf[il])/(Float_t)(ir-il));
        }
        for (Int_t i = 0; i < wf_size; i++) wf[i] = wf1[i];
    }    

    void ChannelEntry::FindDiffWfPars(Short_t &min_diff, Short_t &min_time, Short_t &max_diff, Short_t &max_time, Int_t GATE_BEG, Int_t GATE_END)
    {
        for (Short_t s=GATE_BEG; s < GATE_END; ++s) {
            Short_t v = wf[s];
            if (v < min_diff) {
                min_diff = v;
                min_time = 16*s;
            }
            if (v > max_diff) {
                max_diff = v;
                max_time = 16*s;
            }      
        }
    }


    Int_t ChannelEntry::Get_Zero_Level(Int_t GATE_BEG)
    {
        const Int_t interv_num = 1;
        int zero_lvl = 0;
        int best_spread = -1;
        for (int i=0; i < interv_num; ++i) {
            int vmin = numeric_limits<int>::max();
            int vmax = numeric_limits<int>::min();
            int sum = 0;
            for (int s=GATE_BEG/interv_num * i; s < GATE_BEG/interv_num * (i+1); ++s) {
                int v = wf[s]; 
                sum += v;
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
            int spread = vmax - vmin;
            if (best_spread < 0) best_spread = spread;
            if (spread <= best_spread) {
                best_spread = spread;
                zero_lvl = sum / (GATE_BEG/interv_num);
            }
            zl = zero_lvl;
        }
        return zero_lvl;
    }
    Float_t ChannelEntry::Get_Zero_Level_RMS(Int_t GATE_BEG)
    {
        const Int_t interv_num = 1;
        Float_t best_spread = -1;
        Float_t rms_zl = -1;
        for (Int_t i=0; i < interv_num; ++i) {
            Int_t vmin = numeric_limits<int>::max();
            Int_t vmax = numeric_limits<int>::min();
            Float_t sum = 0; Float_t sumsquare = 0; Float_t sum_counter = 0;
            for (Int_t s=GATE_BEG/interv_num * i; s < GATE_BEG/interv_num * (i+1); ++s) {
                Int_t v = wf[s]; 
                sum += (Float_t)v;
                sum_counter++;
            }
            sum /=sum_counter;
            sumsquare = 0.;
            
            for (Int_t s=GATE_BEG/interv_num * i; s < GATE_BEG/interv_num * (i+1); ++s) {
                sumsquare += (Float_t)(wf[s] - sum)*(wf[s] - sum)/sum_counter;
            }            
            rms_zl = sqrt(sumsquare);
            if (best_spread < 0) best_spread = rms_zl;

        }
        return best_spread;
    }

    Float_t ChannelEntry::Get_Charge(Int_t GATE_BEG, Int_t GATE_END)
    {
        Float_t gateInteg = 0;
        for (int s=GATE_BEG; s < GATE_END+1; ++s) {
            if (((float)zl - (float)wf[s]) < 10 && s > GATE_BEG + 20) {II.signal_length = s - GATE_BEG; II.end_amplitude = (float)zl - (float)wf[s];break;}
            gateInteg +=  (zl - (Float_t)wf[s]) ;
        }
        return gateInteg;
    }
    Short_t ChannelEntry::GetSLength()
    {
        return II.signal_length;
    }
    Int_t ChannelEntry::GetEndAmp()
    {
        return II.end_amplitude;
    }


    Short_t ChannelEntry::Get_time()
    {
        Int_t amp = 0;
        Short_t peakPos = 0;
        for (int s=0; s < wf_size; ++s) {
            int v = wf[s] - zl;
            if (v < amp) {
                amp = v;
                peakPos = s;
            }
        }
        return peakPos;
    }
    Float_t ChannelEntry::Get_time_gauss(Int_t inv_amp, Int_t CH_GATE_END)
    {
        if ( CH_GATE_END ==-1000) CH_GATE_END = wf_size;
        if (wf_size < CH_GATE_END) CH_GATE_END = wf_size;
        if (wf_size == 0) return 0;
        Float_t peak_search = 0.;
        Float_t ampl_sum = 0;
        for (Int_t s= 0; s < CH_GATE_END; ++s) {
            Int_t v = zl - wf[s];
            if (v > inv_amp*0.1)
            {
                ampl_sum += (Float_t) v;
                peak_search+= (Float_t) v*s;
            }
        }
        peak_search /= ampl_sum;
    return 16.0*peak_search;
    }
////////////
    UShort_t ChannelEntry::Get_Amplitude(Int_t GATE_BEG, Int_t CH_GATE_END)
    {
        if ( CH_GATE_END ==-1000) CH_GATE_END = wf_size;
        if (wf_size == 0) {return 0;}
        int s1 = 0;
        Int_t amp = numeric_limits<UShort_t>::min();
        for (int s=GATE_BEG; s < CH_GATE_END; ++s) {
            Int_t v =  (Int_t)zl - (Int_t)wf[s];
            if (v > amp) {amp = v; s1=s;}
        }
        //cout << zl << " " << wf[s1] << " " << amp << endl;
        return (UShort_t)amp;
    }

//##################

    // TString mini_tree_nrg::GetBrName()
    // {
	//     return TString::Format("MiniTree");
    // }
    // TBranch* mini_tree_nrg::CreateBranches(TTree *tree)
    // {
	//     return tree->Branch(GetBrName().Data(), this, "EdepIntermediate/F:EdepScat0/F:EdepScat1/F:EdepDet0/F:EdepDet1/F:DetNum0/S:DetNum1/S");
    // }
    Int_t mini_tree_nrg::Initialize()
    {
        EdepIntermediate0 = 0;
        EdepIntermediate0 = 0;

        EdepScat0 = 0;
        EdepScat1 = 0;
        EdepDet0 = 0;
        EdepDet1 = 0;  
        DetNum0 = 0;
        DetNum1 = 0; 
        EventType = -10;
        return 0;
    }

    // TString mini_tree_time::GetBrName()
    // {
	//     return TString::Format("TimeTree");
    // }

    // TBranch* mini_tree_time::CreateBranches(TTree *tree)
    // {
	//     return tree->Branch(GetBrName().Data(), this");
    // }

    Int_t mini_tree_time::Initialize()
    {
        TimeScat0 = 0;
        TimeScat1 = 0;
        TimeDet0 = 0;
        TimeDet1 = 0;
        TimeIntermediate0 = 0;
        TimeIntermediate1 = 0;

        return 0;
    }
