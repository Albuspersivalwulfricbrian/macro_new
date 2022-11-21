#ifndef CHANNEL_ENTRY_H
#define CHANNEL_ENTRY_H
#include<TTree.h>

const int MAX_N_SAMPLES = 2048;

struct IntegralInfo
{
    Short_t signal_length;
    Short_t npeaks;
    Int_t end_amplitude;
    void Initialize();
};

struct short_energy_ChannelEntry
{
    Float_t charge;
    Float_t time;
    UShort_t amp; 
    Float_t zl;
    Float_t zl_rms;
    IntegralInfo II;

    static TString GetChName(Int_t channel_num);
    TBranch* CreateBranch(TTree *tree, Int_t channel_num);
    Int_t SetBranch(TTree *tree, Int_t channel_num);
    void Initialize();
};

struct diff_short_energy_ChannelEntry
{
    Short_t min_diff;
    Short_t min_diff_time;
    Short_t max_diff;
    Short_t max_diff_time;
    static TString GetChName(Int_t channel_num);
    TBranch* CreateBranch(TTree *tree, Int_t channel_num);
    // Int_t SetBranch(TTree *tree, Int_t channel_num);
    void Initialize();
};

struct ChannelEntry {

    public:
    Short_t wf_size;
    Short_t wf[MAX_N_SAMPLES];

    private:
    Float_t zl;
    IntegralInfo II;
    public:
    static TString GetChName(Int_t channel_num);
    Int_t SetBranch(TTree *tree, Int_t channel_num);
    
    void Initialize();
    void SplineWf();
    void DiffWf();
    //void  CalculateC();
    void FindDiffWfPars(Short_t &min_diff, Short_t &min_time, Short_t &max_diff, Short_t &max_time, Int_t GATE_BEG = 55, Int_t GATE_END = 76);
    Int_t Get_Zero_Level(Int_t GATE_BEG);
    Float_t Get_Zero_Level_RMS(Int_t GATE_BEG);
    Float_t Get_Charge(Int_t GATE_BEG, Int_t GATE_END);
    Short_t Get_time();
    Float_t Get_time_gauss(Int_t inv_amp, Int_t CH_GATE_END = -1000);
    UShort_t Get_Amplitude(Int_t GATE_BEG = 0, Int_t CH_GATE_END = -1000);
    Short_t GetSLength();
    Int_t GetEndAmp();

};

//##################
//#################


struct mini_tree_nrg
{
    Float_t EdepIntermediate0;
    Float_t EdepIntermediate1;
    Float_t EdepScat0;
    Float_t EdepScat1;
    Float_t EdepDet0;
    Float_t EdepDet1;  
    Short_t DetNum0;
    Short_t DetNum1;
    Short_t EventType; 

    // static TString GetBrName();
    // TBranch* CreateBranches(TTree *tree);
    Int_t SetBranch(TTree *tree);
    Int_t Initialize();
};

struct mini_tree_time
{
    Float_t TimeIntermediate0;
    Float_t TimeIntermediate1;

    Float_t TimeScat0;
    Float_t TimeScat1;
    Float_t TimeDet0;
    Float_t TimeDet1;

    // static TString GetBrName();
    // TBranch* CreateBranches(TTree *tree);
    Int_t SetBranch(TTree *tree);
    Int_t Initialize();
};

#endif CHANNEL_ENTRY_H
