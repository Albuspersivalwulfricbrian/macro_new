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
#include "CHSH_class.h"
#include "Double_scattering_cuts.h"

using namespace std;
//using namespace CHSH;
using namespace CUTS;

#define UseManyRoots 0
#define PresentIntermediate 1
#define CalculateRatio 1
#define DrawTime 1
#define UseIntegralCut 1
#define UseNotEntangledPhotons 0
#define DrawToPDF 0
#define calculate_CHSH 0

/////////////////////////////////////////////////////
void try_this_shit()
{
   
	//TString source_path = "/home/doc/entanglement/root_files_data/with_scatterer/big_file/";
	TString source_path = "/home/doc/entanglement/double_GAGG/";

#if UseNotEntangledPhotons
    TString result_path  = source_path + "decoh_mini_tree";
    const TString entangled = "entangled";
#else
    TString result_path = source_path + "entangled_mini_tree";
    const Bool_t entangled = "entangled";
#endif

///////////////////////////////////////////////
    TChain *PMT_tree = new TChain;
	PMT_tree->AddFile( (source_path + "calibrated_time_0.root" + "/adc64_data").Data() );
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\////////////////////////
    const Int_t calculate_Events_Number = PMT_tree->GetEntries();
    Int_t Events_for_cuts = calculate_Events_Number/1;
///////////////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\//////////////////////////
    const Int_t total_channels = 36;

    std::array<short_energy_ChannelEntry*, total_channels> short_channel_info;

    for(Int_t ch = 0; ch < total_channels; ch++)
    {
        short_channel_info[ch] = new short_energy_ChannelEntry();
        short_channel_info[ch]->Initialize();
                
        // cout << &(short_channel_info[ch]) << " \n";

        PMT_tree->SetBranchAddress((TString::Format("channel_%i", ch)).Data(), &(short_channel_info[ch]));
    }


}

