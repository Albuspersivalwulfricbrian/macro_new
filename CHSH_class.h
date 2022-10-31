#ifndef CHSH_CLASS
#define CHSH_CLASS
class CHSH_class
{
    private:
    static const Int_t module_angles = 8;
    static const Int_t all_angles = module_angles*2;
    Float_t N_parallel;
    Float_t N_perpendicular;
    Float_t N_mixed_1;
    Float_t N_mixed_2;
    Int_t d_phi;
    TString z;
    Float_t CHSH[module_angles];
    Float_t CHSH_all[all_angles];
    Float_t CHSH_error[module_angles]; 
    Float_t CHSH_error_all[all_angles];
    Float_t global_CHSH[module_angles]; 
    Float_t global_CHSH_all[all_angles];
    Float_t global_CHSH_error[module_angles]; 
    Float_t global_CHSH_error_all[all_angles];
    Float_t global_Correlation_E[module_angles]; 
    Float_t global_Corr_error[module_angles];
    Float_t global_Correlation_E_all[all_angles];
    Float_t global_Corr_error_all[all_angles];
    Float_t total_Correlation_E[module_angles]; 
    Float_t total_Corr_error[module_angles];
    Float_t total_Correlation_E_all[all_angles]; 
    Float_t total_Corr_error_all[all_angles];
    Int_t NumEvents[all_angles][all_angles];

    public:
    CHSH_class();
    virtual ~CHSH_class();
    Int_t true_number(Int_t a);
    void count_coincidences(Int_t d_a);
    Float_t global_E_coeff(Int_t d_a);
    Float_t global_sqr_E_error(Int_t d_a);
    Float_t global_calculate_CHSH();
    Float_t global_calculate_CHSH_error();
    void Create_global_CHSH_arrays();
/////////////////////////////////////////
/////////////////////////////////////////local_CHSH_for_average
    void local_count_coincidences(Int_t a, Int_t b);
    Float_t E_coeff(Int_t a, Int_t b);
    Float_t sqr_E_error(Int_t a, Int_t b); //calculating square error of correlation coefficient
    Float_t calculate_local_CHSH(Int_t phi);
    Float_t calculate_local_CHSH_error(Int_t phi);
    void Create_average_CHSH_arrays();
    //Setters
    void SetNumEvents(Int_t arr[all_angles][all_angles]);
    //Getters
    Float_t* Get_CHSH();
    Float_t* Get_CHSH_error();
    Float_t* Get_CHSH_all();
    Float_t* Get_CHSH_error_all();
    Float_t* Get_global_CHSH();
    Float_t* Get_global_CHSH_error();
    Float_t* Get_global_CHSH_all();   
    Float_t* Get_global_CHSH_error_all();  
    Float_t* Get_global_Correlation_E();  
    Float_t* Get_global_Correlation_E_all();  
    Float_t* Get_global_Corr_error();    
    Float_t* Get_global_Corr_error_all(); 
    Float_t* Get_total_Correlation_E(); 
    Float_t* Get_total_Corr_error();  
    Float_t* Get_total_Correlation_E_all();
    Float_t* Get_total_Corr_error_all();    
};
#endif CHSH_CLASS
