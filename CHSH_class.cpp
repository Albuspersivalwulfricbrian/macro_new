#include <TString.h>
#include <iostream>

#include "CHSH_class.h"

    CHSH_class::CHSH_class()
    {
        N_parallel = 0;
        N_perpendicular = 0;
        N_mixed_1 = 0;
        N_mixed_2 = 0;
        z = "clockwise";


        for (Int_t i = 0; i < module_angles; i++)
        {
            global_CHSH[i] = 0;
            CHSH_error[i] = 0; 
            CHSH[i] = 0;
            global_CHSH_error[i] = 0; 
            global_Correlation_E[i] = 0; 
            global_Corr_error[i] = 0;
            total_Corr_error[i] = 0;
            total_Correlation_E[i] = 0;             
        }
        for (Int_t i = 0; i < all_angles; i++)
        {
            CHSH_all[i] = 0;
            total_Correlation_E_all[i] = 0; 
            total_Corr_error_all[i] = 0;
            CHSH_error_all[i] = 0;
            global_Correlation_E_all[i] = 0;
            global_CHSH_error_all[i] = 0;
            global_Corr_error_all[i] = 0;
            global_CHSH_all[i] = 0;            
        }
    }
    CHSH_class::~CHSH_class() 
    {
        cout << "CHSH deleted" << endl;
    }

    Int_t CHSH_class::true_number(Int_t a)
    {
        a+=48;
        while(a > 15) a = a-16;
        return a;
    }
    
    void CHSH_class::count_coincidences(Int_t d_a)
    {
        if (z == "counterclockwise") d_a = -d_a;

        N_parallel=0;
        N_perpendicular = 0;
        N_mixed_1 = 0;
        N_mixed_2 = 0;  

        for (Int_t a = 0; a < 16; a++)
        {
            Int_t b = a+d_a+48;
            Int_t a_1 = a+4;
            Int_t b_1 = b+4;
            a=true_number(a);
            b=true_number(b);
            a_1=true_number(a_1);
            b_1=true_number(b_1);   
            N_parallel += NumEvents[a][b];
            N_perpendicular += NumEvents[a_1][b_1];
            N_mixed_1 += NumEvents[a][b_1];
            N_mixed_2 += NumEvents[a_1][b];        
        }
        if (z=="both")
        {
            for (Int_t a = 0; a < 16; a++)
            {
                Int_t b = a-d_a+48;
                Int_t a_1 = a+4;
                Int_t b_1 = b+4;
                a=true_number(a);
                b=true_number(b);
                a_1=true_number(a_1);
                b_1=true_number(b_1);  
                N_parallel += NumEvents[a][b];
                N_perpendicular += NumEvents[a_1][b_1];
                N_mixed_1 += NumEvents[a][b_1];
                N_mixed_2 += NumEvents[a_1][b];        
            }
        }
    }

    Float_t CHSH_class::global_E_coeff(Int_t d_a)
    {
        count_coincidences(d_a);
        return (float)(N_parallel+N_perpendicular-N_mixed_1-N_mixed_2)/(N_parallel+N_perpendicular+N_mixed_1+N_mixed_2);
    }
    Float_t CHSH_class::global_sqr_E_error(Int_t d_a)
    {
        count_coincidences(d_a);

        return (float)4.0/pow((N_parallel+N_perpendicular+N_mixed_1+N_mixed_2),4)*(pow(N_mixed_1+N_mixed_2,2)*(N_parallel+N_perpendicular)+pow(N_parallel+N_perpendicular,2)*(N_mixed_1+N_mixed_2));
    }
    Float_t CHSH_class::global_calculate_CHSH()
    {
        return 3*global_E_coeff(d_phi) - global_E_coeff(3*d_phi);
    }

    Float_t CHSH_class::global_calculate_CHSH_error()
    {
        return global_sqr_E_error(d_phi)*3 + global_sqr_E_error(3*d_phi);
    }

    void CHSH_class::Create_global_CHSH_arrays()
    {
        for (Int_t abc = 1; abc <= module_angles; abc++)
        {
            d_phi = abc;
            z = "clockwise";
            total_Correlation_E_all[abc-1] = global_E_coeff(abc);
            total_Corr_error_all[abc-1] = sqrt(global_sqr_E_error(abc));
            global_CHSH_all[abc-1] = global_calculate_CHSH();
            global_CHSH_error_all[abc-1] = sqrt(global_calculate_CHSH_error());
            z = "counterclockwise";
            total_Correlation_E_all[module_angles +abc-1] = global_E_coeff(abc);
            total_Corr_error_all[module_angles + abc-1] = sqrt(global_sqr_E_error(abc));
            global_CHSH_all[module_angles+abc-1] = global_calculate_CHSH();
            global_CHSH_error_all[module_angles+abc-1] = sqrt(global_calculate_CHSH_error());
            z= "both";
            total_Correlation_E[abc-1] = global_E_coeff(abc);
            total_Corr_error[abc-1] = sqrt(global_sqr_E_error(abc));
            global_CHSH[abc-1] = global_calculate_CHSH();
            global_CHSH_error[abc-1] = sqrt(global_calculate_CHSH_error());
        }
    }

/////////////////////////////////////////
/////////////////////////////////////////local_CHSH_for_average
    void CHSH_class::local_count_coincidences(Int_t a, Int_t b)
    {
        a += 48;
        b += 48;
        Int_t a_1 = a+4;
        Int_t b_1 = b+4;
        a=true_number(a);
        b=true_number(b);
        a_1=true_number(a_1);
        b_1=true_number(b_1);
        N_parallel=NumEvents[a][b];
        N_perpendicular = NumEvents[a_1][b_1];
        N_mixed_1 = NumEvents[a][b_1];
        N_mixed_2 = NumEvents[a_1][b];
    }

    Float_t CHSH_class::E_coeff(Int_t a, Int_t b)
    {
        local_count_coincidences(a,b);
        return (float)(N_parallel+N_perpendicular-N_mixed_1-N_mixed_2)/(float)(N_parallel+N_perpendicular+N_mixed_1+N_mixed_2);
    }

    Float_t CHSH_class::sqr_E_error(Int_t a, Int_t b) //calculating square error of correlation coefficient
    {
        local_count_coincidences(a,b);
        return 
        (float)4.0/pow((N_parallel+N_perpendicular+N_mixed_1+N_mixed_2),4)*(pow(N_mixed_1+N_mixed_2,2)*(N_parallel+N_perpendicular)+pow(N_parallel+N_perpendicular,2)*(N_mixed_1+N_mixed_2));
    }

    Float_t CHSH_class::calculate_local_CHSH(Int_t phi)
    {
        return
        E_coeff(phi,phi+d_phi) -
        E_coeff(phi,phi+3*d_phi) +
        E_coeff(phi+2*d_phi,phi+d_phi) +
        E_coeff(phi+2*d_phi,phi+3*d_phi);
    }

    Float_t CHSH_class::calculate_local_CHSH_error(Int_t phi)
    {
        //if (z == "counterclockwise") d_phi = -d_phi;
        return
        sqr_E_error(phi,phi+d_phi) +
        sqr_E_error(phi,phi+3*d_phi) +
        sqr_E_error(phi+2*d_phi,phi+d_phi) +
        sqr_E_error(phi+2*d_phi,phi+3*d_phi);
    }

    void CHSH_class::Create_average_CHSH_arrays()
    {
        for (Int_t a_angle = 0; a_angle < 16; a_angle++)
        {
            Float_t CHSH_error_local[all_angles] = {0};  Float_t CHSH_error_local_sum[module_angles] = {0};
            Float_t CHSH_local[all_angles] = {0}; Float_t CHSH_local_sum[module_angles] = {0};
            for (Int_t abc = 1; abc <= module_angles; abc++)
            {
                d_phi = abc;
                global_Correlation_E_all[abc-1] += E_coeff(a_angle, a_angle+d_phi);
                global_Corr_error_all[abc-1] += sqr_E_error(a_angle, a_angle+d_phi);                
                CHSH_local[abc-1] = calculate_local_CHSH(a_angle);
                CHSH_error_local[abc-1] = calculate_local_CHSH_error(a_angle);

                d_phi = -abc;
                global_Correlation_E_all[abc+module_angles-1] += E_coeff(a_angle,a_angle+d_phi);                    
                global_Corr_error_all[abc+module_angles-1] += sqr_E_error(a_angle,a_angle+d_phi);
                CHSH_local[module_angles + abc-1] = calculate_local_CHSH(a_angle);
                CHSH_error_local[module_angles + abc-1] = calculate_local_CHSH_error(a_angle);

                CHSH_local_sum[abc-1] = (CHSH_local[abc-1] + CHSH_local[abc+module_angles-1])/2.;
                CHSH_error_local_sum[abc-1] = CHSH_error_local[abc-1] + CHSH_error_local[abc+module_angles-1];
                CHSH[abc-1] += CHSH_local_sum[abc-1];
                CHSH_error[abc-1] += CHSH_error_local_sum[abc-1];       
                CHSH_error_local_sum[abc-1] = sqrt(CHSH_error_local_sum[abc-1])/2;  
                CHSH_all[abc-1] += CHSH_local[abc-1];
                CHSH_all[abc + module_angles-1] += CHSH_local[abc + module_angles-1];
                CHSH_error_all[abc-1] += CHSH_error_local[abc-1];  
                CHSH_error_all[abc + module_angles-1] += CHSH_error_local[abc + module_angles-1];  
                CHSH_error_local[abc-1] = sqrt(CHSH_error_local[abc-1]); 
                CHSH_error_local[abc + module_angles-1] = sqrt(CHSH_error_local[abc + module_angles-1]); 
            }

                // TCanvas *canvas_CHSH = new TCanvas ( Form("canvas_CHSH_a=%i",a_angle), Form("canvas_CHSH_a=%i",a_angle));
                // TGraphErrors * gr_CHSH = new TGraphErrors(module_angles,phi_angle,CHSH_local_sum, angle_arr_err, CHSH_error_local_sum);
                // //TF1 *my_fit_2 = new TF1 ("my_fit_2","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",10,190);
                // canvas_CHSH->Divide(2);

                // gr_CHSH->Fit("my_fit","R");
                // Float_t a_par_CHSH = abs(my_fit->GetParameter(0));
                // Float_t a_error_CHSH = abs(my_fit->GetParError(0));
                // graph_like_ivashkin_wants_it(gr_CHSH,"angle [degrees]","S", 
                // Form("canvas_CHSH_a=%i <> E(90) = %4.3f+-%4.3f",a_angle,
                // (E_coeff(NumEvents,a_angle,true_number(a_angle+4))+E_coeff(NumEvents,a_angle,true_number(a_angle-4)))/2,
                // sqrt(sqr_E_error(NumEvents,a_angle,true_number(a_angle+4))+sqr_E_error(NumEvents,a_angle,true_number(a_angle-4)))/2),1);
                // canvas_CHSH->cd(1);
                // gr_CHSH->Draw("AP");

                // //TF1 *my_fit = new TF1 ("my_fit","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",-190,190);
                // TGraphErrors * gr_CHSH_all = new TGraphErrors(all_angles-1,phi_angle,CHSH_local, angle_arr_err, CHSH_error_local);
                // gr_CHSH_all->Fit("my_fit","R");
                // a_par_CHSH = abs(my_fit->GetParameter(0));
                // a_error_CHSH = abs(my_fit->GetParError(0));
                // graph_like_ivashkin_wants_it(gr_CHSH_all,"angle [degrees]","S", Form("canvas_CHSH_a=%i <> E(90) = %4.3f+-%4.3f",a_angle,E_coeff(NumEvents,a_angle,true_number(a_angle+4)),sqrt(sqr_E_error(NumEvents,a_angle,true_number(a_angle+4)))),1);
                // canvas_CHSH->cd(2);
                // gr_CHSH_all->Draw("AP");
                // gr_CHSH_all->Write();
                // canvas_CHSH->SaveAs(result_path+".pdf",".pdf");
                // delete canvas_CHSH;
        }

        for (Int_t i = 0; i < all_angles; i++) 
        {

            if ( i < module_angles)
            {
                CHSH[i] = CHSH[i]/16;
                cout << CHSH[i] << " ";
                CHSH_error[i] = sqrt(CHSH_error[i])/32;
                global_Correlation_E[i] = (global_Correlation_E_all[i] + global_Correlation_E_all[i+module_angles])/32;
                global_Corr_error[i] = sqrt(global_Corr_error_all[i]+ global_Corr_error_all[i+module_angles])/32;
            }
            cout << endl;
            CHSH_all[i] = CHSH_all[i]/16;
            CHSH_error_all[i] = sqrt(CHSH_error_all[i])/16;
            global_Corr_error_all[i] = sqrt(global_Corr_error_all[i])/16;
            global_Correlation_E_all[i] = global_Correlation_E_all[i]/16;

        }        
    }


    //Setters
    void CHSH_class::SetNumEvents(Int_t arr[16][16])
    {
        for (Int_t i = 0; i < 16; i++) for (Int_t j = 0; j < 16; j++) NumEvents[i][j] = arr[i][j];
    }

    //Getters

    Float_t* CHSH_class::Get_CHSH()
    {
        return CHSH;
    }
    Float_t* CHSH_class::Get_CHSH_error()
    {
        return CHSH_error;
    }
    Float_t* CHSH_class::Get_CHSH_all()
    {
        return CHSH_all;
    }
    Float_t* CHSH_class::Get_CHSH_error_all()
    {
        return CHSH_error_all;
    }
    Float_t* CHSH_class::Get_global_CHSH()
    {
        return global_CHSH;
    }
    Float_t* CHSH_class::Get_global_CHSH_error()
    {
        return global_CHSH_error;
    }
    Float_t* CHSH_class::Get_global_CHSH_all()
    {
        return global_CHSH_all;
    }    
    Float_t* CHSH_class::Get_global_CHSH_error_all()
    {
        return global_CHSH_error_all;
    }    
    Float_t* CHSH_class::Get_global_Correlation_E()
    {
        return global_Correlation_E;
    }    
    Float_t* CHSH_class::Get_global_Correlation_E_all()
    {
        return global_Correlation_E_all;
    }    
    Float_t* CHSH_class::Get_global_Corr_error()
    {
        return global_Corr_error;
    }    
    Float_t* CHSH_class::Get_global_Corr_error_all()
    {
        return global_Corr_error_all;
    }    
    Float_t* CHSH_class::Get_total_Correlation_E()
    {
        return total_Correlation_E;
    }    
    Float_t* CHSH_class::Get_total_Corr_error()
    {
        return total_Corr_error;
    }    
    Float_t* CHSH_class::Get_total_Correlation_E_all()
    {
        return total_Correlation_E_all;
    }
    Float_t* CHSH_class::Get_total_Corr_error_all()
    {
        return total_Corr_error_all;
    } 