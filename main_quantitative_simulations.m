function main_quantitative_simulations(cohort_catergory,hmo,mucin,fiber,rs)
diet_baline=[hmo(1) mucin(1) fiber(1) rs(1)];
diet_intervention=[hmo(2) mucin(2) fiber(2) rs(2)];

diet_ref=[14     1      0     0
           0     1.5   12.5   10
           0     2     33.2   23.8
           0     2     19.7   4];      


switch cohort_catergory

    case 'Infant'
        if isequal(diet_baseline,diet_ref(1,:) && isequal(diet_intervention,diet_ref(2,:)
            result_file='';    
        else   %% perform new simulations
            
            nspc=1:5;
            Result_BF=main_Infant_BFeeding(hmo,mucin,fiber,rs,nspc);

            nspc=1:8;
            Result_SF=main_Infant_SFeeding(hmo,mucin,fiber,rs,nspc,Result_BF);

            out_filename='Community_10tanks_infant_BF_SF_new.xls';
            Result_total=concate_BF_SF(Result_BF,Result_SF,out_filename);   %% out_filename is the new file name stored for infant data
            
        end
        
    case 'Adult'
        if isequal(diet_baseline,diet_ref(3,:) && isequal(diet_intervention,diet_ref(4,:)
            result_file='';    
        else   %% perform new simulations
            
            nspc=2:8;
            Result_BSL=main_adult_Baseline(hmo,mucin,fiber,rs,nspc);

            nspc=2:8;
            Result_Intv=main_adult_Intv(hmo,mucin,fiber,rs,nspc,Result_BSL);

            out_filename='community_10tank_adult_concate_profiles_new.xls';
            Result_total=concate_BSL_Intv(Result_BSL,Result_Intv,out_filename);   %% out_filename is the new file name stored for infant data
           
        end