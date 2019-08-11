function Result_total=run_Simulation(diet,cohort) 
diet_BF=diet(1,:);
diet_SF=diet(2,:);
Tanks={'Lumen_1','Lumen_2','Lumen_3','Lumen_4','Mucosa_1','Mucosa_2','Mucosa_3','Mucosa_4','Blood','Feces'};

if strcmpi(cohort,'Infant')
    if exist('Preidiction_RA_exp_infant.csv', 'file')==2
            delete('Preidiction_RA_exp_infant.csv');
        end
    for i=1:length(Tanks)
        outfile=strcat('Community_10tanks_infant_BF_SF_new_',Tanks{i},'.csv');
        if exist(outfile, 'file')==2
            delete(outfile);
        end
    end
    nspc=1:5;
    Result_BF=main_Infant_BFeeding(diet_BF,nspc);  
    nspc=1:8;
    Result_SF=main_Infant_SFeeding(diet_SF,nspc,Result_BF);
    out_filename='Community_10tanks_infant_BF_SF_new';
    Result_total=concate_BF_SF(diet,Result_BF,Result_SF,out_filename);
                       
elseif strcmpi(cohort,'Adult')
     if exist('prediction_format.csv', 'file')==2
            delete('prediction_format.csv');
        end
    for i=1:length(Tanks)
        outfile=strcat('community_10tank_adult_concate_profiles_',Tanks{i},'.csv');
        if exist(outfile, 'file')==2
            delete(outfile);
        end
    end
    nspc=2:8;
    Result_BF=main_CM_adult_Baseline(diet_BF,nspc);
    nspc=2:8;
    Result_SF=main_CM_adult_Intervention(diet_SF,nspc,Result_BF);
    out_filename='community_10tank_adult_concate_profiles';
    Result_total=concate_BSL_Intv(diet,Result_BF,Result_SF,out_filename);
                        
end
fprintf('\n\n');
disp('Thanks for waiting: Simulation Finished!');
fprintf('\n');

