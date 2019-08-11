function Result_total=run_Simulation_3d(cohort,index)

index_matlab=index;            
                 
if strcmpi(cohort,'Infant')
    Result_BF='20181024_4m_diet_302h.mat';
    Result_SF='20181029_12m_diet_732h.mat';
    Result_total=concate_4m_12m_mets_3d(Result_BF,Result_SF,index_matlab);
                
    
elseif strcmpi(cohort,'Adult')

    Result_BF='20190312_CM_diet_baseline.mat';
    Result_SF='20190319_CM_diet_intervention.mat';
    Result_total=concate_4m_12m_mets_3d_adult(Result_BF,Result_SF,index_matlab);
                 
                 
end
