function write_tanks_cmp(Res)
Y=Res.SimulationResult;      
Ymet=Res.SimulationMetabolites;
Ymodel=Res.Result_cell;        
ylg=Res.length_tank;
lxini=length(Ymet{1}(1,:));    
nspc=Res.SpeciesNumber;
met_co=Res.met_udf_co;
Time=Res.SimulationTime;
rct=Res.RectorNumber;

header={'Species' 'Ascending' 'Transverse' 'Descending' 'Sigmoid' 'Ascending_M' 'Transverse_M' 'Descending_M' 'Sigmoid_M' 'Rectum'}';
header_mets={'Metabolites' 'Ascending' 'Transverse' 'Descending' 'Sigmoid' 'Ascending_M' 'Transverse_M' 'Descending_M' 'Sigmoid_M' 'Rectum' 'Feces' 'Blood'}';

colnames_bac=met_co(1:nspc)';
colnames_bac=strrep(colnames_bac,'BIOMASS_','');  % remove prefix of Biomass_, as well as the blank
colnames_bac=strrep(colnames_bac,' ','_');
colnames_met=met_co(nspc+1:lxini)';
% conc_bac=colnames_bac;
% conc_mets=colnames_met;
conc_bac=[];
conc_mets=[];

for i=1:rct-1
    conc_bac=[conc_bac;Ymet{i}(end,1:nspc)];
    conc_mets=[conc_mets;Ymet{i}(end,nspc+1:lxini)];   
end
i=9;
conc_mets_feces=Y(end,(i-1)*ylg+ylg+1+1:(i-1).*ylg+ylg+1+lxini);   
if nspc==5   
    coef=20;
else
    coef=6;
end
conc_mets_feces=conc_mets_feces(nspc+1:lxini).*coef;              
i=10;  % add blood metabolites
conc_mets_blood=Ymodel{i}(end,:);                
conc_mets=[conc_mets;conc_mets_feces;conc_mets_blood];
% conc_mets=[conc_mets;num2cell(Ymet{i}(end,:))];    
conc_mets=[colnames_met;num2cell(conc_mets)];

conc_bac=[colnames_bac;num2cell(conc_bac)];
output_bac=[header conc_bac];
output_mets=[header_mets conc_mets];
% cur_dir=pwd;
% dir=[cur_dir '/tank_cmp_data'];

if nspc==5
%     filename_bac='/bacteria_4m_1102.xlsx';
%     filename_mets='/mets_4m_1102.xlsx';
    filename_bac='/Microbial_BF.xlsx';
    filename_mets='/Metabolite_BF.xlsx';
else
%     filename_bac='/bacteria_12m_1102.xlsx';
%     filename_mets='/mets_12m_1102.xlsx';
    filename_bac='/Microbial_SF.xlsx';
    filename_mets='/Metabolite_SF.xlsx';

end

writetable(cell2table(output_bac),[dir filename_bac],'WriteVariableNames',false);  
writetable(cell2table(output_mets),[dir filename_mets],'WriteVariableNames',false);  

end