function  Result_total=concate_BSL_Intv(diet,Result_BF,Result_SF,out_filename)
% read 4m results
% Res_4m=load('20190312_CM_diet_baseline.mat');
% Res_4m=Res_4m.ReferenceStates;
Res_4m=Result_BF;
% out_filename='community_10tank_adult_concate_profiles';

% Res_4m=load('20190312_CM_diet_baseline.mat');
% Res_4m=Res_4m.ReferenceStates;
Time_4m=Res_4m.SimulationTime;
Y_4m=Res_4m.SimulationResult;      %% vector
Ymet_4m=Res_4m.SimulationMetabolites;%%
Ymodel_4m=Res_4m.Result_cell;    

% Res_12m=load('20190319_CM_diet_intervention.mat');
% Res_12m=Res_12m.ReferenceStates;
Res_12m=Result_SF;
Time_12m=Res_12m.SimulationTime;
Y_12m=Res_12m.SimulationResult;      %% vector
Ymet_12m=Res_12m.SimulationMetabolites;%% 
Ymodel_12m=Res_12m.Result_cell;  
met_co_12m=Res_12m.met_udf_co;


ylg=Res_12m.ylg;
nspc=Res_12m.SpeciesNumber;
met_udf_co=Res_12m.met_udf_co;
lxini=length(met_udf_co);
n_rct=Res_12m.RectorNumber;

SmpPnt=4;

Time_str=find((Time_4m-600)>0.01); 
Time_str_index=Time_str(1);
Time_4m(Time_str_index);
Tn_4=Time_str_index;
Time_4m_new=Time_4m;
New_Time_4m=Time_4m_new(1:SmpPnt:Tn_4);


T=[Time_4m(1:Tn_4-1,:);Time_12m+600];
Y=[Y_4m(1:Tn_4-1,:);Y_12m];

Tn=length(T);
NewT=T(1:SmpPnt:Tn);
newY=Y(1:SmpPnt:Tn,:);
Concate_Result=[NewT newY];
T=Concate_Result(:,1);
Y=Concate_Result(:,2:end);
Result_total=[T,Y];

Y(Y<0)=eps;

for k=1:n_rct-2
    Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+(k-1).*ylg); 
    Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg); 
end
k=9;
Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+1+lxini+(k-1).*ylg);
Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg);
V5_rctm=Y(:,(k-1)*ylg+ylg+1);  
Ymodel_mass_rctm=Y(:,(k-1)*ylg+ylg+1+1:(k-1).*ylg+ylg+1+lxini);  
bls_index=find(T==New_Time_4m(length(New_Time_4m)-1));

AC_baseline=Ymodel_mass_rctm(bls_index,1:nspc);
AC_intervention=Ymodel_mass_rctm(end,1:nspc);

k=10;
Ymodel{k}=Y(:,(k-1).*ylg+1+lxini+1:(k-1).*ylg+1+lxini+lxini-nspc);   

lg=ylg;
met_udf_co_new=met_udf_co;
met_udf_co_new(1:nspc)=strrep(met_udf_co(1:nspc),'BIOMASS_','');  
met_udf_co_new(1:nspc)=strrep(met_udf_co_new(1:nspc),' ','_');
write_tank_table_adult_csv([T Y],Ymodel,Ymodel_met,met_udf_co_new,lg,lxini,n_rct,nspc,out_filename);


Result_total=[T Y];
RA_file='Mean_RA_species.txt';
RA_data=readtable(RA_file);
Fpr_dup=find(contains(RA_data.Var1,'Faecalibacterium'));
Fpr_dup_data=RA_data(Fpr_dup,:);
Fpr_RA=sum(Fpr_dup_data{:,2:end},1);
RA_data(Fpr_dup,:)=[];
Fpr_new=RA_data(Fpr_dup(1),:);
Fpr_new.Var1='Faecalibacterium prausnitzii';
Fpr_new(1,2:end)=num2cell(Fpr_RA);
RA_data=[RA_data;Fpr_new];
RA_data= [RA_data(:,1) RA_data(:,end) RA_data(:,2:end-1)];
header=table2cell(RA_data(1,2:end));
header=cellfun(@string,header,'un',0);
header=cellfun(@char,header,'un',0);
header=strcat('day_',header);
RA_data.Properties.VariableNames=['Species' header];
RA_data(1,:)=[];
species={'Bacteroides fragilis','Bifidobacterium longum','Bifidobacterium breve', 'Bifidobacterium adolescentis' 'Eubacterium hallii' 'Faecalibacterium prausnitzii' 'Roseburia intestinalis'}';
index=[];
for i=1:length(species)
    index(i)=find(contains(RA_data.Species,species(i)));
end
RA_data=RA_data(index,:);
Exp_baseline=RA_data.day_0;
Exp_intervention=RA_data.day_14;



% RA_4m_pre=AC_baseline./sum(AC_baseline).*sum(Exp_baseline);
% RA_12m_pre=AC_intervention./sum(AC_intervention).*sum(Exp_intervention);



diet_sum=sum(diet,2);
RA_4m_pre=AC_baseline./sum(AC_baseline).*sum(Exp_baseline)*(sqrt(diet_sum(1)/59));
RA_12m_pre=AC_intervention./sum(AC_intervention).*sum(Exp_intervention)*(sqrt(diet_sum(2)/25.7));


k=10;
Ymodel{k}=Y(:,(k-1).*ylg+1+lxini+1:(k-1).*ylg+1+lxini+lxini-nspc);    

colnames_bac=strrep(met_udf_co(1:nspc),'BIOMASS_','');  
colnames_bac=strrep(colnames_bac,' ','_');
output=[colnames_bac' ;num2cell(AC_baseline);num2cell(RA_4m_pre);num2cell(Exp_baseline');num2cell(AC_intervention);num2cell(RA_12m_pre);num2cell(Exp_intervention');];
rowname={'Time' 'Baseline AA' 'Baseline RA' 'Baseline RA_exp' 'Intervention AA' 'Intervention RA' 'Intervention RA_exp'}';
output=[rowname output];
writetable(cell2table(output),'Preidiction_RA_exp_adult.csv','WriteVariableNames',false,'Delimiter',';'); % WriteVariableNames means write table without column header


output=[colnames_bac' ;num2cell(RA_4m_pre);num2cell(RA_12m_pre)];
rowname={'Time'  'Baseline RA' 'Intervention RA' }';
output=[rowname output];
writetable(cell2table(output),'prediction_format.csv','WriteVariableNames',false); % WriteVariableNames means write table without column header

