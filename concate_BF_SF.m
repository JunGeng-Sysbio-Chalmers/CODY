
function  Result_total=concate_BF_SF(diet,Result_BF,Result_SF,out_filename)
% read 4m results
% Res_4m=load('20181024_4m_diet_302h.mat');
% Res_4m=Res_4m.Result_BF;

Res_4m=Result_BF;
% Res_4m=Res_4m.ReferenceStates;

Y_4m=Res_4m.SimulationResult;      
Ymet_4m=Res_4m.SimulationMetabolites;
Ymodel_4m=Res_4m.Result_cell;        
ylg_4m=Res_4m.length_tank;
lxini_4m=length(Ymet_4m{1}(1,:));    
nspc_4m=Res_4m.SpeciesNumber;
met_co_4m=Res_4m.met_udf_co;
Time_4m=Res_4m.SimulationTime;
Time_str=find((Time_4m-301.5)>0.01);
Time_str_index=Time_str(1);
Time_4m(Time_str_index);
% final state of 4month
Final_State_4m=Y_4m(Time_str_index,:)';  

% read 12m results
% Res_12m=load('20181029_12m_diet_732h.mat');
% Res_12m=Res_12m.Result_SF;

Res_12m=Result_SF;
% Res_12m=Res_12m.ReferenceStates;

Y_12m=Res_12m.SimulationResult;      
Ymet_12m=Res_12m.SimulationMetabolites;
Ymodel_12m=Res_12m.Result_cell;        
ylg_12m=Res_12m.length_tank;
lxini_12m=length(Ymet_12m{1}(1,:));
nspc_12m=Res_12m.SpeciesNumber;
met_co_12m=Res_12m.met_udf_co;
Time_12m=Res_12m.SimulationTime;
rct=Res_12m.RectorNumber;
Ini_State_12m=Res_12m.InitCondition;
Ini_State_12m=Res_12m.SimulationResult(1,:);


index_12m_in_4m=find(ismember(met_co_12m,met_co_4m)); 
full_range=1:length(met_co_12m);
diff_index=setdiff(full_range,index_12m_in_4m);
enz_4m=lxini_4m+1:ylg_4m;
enz_12m=lxini_12m+1:ylg_12m;
index_enz_12m_in_4m=enz_12m(1:length(enz_4m));    
diff_enz=enz_12m(length(enz_4m)+1:end);            

Ini_State=[];
for i=1:rct-1
    Ini_met=Ini_State_12m(index_12m_in_4m+(i-1).*ylg_12m);
    Ini_enz=Ini_State_12m(index_enz_12m_in_4m+(i-1).*ylg_12m);
    Ini_State=[Ini_State;Ini_met';Ini_enz'];
end  
k=9;
V9_12m=Ini_State_12m((k-1)*ylg_12m+ylg_12m+1);
rectum_mass=Ini_State_12m(index_12m_in_4m+k.*ylg_12m)';
Ini_State_V9_mass=[V9_12m;rectum_mass];
k=10;
blood_index=(k-1)*ylg_12m+1;
Ini_State_Blood=Ini_State_12m((k-1).*ylg_12m+1+lxini_12m+1:(k-1).*ylg_12m+1+lxini_12m+lxini_12m-nspc_12m-2);
Ini_State=[Ini_State;Ini_State_V9_mass;Ini_State_Blood'];
cmp=[Ini_State Final_State_4m];    


met_index=[];
enz_index=[];
diff_index_4m=[];
index_met=diff_index;
index_enz=diff_enz;
for i=1:rct-1
    index_met=diff_index+(i-1).*ylg_12m;
    index_enz=diff_enz+(i-1).*ylg_12m;
    met_index=[met_index;index_met'];
    enz_index=[enz_index;index_enz'];
    diff_index_4m=[diff_index_4m;index_met';index_enz'];
end
k=9;
% first is volumn, same
index_mass_rctm=(k-1)*ylg_12m+ylg_12m+1+1:(k-1).*ylg_12m+ylg_12m+1+lxini_12m;    
index_mass_rctm=index_mass_rctm(diff_index);

k=10;
met_index_fill_blood=(diff_index(end-1:end)-nspc_12m);
index_mass_blood=(k-1).*ylg_12m+1+lxini_12m+1:(k-1).*ylg_12m+1+lxini_12m+lxini_12m-nspc_12m;
index_mass_blood=index_mass_blood(met_index_fill_blood);
diff_index_4m=[diff_index_4m;index_mass_rctm';index_mass_blood'];   

new_Y_4m=zeros(length(Time_4m),length(Ini_State_12m));
full_range_12m=(1:length(Ini_State_12m))';
org_4m_index=setdiff(full_range_12m,diff_index_4m);
new_Y_4m(:,org_4m_index')=Y_4m;  % give the original 4month data to the new 4month data
% test whether new_Y_4m contains correct of zero columns
% b=find(~all(new_Y_4m))'
% % b=[intersect(b,diff_index_4m) diff_index_4m]
% c=setdiff(b,diff_index_4m)
% confirm the the final state of 4month is the same as the initial state of 12 month
temp_4m_3015=new_Y_4m(Time_str_index,:)';
temp_cmp=[temp_4m_3015 Ini_State_12m'];     
time_4m_end=Time_4m(Time_str_index);

new_Y_4m;
%% Step4: concate two results together
% Preprocess of 4m data:  sample 4m data, every 5 points
SmpPnt=1;
% Tn=length(Time_4m);
Tn=Time_str_index;   %% end with the 301.5h
NewTime_4m=Time_4m(1:SmpPnt:Tn);
NewY_4m=new_Y_4m(1:SmpPnt:Tn,:);
Concate_NewY_4m=[NewTime_4m NewY_4m];

SmpPnt=1;
Time_str=find((Time_12m-731.96)>=0.01); %
Time_str_index=Time_str(1);
Time_12m(Time_str_index);
Tn_12=Time_str_index;
Time_12m_new=Time_12m+time_4m_end;
Time_12m_new=Time_12m_new(1:SmpPnt:Tn_12);
Y_12m=Y_12m(1:SmpPnt:Tn_12,:);

Concate_NewY_12m=[Time_12m_new Y_12m];

Concate_0to12m_Y=[Concate_NewY_4m;Concate_NewY_12m]; 


RA_data=readtable('68_BF_sample_8spc_RA.csv');
NewBorn_RA_exp=RA_data.B(1:nspc_12m)'.*100;
RAexp_4m=RA_data.x4M(1:nspc_12m)'.*100;%% 
RAexp_12m=RA_data.x12M(1:nspc_12m)'.*100;%% 

Y=Concate_0to12m_Y(:,2:end);
Y(Y<0)=eps;
for k=1:rct-2
    Ymodel{k}=Y(:,1+(k-1).*ylg_12m:ylg_12m+(k-1).*ylg_12m); 
    Ymodel_met{k}=Y(:,1+(k-1).*ylg_12m:lxini_12m+(k-1).*ylg_12m); 
end
% ylg=lg_standard;
k=9;
Ymodel{k}=Y(:,1+(k-1).*ylg_12m:ylg_12m+1+lxini_12m+(k-1).*ylg_12m);
Ymodel_met{k}=Y(:,1+(k-1).*ylg_12m:lxini_12m+(k-1).*ylg_12m);
V5_rctm=Y(:,(k-1)*ylg_12m+ylg_12m+1);  
Ymodel_mass_rctm=Y(:,(k-1)*ylg_12m+ylg_12m+1+1:(k-1).*ylg_12m+ylg_12m+1+lxini_12m);    
AC_4m=Ymodel_mass_rctm(length(NewTime_4m),1:nspc_12m);
AC_12m=Ymodel_mass_rctm(end,1:nspc_12m);
%% RA:
% RA_4m_pre=AC_4m./sum(AC_4m).*sum(RAexp_4m);   
% RA_12m_pre=AC_12m./sum(AC_12m).*sum(RAexp_12m);


diet_sum=sum(diet,2);
RA_4m_pre=AC_4m./sum(AC_4m).*sum(RAexp_4m)*(sqrt(diet_sum(1)/15));   
RA_12m_pre=AC_12m./sum(AC_12m).*sum(RAexp_12m)*(sqrt(diet_sum(2)/24)); 


k=10;
Ymodel{k}=Y(:,(k-1).*ylg_12m+1+lxini_12m+1:(k-1).*ylg_12m+1+lxini_12m+lxini_12m-nspc_12m);    

colnames_bac=strrep(met_co_12m(1:nspc_12m),'BIOMASS_','');  
colnames_bac=strrep(colnames_bac,' ','_');
output=[colnames_bac' ;num2cell(NewBorn_RA_exp); num2cell(AC_4m);num2cell(RA_4m_pre);num2cell(RAexp_4m);num2cell(AC_12m);num2cell(RA_12m_pre);num2cell(RAexp_12m);];
rowname={'Time' 'Newborn_RA_exp' 'Forth month AA' 'Forth month RA' 'Forth month RA_exp' 'Twelvth month AA' 'Twelvth month RA' 'Twelvth month RA_exp'}';
output=[rowname output];
writetable(cell2table(output),'Preidiction_RA_exp_infant.csv','WriteVariableNames',false,'Delimiter',';'); 
 

lg=ylg_12m;
met_udf_co_new=met_co_12m;
met_udf_co_new(1:nspc_12m)=strrep(met_co_12m(1:nspc_12m),'BIOMASS_','');  % remove prefix of Biomass_, as well as the blank
met_udf_co_new(1:nspc_12m)=strrep(met_udf_co_new(1:nspc_12m),' ','_');
write_tank_table_infant_csv(Concate_0to12m_Y,Ymodel,Ymodel_met,met_udf_co_new,lg,lxini_12m,rct,nspc_12m,out_filename);
Result_total=Concate_0to12m_Y;



