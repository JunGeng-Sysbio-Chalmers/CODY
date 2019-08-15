
function Result_total=concate_4m_12m_mets_3d(Result_BF,Result_SF,index)
%% index is the bacterial/metabolites that want to vitulize
addpath('./export_fig-master');
% Res_4m=load('20181024_4m_diet_302h.mat');
Res_4m=load(Result_BF);
Res_4m=Res_4m.Result_BF;

Y_4m=Res_4m.SimulationResult;      %
Ymet_4m=Res_4m.SimulationMetabolites;%% 
Ymodel_4m=Res_4m.Result_cell;        
ylg_4m=Res_4m.length_tank;
lxini_4m=length(Ymet_4m{1}(1,:));    
nspc_4m=Res_4m.SpeciesNumber;
met_co_4m=Res_4m.met_udf_co;
Time_4m=Res_4m.SimulationTime;
Time_str=find((Time_4m-301.5)>0.01); 
Time_str_index=Time_str(1);
Time_4m(Time_str_index);
Final_State_4m=Y_4m(Time_str_index,:)';  

Res_12m=load(Result_SF);
Res_12m=Res_12m.Result_SF;

Y_12m=Res_12m.SimulationResult;      %%
Ymet_12m=Res_12m.SimulationMetabolites;%%
Ymodel_12m=Res_12m.Result_cell;        
ylg_12m=Res_12m.length_tank;
lxini_12m=length(Ymet_12m{1}(1,:));
nspc_12m=Res_12m.SpeciesNumber;
met_co_12m=Res_12m.met_udf_co;
Time_12m=Res_12m.SimulationTime;
rct=Res_12m.RectorNumber;
Ini_State_12m=Res_12m.InitCondition;
Ini_State_12m=Res_12m.SimulationResult(1,:);


index_12m_in_4m=find(ismember(met_co_12m,met_co_4m)); % bacterial+metabolites index
full_range=1:length(met_co_12m);
diff_index=setdiff(full_range,index_12m_in_4m); 
% enzyme of 12m in 4m index
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
new_Y_4m(:,org_4m_index')=Y_4m;  %


temp_4m_3015=new_Y_4m(Time_str_index,:)';
temp_cmp=[temp_4m_3015 Ini_State_12m'];     
time_4m_end=Time_4m(Time_str_index);

new_Y_4m;
%% Step4: concate two results together
SmpPnt=10;
Tn=Time_str_index;   
NewTime_4m=Time_4m(1:SmpPnt:Tn);
NewY_4m=new_Y_4m(1:SmpPnt:Tn,:);
Concate_NewY_4m=[NewTime_4m NewY_4m];

SmpPnt=1;

Time_str=find((Time_12m-731.96)>=0.01); 
Time_str_index=Time_str(1);
Time_12m(Time_str_index);
Tn_12=Time_str_index;
Time_12m_new=Time_12m+time_4m_end;
Time_12m_new=Time_12m_new(1:SmpPnt:Tn_12);
Y_12m=Y_12m(1:SmpPnt:Tn_12,:);

Concate_NewY_12m=[Time_12m_new Y_12m];

Concate_0to12m_Y=[Concate_NewY_4m;Concate_NewY_12m]; 

Result_total=Concate_0to12m_Y;
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

k=10;
Ymodel{k}=Y(:,(k-1).*ylg_12m+1+lxini_12m+1:(k-1).*ylg_12m+1+lxini_12m+lxini_12m-nspc_12m);  

T_feces=Concate_0to12m_Y(:,1);
str_index=1;
defe_index_4m=[30.7,53.3,75.6,90,...
            112.3,134.3,156.3,178.3,...
            200, 222, 243.9,265.6,287.2, 301.5];
defe_index_12m=[359.1];
i=1;
while i<29
    i=i+1;
    defe_index_12m=[defe_index_12m defe_index_12m(i-1)+24];
end
defe_index=[defe_index_4m defe_index_12m];
defe_time_index=cell2mat(arrayfun(@(x) find(abs(T_feces-x)<0.23,1), defe_index,'un',0));

Fece_bacterial=Ymodel_mass_rctm(defe_time_index,2:end);  

smpt=1;
Fece_bacterial_plot=[];
    Fece_bacterial_plot=[Fece_bacterial_plot Fece_bacterial(1:smpt:end,1:8)];
T_feces_plot=[];
T_feces_plot=[T_feces_plot T_feces(defe_time_index)];
Spatial_Blg=1:9;
smpt=500;

Tmodel=Concate_0to12m_Y(1:smpt:end,1);
Blg_SpTm=[];    
if index<9    %% microbial
    for i=1:9
    Blg_SpTm=[Blg_SpTm Ymodel{i}(1:smpt:end,index)];
    end
else          %% metabolites
    for i=1:9
    Blg_SpTm=[Blg_SpTm Ymodel{i}(1:smpt:end,index)];
    end
    Spatial_Blg=1:10;
    Blg_SpTm=[Blg_SpTm Ymodel{10}(1:smpt:end,index-8)];
end
[Blg_X, Blg_Y]=meshgrid(Tmodel,Spatial_Blg);


% set(gcf,'unit','centimeters','position',[17 27 25 30]);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
h1=figure('Renderer', 'painters', 'Position', [10 10 600 1000]);

surf(Blg_X',Blg_Y',Blg_SpTm,'FaceAlpha',0.75,'edgeColor','none')
colorbar



% set(0,'defaultfigurecolor',[1 1 1])
grid off
% colorbar
set(gca,'FontSize',18,'FontName','Helvetica','xtick',[],'ytick',[],'ztick',[],'color',[190, 222, 247]./255,'linewidth',2);
set(get(gca,'XLabel'),'String','Time');
set(get(gca,'YLabel'),'String','Colon Site');

xlim([0 1000]);
zlim([0 30]); 

% a='Metabolite_3d_12.jpg'

view(-35,30)
% axis off
output_file=strcat('Metabolite_3d_',num2str(index),'.jpg');
% export_fig Metabolites_3d.png -transparent
% export_fig Metabolites_3d.jpg -transparent
export_fig(sprintf(output_file), h1);

fig_3d_dir='./www/';
% copyfile('Metabolites_3d.jpg',fig_3d_dir);
copyfile(output_file,fig_3d_dir);

end

