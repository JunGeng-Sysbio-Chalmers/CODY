function Result_Intv=main_adult_Intervention(diet,nspc,Result_BSL)
global food n_species Aim  lg_standard
global store_flag
global SPC e_ini  e_ini_X spc SxZ kmax e_rel ke alpha beta Y_BM maxmue maxEnzyme e0 met_udf met_length met_co_hmo
global RA_data Km_lumen Pn_mucosa  Mass BL_Data Carb_interv total_aim
hmo=double(diet(1))
mucin=double(diet(2))
fiber=double(diet(3))
rs=double(diet(4))
addpath('./data');
load('8SPC_0405.mat');
SPC=SPC(2:8);
n_species=length(SPC);

e_ini=[];
e_ini_X=[];
 for s=1:n_species
     spc{s}=SPC{s};
     SxZ{s}=spc{s}.S;
     kmax{s}=spc{s}.kmax;
     e_rel{s}=spc{s}.e_rel;
     ke{s}=spc{s}.ke;
     alpha{s}=spc{s}.alpha;
     beta{s}=spc{s}.beta;
     Y_BM{s}=SxZ{s}(1,:)';
    maxmue{s}=kmax{s}.*Y_BM{s};
    maxEnzyme{s}=(ke{s}+alpha{s})./(beta{s}+maxmue{s});
    e0{s}=e_rel{s}.*maxEnzyme{s};
     met_udf{s}=spc{s}.met_udf;
     e_ini=[e_ini;e0{s}];
     e_ini_X=[e_ini_X;e0{s}];
 end

met_co={};
for s=1:n_species
        met_co=union(met_co,met_udf{s},'stable');
%         bm_co=union(bm_co,met_udf{s}(1),'stable');
end
met_co(1)=[];    %% remove biomass from metabolite list
% met_length=length(met_co);

Metabolites_list{1}={'HMO'  'hexose' 'succinate' 'acetate' 'lactate' 'formate' 'butyrate' 'ethanol' 'H2' 'propionate'};
Metabolites_list{2}={1000  180.16   118.09       59.04      90.08      46.03    87.098     46.07    2    73.07};  
Pn_list={power(10,-6); 2*power(10,-6);power(10,-6);5*power(10,-6);power(10,-6);5*power(10,-6);5*power(10,-6);power(10,-6);0; 5*power(10,-6)};    
Km_list={0.3;         0.2;          0.1;        0.1;           0.1;         0.1;           0.2;          0.1;        0;   0.1      };
Pn_list={0.5*power(10,-6); 10*power(10,-6);power(10,-6); 10*power(10,-5);power(10,-6);5*power(10,-6);10*power(10,-5);power(10,-6);0; 10*power(10,-5)};

met_co_hmo=['HMO';met_co];
met_length=length(met_co);

for i=1:length(met_co_hmo)
    met_index(i)=find(strcmpi(met_co_hmo(i),Metabolites_list{1}));
end
M_biomass=113;   
Mass=cell2mat(Metabolites_list{2}(met_index));
Pn=cell2mat(Pn_list(met_index));
Km_lumen=cell2mat(Km_list(met_index));
Pn=Pn.*3600/100; 
Pn=Pn.*10;
av=1/0.05;
Pn_mucosa=Pn.*av;

%% deal with RA_data
RA_file='Mean_RA_species.txt';
RA_data=readtable(RA_file);
% RA_data=importdata(RA_file);
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


baseline=RA_data{:,3};  
total_RA_bl=sum(baseline);
ratio=baseline./total_RA_bl;

ResulsBL=Result_BSL;
T=ResulsBL.SimulationTime;
Y=ResulsBL.SimulationResult;
Y_model=ResulsBL.Result_cell;
Y_mets=ResulsBL.SimulationMetabolites;

lg_standard=90;
    
Y_bl=Y;
Y_model_bl=Y_model;
Y_met_bl=Y_mets;
ylg_bl=lg_standard;
lxini_bl=length(Y_met_bl{1}(1,:));
n_spc=n_species;
met_co=met_co_hmo;
Time_bl=T;
% Time_str=find((Time_bl-611.96)>0.01); %% this is the time for 4th meal of 4month infant
Time_str=find((Time_bl-600)>0.01); %% this is the time for 4th meal of 4month infant

Time_str_index=Time_str(1);
    for i=1:length(met_co)
        met_int(i)=find(strcmpi(met_co(i),met_co_hmo));
    end
length_enz=length(Y_model_bl{1}(1,:));
for i=1:9
    ini_bl_x0{i}=Y_met_bl{i}(Time_str_index,1:n_species); %% 5species initial value
    ini_bl_met{i}=Y_met_bl{i}(Time_str_index,n_species+1:end);
%     ini_4m_x0{i}=Y_met_4m{i}(end,1:n_spc_4m); %% 5species initial value
%     ini_4m_met{i}=Y_met_4m{i}(end,n_spc_4m+1:end);
    ini_bl_met{i}=[ini_bl_met{i}'];
    ini_bl_enz{i}=Y_model_bl{i}(Time_str_index,lxini_bl+1:length_enz);
%         ini_4m_enz{i}=Y_model_4m{i}(end,lxini_4m+1:length_enz);

end

ini_bl_met{10}=Y_model_bl{10}(end,:)';

x10=[ini_bl_x0{1} ]';
x20=[ini_bl_x0{2} ]';
x30=[ini_bl_x0{3} ]';
x40=[ini_bl_x0{4} ]';
x50=0.00001.*ones(n_species,1).*ratio; 
XB10=[ini_bl_x0{5} ]';
XB20=[ini_bl_x0{6} ]';
XB30=[ini_bl_x0{7} ]';
XB40=[ini_bl_x0{8} ]';

lenz_bl=length(ini_bl_enz{1});
for i=1:4
    ini_bl_enz{i}=[ini_bl_enz{i} e_ini(lenz_bl+1:end)']';
end
for i=5:8
    ini_bl_enz{i}=[ini_bl_enz{i} e_ini_X(lenz_bl+1:end)']';
end
ini_bl_enz{9}=[ini_bl_enz{9} e_ini(lenz_bl+1:length(e_ini))']';  %% cannot go to the end!! due to longer length in the rectum!!

BL_Data.ini_stat=Y_bl(end,:);
k=9;
ylg=ylg_bl;
lxini=lxini_bl;

rectum_mass_index=(k-1)*ylg+ylg+1+1:(k-1).*ylg+ylg+1+lxini;
BL_Data.ini_stat(rectum_mass_index)=0.0000001.*ones(1,lxini);


store_flag=0;

Result_Intv=main_10tank(hmo,mucin,fiber,rs,nspc,Result_BSL);

clear food n_species RA_intv  lg_standard
clear store_flag
clear SPC e_ini  e_ini_X spc SxZ kmax e_rel ke alpha beta Y_BM maxmue maxEnzyme e0 met_udf met_length met_co_hmo
clear RA_data Km_lumen Pn_mucosa  Mass BL_Data Carb_interv total_aim
end

function Result_Intv=main_10tank(hmo,mucin,fiber,rs,nspc,Result_BSL)
global food n_species RA_intv  lg_standard
global store_flag
global SPC e_ini  e_ini_X spc SxZ kmax e_rel ke alpha beta Y_BM maxmue maxEnzyme e0 met_udf met_length met_co_hmo
global RA_data Km_lumen Pn_mucosa  Mass BL_Data Carb_interv total_aim

RA_file='Mean_RA_species.txt';
RA_data=readtable(RA_file);
% RA_data=importdata(RA_file);
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
baseline=RA_data{:,3};  
RA_intv=RA_data{:,2}';  
TT_RA=sum(RA_data{:,2:end});
Total_RA=mean(TT_RA([1,3:end]));
n_rct=10;
V_colon_t=4.5; 
V_colon=[0.6 0.6 0.6 0.6 0.05 0.05 0.05 0.05 0.3 2]';  
x=[0.00943    0.0063    0.0037     0.0042
    0.4839    2.2859    0.4272    1.0630
    0.8176    1.6043    0.0468    2.1797
    0.0107    0.0003    0.0001    0.0001
    0.1900    0.5406    0.0148    0.8614
    0.0898    3.2433    0.4290    0.1206
    0.8968    0.9452    0.5925    0.0361
    0.1161    2.8742    0.0065    0.5014];
diet_matrix=[hmo mucin fiber rs]'./(hmo+mucin+fiber+rs);
deg=x*diet_matrix;
HMO=(hmo+mucin+fiber+rs)./0.45;  
GLU_in=150;     
F0=[0.15  0.15  0  0.5]';   
f=[0  0.01  0.01  0.01]';  
kd=[0.1 0.05 0.025 0.01
    16 16 25 25];
sim_time=612;              
k_hyd=1.2*power(10,3)/2400;
K_XZ=10;
K_XZ_5spc=[15 10 1000 10 25 14 5]';
k_hyd=deg(nspc).*K_XZ_5spc;
k_hyd_5spc=k_hyd.*[0.9524 1 1 2 0.5 1.229 0.4]';
Y_SZ=5.*ones(n_species,1);
hyd_spc=ones(n_species,1);
k_hyd=k_hyd_5spc;
K_XZ=K_XZ_5spc;
hmo_para=[HMO;hyd_spc; k_hyd; K_XZ; Y_SZ];
HMO_ini=0.5;
ini_stat=BL_Data.ini_stat;
Km={1/9 1/15 1/20 1/30};
density=0.7*1000;
fprintf('\n\n');
disp('Simulation for 2nd period: Please keep on waiting......');
fprintf('\n');
[T,Y,Y_model,Y_mets,met_udf_co]=dynamic_simulation_CM_adult_Intervention(SPC,n_species,n_rct,GLU_in,F0,f,kd,V_colon,ini_stat,sim_time,lg_standard,Km_lumen,Pn_mucosa,Km,density,Mass,hmo_para);
ResulsBL.T=T+612;
ResulsBL.Y=Y;
ResulsBL.Y_model=Y_model;
ResulsBL.Y_mets=Y_mets;
% ylg=lg_8tank;   
lxini=length(met_co_hmo)+n_species;
ylg=lg_standard;
Tmodel=T;
[ys,yl]=size(Y);
Y(Y<0)=0;

for k=1:n_rct-2
    Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+(k-1).*ylg); 
    Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg); 
end
k=9;
Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+1+lxini+(k-1).*ylg);
Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg);
V5_rctm=Y(:,(k-1)*ylg+ylg+1);  
Ymodel_mass_rctm=Y(:,(k-1)*ylg+ylg+1+1:(k-1).*ylg+ylg+1+lxini);   
AC=Ymodel_mass_rctm(end,1:n_species);
RA=AC./sum(AC)*sum(RA_intv);      
% RAmodel=RA.*100.*(sqrt((hmo+mucin+fiber+rs)/25.7))
% RAexp=RA_intv.*100   

    ReferenceStates=StoreResults_10tank(SPC,n_species,n_rct,GLU_in,F0,f,kd,V_colon,ini_stat,sim_time,lg_standard,Km_lumen,Pn_mucosa,Km,density,T,Y,Y_model,Y_mets,met_udf_co,lg_standard);
    refdata_name='20190319_CM_diet_intervention.mat';
        Result_Intv=ReferenceStates;

    save(refdata_name,'ReferenceStates','-v7.3');
%     store_data(T,Y,Y_model,Y_mets,met_udf_co);
end