function Result_SF=main_Infant_SFeeding(diet,nspc,Result_BF)
global food n_species RA_m12  lg_standard
global store_flag
global SPC e_ini  e_ini_X spc SxZ kmax e_rel ke alpha beta Y_BM maxmue maxEnzyme e0 met_udf met_length met_co_hmo
global RA_data Km_lumen Pn_mucosa  Mass
global ini_4m_x0  ini_4m_met ini_4m_enz
global x10 x20 x30 x40 x50 XB10 XB20 XB30 XB40
hmo=diet(1)
mucin=diet(2)
fiber=diet(3)
rs=diet(4)
Pn_mucosa={}; 
addpath('./data');
    load 8SPC_0405
    n_species=size(SPC,2);
    RA_data=readtable('68_BF_sample_8spc_RA.csv');
    RA_m12=RA_data.x12M(1:n_species)'.*100;
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
    Metabolites_list{2}={1000  180.16   118.09       59.04      90.08      46.03    87.098     46.07    2    73.07};  %% 
    Pn_list={power(10,-6); 2*power(10,-6);power(10,-6);5*power(10,-6);power(10,-6);5*power(10,-6);5*power(10,-6);power(10,-6);0; 5*power(10,-6)};    %%
    Km_list={0.3;         0.2;          0.1;        0.1;           0.1;         0.1;           0.2;          0.1;        0;   0.1      };
    Pn_list={power(10,-6); 2*power(10,-6);power(10,-6); 10*power(10,-5);power(10,-6);5*power(10,-6);10*power(10,-5);power(10,-6);0; 10*power(10,-5)};   
    
 
    met_co_hmo=['HMO';met_co];
    met_length=length(met_co);

    for i=1:length(met_co_hmo)
        met_index(i)=find(strcmpi(met_co_hmo(i),Metabolites_list{1}));
    end
    Mass=cell2mat(Metabolites_list{2}(met_index));
    Pn=cell2mat(Pn_list(met_index));
    Km_lumen=cell2mat(Km_list(met_index));
    Pn=Pn.*3600/100; %% unit transfer to dm/h
    Pn=Pn.*10;
    av=1/0.05;
    Pn_mucosa=Pn.*av;
ratio=RA_data.x4M(1:n_species).*100;
ratio(6:8)=[0.0005  0.008  0.009];
total_ratio=sum(ratio);
ratio=ratio./total_ratio;


%% old setting, no relate with 4month results
x10=3.65.*ones(n_species,1).*ratio;     
x20=10.*ones(n_species,1).*ratio;
x30=15.*ones(n_species,1).*ratio;
x40=18.*ones(n_species,1).*ratio;
XB10=15.*ones(n_species,1).*ratio;    
XB20=10.*ones(n_species,1).*ratio;
XB30=4.*ones(n_species,1).*ratio;
XB40=2.*ones(n_species,1).*ratio;
x50=0.00001.*ones(n_species,1).*ratio;

    
%% load 4month data 
% load 20181024_4m_diet.mat
% load 20181024_4m_diet_302h.mat
ReferenceStates=Result_BF;


Y_4m=ReferenceStates.SimulationResult;
Y_model_4m=ReferenceStates.Result_cell;
Y_met_4m=ReferenceStates.SimulationMetabolites;
ylg_4m=ReferenceStates.length_tank;
lxini_4m=length(Y_met_4m{1}(1,:));
n_spc_4m=ReferenceStates.SpeciesNumber;
met_co_4m=ReferenceStates.met_udf_co(n_spc_4m+1:end);
Time_4m=ReferenceStates.SimulationTime;
Time_str=find((Time_4m-301.5)>0.01); 
Time_str_index=Time_str(1);
    for i=1:length(met_co_4m)
        met_4to12(i)=find(strcmpi(met_co_4m(i),met_co_hmo));
    end
length_enz=length(Y_model_4m{1}(1,:));
for i=1:9
    ini_4m_x0{i}=Y_met_4m{i}(Time_str_index,1:n_spc_4m); %% 5species initial value
    ini_4m_met{i}=Y_met_4m{i}(Time_str_index,n_spc_4m+1:end);
%     ini_4m_met{i}=Y_met_4m{i}(end,n_spc_4m+1:end);
    ini_4m_met{i}=[ini_4m_met{i}';ones(2,1).*0.1];
    ini_4m_enz{i}=Y_model_4m{i}(Time_str_index,lxini_4m+1:length_enz);
%         ini_4m_enz{i}=Y_model_4m{i}(end,lxini_4m+1:length_enz);

end

ini_4m_met{10}=Y_model_4m{10}(Time_str_index,:)';
ini_4m_met{10}=[ini_4m_met{10};ones(2,1).*0.1];


x10=[ini_4m_x0{1} x10(6:8)']';
x20=[ini_4m_x0{2} x20(6:8)']';
x30=[ini_4m_x0{3} x30(6:8)']';
x40=[ini_4m_x0{4} x40(6:8)']';
x50=0.00001.*ones(n_species,1).*ratio; %% 
XB10=[ini_4m_x0{5} XB10(6:8)']';
XB20=[ini_4m_x0{6} XB20(6:8)']';
XB30=[ini_4m_x0{7} XB30(6:8)']';
XB40=[ini_4m_x0{8} XB40(6:8)']';

lenz_4m=length(ini_4m_enz{1});
for i=1:4
    ini_4m_enz{i}=[ini_4m_enz{i} e_ini(lenz_4m+1:end)']';
end
for i=5:8
    ini_4m_enz{i}=[ini_4m_enz{i} e_ini_X(lenz_4m+1:end)']';
end
ini_4m_enz{9}=[ini_4m_enz{9} e_ini(lenz_4m+1:length(e_ini))']';  
n_rct=10;
V_colon_t=4.5; 
V_colon=[0.15 0.15 0.15 0.15 0.03 0.03 0.03 0.03 0.15 1]';  
GLU_in=150;    
K_XZ_5spc=[55.54 10.65 6 95 10.37 60 27.34 56.56];
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
HMO=(hmo+mucin+fiber+rs)./0.3;  
F0=[0.1  0.1  0  0.5]';   
f=[0  0.01  0.01  0.01]';   
kd=[0.075 0.03 0.015 0.008
    16 16 16 16];   
sim_time=600;            
M_biomass=113; 
Y_SZ=5.*ones(n_species,1);
hyd_spc=ones(n_species,1);  
k_hyd=deg.*K_XZ_5spc';
HMO_ini=0.5;
k_hyd=k_hyd.*[1 2/1.1 4/3 1 4/3 1/30 1/1.093 1/3.78]';
K_XZ_5spc=K_XZ_5spc';
hmo_para=[HMO;hyd_spc; k_hyd; K_XZ_5spc; Y_SZ];
ini_stat_Lumen1=[x10;ini_4m_met{1};ini_4m_enz{1}];
ini_stat_Lumen2=[x20;ini_4m_met{2};ini_4m_enz{2}];
ini_stat_Lumen3=[x30;ini_4m_met{3};ini_4m_enz{3}];
ini_stat_Lumen4=[x40;ini_4m_met{4};ini_4m_enz{4}];
ini_stat_Mucosa1=[XB10;ini_4m_met{5};ini_4m_enz{5}];
ini_stat_Mucosa2=[XB20;ini_4m_met{6};ini_4m_enz{6}];
ini_stat_Mucosa3=[XB30;ini_4m_met{7};ini_4m_enz{7}];
ini_stat_Mucosa4=[XB40;ini_4m_met{8};ini_4m_enz{8}];
lg_standard=length(ini_stat_Lumen1); 
ini_conc_V5=zeros(met_length+1,1);        
% ini_stat_V5=[x50;ini_4m_met{9};ini_4m_enz{9};V_colon(9);x50.*V_colon(9);ini_4m_met{9}.*V_colon(9)];  
ini_stat_V5=[x50;ini_conc_V5;e_ini;V_colon(9);x50.*V_colon(9);ini_conc_V5.*V_colon(9)];  %% 
ini_stat_Blood=[ini_4m_met{10}];
ini_stat=[ini_stat_Lumen1;ini_stat_Lumen2;ini_stat_Lumen3;ini_stat_Lumen4;ini_stat_Mucosa1;ini_stat_Mucosa2;ini_stat_Mucosa3;ini_stat_Mucosa4;ini_stat_V5;ini_stat_Blood];
Km={1/2 1/6 1/40 1/70};
density=0.7*1000;
fprintf('\n\n');
disp('Simulation for 2nd period: Please keep on waiting......');
fprintf('\n');

[T,Y,Y_model,Y_mets,met_udf_co]=dynamic_simulation_Infant_SF(SPC,n_species,n_rct,GLU_in,F0,f,kd,V_colon,ini_stat,sim_time,lg_standard,Km_lumen,Pn_mucosa,Km,density,Mass,hmo_para);

lxini=length(met_co_hmo)+n_species;
ylg=lg_standard;
Tmodel=T;
[ys,yl]=size(Y);
Y(Y<0)=eps;

for k=1:n_rct-2
    Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+(k-1).*ylg); 
    Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg); 
end
% ylg=lg_standard;
k=9;
Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+1+lxini+(k-1).*ylg);
Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg);
V5_rctm=Y(:,(k-1)*ylg+ylg+1); 
Ymodel_mass_rctm=Y(:,(k-1)*ylg+ylg+1+1:(k-1).*ylg+ylg+1+lxini);  
AC=Ymodel_mass_rctm(end,1:n_species);
RA=AC./sum(AC).*sum(RA_m12)./100;

    Result_SF=StoreResults_10tank(SPC,n_species,n_rct,GLU_in,F0,f,kd,V_colon,ini_stat,sim_time,lg_standard,Km_lumen,Pn_mucosa,Km,density,T,Y,Y_model,Y_mets,met_udf_co,lg_standard);
    refdata_name='20181029_12m_diet_732h.mat';
    save(refdata_name,'Result_SF','-v7.3');

clear food n_species RA_m12  lg_standard
clear store_flag
clear SPC e_ini  e_ini_X spc SxZ kmax e_rel ke alpha beta Y_BM maxmue maxEnzyme e0 met_udf met_length met_co_hmo
clear RA_data Km_lumen Pn_mucosa  Mass
clear ini_4m_x0  ini_4m_met ini_4m_enz
clear x10 x20 x30 x40 x50 XB10 XB20 XB30 XB40 SxZ
end

