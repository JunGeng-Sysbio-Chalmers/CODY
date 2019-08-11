function [Tmodel,Y,Ymodel,Ymodel_met,met_udf_co]=dynamic_simulation_Infant_BF_3feces_6fed(community,n_spc,num_rct,FSG,Feed,flow_back,k,V,ini_stat,time,length_each_tank,km_lumen,Pn_mucosa,km_water,density,Met_Mass,hmo_para,vld,uv)

global lxini spc 
global SxZ_tot kmax KM ke alpha beta subs_cn biom_inx met_udf met_udf_co n_carbon met_udf_co_new
global Y_BM Y_SUB maxEnzyme kl rct obsv V_obsv t_obsv F5 length_additional kd C_Kd GLU_in GLU_obsv HMO_obsv
global lg_8tank  f_back   n_rct  km   km_w   Pn   rho mol_mass
global flag_meal1 flag_meal2 flag_meal3 flag_fece 
global HMO_in spc_hyd kmax_hmo KM_hmo Y_hmo HMO_const F5_obsv


n_rct=num_rct;
rct=n_rct;
if nargin<19
    uv=false;
end
 
if nargin<18
    vld=false;
end 

f_back=flow_back;
n_s=n_spc;
kl={};
kd=k(1,:)';
C_Kd=k(2,:)';  
km=km_lumen;
km_w=km_water;
Pn=Pn_mucosa;
rho=density;
mol_mass=Met_Mass./1000;

HMO_const=hmo_para(1);
spc_hyd=hmo_para(2:n_s+1);
kmax_hmo=hmo_para(n_s+2:2*n_s+1);
KM_hmo=hmo_para(2*n_s+2:3*n_s+1);
Y_hmo=hmo_para(3*n_s+2:end);

 for s=1:n_s
     spc{s}=community{s};

     name{s}=spc{s}.id;
     SxZ{s}=spc{s}.S;
     kmax{s}=spc{s}.kmax;
     K_MM{s}=spc{s}.KM;
     ke{s}=spc{s}.ke;
     alpha{s}=spc{s}.alpha;
     beta{s}=spc{s}.beta;
     met_udf{s}=spc{s}.met_udf;
     e_rel{s}=spc{s}.e_rel;
     subs_cn{s}=spc{s}.subs_cn;
     
     if exist('bm_coef')
         biom_coef{s}=bm_coef{s};
     else
         biom_coef{s}=spc{s}.biom_coef;
     end
     
     n_carbon{s}=spc{s}.n_carbon;  
     [cr,cl]=size(n_carbon{s}); 
     if (cr>1)&&(cl>1)
        disp('this species use more than one cl species');
     end
     biom_inx{s}=strmatch('BIOM',met_udf{s});

        
    [sr,sl]=size(SxZ{s});
    n_ezm=sl;
    kl{s}=length(kmax{s});
    KM{s}=[];
    for i=1:length(subs_cn{s})
        KM{s}(:,i)=K_MM{s}(:,i).*ones(kl{s},1);  
    end
    alpha{s}=alpha{s}.*ones(kl{s},1);   
    beta{s}=beta{s}.*ones(kl{s},1);
    ke{s}=ke{s}.*ones(kl{s},1); 
    e_rel{s}=e_rel{s}.*ones(kl{s},1); 

    Y_BM{s}=SxZ{s}(1,:)';
    Y_SUB{s}=SxZ{s}(subs_cn{s},:)';    
    maxmue{s}=kmax{s}.*Y_BM{s};
    maxEnzyme{s}=(ke{s}+alpha{s})./(beta{s}+maxmue{s}); 
    e0{s}=e_rel{s}.*maxEnzyme{s};     
    para{s}=kmax{s};
    options_ode=[];
 end
 
met_co={};
for s=1:n_s
        met_co=union(met_co,met_udf{s},'stable');
end
met_co(1)=[]; 
for s=1:n_s
    met_udf{s}{1}=[met_udf{s}{1} '_' name{s}];
end

met_udf_co={};
bm_co={};
for s=1:n_s
        met_udf_co=union(met_udf_co,met_udf{s},'stable');
        bm_co=union(bm_co,met_udf{s}(1),'stable');
end
clear_index=[];
for i=1:length(met_udf_co)   
    if ~isempty(strmatch('BIOMASS',met_udf_co{i}))
        clear_index=[clear_index;i];
    end
end
met_udf_co(clear_index)=[]; 
met_udf_co=union(bm_co,met_udf_co,'stable');

cm_met={};  
for i=1:n_s
    for j=1:length(met_udf_co)
        temp=strmatch(met_udf_co{j},met_udf{i},'exact');
        if (isempty(temp))   
            cm_met{i}{j}=1;   
        else
            cm_met{i}{j}=[];  
        end
    end                   
end
row_index={};
for i=1:n_s
    temp=[];
    for j=1:length(cm_met{i})
        if find(cm_met{i}{j})            
            temp=[temp;j];  
        end
    end
    row_index{i}=temp;  
end
SxZ_tot={};
%%
for i=1:n_s    
    temp_S=SxZ{i};
    for j=1:length(row_index{i})  
        [rlos,clos]=size(temp_S);
        if row_index{i}(j)<=rlos    
                ij=row_index{i}(j);  
                temp_S=[temp_S(1:ij-1,:);zeros(1,length(temp_S(1,:)));temp_S(ij:end,:)];   
        elseif row_index{i}(j)>rlos
            ij=row_index{i}(j);  
            temp_S=[temp_S(1:ij-1,:);zeros(1,length(temp_S(1,:)))];
        end
    end
    SxZ_tot{i}=temp_S;
end
%%
met_udf_co_new={};
bm_new=[met_udf_co(1:n_s); met_udf_co(1:n_s)];
met_udf_co_new(1:2*n_s)=bm_new;
met_udf_co_new(2*n_s+1:2*n_s+length(met_udf_co)-n_s)=met_udf_co(n_s+1:end);
met_udf_co_new=met_udf_co_new';

met_udf_co(n_s+2:end+1)=met_udf_co(n_s+1:end);
met_udf_co(n_s+1)={'HMO'};

    sub_index_temp=[];
    for i=1:n_s
        for j=1:length(subs_cn{i})
           sub_index_temp(i,j)=strmatch(met_udf{i}(subs_cn{i}(j)),met_udf_co);   
        end
    end
    hmo_ix=find(strcmp('HMO',met_udf_co));
for k=1:n_rct-1
    HMO_index(k)=hmo_ix+(k-1).*length_each_tank;  
end
k=n_rct;  
    HMO_index(k)=(k-1).*length_each_tank+length(met_udf_co)+1+hmo_ix-n_s;
   [sr,sc]=size(sub_index_temp);
   sub_shape=reshape(sub_index_temp,sr*sc,1);
   t_inx=length(sub_shape);
   subset_t=1:t_inx;
   sub_set_inx=find(sub_shape>0); 
   nan_set_inx=setxor(subset_t,sub_set_inx);
sub_index_new=zeros(t_inx,1);
for k=1:n_rct-1
    sub_index_new=zeros(t_inx,1);
    sub_index_new(sub_set_inx)=sub_shape(sub_set_inx)+(k-1).*length_each_tank;
    sub_index_new(nan_set_inx)=nan;
    temp_sub_cell=reshape(sub_index_new,sr,sc);
    pre_sub_cell=num2cell(temp_sub_cell',[1,n_s]); 
    for ss=1:n_s
        pre_sub_cell{ss}(isnan(pre_sub_cell{ss}))=[];   %%
    end
    substrate_index{k}=pre_sub_cell;   
end
    
lg_met_co=length(met_udf_co);
lg_enzyme=length_each_tank-lg_met_co;
ylg=length_each_tank;
for k=1:rct-1
    for i=1:n_s
        enzyme_index{k}{i}=((lg_met_co+1+(k-1).*ylg):(lg_met_co+length(kmax{i})+(k-1).*ylg))';    
        lg_met_co=lg_met_co+length(kmax{i});
    end
    lg_met_co=length(met_udf_co);
end 
    for i=1:n_s
        if length(subs_cn{i}>1)
            [rkm,ckm]=size(KM{i});
            if ckm==1
                for j=1:length(subs_cn{i})
                    KM{i}(:,j)=KM{i}(:,1);
                end
            end
        end
    end
%%

lxini=length(met_udf_co);
e_ini=[];
e_ini_X=[];
for s=1:n_s
    e_ini=[e_ini;e0{s}];
    e_ini_X=[e_ini_X;e0{s}];   
end

V_rct=0.00001;  
y0=ini_stat;
lg_8tank=length_each_tank;
y0(9*lg_8tank+1)=V_rct;
options_ode=odeset('maxstep',0.01);
obsv=1;   
V_obsv=V_rct;
t_obsv=0;
F5=0; 
GLU_in=0;%
HMO_in=0;  %%
GLU_obsv=[];HMO_obsv=[];F5_obsv=F5;

tspan=[0;302;0.05];

[T Y]=my_runge_kutta4_smp_wait(@qssm_comm_Infant_BF_3feces_6fed,y0,tspan,kmax,n_s,FSG,Feed,V,substrate_index,HMO_index,enzyme_index);


ylg=lg_8tank;   

Tmodel=T;
[ys,yl]=size(Y);
Y(Y<0)=eps;
for k=1:rct-2
    Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+(k-1).*ylg); %%
    Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg); %%
end
k=9;
Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+1+lxini+(k-1).*ylg);
Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg);
V5_rctm=Y(:,(k-1)*ylg+ylg+1);  %%
Ymodel_mass_rctm=Y(:,(k-1)*ylg+ylg+1+1:(k-1).*ylg+ylg+1+lxini);  

k=10;
Ymodel{k}=Y(:,(k-1).*ylg+1+lxini+1:(k-1).*ylg+1+lxini+lxini-n_s);    %% no bacteria, only metabolites




