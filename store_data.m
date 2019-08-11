function store_data(T,Y,Y_model,Y_mets,met_udf_co,file_out_name)
global lxini spc  n_species Aim lg_standard
global SxZ_tot kmax KM ke alpha beta subs_cn biom_inx met_udf met_udf_co n_carbon met_udf_co_new
global Y_BM Y_SUB maxEnzyme kl rct obsv V_obsv t_obsv F5 length_additional kd C_Kd GLU_in GLU_obsv HMO_obsv
global lg_8tank  f_back   n_rct  km   km_w   Pn   rho mol_mass
global flag_meal1 flag_meal2 flag_meal3 flag_fece 
global HMO_in spc_hyd kmax_hmo KM_hmo Y_hmo HMO_const F5_obsv


javaaddpath('../poi_library/poi-3.8-20120326.jar');
javaaddpath('../poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('../poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('../poi_library/xmlbeans-2.3.0.jar');
javaaddpath('../poi_library/dom4j-1.6.1.jar');
javaaddpath('../poi_library/stax-api-1.0.1.jar');

ylg=lg_8tank;   
lg=ylg;
ns=n_species;
met_co=met_udf_co;
n_s=ns;
rct_num=rct;
Tmodel=T;
[ys,yl]=size(Y);
for p=1:ys
    for q=1:yl
        if Y(p,q)<0
            Y(p,q)=0;
        end
    end
end
for k=1:rct-2
    Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+(k-1).*ylg); %% 
    Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg); %% 
end
% ylg=lg_standard;
k=9;
Ymodel{k}=Y(:,1+(k-1).*ylg:ylg+1+lxini+(k-1).*ylg);
Ymodel_met{k}=Y(:,1+(k-1).*ylg:lxini+(k-1).*ylg);
V5_rctm=Y(:,(k-1)*ylg+ylg+1);  %%
Ymodel_mass_rctm=Y(:,(k-1)*ylg+ylg+1+1:(k-1).*ylg+ylg+1+lxini);    %% 

k=10;
Ymodel{k}=Y(:,(k-1).*ylg+1+lxini+1:(k-1).*ylg+1+lxini+lxini-ns);    %% no bacteria, only metabolites


for ij=1:9
    Conc_Ymodel{ij}=Ymodel{ij};
    Conc_Ymodel_met{ij}=Ymodel_met{ij};
end
ij=10;
Conc_Ymodel{ij}=Ymodel{ij};
SmpPnt=1;
Tn=length(Tmodel);
% nspc=5; 
% nspc=8;
type='lumen';
% for i=1:SmpPnt:Tn
    NewTime=Tmodel(1:SmpPnt:Tn);
    NewY=Y(1:SmpPnt:Tn,:);
    for ij=1:9
        New_Ymodel{ij}=Conc_Ymodel{ij}(1:SmpPnt:Tn,:);
        New_Ymodel_met{ij}=Conc_Ymodel_met{ij}(1:SmpPnt:Tn,:);
    end
ij=10;
New_Ymodel{ij}=Ymodel{ij}(1:SmpPnt:Tn,:);


for i=1:n_species
    met_co{i}=strrep(met_co{i}, ' ', '_');
end

col_header=['Time',met_co']
met_co1=met_co(ns+1:end);
col_header2=['Time',met_co1']
Tanks={'Lumen_1','Lumen_2','Lumen_3','Lumen_4','Mucosa_1','Mucosa_2','Mucosa_3','Mucosa_4','Rectum','Blood','Feces','HMO'};
for k=1:rct_num-1   %% blood????????bacteria??????Ymodel??????????????????????????
    Community_Result{k}.result=num2cell([NewTime,New_Ymodel_met{k}]);
    output=[col_header;Community_Result{k}.result];

    xlwrite(file_out_name,output,Tanks{k});   %%????????????????????????????
%     T=cell2table(Community_Result{k}.result,'VariableNames',col_header);
%     T.Properties.VariableNames=col_header;
%     writetable(T,file_out_name,'Sheet',Tanks{k});
end

k=9;
    Ymass_fece=NewY(:,(k-1)*lg+lg+1+1:(k-1).*lg+lg+1+lxini);
    Ymets{k+2}=num2cell([NewTime,Ymass_fece]);
    output=[col_header;Ymets{k+2}];
    xlwrite(file_out_name,output,Tanks{k+2});
%     T=Ymets{k+2};
%     T=cell2table(Ymets{k+2},'VariableNames',col_header);
%     writetable(T,file_out_name,'Sheet',Tanks{k+2});
    
    V5_rctm=NewY(:,(k-1)*lg+lg+1);  %%????????????????V5??????????????????????????????????????????????
    Volume=num2cell([NewTime,V5_rctm]);
    col_header={'Time','Fece_Vol'};
    output=[col_header;Volume];
%     T=Volume;
%     T=cell2table(Volume,'VariableNames',col_header);
%     writetable(T,file_out_name,'Sheet','Feces_Volume');
    xlwrite(file_out_name,output,'Feces_Volume');
k=10;
    Community_Result{k}.result=num2cell([NewTime,New_Ymodel{k}]);   %% ??10????????Ymets
    output=[col_header2;Community_Result{k}.result];
%     T=cell2table(Community_Result{k}.result,'VariableNames',col_header2);
%     writetable(T,file_out_name,'Sheet',Tanks{k});
    xlwrite(file_out_name,output,Tanks{k});  
k=11;



    
