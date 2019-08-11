function write_tank_table(Y,Ymodel,Ymets,met_co,lg,lxini,rct_num,ns,file_out_name)
curr_dir=pwd;
cd(curr_dir)
addpath('/Users/gejun/Desktop/Projects/Dynamic_Infant_Microbiome/Wuhao_CM_data/CM_diet_4cluster_20190301_return/poi_library');

javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');

% T=ReferenceStates.SimulationTime;
% Y=ReferenceStates.SimulationResult;
% Ymodel=ReferenceStates.Result_cell;
% Ymets=ReferenceStates.SimulationMetabolites;
% lg=ReferenceStates.ylg;
% rct_num=ReferenceStates.RectorNumber;
% met_co=ReferenceStates.met_udf_co;
% ns=ReferenceStates.SpeciesNumber;
% lxini=length(met_co);
% % HMO=ReferenceStates.PS_in;
% 
% %% sampling to 5000 points for drawing and PCA %%
% A=find(T>=0);
% B=find(T<=800);
% C=intersect(A,B);
% Conc_Time=T(C);
% Conc_Y=Y(C,:);
% for ij=1:9
%     Conc_Ymodel{ij}=Ymodel{ij}(C,:);
%     Conc_Ymodel_met{ij}=Ymets{ij}(C,:);
% end
% ij=10;
% Conc_Ymodel{ij}=Ymodel{ij}(C,:);
% 
% Tn=length(Conc_Time);
% nspc=ns;
% if nspc==5
%     SmpPnt=5;
% else          % nspc=8;
%     SmpPnt=1;
% end
% type='lumen';
% % for i=1:SmpPnt:Tn
%     NewTime=Conc_Time(1:SmpPnt:Tn);
%     NewY=Conc_Y(1:SmpPnt:Tn,:);
%     for ij=1:9
%         New_Ymodel{ij}=Conc_Ymodel{ij}(1:SmpPnt:Tn,:);
%         New_Ymodel_met{ij}=Conc_Ymodel_met{ij}(1:SmpPnt:Tn,:);
%     end
% ij=10;
% New_Ymodel{ij}=Ymodel{ij}(1:SmpPnt:Tn,:);


Tanks={'Lumen_1','Lumen_2','Lumen_3','Lumen_4','Mucosa_1','Mucosa_2','Mucosa_3','Mucosa_4','Rectum','Blood','Feces','HMO'};

Time=Y(:,1);
Y=Y(:,2:end);

col_header=['Time',met_co'];
met_co1=met_co(ns+1:end);
col_header2=['Time',met_co1'];

for k=1:rct_num-1   %% blood????????bacteria??????Ymodel??????????????????????????
    Community_Result{k}.result=num2cell([Time,Ymets{k}]);
    output=[col_header;Community_Result{k}.result];
    xlwrite(file_out_name,output,Tanks{k});   %%????????????????????????????
%     T=cell2table(Community_Result{k}.result,'VariableNames',col_header);
%     writetable(T,file_out_name,'Sheet',Tanks{k});
end

k=9;
    Ymass_fece=Y(:,(k-1)*lg+lg+1+1:(k-1).*lg+lg+1+lxini);
    Ymets{k+2}=num2cell([Time,Ymass_fece]);
    output=[col_header;Ymets{k+2}];
    xlwrite(file_out_name,output,Tanks{k+2});
%     T=cell2table(Community_Result{k+2}.result,'VariableNames',col_header);
%     writetable(T,file_out_name,'Sheet',Tanks{k+2});
    V5_rctm=Y(:,(k-1)*lg+lg+1);  %%????????????????V5??????????????????????????????????????????????
    Volume=num2cell([Time,V5_rctm]);
    col_header={'Time','Fece_Vol'};
    output=[col_header;Volume];
    xlwrite(file_out_name,output,'Feces_Volume');
%     T=cell2table(Volume,'VariableNames',col_header);
%     writetable(T,file_out_name,'Sheet','Feces_Volume');
k=10;
    Community_Result{k}.result=num2cell([Time,Ymodel{k}]);   %% ??10????????Ymets
    output=[col_header2;Community_Result{k}.result];
    xlwrite(file_out_name,output,Tanks{k});  
%     T=cell2table(Community_Result{k}.result,'VariableNames',col_header2);
%     writetable(T,file_out_name,'Sheet',Tanks{k});
k=11;
% xlswrite(file_out_name,HMO,Tanks{k+1})

