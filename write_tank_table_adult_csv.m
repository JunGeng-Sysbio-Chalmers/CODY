function write_tank_table_adult_csv(Y,Ymodel,Ymets,met_co,lg,lxini,rct_num,ns,file_out_name)
% curr_dir=pwd;
% cd(curr_dir)

Tanks={'Lumen_1','Lumen_2','Lumen_3','Lumen_4','Mucosa_1','Mucosa_2','Mucosa_3','Mucosa_4','Rectum','Blood','Feces','HMO'};

Time=Y(:,1);
Y=Y(:,2:end);

col_header=['Time',met_co'];
met_co1=met_co(ns+1:end);
col_header2=['Time',met_co1'];

for k=1:rct_num-1   %% blood????????bacteria??????Ymodel??????????????????????????
    Community_Result{k}.result=num2cell([Time,Ymets{k}]);
    output=[col_header;Community_Result{k}.result];
%     xlwrite(file_out_name,output,Tanks{k});   %%????????????????????????????
    T=cell2table(Community_Result{k}.result,'VariableNames',col_header);
    file_out=strcat(file_out_name,'_',Tanks{k},'.csv');
%     writetable(T,file_out_name,'Sheet',Tanks{k});
    writetable(T,file_out,'Delimiter',',');

end

k=9;
    Ymass_fece=Y(:,(k-1)*lg+lg+1+1:(k-1).*lg+lg+1+lxini);
    Ymets{k+2}=num2cell([Time,Ymass_fece]);
    output=[col_header;Ymets{k+2}];
%     xlwrite(file_out_name,output,Tanks{k+2});
    T=cell2table(Ymets{k+2},'VariableNames',col_header);
        file_out=strcat(file_out_name,'_',Tanks{k+2},'.csv');
    writetable(T,file_out,'Delimiter',',');
    V5_rctm=Y(:,(k-1)*lg+lg+1);  %%????????????????V5??????????????????????????????????????????????
    Volume=num2cell([Time,V5_rctm]);
    col_header={'Time','Fece_Vol'};
    output=[col_header;Volume];
%     xlwrite(file_out_name,output,'Feces_Volume');
    T=cell2table(Volume,'VariableNames',col_header);
%     writetable(T,file_out_name,'Sheet','Feces_Volume');
    file_out=strcat(file_out_name,'_','Feces_Volume','.csv');
        writetable(T,file_out,'Delimiter',',');


k=10;
    Community_Result{k}.result=num2cell([Time,Ymodel{k}]);   %% ??10????????Ymets
    output=[col_header2;Community_Result{k}.result];
%     xlwrite(file_out_name,output,Tanks{k});  
    T=cell2table(Community_Result{k}.result,'VariableNames',col_header2);
%     writetable(T,file_out_name,'Sheet',Tanks{k});
    file_out=strcat(file_out_name,'_',Tanks{k},'.csv');
    writetable(T,file_out,'Delimiter',',');

k=11;
% xlswrite(file_out_name,HMO,Tanks{k+1})

