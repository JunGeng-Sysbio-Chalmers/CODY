function write_tank_table_infant_csv(Y,Ymodel,Ymets,met_co,lg,lxini,rct_num,ns,file_out_name)
% curr_dir=pwd;
% cd(curr_dir)


Tanks={'Lumen_1','Lumen_2','Lumen_3','Lumen_4','Mucosa_1','Mucosa_2','Mucosa_3','Mucosa_4','Rectum','Blood','Feces','HMO'};
SmpPnt=2;

Time=Y(1:SmpPnt:end,1);
Y=Y(1:SmpPnt:end,2:end);

col_header=['Time',met_co'];
met_co1=met_co(ns+1:end);
col_header2=['Time',met_co1'];

for k=1:rct_num-1   %% blood????????bacteria??????Ymodel??????????????????????????
    Ymets_new=Ymets{k}(1:SmpPnt:end,:);
    Community_Result{k}.result=num2cell([Time,Ymets_new]);
    output=[col_header;Community_Result{k}.result];
%     xlwrite(file_out_name,output,Tanks{k});   %%????????????????????????????
    T=cell2table(Community_Result{k}.result,'VariableNames',col_header);
        file_out=strcat(file_out_name,'_',Tanks{k},'.csv');

    writetable(T,file_out,'Delimiter',',');
end

k=9;
    Ymass_fece=Y(:,(k-1)*lg+lg+1+1:(k-1).*lg+lg+1+lxini);
    Ymets_new2=num2cell([Time,Ymass_fece]);
    output=[col_header;Ymets_new2];
%     xlwrite(file_out_name,output,Tanks{k+2});
    T=cell2table(Ymets_new2,'VariableNames',col_header);
            file_out=strcat(file_out_name,'_',Tanks{k+2},'.csv');
    writetable(T,file_out,'Delimiter',',');
    V5_rctm=Y(:,(k-1)*lg+lg+1);  %%????????????????V5??????????????????????????????????????????????
    Volume=num2cell([Time,V5_rctm]);
    col_header={'Time','Fece_Vol'};
    output=[col_header;Volume];
%     xlwrite(file_out_name,output,'Feces_Volume');
    T=cell2table(Volume,'VariableNames',col_header);
        file_out=strcat(file_out_name,'_','Feces_Volume','.csv');

    writetable(T,file_out,'Delimiter',',');
k=10;
    Ymodel_new=Ymodel{k}(1:SmpPnt:end,:);
    Community_Result{k}.result=num2cell([Time,Ymodel_new]);   %% ??10????????Ymets
    output=[col_header2;Community_Result{k}.result];
%     xlwrite(file_out_name,output,Tanks{k});  
    T=cell2table(Community_Result{k}.result,'VariableNames',col_header2);
        file_out=strcat(file_out_name,'_',Tanks{k},'.csv');

    writetable(T,file_out,'Delimiter',',');
k=11;
% xlswrite(file_out_name,HMO,Tanks{k+1})

