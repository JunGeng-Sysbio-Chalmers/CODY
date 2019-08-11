function [cm_met,SxZ_tot]=calc_SxZ_4m(met_udf_co,met_udf,SxZ,n_s)

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