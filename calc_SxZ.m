function [cm_met,SxZ_tot]=calc_SxZ(met_udf_co,met_udf,SxZ,n_s)

for i=1:n_s
    for j=1:length(met_udf_co)
        temp=strmatch(met_udf_co{j},met_udf{i},'exact');
        if (isempty(temp))   
            cm_met{i}{j}=0;   
            SxZ_tot{i}(j,:)=zeros(1,length(SxZ{i}(1,:)));
        else
            cm_met{i}{j}=temp;  
            SxZ_tot{i}(j,:)=SxZ{i}(temp,:);
        end
    end
end
