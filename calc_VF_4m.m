function [V5,F5]=calc_VF_4m(V5,F5,t,T,V_colon)
       if F5==0
            if (V5-V_colon(9))>=0
                F5=V5./T;
            end            

       else
            if (V5-0.0001)<=0
                V5=0.0001;
                F5=0;
            end
       end
        
        if (V5-0.0001)<=0
                V5=0.0001;
        end
