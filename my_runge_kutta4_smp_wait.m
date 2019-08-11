function [t,y]=my_runge_kutta4(ufunc,y0,tspan,k,ns,FS,FV,V_colon,substrate_index,HMO_index,enzyme_index)    
t_ini=tspan(1);
t_final=tspan(2);
step=tspan(3);
N=floor((t_final-t_ini)/step);                                  
t(1)=t_ini;                                      
y(1,:)=y0';   
% h = waitbar(0,'Please wait...');
steps = N;

for ii=1:N
   t(ii+1)=t(ii)+step;
   k1=ufunc(t(ii),y(ii,:)',step,k,ns,FS,FV,V_colon,substrate_index,HMO_index,enzyme_index);
   k2=ufunc(t(ii)+step/2,y(ii,:)'+step*k1/2,step/2,k,ns,FS,FV,V_colon,substrate_index,HMO_index,enzyme_index);
   k3=ufunc(t(ii)+step/2,y(ii,:)'+step*k2/2,step/2,k,ns,FS,FV,V_colon,substrate_index,HMO_index,enzyme_index);
   k4=ufunc(t(ii)+step,y(ii,:)'+step*k3,step,k,ns,FS,FV,V_colon,substrate_index,HMO_index,enzyme_index);
   y(ii+1,:)=(y(ii,:)'+step*(k1+2*k2+2*k3+k4)/6)';      
%    waitbar(ii / N)
   progressbar(ii / N)
end
% close(h)
t=t';  