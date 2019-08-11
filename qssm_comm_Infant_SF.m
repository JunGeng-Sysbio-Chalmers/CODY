function dy=qssm_comm_Infant_SF(t,y,deltat,para,ns,FS,FV,V_colon,sub_index,hmo_in_y,enz_inx_y)
global  spc 
global SxZ_tot kmax KM ke alpha beta subs_cn biom_inx  met_udf met_udf_co n_carbon met_udf_co_new
global Y_BM Y_SUB maxEnzyme kl rct obsv V_obsv t_obsv    F5    kd  C_Kd   GLU_in     GLU_obsv   HMO_obsv
global lg_8tank  f_back   n_rct      km   km_w   Pn   rho  mol_mass
global HMO_in spc_hyd  kmax_hmo  KM_hmo  Y_hmo  HMO_const V5 
for i=1:length(y)
    if y(i)<0,y(i)=eps;end
end

ind_hex=find(strcmpi('Hexose',met_udf_co));
hmo_ix=find(strcmp('HMO',met_udf_co));
ylg=lg_8tank; 
lx=length(met_udf_co);
substrate_inx=sub_index;
enzyme_inx=enz_inx_y;
HMO_inx=hmo_in_y;
         V5=y(ylg*9+1);  %%
        if (mod(t,24)-0-0.01)>=0&&(mod(t,24)-1-0.01)<=0
            HMO_in=HMO_const;
        elseif (mod(t,24)-6-0.01)>=0&&(mod(t,24)-7-0.01)<=0
            HMO_in=HMO_const;
        elseif (mod(t,24)-13-0.01)>=0&&(mod(t,24)-14-0.01)<=0
            HMO_in=HMO_const;
        else
            HMO_in=0;
        end


HMO=y(HMO_inx);

BIOM_k=[];
rhyd=[];rhyd_hex_spc=[];
r_hyd_hmo={};r_hyd_hex={};
sub_mdf={};
substrate={};
for k=1:n_rct-1
    BIOM_k=[];
    BIOM_k=y((k-1).*ylg+1:(k-1).*ylg+ns);   
    rhyd(:,k)=kmax_hmo.*HMO(k).*BIOM_k./(KM_hmo.*BIOM_k+HMO(k)).*spc_hyd;    

    rhdrs=max(rhyd(:,k),zeros(ns,1));
    pu_hyd=max(rhdrs,zeros(ns,1)); pv_hyd=max(rhdrs,zeros(ns,1));
    sumpu_hyd=sum(pu_hyd); maxpv_hyd=max(pv_hyd);
    if sumpu_hyd>0, u_hyd{k}=pu_hyd/sumpu_hyd; else u_hyd{k}=zeros(ns,1); end
    if maxpv_hyd>0, v_hyd{k}=pv_hyd/maxpv_hyd; else v_hyd{k}=zeros(ns,1); end    
    
    r_hyd_hmo{k}=sum(rhyd(:,k));        
    rhyd_hex_spc(:,k)=Y_hmo.*rhyd(:,k); 
    r_hyd_hex{k}=sum(rhyd_hex_spc(:,k));   

   for s=1:ns
       substrate{k}{s}=y(substrate_inx{k}{s});
   end
   sub_add_total{k}=r_hyd_hex{k}.* deltat;

end
sub_mdf=substrate;
sub_add=zeros(ns,n_rct);
    
        sub_add=rhyd_hex_spc.* deltat;

    sub_add(find(sub_add<0))=0;


k=10;  
    r_hyd_hmo{k}=0;
    r_hyd_hex{k}=0;
    
    FT0=FV(1);FL0=FV(2);FM0=FV(3);FB0=FV(4); 
    f1=f_back(1);f2=f_back(2);f3=f_back(3);f4=f_back(4);
    VL1=V_colon(1);VL2=V_colon(2);VL3=V_colon(3);VL4=V_colon(4);
    VM1=V_colon(5);VM2=V_colon(6);VM3=V_colon(7);VM4=V_colon(8);
    V_rtm=V_colon(9); V_BLD=V_colon(10);    
    rho=0.7*1000;    
    JLT={};
    JL_temp=0;
    Mass_temp=0;  
    for k=1:4
        for i=ns+1:lx
            Mass_temp=Mass_temp+mol_mass(i-ns).*(y(i+(k-1)*ylg));       
            JL_temp=JL_temp+km(i-ns)*mol_mass(i-ns).*(y(i+(k-1)*ylg)-y(i+(k-1+4)*ylg));  
        end
        JLT{k}=JL_temp; 
        WL{k}=rho-Mass_temp; 
        FL_water{k}=km_w{k}*(WL{k}-0.1*WL{k});  

        Flow_L{k}=(JLT{k}+FL_water{k})*V_colon(k)/rho;   
        JL_temp=0;
        Mass_temp=0;
    end

    JMT={};
    JM_temp=[0];
    Mass_mucosa_temp=0;  
    icr_blood=9*ylg+1+lx;  

    for k=5:8
        for i=ns+1:lx
            Mass_mucosa_temp=Mass_mucosa_temp+mol_mass(i-ns).*(y(i+(k-1)*ylg));   
            JM_temp=JM_temp+Pn(i-ns)*mol_mass(i-ns)*(y(i+(k-1)*ylg)-y(i-ns+icr_blood));
        end
        JMT{k-4}=JM_temp;
        WML{k-4}=rho-Mass_mucosa_temp;
        FM_water{k-4}=Flow_L{k-4}-JMT{k-4}*V_colon(k)/rho;
        Flow_M{k-4}=FM_water{k-4}+JMT{k-4}*V_colon(k)/rho;
        JM_temp=[0];
        Mass_mucosa_temp=0;

    end

    FL1=FL0+f2-Flow_L{1}-f1;
    FL2=FL1-Flow_L{2}+f3-f2;
    FL3=FL2-Flow_L{3}+f4-f3;
    FL4=FL3-Flow_L{4}-f4;
    FM1=FM0+Flow_L{1}-Flow_M{1};
    FM2=FM1+Flow_L{2}-Flow_M{2};
    FM3=FM2+Flow_L{3}-Flow_M{3};
    FM4=FM3+Flow_L{4}-Flow_M{4};
    FL=[FL0;FL1;FL2;FL3;FL4];FM=[FM0;FM1;FM2;FM3;FM4];
    FM(find(FM<0))=0;    
    
 tc=24;   
 T=5/60;     
       if F5==0
            if (V5-0.2*V_colon(9))<0
                F5=0;
            else
                if ((mod(t,24)-12-0.0001)>=0)&&((mod(t,24)-12-T-0.0001)<=0)    
                   F5=V5./T;
                end
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

lxe=lg_8tank-lx; 
lx=length(met_udf_co);

u_co={};
v_co={};
rg_co={};
dxdt={};
re_co={};
for k=1:rct-1   
  
  for s=1:ns
    e{k}{s}=y(enzyme_inx{k}{s});
    e_rel=e{k}{s}./maxEnzyme{s};   
    rkin=kmax{s};
    re_s=ke{s};
    sub=sub_mdf{k}{s};   
    KS=KM{s};
    BIOM=y(s+(k-1).*ylg);  
    SxZ=SxZ_tot{s};
    klz=length(kmax{s});
    S=[];   
    for i=1:length(sub)  
      S(:,i)=sub(i).*ones(klz,1); 
    end      
    for ss=1:length(sub)  
        rkin=rkin.* S(:,ss) ./(KS(:,ss) + S(:,ss));  
        re_s=re_s.*S(:,ss) ./(KS(:,ss) + S(:,ss));
    end     
    r=rkin.*e_rel;
    rc=n_carbon{s};  
    rc=-rc.*Y_SUB{s};
    rc=sum(rc,2);  
    r_c=r.*rc;
    roi=e_rel.*rkin;
    roi=rc.*roi;
    roi=max(roi,zeros(kl{s},1));
    pu=max(roi,zeros(kl{s},1)); pv=max(roi,zeros(kl{s},1));
    sumpu=sum(pu); maxpv=max(pv);
    if sumpu>0, u=pu/sumpu; else u=zeros(kl{s},1); end
    if maxpv>0, v=pv/maxpv; else v=zeros(kl{s},1); end
    u_co{k}{s}=u;   
    v_co{k}{s}=v;   
    rM=v.*e_rel.*rkin;
    r=r.*(0.15+u_hyd{k}(s));
    re_s=re_s.*(0.15+u_hyd{k}(s));
    rg_co{k}{s}=sum(Y_BM{s}.*v.*r); 
    re_co{k}{s}=re_s;   
    diagV=diag(v);
    dxdt{k}{s}=SxZ*diagV*r*BIOM;  

  end
end 
    dy=zeros(length(y),1);
    lx=length(met_udf_co);
    dxdt_fnl=zeros(lx,1); %%
    for k=1:rct-1
        for s=1:ns
            mmcc=dxdt{k}{s}(ns+1:lx-1);
            dxdt{k}{s}(ns+2:lx)=mmcc;
            dxdt{k}{s}(ns+1)=0;
        end
    end
    for k=1:rct-1  
        for s=1:ns   
            dxdt_fnl(1:lx)=dxdt_fnl(1:lx)+dxdt{k}{s}(1:lx);  
        end
        dy(1+(k-1).*ylg:lx+(k-1).*ylg)=dy(1+(k-1).*ylg:lx+(k-1).*ylg)+dxdt_fnl(1:lx); 
        dxdt_fnl=zeros(lx,1);
        lx=length(met_udf_co);
    end

k=1;   
     icr=(k-1).*ylg;
     icr_pre=(k-2).*ylg;
     icr_next=(k).*ylg; 
     icr_m=(k-1+4)*ylg;
     detachment(:,k)=kd(k).*y(1+icr_m:ns+icr_m)./(C_Kd(k)+y(1+icr_m:ns+icr_m)).*y(1+icr_m:ns+icr_m).*V_colon(k+4)./V_colon(k);
     dy(1+icr:ns+icr)=dy(1+icr:ns+icr)-FL(k+1)/V_colon(k)*y(1+icr:ns+icr)+f_back(k+1)/V_colon(k).*y((1+icr_next):(ns+icr_next))+detachment(:,k);       
%      dy(ind_hex)=dy(ind_hex)+FL0/V_colon(k).*GLU_in-FL(k+1)/V_colon(k)*y(ind_hex)+f_back(k+1)/V_colon(k).*y(ind_hex+icr)-km(1)*(y(ind_hex)-y(ind_hex+icr_m));   
     dy(hmo_ix+icr)=dy(hmo_ix+icr)+FL0/V_colon(k).*HMO_in-FL(k+1)/V_colon(k)*y(hmo_ix+icr)+f_back(k+1)/V_colon(k).*y(hmo_ix+icr_next)-km(1)*(y(hmo_ix+icr)-y(hmo_ix+icr_m))-r_hyd_hmo{k};   
     dy(ind_hex+icr)=dy(ind_hex+icr)-FL(k+1)/V_colon(k)*y(ind_hex+icr)+f_back(k+1)/V_colon(k).*y(ind_hex+icr_next)-km(2)*(y(ind_hex+icr)-y(ind_hex+icr_m))+r_hyd_hex{k};  
     dy(ind_hex+1+icr:lx+icr)=dy(ind_hex+1+icr:lx+icr)-FL(k+1)/V_colon(k)*y(ind_hex+1+icr:lx+icr)+f_back(k+1)/V_colon(k).*y(ind_hex+1+icr_next:lx+icr_next)-km(3:end).*(y(ind_hex+1+icr:lx+icr)-y(ind_hex+1+icr_m:lx+icr_m));
 
   for k=2:3
     icr=(k-1).*ylg;
     icr_pre=(k-2).*ylg;
     icr_next=(k).*ylg;  
     icr_m=(k-1+4)*ylg;
     detachment(:,k)=kd(k).*y(1+icr_m:ns+icr_m)./(C_Kd(k)+y(1+icr_m:ns+icr_m)).*y(1+icr_m:ns+icr_m).*V_colon(k+4)./V_colon(k);
     dy(1+icr:ns+icr)=dy(1+icr:ns+icr)+FL(k)/V_colon(k).*y(1+icr_pre:ns+icr_pre)-FL(k+1)/V_colon(k)*y(1+icr:ns+icr)-f_back(k)/V_colon(k).*y(1+icr:ns+icr)+f_back(k+1)/V_colon(k).*y(1+icr_next:ns+icr_next)+detachment(:,k);        %% biomass of each species in the fluid  
     dy(hmo_ix+icr)=dy(hmo_ix+icr)+FL(k)/V_colon(k).*y(hmo_ix+icr_pre)-FL(k+1)/V_colon(k)*y(hmo_ix+icr)-f_back(k)/V_colon(k).*y(ind_hex+icr)+f_back(k+1)/V_colon(k).*y(hmo_ix+icr_next)-km(1)*(y(hmo_ix+icr)-y(hmo_ix+icr_m))-r_hyd_hmo{k};
     dy(ind_hex+icr)=dy(ind_hex+icr)+FL(k)/V_colon(k).*y(ind_hex+icr_pre)-FL(k+1)/V_colon(k)*y(ind_hex+icr)-f_back(k)/V_colon(k).*y(ind_hex+icr)+f_back(k+1)/V_colon(k).*y(ind_hex+icr_next)-km(2)*(y(ind_hex+icr)-y(ind_hex+icr_m))+r_hyd_hex{k};
     dy(ind_hex+1+icr:lx+icr)=dy(ind_hex+1+icr:lx+icr)+FL(k)/V_colon(k)*y(ind_hex+1+icr_pre:lx+icr_pre)-FL(k+1)/V_colon(k)*y(ind_hex+1+icr:lx+icr)-f_back(k)/V_colon(k).*y(ind_hex+1+icr:lx+icr)...
                              +f_back(k+1)/V_colon(k).*y(ind_hex+1+icr_next:lx+icr_next)-km(3:end).*(y(ind_hex+1+icr:lx+icr)-y(ind_hex+1+icr_m:lx+icr_m));
   end
k=4;    
   icr=(k-1).*ylg;
     icr_pre=(k-2).*ylg;
     icr_next=(k).*ylg;  
     icr_m=(k-1+4)*ylg;
     detachment(:,k)=kd(k).*y(1+icr_m:ns+icr_m)./(C_Kd(k)+y(1+icr_m:ns+icr_m)).*y(1+icr_m:ns+icr_m).*V_colon(k+4)./V_colon(k);

     dy(1+icr:ns+icr)=dy(1+icr:ns+icr)+FL(k)/V_colon(k).*y(1+icr_pre:ns+icr_pre)-FL(k+1)/V_colon(k)*y(1+icr:ns+icr)-f_back(k)/V_colon(k).*y(1+icr:ns+icr)+detachment(:,k);        
     dy(hmo_ix+icr)=dy(hmo_ix+icr)+FL(k)/V_colon(k).*y(hmo_ix+icr_pre)-FL(k+1)/V_colon(k)*y(hmo_ix+icr)-f_back(k)/V_colon(k).*y(ind_hex+icr)-km(1)*(y(hmo_ix+icr)-y(hmo_ix+icr_m))-r_hyd_hmo{k};
     dy(ind_hex+icr)=dy(ind_hex+icr)+FL(k)/V_colon(k).*y(ind_hex+icr_pre)-FL(k+1)/V_colon(k)*y(ind_hex+icr)-f_back(k)/V_colon(k).*y(ind_hex+icr)-km(2)*(y(ind_hex+icr)-y(ind_hex+icr_m))+r_hyd_hex{k};
     dy(ind_hex+1+icr:lx+icr)=dy(ind_hex+1+icr:lx+icr)+FL(k)/V_colon(k)*y(ind_hex+1+icr_pre:lx+icr_pre)-FL(k+1)/V_colon(k)*y(ind_hex+1+icr:lx+icr)-f_back(k)/V_colon(k).*y(ind_hex+1+icr:lx+icr)...
                              -km(3:end).*(y(ind_hex+1+icr:lx+icr)-y(ind_hex+1+icr_m:lx+icr_m));
 k=5;   %% VM1
 icr=(k-1).*ylg;
 icr_pre=(k-2).*ylg;
 icr_next=(k).*ylg; 
 icr_lm=(k-1-4).*ylg;
 icr_bld=9*ylg+1+lx;
     detachment(:,k)=kd(k-4).*y(1+icr:ns+icr)./(C_Kd(k-4)+y(1+icr:ns+icr)).*y(1+icr:ns+icr);
%      dy(1+icr:ns+icr)=dy(1+icr:ns+icr)-detachment(:,k-4).*V_colon(k-4)./V_colon(k);  
     dy(1+icr:ns+icr)=dy(1+icr:ns+icr)-detachment(:,k);   
     dy(hmo_ix+icr)=dy(hmo_ix+icr)+FM(k-4)/V_colon(k).*HMO_in-FM(k-4+1)/V_colon(k)*y(hmo_ix+icr)+V_colon(k-4)/V_colon(k).*km(1)*(y(hmo_ix+icr_lm)-y(hmo_ix+icr))-Pn(1)*(y(hmo_ix+icr)-y(hmo_ix-ns+icr_bld))-r_hyd_hmo{k};   %%????????????????????????????hmo??????????????????????????????hmo??????????????????g/L??????????????????????????????????????????????????????
     dy(ind_hex+icr)=dy(ind_hex+icr)-FM(k-4+1)/V_colon(k)*y(ind_hex+icr)+V_colon(k-4)/V_colon(k).*km(2)*(y(ind_hex+icr_lm)-y(ind_hex+icr))-Pn(2)*(y(ind_hex+icr)-y(ind_hex-ns+icr_bld))+r_hyd_hex{k};   %%????????????????????????????????FS(1)
     dy(ind_hex+1+icr:lx+icr)=dy(ind_hex+1+icr:lx+icr)-FM(k-4+1)/V_colon(k)*y(ind_hex+1+icr:lx+icr)+V_colon(k-4)/V_colon(k).*km(3:end).*(y(ind_hex+1+icr_lm:lx+icr_lm)-y(ind_hex+1+icr:lx+icr))-Pn(3:end).*(y(ind_hex+1+icr:lx+icr)-y(3+icr_bld:lx-ns+icr_bld));
  for k=6:8  
    icr=(k-1).*ylg;
    icr_pre=(k-2).*ylg;
    icr_next=(k).*ylg;
    icr_lm=(k-1-4).*ylg;
    icr_bld=9*ylg+1+lx;
%     dy(1+icr:ns+icr)=dy(1+icr:ns+icr)-kd(k-4).*y(1+icr:ns+icr);  
%     dy(1+icr:ns+icr)=dy(1+icr:ns+icr)-detachment(:,k-4).*V_colon(k-4)./V_colon(k);   
     detachment(:,k)=kd(k-4).*y(1+icr:ns+icr)./(C_Kd(k-4)+y(1+icr:ns+icr)).*y(1+icr:ns+icr);
     dy(1+icr:ns+icr)=dy(1+icr:ns+icr)-detachment(:,k);  
    dy(hmo_ix+icr)=dy(hmo_ix+icr)+FM(k-4)/V_colon(k).*y(hmo_ix+icr_pre)-FM(k-4+1)/V_colon(k)*y(hmo_ix+icr)+V_colon(k-4)/V_colon(k).*km(1)*(y(hmo_ix+icr_lm)-y(hmo_ix+icr))-Pn(1)*(y(hmo_ix+icr)-y(hmo_ix-ns+icr_bld))-r_hyd_hmo{k};   %%????????????????????????????hmo??????????????????????????????hmo??????????????????g/L??????????????????????????????????????????????????????

    dy(ind_hex+icr)=dy(ind_hex+icr)+FM(k-4)/V_colon(k).*y(ind_hex+icr_pre)-FM(k-4+1)/V_colon(k)*y(ind_hex+icr)+V_colon(k-4)/V_colon(k).*km(2)*(y(ind_hex+icr_lm)-y(ind_hex+icr))-Pn(2)*(y(ind_hex+icr)-y(1+icr_bld))+r_hyd_hex{k};   %%????????????????????????????????FS(1)
    dy(ind_hex+1+icr:lx+icr)=dy(ind_hex+1+icr:lx+icr)+FM(k-4)/V_colon(k).*y(ind_hex+1+icr_pre:lx+icr_pre)-FM(k-4+1)/V_colon(k)*y(ind_hex+1+icr:lx+icr)+V_colon(k-4)/V_colon(k).*km(3:end).*(y(ind_hex+1+icr_lm:lx+icr_lm)-y(ind_hex+1+icr:lx+icr))-Pn(3:end).*(y(ind_hex+1+icr:lx+icr)-y(3+icr_bld:lx-ns+icr_bld));
  end
 
 k=9;    
     icr=(k-1).*ylg;
     icr_pre=(k-1-1).*ylg;   
     icr_lm=(k-1-1-4).*ylg; 
     icr_next=(k).*ylg;

            
               dx_dt_V5_temp=dy(1+icr:ns+icr);dSP_dt_V5_temp=dy(hmo_ix+icr:lx+icr); 
               dy(1+icr:ns+icr)=dy(1+icr:ns+icr)-FM(k-4)/V5*y(1+icr:ns+icr)+FL(k-4)/V5*(y(1+icr_lm:ns+icr_lm)-y(1+icr:ns+icr));       
               dy(hmo_ix+icr)=dy(hmo_ix+icr)+FM(k-4)/V5.*(y(hmo_ix+icr_pre)-y(hmo_ix+icr))+FL(k-4)/V5*(y(hmo_ix+icr_lm)-y(hmo_ix+icr))-r_hyd_hmo{k};

               dy(ind_hex+icr)=dy(ind_hex+icr)+FM(k-4)/V5.*(y(ind_hex+icr_pre)-y(ind_hex+icr))+FL(k-4)/V5*(y(ind_hex+icr_lm)-y(ind_hex+icr))+r_hyd_hex{k};
               dy(ind_hex+1+icr:lx+icr)=dy(ind_hex+1+icr:lx+icr)+FM(k-4)/V5*(y(ind_hex+1+icr_pre:lx+icr_pre)-y(ind_hex+1+icr:lx+icr))+FL(k-4)/V5*(y(ind_hex+1+icr_lm:lx+icr_lm)-y(ind_hex+1+icr:lx+icr));
               dy(icr+ylg+1)=FL(k-4)+FM(k-4)-F5; 
               dy(icr+ylg+1+1:icr+ylg+1+ns)=V5*dx_dt_V5_temp+FL(k-4)*y(1+icr_lm:ns+icr_lm)-F5*y(1+icr:ns+icr);       
               dSP_dt_V5_temp(1)=dSP_dt_V5_temp(1)-r_hyd_hmo{k};  
               dSP_dt_V5_temp(2)=dSP_dt_V5_temp(2)+r_hyd_hex{k};
               dy(icr+ylg+1+ns+1:icr+ylg+1+lx)=V5*dSP_dt_V5_temp+FL(k-4)*y(ns+1+icr_lm:lx+icr_lm)+FM(k-4)*y(ns+1+icr_pre:lx+icr_pre)-F5*y(ns+1+icr:lx+icr);      
            

k=10;  
icr=(k-1).*ylg;
icr_M1=(k-4-1-1).*ylg;
icr_M2=(k-4-1-1+1).*ylg;
icr_M3=(k-4-1-1+2).*ylg;
icr_M4=(k-4-1-1+3).*ylg;
lx_bld=icr+1+lx;   
CB0=0.5.*ones(lx-ns,1);   
CB0(1)=0.1;   
% dy(lx_bld+1:lx_bld+lx-ns)=FB0/V_colon(k).*(CB0-y(lx_bld+1:lx_bld+lx-ns))+VM1/V_BLD.*Pn.*(y(ind_hex+icr_M1:lx+icr_M1)-CB0)+VM2/V_BLD.*Pn.*(y(ind_hex+icr_M2:lx+icr_M2)-CB0)...
%                                                                         +VM3/V_BLD.*Pn.*(y(ind_hex+icr_M3:lx+icr_M3)-CB0)+VM4/V_BLD.*Pn.*(y(ind_hex+icr_M4:lx+icr_M4)-CB0);
% dy(lx_bld+1:lx_bld+lx-ns)=FB0/V_colon(k).*(CB0-y(lx_bld+1:lx_bld+lx-ns))+VM1/V_BLD.*Pn(1).*(y(ind_hex+icr_M1:lx+icr_M1)-y(lx_bld+1:lx_bld+lx-ns))+VM2/V_BLD.*Pn(2).*(y(ind_hex+icr_M2:lx+icr_M2)-y(lx_bld+1:lx_bld+lx-ns))...
%                                                                         +VM3/V_BLD.*Pn(3).*(y(ind_hex+icr_M3:lx+icr_M3)-y(lx_bld+1:lx_bld+lx-ns))+VM4/V_BLD.*Pn(4).*(y(ind_hex+icr_M4:lx+icr_M4)-y(lx_bld+1:lx_bld+lx-ns));
dy(lx_bld+1:lx_bld+lx-ns)=FB0/V_colon(k).*(CB0-y(lx_bld+1:lx_bld+lx-ns))+VM1/V_BLD.*Pn.*(y(hmo_ix+icr_M1:lx+icr_M1)-y(lx_bld+1:lx_bld+lx-ns))+VM2/V_BLD.*Pn.*(y(hmo_ix+icr_M2:lx+icr_M2)-y(lx_bld+1:lx_bld+lx-ns))...
                                                                        +VM3/V_BLD.*Pn.*(y(hmo_ix+icr_M3:lx+icr_M3)-y(lx_bld+1:lx_bld+lx-ns))+VM4/V_BLD.*Pn.*(y(hmo_ix+icr_M4:lx+icr_M4)-y(lx_bld+1:lx_bld+lx-ns));


  kmax_new={};
  alpha_new={};
  beta_new={};
  for i=1:ns
      kmax_new{i}=kmax{i};
      kmax_new{i+ns}=kmax{i};
      alpha_new{i}=alpha{i};
      alpha_new{i+ns}=alpha{i};
      beta_new{i}=beta{i};
      beta_new{i+ns}=beta{i};
  end
            
            
  for k=1:rct-1
      icr=(k-1).*ylg;
      for s=1:ns        
        dy(lx+1+icr:lx+length(kmax{s})+icr)=alpha{s}+re_co{k}{s}.*u_co{k}{s}-(beta{s}+rg_co{k}{s}).*e{k}{s};
        lx=lx+length(kmax{s});
      end
     lx=length(met_udf_co);
  end