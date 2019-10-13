function sub_empirical_equation(exData,fid,name,isliner,theta)
opt=struct('u1',1,'d1',55,'u2',-50,'d2',50,'u3',-50,'d3',50,...
           'u4',-50,'d4',50,'u5',1,'d5',10,'u6',1,'d6',7,'u7',1,'d7',5,'liner',1);
opt.liner=isliner;
cof=Gencofmatrix(opt);
num1=find(exData(:,2)==90);
num2=find(exData(:,2)==0);
[~,num3]=min(exData(:,3));
sigmat90=(exData(num1(1),3)-exData(num1(1),1)*((1+sin(cof(:,1)*pi/180))./(1-sin(cof(:,1)*pi/180))))+cof(:,3);
sigmat0=(exData(num2(1),3)-exData(num2(1),1)*((1+sin(cof(:,1)*pi/180))./ (1-sin(cof(:,1)*pi/180))))+cof(:,2);
sigmaMin=(exData(num3(1),3)-exData(num3(1),1)*((1+sin(cof(:,1)*pi/180))./(1-sin(cof(:,1)*pi/180))))+cof(:,4);
use_cof=cof;
use_cof(:,2:4)=[sigmat0,sigmat90,sigmaMin];
use_cof(sigmaMin<0,:)=[];
nums=[40000,20000,10000,5000,2500,1250,625,500,400,200];
for cilcle=1:10
    if cilcle==10
       Out_Flag=1;
    else
       Out_Flag=0;
    end
    [RMSE,id,R_square,MAPE]=Find_value_id(exData,use_cof,theta,opt,Out_Flag,name);
            %RMSE(1)
     my=use_cof(id(1:nums(cilcle)),:);
     for temp1=1:6
          str=num2str(temp1);
          eval(['opt.u',str,'=min(my(:,',str,'));opt.d',str,'=max(my(:,',str,'));'])
     end
     if cilcle<10
         use_cof=Gencofmatrix(opt);
     end
end
id=id(1);
sigma3=unique(exData(:,1));
 Colors = linspecer(length(sigma3));
strength1=calstrength1(0,theta,use_cof,id,opt);
SAR=max(strength1)/min(strength1);
for num2=1:length(sigma3)
    strength3=sigma3(num2);
    if strength3>0
       strength1=calstrength1(strength3,theta,use_cof,id,opt);
    end
    plot(0:2:90,strength1,'color',Colors(num2,:));
    hold on;
    index=find(exData(:,1)==strength3);
    %scatter(datax',data(ii,:),60,char(pointStyle{2}),'MarkerEdgeColor',Colors(1,:));
    scatter(exData(index,2),exData(index,3),5,'o','MarkerEdgeColor',Colors(num2,:));
end
fprintf(fid,'%s\n',name);
if opt.liner==1
    title(['liner-',strrep(name,'_','-')]);
    fprintf(fid,'%f  %f  %f  %f  %f  %f   %f  %f   %f  %f\n',[use_cof(id,:),R_square,MAPE,RMSE(1),SAR]);
    saveas(gcf,[cd,'\png\','liner_',name([end-1,end]),'.png'])
else
    title(['non-',strrep(name,'_','-')]);
    fprintf(fid,'%f  %f  %f  %f  %f  %f   %f  %f   %f   %f  %f\n',[use_cof(id,:),R_square,MAPE,RMSE(1),SAR]);
    saveas(gcf,[cd,'\png\','nonliner_',name([end-1,end]),'.png'])
end
hold off;
end
%clear all;
function cof=Gencofmatrix(opt)
   FrictionAngle=mylinspace(opt.u1,opt.d1) ;%9
   sigmaC0=mylinspace(opt.u2,opt.d2) ;%11
   sigmaC90=mylinspace(opt.u3,opt.d3);%11
   sigmaCmin=mylinspace(opt.u4,opt.d4);%11
   m=mylinspace(opt.u5,opt.d5);%
   n=mylinspace(opt.u6,opt.d6);%
   if opt.liner==1
       num=length(FrictionAngle)*length(sigmaC0)*length(sigmaC90)*length(sigmaCmin)*length(m)*length(n);
       [x,y,z,w,g,h]=ndgrid(FrictionAngle,sigmaC0,sigmaC90,sigmaCmin,m,n);
       cof= [reshape(x,num,1),reshape(y,num,1),reshape(z,num,1),reshape(w,num,1),...
       reshape(g,num,1),reshape(h,num,1)];
   else
        k=mylinspace(opt.u7,opt.d7);
        num=length(FrictionAngle)*length(sigmaC0)*length(sigmaC90)*length(sigmaCmin)*length(m)*length(n)*length(k);
        [x,y,z,w,g,h,I]=ndgrid(FrictionAngle,sigmaC0,sigmaC90,sigmaCmin,m,n,k);
        cof= [reshape(x,num,1),reshape(y,num,1),reshape(z,num,1),reshape(w,num,1),...
              reshape(g,num,1),reshape(h,num,1),reshape(I,num,1)];
   end
end
function out=mylinspace(a,b)
   if a==b
       out=a;
   else 
       out=linspace(a,b,7);
   end
end
function [RMSE,id,R_square,MAPE]=Find_value_id(exData,use_cof,theta,opt,Out_Flag,name)
   num=length(use_cof);
   allreslut=ones(num,length(exData));
   for cc=1:length(exData)
       if opt.liner==1
          if(exData(cc,2)<theta)
             allreslut(:,cc)=exData(cc,1)*((1+sin(use_cof(:,1)*pi/180))./(1-sin(use_cof(:,1)*pi/180)))+use_cof(:,2)-(use_cof(:,2)-use_cof(:,4)).*sin((exData(cc,2)*90/theta)*pi/180).^use_cof(:,5);
          else
             allreslut(:,cc)=exData(cc,1)*((1+sin(use_cof(:,1)*pi/180))./(1-sin(use_cof(:,1)*pi/180)))+use_cof(:,3)-(use_cof(:,3)-use_cof(:,4)).*cos(((exData(cc,2)-theta)*90/(90-theta))*pi/180).^use_cof(:,6);
          end
       else
            if(exData(cc,2)<theta)
                UCS=use_cof(:,2)-(use_cof(:,2)-use_cof(:,4)).*sin((exData(cc,2)*90/theta)*pi/180).^use_cof(:,5);
            else
                UCS=use_cof(:,3)-(use_cof(:,3)-use_cof(:,4)).*cos(((exData(cc,2)-theta)*90/(90-theta))*pi/180).^use_cof(:,6);
            end
            max_UCS=max(use_cof(:,2),use_cof(:,3));
            detaA = exData(cc,1)*((1+sin(use_cof(:,1)*pi/180))./(1-sin(use_cof(:,1)*pi/180)))-exData(cc,1)^2*(sin(use_cof(:,1)*pi/180))./((1-sin(use_cof(:,1)*pi/180)).*use_cof(:,7).*max_UCS);
            detaB = exData(cc,1)+((use_cof(:,7).*max_UCS.*sin(use_cof(:,1)*pi/180))./(1-sin(use_cof(:,1)*pi/180)));
            allreslut(:,cc)=(exData(cc,1)<=use_cof(:,7).*(max_UCS)).*detaA+(exData(cc,1)>use_cof(:,7).*(max_UCS)).*detaB+UCS;
       end
   end
   %RMSE 
   
   RMSE=allreslut-exData(:,3)';
   RMSE=RMSE.^2;
   [RMSE,id]=sort((mean(RMSE,2)).^0.5);
   % 
    R_square=1-sum((allreslut(id(1),:)-exData(:,3)').^2)/sum((exData(:,3)-mean(exData(:,3))).^2);
   % 
   MAPE=mean(abs(allreslut(id(1),:)-exData(:,3)')./exData(:,3)')*100;
   if Out_Flag==1
       outdata=[allreslut(id(1),:);exData(:,3)']';
       if opt.liner==1
       save([cd,'\mat\','liner_',name([end-1,end]),'.mat'],'outdata');
       else
       save([cd,'\mat\','nonliner_',name([end-1,end]),'.mat'],'outdata');
       end
   end
end
function strength1=calstrength1(strength3,theta,use_cof,id,opt)
    tt=0:2:90;
    if opt.liner==1
         strength1= (tt<theta).*((1+sin(use_cof(id,1)*pi/180))*strength3/(1-sin(use_cof(id,1)*pi/180))+use_cof(id,2)-(use_cof(id,2)-use_cof(id,4)).*sin((tt*90/theta)*pi/180).^use_cof(id,5))+...
                   (tt>=theta).*((1+sin(use_cof(id,1)*pi/180))*strength3/(1-sin(use_cof(id,1)*pi/180))+use_cof(id,3)-(use_cof(id,3)-use_cof(id,4)).*cos(((tt-theta)*90/(90-theta))*pi/180).^use_cof(id,6));
    else
        max_UCS=max(use_cof(id,2),use_cof(id,3));
        detaA = strength3*((1+sin(use_cof(id,1)*pi/180))./(1-sin(use_cof(id,1)*pi/180)))-strength3^2*(sin(use_cof(id,1)*pi/180))./((1-sin(use_cof(id,1)*pi/180)).*use_cof(id,7).*max_UCS);
        detaB = strength3+((use_cof(id,7).*max_UCS.*sin(use_cof(id,1)*pi/180))./(1-sin(use_cof(id,1)*pi/180)));
        UCS= (tt<theta).*(use_cof(id,2)-(use_cof(id,2)-use_cof(id,4)).*sin((tt*90/theta)*pi/180).^use_cof(id,5))+...
            (tt>=theta).*(use_cof(id,3)-(use_cof(id,3)-use_cof(id,4)).*cos(((tt-theta)*90/(90-theta))*pi/180).^use_cof(id,6));
        strength1=(strength3<=use_cof(id,7).*(max_UCS)).*detaA+(strength3>use_cof(id,7).*(max_UCS)).*detaB+UCS;
    end
end