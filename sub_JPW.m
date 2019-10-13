function sub_JPW(exData,fid,name)
opt=struct('u1',1,'d1',100,'u2',-10,'d2',10,'u3',10,'d3',55,...
           'u4',-20,'d4',10);
cof=Gencofmatrix(opt);
num1=find(exData(:,2)==90);
num2=find(exData(:,2)==0);
[~,num3]=min(exData(:,3));
sigmat90=(exData(num1(1),3)-exData(num1(1),1)*((1+sin(cof(:,3)*pi/180))./(1-sin(cof(:,3)*pi/180))));
sigmat0=(exData(num2(1),3)-exData(num2(1),1)*((1+sin(cof(:,3)*pi/180))./ (1-sin(cof(:,3)*pi/180))));
sigmaMin=(exData(num3(1),3)-exData(num3(1),1)*((1+sin(cof(:,3)*pi/180))./(1-sin(cof(:,3)*pi/180))));
cohensive_joint_ad=0.2*sigmaMin+cof(:,2);
cohensive_rock_ad=(sigmat0+sigmat90)/4+cof(:,4);
delete=(cohensive_joint_ad<=0|cohensive_rock_ad<=0);
use_cof=cof;
use_cof(:,[2,4])=[cohensive_joint_ad,cohensive_rock_ad];
use_cof(delete>0,:)=[];
for cilcle=1:10
   if cilcle==10
       Out_Flag=1;
   else
       Out_Flag=0;
   end
   [RMSE,id,R_square,MAPE]=Find_value_id(exData,use_cof,Out_Flag,name);
    my=use_cof(id(1:2000),:);
    for temp1=1:4
        str=num2str(temp1);
        eval(['opt.u',str,'=min(my(:,',str,'));opt.d',str,'=max(my(:,',str,'));'])
    end
    if cilcle<10
       use_cof=Gencofmatrix(opt);
    end
end
  disp(num2str(MAPE));
id=id(1);
strength1=calstrength1(0,use_cof,id);
SAR=max(strength1)/min(strength1);
sigma3=unique(exData(:,1));
 Colors = linspecer(length(sigma3));
for num2=1:length(sigma3)
    strength3=sigma3(num2);
    strength1=calstrength1(strength3,use_cof,id);
    plot(0:2:90,strength1,'color',Colors(num2,:));
    hold on;
    index=find(exData(:,1)==strength3);
      scatter(exData(index,2),exData(index,3),5,'o','MarkerEdgeColor',Colors(num2,:));
    %plot(exData(index,2),exData(index,3),'o','color',Colors(num2,:));
end
   title(['JPW-',strrep(name,'_','-')]);
   fprintf(fid,'%s\n',name);
   fprintf(fid,'%f  %f  %f  %f  %f  %f   %f  %f\n',[use_cof(id,:),R_square,MAPE,RMSE(1),SAR]);
   saveas(gcf,[cd,'\png\','JPW',name([end-1,end]),'.png']);
   hold off;
end
%clear all;
function cof=Gencofmatrix(opt)
   FrictionAngle_joint=mylinspace(opt.u1,opt.d1) ;%9
   cohensive_joint=mylinspace(opt.u2,opt.d2) ;%11
   FrictionAngle_rock=mylinspace(opt.u3,opt.d3);%11
   cohensive_rock=mylinspace(opt.u4,opt.d4);%11
   num=length(FrictionAngle_joint)*length(cohensive_joint)*length(FrictionAngle_rock)*length(cohensive_rock);
   [x,y,z,w]=ndgrid(FrictionAngle_joint,cohensive_joint,FrictionAngle_rock,cohensive_rock);
   cof= [reshape(x,num,1),reshape(y,num,1),reshape(z,num,1),reshape(w,num,1)];
end
function out=mylinspace(a,b)
   if a==b
       out=a;
   else 
       out=linspace(a,b,20);
   end
end
function [RMSE,id,R_square,MAPE]=Find_value_id(exData,use_cof,Out_Flag,name)
   num=length(use_cof);
   allreslut=ones(num,length(exData));
   for cc=1:length(exData)
             RockMC=exData(cc,1)*(1+sin(use_cof(:,3)*pi/180))./(1-sin(use_cof(:,3)*pi/180))+(2*use_cof(:,4).*cos(use_cof(:,3)*pi/180))./(1-sin(use_cof(:,3)*pi/180));
             JointMC=exData(cc,1)+2*(use_cof(:,2)+exData(cc,1)*tan(use_cof(:,1)*pi/180))./(sin(2*(exData(cc,2)*pi/180)).*(1-tan(use_cof(:,1)*pi/180).*cot(exData(cc,2)*pi/180)));
             JointMC=(JointMC<=0|isnan(JointMC))*10^50+(JointMC>0).*JointMC;
             allreslut(:,cc)=min(RockMC,JointMC);
   end
   RMSE=allreslut-exData(:,3)';
   RMSE=RMSE.^2;
   [RMSE,id]=sort((mean(RMSE,2)).^0.5);
   R_square=sum((allreslut(id(1),:)-mean(exData(:,3))).^2)/sum((exData(:,3)-mean(exData(:,3))).^2);
   MAPE=mean(abs(allreslut(id(1),:)-exData(:,3)')./exData(:,3)')*100;
   if Out_Flag==1
       outdata=[allreslut(id(1),:);exData(:,3)']';
       save([cd,'\mat\','JPW_',name([end-1,end]),'.mat'],'outdata');
   end
end
function strength1=calstrength1(strength3,use_cof,id)
  tt=0:2:90;
  disp('JPW');
  use_cof(id,:)
  RockMC=strength3*(1+sin(use_cof(id,3)*pi/180))./(1-sin(use_cof(id,3)*pi/180))+(2*use_cof(id,4).*cos(use_cof(id,3)*pi/180))./(1-sin(use_cof(id,3)*pi/180));
  JointMC=strength3+2*(use_cof(id,2)+strength3*tan(use_cof(id,1)*pi/180))./(sin(2*(tt*pi/180)).*(1-tan(use_cof(id,1)*pi/180)*cot(tt*pi/180)));
  JointMC=(JointMC<=0)*10^50+(JointMC>0).*JointMC;
  JointMC(isnan(JointMC))=10^50;
  strength1= (JointMC<RockMC).*JointMC+(JointMC>=RockMC)*RockMC;
end