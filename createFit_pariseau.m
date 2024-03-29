function [fitresult, gof] = createFit_pariseau(exData)
%CREATEFIT(X,Y,Z)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Input : y
%      Z Output: z
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 09-Apr-2019 22:57:21 自动生成
%% Fit: 'untitled fit 1'.
x=exData(:,1);y=exData(:,2);z=exData(:,3);
[xData, yData, zData] = prepareSurfaceData( x, y, z );
% Set up fittype and options.
ft = fittype( '(1+x*((1/u)+2*(1/v)))./((((1/f^2)*(sin(y*pi/180))^4)+(1/g^2)*((cos(y*pi/180))^4+((cos(2*y*pi/180))^2))+0.25*(1/m^2)*((sin(2*y*pi/180)))^2)^0.5-((1/u)*(cos(y*pi/180))^2+(1/v)*(sin(y*pi/180))^2))+x', 'independent', {'x', 'y'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [5 5 5 5 5];
opts.MaxIter = 400;
opts.Robust = 'LAR';
opts.StartPoint =[100,100,100,100,100];
opts.Upper = [5000 5000 5000 5000 5000];

% Fit model to data.
[fitresult, gof] = fit( [xData, yData], zData, ft, opts );
if (gof.rmse>100)
    opt=struct('u1',5,'d1',2000,'u2',5,'d2',500,'u3',5,'d3',500,...
               'u4',50,'d4',500,'u5',50,'d5',2000,'liner',1);
    use_cof=Gencofmatrix(opt);
    inital=Find_value_id(exData,use_cof);
    opts.StartPoint =inital;
    [fitresult, gof] = fit( [xData, yData], zData, ft, opts );
end
while (gof.rmse>100)
    inital=inital+rand(1,5);
    opts.StartPoint =inital;
    [fitresult, gof] = fit( [xData, yData], zData, ft, opts );
end
end

function out=Find_value_id(exData,use_cof)
   num=length(use_cof);
   allreslut=ones(num,length(exData)); 
   F=1./(use_cof(:,1).^2);G=1./(use_cof(:,2).^2);M=1./(use_cof(:,3).^2);
   U=1./use_cof(:,4);V=1./use_cof(:,5);
   remove=zeros(num,1);
   for cc=1:length(exData)
        t=exData(cc,2);
        xPos=((F*(sin(t*pi/180)).^4)+G*((cos(t*pi/180)).^4+((cos(2*t*pi/180)).^2))+0.25*M*((sin(2*t*pi/180)).^2)).^0.5-(U*(cos(t*pi/180)).^2+V*(sin(t*pi/180)).^2); 
        remove(xPos<0)=1;
        allreslut(:,cc)=(1+exData(cc,1)*(U+2*V))./xPos;
   end
   allreslut(remove==1,:)=[];
   use_cof(remove==1,:)=[];
   RMSE=allreslut-exData(:,3)';
   RMSE=RMSE.^2;
   [~,id]=sort((mean(RMSE,2)).^0.5);
   out=use_cof(id(1),:);
    %R_square=sum((allreslut(id(1),:)-mean(exData(:,3))).^2)/sum((exData(:,3)-mean(exData(:,3))).^2);
end
function cof=Gencofmatrix(opt)
   f=mylinspace(opt.u1,opt.d1) ;%9
   g=mylinspace(opt.u2,opt.d2) ;%11
   m=mylinspace(opt.u3,opt.d3);%11
   u=mylinspace(opt.u4,opt.d4);%11
   v=mylinspace(opt.u5,opt.d5);%11
   
   num=length(f)*length(g)*length(u)*length(v)*length(m);
   [x,y,z,w,h]=ndgrid(f,g,m,u,v);
   cof= [reshape(x,num,1),reshape(y,num,1),reshape(z,num,1),reshape(w,num,1),reshape(h,num,1)];
end
function out=mylinspace(a,b)
   if a==b
       out=a;
   else 
       out=linspace(a,b,20);
   end
end