     clc;clear;fclose all;
    [data,name]=xlsread('ESM.xls');
    [num,~]=find(~isnan(data(:,1)));
    opt=struct('Liner',1,'nonLiner',1,'Pariseau',1,'JPW',0,'Hoffman',1);
    [fid1,fid2,fid3,fid4,fid5]=wirte_txt_head(opt);
    for ii=1:43
        rock_ty=char(name(num(ii)+1,4));
        disp(rock_ty);
        if ii<length(num)
            exData=data(num(ii)+1:num(ii+1)-1,2:4);
            exData(isnan(exData(:,1))&isnan(exData(:,2))&isnan(exData(:,3))==1,:)=[];
        else
            exData=data(num(ii)+1:end,2:4);
            exData(isnan(exData(:,1))&isnan(exData(:,2))&isnan(exData(:,3))==1,:)=[];
        end
        exData=standard_type_change(exData);
        if opt.Liner==1
           theta=mugrpstats(exData);
           sub_empirical_equation(exData,fid1,rock_ty,1,theta)
           set(gcf,'unit','centimeters','position',[3,5,8,6]);
           set(gca,'FontSize',12);
           xlabel('\beta (бу)');
           ylabel('\sigma_{1(\beta)}(MPa)');
                title '';
           axis([0,90,0,260]);
           xticks([0 30 60 90])
        end
        if opt.nonLiner==1
           theta=mugrpstats(exData);
           sub_empirical_equation(exData,fid2,rock_ty,0,theta)
           set(gcf,'unit','centimeters','position',[3,5,8,6]);
           set(gca,'FontSize',12);
           xlabel('\beta (бу)')
           ylabel('\sigma_{1(\beta)}(MPa)')
                title '';
           axis([0,90,0,260]);
             xticks([0 30 60 90])
        end
        if opt.Pariseau==1
           sub_pariseau(exData,fid3,rock_ty);
           set(gcf,'unit','centimeters','position',[3,5,8,6]);
           set(gca,'FontSize',12);
           xlabel('\beta (бу)')
           ylabel('\sigma_{1(\beta)}(MPa)')
                title '';
           axis([0,90,0,260]);
              xticks([0 30 60 90])
        end
        if opt.JPW==1
           sub_JPW(exData,fid4,rock_ty);
            set(gcf,'unit','centimeters','position',[3,5,8,6]);
           set(gca,'FontSize',12);
           xlabel('\beta (бу)')
           ylabel('\sigma_{1(\beta)}(MPa)')
                title '';
           axis([0,90,0,260]);
             xticks([0 30 60 90])
        end
        if opt.Hoffman==1
            outdata=ones(10^6,3);
            for jj=1:10^6
                jj
                outdata(jj,:)=sub_hoffman(exData,fid5,rock_ty);
            end
%             set(gcf,'unit','centimeters','position',[3,5,8,6]);
%            set(gca,'FontSize',12);
%            xlabel('\beta (бу)')
%            ylabel('\sigma_{1(\beta)}(MPa)')
%            title '';
%            axis([0,90,0,260]);
%              xticks([0 30 60 90])
        end
    end
    fclose all;


function out=standard_type_change(exData)
    [angle_re_nan,~]=find(~isnan(exData(:,1)));
    for ii=1:length(angle_re_nan)
        if (ii<length(angle_re_nan))&&(angle_re_nan(ii+1)-angle_re_nan(ii))>1
             exData(angle_re_nan(ii)+1:angle_re_nan(ii+1)-1,1)=exData(angle_re_nan(ii),1);
        elseif (ii==length(angle_re_nan))&&(length(exData)-angle_re_nan(ii))>1
            exData(angle_re_nan(ii)+1:length(exData),1)=exData(angle_re_nan(ii),1);
        end
    end
    exData=exData(:,[2,1,3]);
    out=sortrows(exData,1);
end

function  [fid1,fid2,fid3,fid4,fid5]=wirte_txt_head(opt)
    fid1=0; 
    fid2=0; 
    fid3=0; 
    fid4=0; 
    fid5=0;
    if opt.Pariseau==1
        fid3=fopen([cd,'\statics\Pariseau.txt'],'w');
        fprintf(fid3,'%s\n','F,  G,  M,  U,  V ,  R_square, MAPE, RMSE, SAR');
    end
    if opt.Liner==1
        fid1=fopen([cd,'\statics\Liner.txt'],'w');
        fprintf(fid1,'%s\n','FrictionAngle, sigmaC0, sigmaC90, sigmaCmin, m,  n,  R_square, MAPE, RMSE, SAR');
    end
    if opt.nonLiner==1
        fid2=fopen([cd,'\statics\nonLiner.txt'],'w');
        fprintf(fid2,'%s\n','FrictionAngle, sigmaC0, sigmaC90, sigmaCmin, m,  n, k,  R_square, MAPE, RMSE, SAR');
    end
    if opt.JPW==1
        fid4=fopen([cd,'\statics\JPW.txt'],'w');
        fprintf(fid4,'%s\n','FrictionAngle_joint, cohensive_joint, FrictionAngle_rock, cohensive_rock,  R_square, MAPE, RMSE, SAR');
    end
    if opt.Hoffman==1
        fid5=fopen([cd,'\statics\Hoffman.txt'],'w');
        fprintf(fid5,'%s\n','F, G, M, U, V, R_square, MAPE, RMSE, SAR');
    end
end     

function  theta=mugrpstats(exData)
     nums=unique(exData(:,2));
     av=zeros(1,length(nums));
     for i=1:length(nums)
         av(i)=mean(exData((exData(:,2)==nums(i)),3));
     end
     [~,id]=min(av);
     if nums(id)>=45&&nums(id)<=75
         theta=nums(id);
     else
         theta=60;
     end
end

