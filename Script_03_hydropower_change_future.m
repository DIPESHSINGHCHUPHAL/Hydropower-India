%%% for mahi basin dams

clear all;
clc;
dam_data = importdata('/media/wcl/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/dam_params/dam_params/damparam_mahi_06min_2020.csv');
dam_lalo_sto = dam_data.data;
dam_name = dam_data.textdata;
per={'hist';'near_ssp126';'near_ssp585';'mid_ssp126';'mid_ssp585';'end_ssp126';'end_ssp585'};
sce={'historical';'ssp126';'ssp585';'ssp126';'ssp585';'ssp126';'ssp585'};
for sss=1:7
    %% dam
    if sss==1
        year=1995:2014;
    elseif sss==2
        year = 2021:2040;
    elseif sss==3
        year = 2021:2040;
    elseif sss==4
        year = 2041:2060;
    elseif sss==5
        year = 2041:2060;
    elseif sss==6
        year = 2081:2100;
    elseif sss==7
        year = 2081:2100;
    end
    
    Thp=[];monthly_all2=[];
    for mm=1:5
        f2={'INM-CM5-0','MIROC6','NorESM2-MM','TaiESM1','BCC-CSM2-MR','hist'};
        f1='mahi_';
        folder=horzcat(f1,f2{mm});
        for k=[2]
            final1=[];
            for j = 1:size(year,2)
                disp(year(1,j))
                final=[];
                lat = ncread(['/media/wcl/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'lat');
                lon = ncread(['/media/wcl/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'lon');
                lat=double(lat); lat=round(lat,3);
                lon=double(lon); lon=round(lon,3);
                sto = ncread(['/media/wcl/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'damsto');
                
                t = datetime(year(1,j),01,01):datetime(year(1,j),12,31);
                t=t';
                [y, m, d] =ymd(t);
                time = [y m d];
                
                out = sto(dam_lalo_sto(k,4),dam_lalo_sto(k,5),:);
                final=[time squeeze(out)];
                final1=vertcat(final1,final);
            end
            dam_sto=final1;  %dam_sto total
            dam_sto(:,4)=final1(:,4)-(dam_lalo_sto(k,7)*10^6/2);  %dam_sto 
            
            %% inflow
            data=[];
            t = datetime(year(1,1),01,01):datetime(year(1,end),12,31);
            t=t';
            [y, m, d] =ymd(t);
            time = [y m d];
            
            for j = 1:size(year,2)
                disp(year(1,j))
                inflow = ncread(['/media/wcl/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_daminf',num2str(year(1,j)),'.nc'],'daminf');
                
                out = squeeze(inflow(dam_lalo_sto(k,4),dam_lalo_sto(k,5),:));
                data=vertcat(data,out);
            end
            inflw_data=[time data];
            
            %% %% dam sto monthly mean
            S=[];
            fd=[];
            for i=1:size(year,2)
                for mm=1:12
                    if mm<=9
                        fd=find(dam_sto(:,1)==year(1,i)+1 & dam_sto(:,2)==mm);
                        S(mm,i)= mean(dam_sto(fd,4));
                    else
                        fd=find(dam_sto(:,1)==year(1,i) & dam_sto(:,2)==mm);
                        S(mm,i)= mean(dam_sto(fd,4));
                        
                    end
                end
            end
            S=vertcat(S(10:12,1:end-1),S(1:9,1:end-1)) ;

            %storage at beggining of year
            S_beg=[];
            fd=[];
            for yy=1:size(year,2)
                fd=find(dam_sto(:,1)==year(1,yy) & dam_sto(:,2)==10 &dam_sto(:,3)==1);
                S_beg(yy,1)= (dam_sto(fd,4));
            end
            S_beg=S_beg(1:end-1,1);
            
            %mean_annual_inflow
            monthly_inflw=[];
            fd=[];
            for i=1:size(year,2)
                for mm=1:12
                    if mm<=9
                        fd=find(inflw_data(:,1)==year(1,i)+1 & inflw_data(:,2)==mm);
                        monthly_inflw(mm,i)= mean(inflw_data(fd,4));
                    else
                        fd=find(inflw_data(:,1)==year(1,i) & inflw_data(:,2)==mm);
                        monthly_inflw(mm,i)= mean(inflw_data(fd,4));
                        
                    end
                end
            end
            im=vertcat(monthly_inflw(10:12,1:end-1),monthly_inflw(1:9,1:end-1)) ;im_mean=mean(im,2);im_mean2=vertcat(im_mean(9:12,:),im_mean(1:8,:));
            
            monthly_all(:,k)=im_mean2;
            ia=mean(monthly_inflw,1);
            ia=ia(1,1:end-1);
            Ia= ia*(24*3600*365);
            
            max_Hp_all=[140 240];
            alpha=0.85;
            C_all=dam_lalo_sto(:,8);
            H_all = [74.5 66];
            C=C_all(k,1)*10^6-(dam_lalo_sto(k,7)*10^6/2);                %m3
            max_Hp=max_Hp_all(1,k);
            H=H_all(1,k);
            A=C/H;
            h=(S/A);
            K=S_beg/(alpha*C);
            
            R=[];
            for yy = 1:size(year,2)-1
                c=C/Ia(1,yy);
                c_cal(yy,1)=c;
                if c>=0.5

                    for mm=1:12
                        R(mm,yy)= K(yy,1)*ia(1,yy);
                    end

                elseif c==0

                    for mm=1:12
                        R(mm,yy)=im(mm,yy);
                    end
                else

                    for mm=1:12
                        R(mm,yy)=[(c/0.5)^2]*[K(yy,1)*ia(1,yy)] + [1-(c/0.5)^2]*im(mm,yy);
                    end
                end
            end
            Hp=[];
            for i=1:size(year,2)-1
                for mm=1:12
                    Hp(mm,i)=[R(mm,i)*h(mm,i)*9.81]/1000;
                end
            end

            for l=1:size(Hp,2)
                for i=1:12
                    if Hp(i,l)>max_Hp
                        Hp(i,l)=max_Hp;
                    else
                        Hp(i,l)=Hp(i,l);
                    end
                end
            end
            Hp=mean(Hp,2);
            
            Hydropower(:,k)=Hp;
        end
        Thp=vertcat(Thp,Hydropower);monthly_all2=vertcat(monthly_all2,monthly_all);
        
        Hp_J_M=vertcat(Hydropower(9:12,:),Hydropower(1:8,:));
        %     figure
        %     subplot(3,3,1)
        %     plot(Hp_J_M(:,1));
        %     xlim([1 12]);
        %     ylim([0 150])
        %     title('Mahi')
        %     xlabel('Months')
        %     ylabel('Hydropower MW')
        %
        %     subplot(3,3,2)
        %     plot(Hp_J_M(:,2));
        %     xlim([1 12]);
        %     ylim([0 300])
        %     title('Kadana')
        %     xlabel('Months')
        %     ylabel('Hydropower MW')
    end
    Thp_all=reshape(Thp,12,5*k);monthly_all3=reshape(monthly_all2,12,5*k);
    
    
    mean_hp=[];
    for i=1:2
        for mm=1:12
            mean_hp(mm,i)=[sum(Thp(mm,i)+Thp(12+mm,i)+Thp(24+mm,i)+Thp(36+mm,i)+Thp(48+mm,i))]/5;
        end
    end
    for i=1:2
        for mm=1:12
            mean_sd(mm,(i-1)*2+1:(i-1)*2+2) = [mean(Thp_all(mm,(i-1)*5+1:(i-1)*5+5),2) std(Thp_all(mm,(i-1)*5+1:(i-1)*5+5))];
        end
    end
    mean_sd2=vertcat(mean_sd(9:12,:),mean_sd(1:8,:));
    
    for i=1:2
        for mm=1:12
            mean_sd4(mm,(i-1)*2+1:(i-1)*2+2) = [mean(monthly_all3(mm,(i-1)*5+1:(i-1)*5+5),2) std(monthly_all3(mm,(i-1)*5+1:(i-1)*5+5))];
        end
    end
    
    inflow_dam=[];
    for i=1:2
        for mm=1:12
            inflow_dam(mm,i)=[sum(monthly_all2(mm,i)+monthly_all2(12+mm,i)+monthly_all2(24+mm,i)+monthly_all2(36+mm,i)+monthly_all2(48+mm,i))]/5;
        end
    end
    
    inflow_final=mean(inflow_dam,1);inflow_dam=[[1:12]' inflow_dam];
    Hp_June_May=vertcat(mean_hp(9:12,:),mean_hp(1:8,:))
    Hp_final=mean(Hp_June_May,1); Hp_final1=Hp_final';
    Hp_all_final(:,sss)=Hp_final1;
    inflow_final_all(:,sss)=inflow_final';
    
 
    
    mean_sd3=[[1:12]' mean_sd2];mean_sd5=[[1:12]' mean_sd4];
    for i=1:2
        sd_minus(:,i)=mean_sd2(:,(i-1)*2+1)-mean_sd2(:,(i-1)*2+2);
    end
    for i=1:2
        sd_plus(:,i)=mean_sd2(:,(i-1)*2+1)+mean_sd2(:,(i-1)*2+2);
    end
    
    for j=1:2
        for i=1:12
            if sd_plus(i,j)>max_Hp_all(1,j)
                sd_plus(i,j)=max_Hp_all(1,j);
            end
        end
    end
    sd_minus=[[1:12]' sd_minus];
    sd_plus=[[1:12]' sd_plus];  
  
end
for j=1:size(Hp_all_final,1)
    for i=1:6
        Hp_final_all_percentage(j,i)=[[Hp_all_final(j,i+1)-Hp_all_final(j,1)]/Hp_all_final(j,1)]*100;
        inflow_final_all_percentage(j,i)=[[inflow_final_all(j,i+1)-inflow_final_all(j,1)]/inflow_final_all(j,1)]*100;
    end
end
% dlmwrite(['/media/wcl/Elements/data/Results/mahi/mahi_hp_all'],Hp_all_final,' ');
% dlmwrite(['/media/wcl/Elements/data/Results/mahi/mahi_daminflow_all'],inflow_final_all,' ');
% dlmwrite(['/media/wcl/Elements/data/Results/mahi/mahi_hp_all_percentage'],Hp_final_all_percentage,' ');
% dlmwrite(['/media/wcl/Elements/data/Results/mahi/mahi_daminflow_all_percentage'],inflow_final_all_percentage,' ');
