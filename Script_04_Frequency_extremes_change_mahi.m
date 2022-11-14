%% Change in frequency of days when reservoir storage and inflow simultaneously exceeds 95th percentile threshold

clear all;
clc;
dam_data = importdata('/media/rahul/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/dam_params/dam_params/damparam_mahi_06min_2020.csv');
dam_lalo_sto = dam_data.data;
dam_name = dam_data.textdata;
per={'hist';'near_ssp126';'near_ssp585';'mid_ssp126';'mid_ssp585';'end_ssp126';'end_ssp585'};
sce={'historical';'ssp126';'ssp585';'ssp126';'ssp585';'ssp126';'ssp585'};
per=95;
for sss=1:7
    disp(sss)
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
    hist2=[];
  for mm=1:5
        f2={'INM-CM5-0','MIROC6','NorESM2-MM','TaiESM1','BCC-CSM2-MR','hist'};
        f1='mahi_';
        folder=horzcat(f1,f2{mm});
        hist1=[];
        
        for k=[2]
            final1=[];
            for j = 1:size(year,2)
%                 disp(year(1,j))
                final=[];
                lat = ncread(['/media/rahul/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'lat');
                lon = ncread(['/media/rahul/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'lon');
                lat=double(lat); lat=round(lat,3);
                lon=double(lon); lon=round(lon,3);
                sto = ncread(['/media/rahul/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'damsto');
                
                t = datetime(year(1,j),01,01):datetime(year(1,j),12,31);
                t=t';
                [y, m, d] =ymd(t);
                time = [y m d];
                
                out = sto(dam_lalo_sto(k,4),dam_lalo_sto(k,5),:);
                final=[time squeeze(out)];
                final1=vertcat(final1,final);
            end
            dam_sto=final1;  %dam_sto total
            dam_sto(:,4)=final1(:,4)-(dam_lalo_sto(k,7)*10^6/2);  %dam_sto total
            
            %% inflow
            data=[];
            t = datetime(year(1,1),01,01):datetime(year(1,end),12,31);
            t=t';
            [y, m, d] =ymd(t);
            time = [y m d];
            
            for j = 1:size(year,2)
%                 disp(year(1,j))
                inflow = ncread(['/media/rahul/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{sss},'/',folder,'/o_daminf',num2str(year(1,j)),'.nc'],'daminf');
                
                out = squeeze(inflow(dam_lalo_sto(k,4),dam_lalo_sto(k,5),:));
                data=vertcat(data,out);
            end
            inflw_data=[time data];
            
             %% dam sto monthly mean
            S=[];
            fd=[];
            for i=1:size(year,2)
                fd=find(dam_sto(:,1)==year(1,i) & dam_sto(:,2)>=6 & dam_sto(:,2)<=9);
                S1=dam_sto(fd,:);
                S=vertcat(S,S1);
            end
            
            %mean_annual_inflow
            I=[];
            fd=[];
            for i=1:size(year,2)
                fd=find(inflw_data(:,1)==year(1,i) & inflw_data(:,2)>=6 & inflw_data(:,2)<=9);
                I1=inflw_data(fd,:);
                I=vertcat(I,I1);
            end
            %%
            if sss==1
                S99(mm,k)=prctile(S(:,4),per);
                I99(mm,k)=prctile(I(:,4),per);
            end
            
            fd=[];fd2=[];fd3=[];
            for i =1:length(I)
                %                 disp(i)
                
                fd22=find(I(i,4)>=I99(mm,k));
                fd33=find(S(i,4)>=S99(mm,k));
                fd1=find(S(i,4)>=S99(mm,k) & I(i,4)>=I99(mm,k));
                if length(fd1)>0
                    fd=vertcat(fd,fd1);
                end
                if length(fd22)>0
                    fd2=vertcat(fd2,fd22);
                end
                if length(fd33)>0
                    fd3=vertcat(fd3,fd33);
                end
            end
            hist1(:,k)=length(fd);
            hist12(:,k)=length(fd2);
            hist13(:,k)=length(fd3);
        end
        hist2(mm,:)=hist1;
        hist22(mm,:)=hist12;
        hist23(mm,:)=hist13;
        
    end
    hist3(sss,:)=mean(hist2,1);
    hist32(sss,:)=mean(hist22,1);
    hist33(sss,:)=mean(hist23,1);
    
end

for m=1:size(hist3,1)-1
    change_days(m,:) = [hist3(m+1,:)-hist3(1,:)];
end
for m=1:size(hist3,1)-1
    percent_change_days(m,:) =[ [hist3(m+1,:)-hist3(1,:)]./hist3(1,:)]*100;
end
for m=1:size(hist3,1)
    percent_relative_frequency(m,:) =[hist3(m,:)/2440]*100;
end
for m=1:size(hist3,1)-1
    change_percent_relative_frequency(m,:) =percent_relative_frequency(m+1,:)-percent_relative_frequency(1,:);
end
%%

for m=1:size(hist32,1)
    percent_relative_frequency_inflow(m,:) =[hist32(m,:)/2440]*100;
end
for m=1:size(hist32,1)-1
    change_percent_relative_frequency_inflow(m,:) =percent_relative_frequency_inflow(m+1,:)-percent_relative_frequency_inflow(1,:);
end
%%

for m=1:size(hist33,1)
    percent_relative_frequency_sto(m,:) =[hist33(m,:)/2440]*100;
end
for m=1:size(hist33,1)-1
    change_percent_relative_frequency_sto(m,:) =percent_relative_frequency_sto(m+1,:)-percent_relative_frequency_sto(1,:);
end
% fd2=find(S(:,4)>=S99)
% S1=S(fd2,:);
% fd3=find(I(:,4)>=I99)
% I1=I(fd3,:);
%
mkdir('/media/rahul/Elements/data/Results_final/mahi/Frequency');
dlmwrite(['/media/rahul/Elements/data/Results_final/mahi/Frequency/change_days'],change_days,' ');
dlmwrite(['/media/rahul/Elements/data/Results_final/mahi/Frequency/percent_change_days'],percent_change_days,' ');
dlmwrite(['/media/rahul/Elements/data/Results_final/mahi/Frequency/percent_relative_frequency'],percent_relative_frequency,' ');
dlmwrite(['/media/rahul/Elements/data/Results_final/mahi/Frequency/change_percent_relative_frequency'],change_percent_relative_frequency,' ');
dlmwrite(['/media/rahul/Elements/data/Results_final/mahi/Frequency/percent_relative_frequency_inflow'],change_percent_relative_frequency_inflow,' ');
dlmwrite(['/media/rahul/Elements/data/Results_final/mahi/Frequency/change_percent_relative_frequency_sto'],change_percent_relative_frequency_sto,' ');
