%%% For 3 degree warming
%%% Written by Dipes Chuphal

clear all;
clc;

dam_data = importdata('/media/iitgn/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/dam_params/dam_params/damparam_mahi_06min_2020.csv');
dam_lalo_sto = dam_data.data;
dam_name = dam_data.textdata;
sce={'historical';'historical';'historical';'historical';'historical';'ssp585';'ssp585';'ssp585';'ssp585';'ssp585'};
warming={'1';'1_1.5';'1.5_2';'2_2.5';'2.5_3'};

Thp=[];monthly_all2=[];
for mmm=1:10
if mmm<6
        year=1995:2014
    elseif mmm==6
        year=2069:2078;
    elseif mmm==7
        year = 2085:2094;
    elseif mmm==8
        year = 2079:2088;
    elseif mmm==9
        year = 2076:2085;
    elseif mmm==10
        year = 2050:2059;
    end
    f2={'BCC-CSM2-MR','INM-CM5-0','MIROC6','NorESM2-MM','TaiESM1','BCC-CSM2-MR','INM-CM5-0','MIROC6','NorESM2-MM','TaiESM1'};    
    f1='mahi_';
    folder=horzcat(f1,f2{mmm});
    for k=[2]
        final1=[];
        for j = 1:size(year,2)
            disp(year(1,j))
            final=[];
            lat = ncread(['/media/iitgn/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{mmm},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'lat');
            lon = ncread(['/media/iitgn/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{mmm},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'lon');
            lat=double(lat); lat=round(lat,3);
            lon=double(lon); lon=round(lon,3);
            sto = ncread(['/media/iitgn/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{mmm},'/',folder,'/o_damsto',num2str(year(1,j)),'.nc'],'damsto');
            
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
            inflow = ncread(['/media/iitgn/Elements/cmf_v401_pkg_20210331/etc/reservoir_operation_mahi/out/',sce{mmm},'/',folder,'/o_daminf',num2str(year(1,j)),'.nc'],'daminf');
            
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
        %% dam sto begg
        %storage at beggining of year
        S_beg=[];
        fd=[];
        for yy=1:size(year,2)
            fd=find(dam_sto(:,1)==year(1,yy) & dam_sto(:,2)==10 &dam_sto(:,3)==1);
            S_beg(yy,1)= (dam_sto(fd,4));
        end
        S_beg=S_beg(1:end-1,1);
        %%
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
        
        %% dam parameter for all
        max_Hp_all=[140 240];
        alpha=0.85;
        C_all=dam_lalo_sto(:,8);
        H_all = [74.5 66];
        %%
        C=C_all(k,1)*10^6-(dam_lalo_sto(k,7)*10^6/2);
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
                %     for i=1:size(year,2)
                for mm=1:12
                    R(mm,yy)= K(yy,1)*ia(1,yy);
                end
                %     end
            elseif c==0
                %      for i=1:size(year,2)
                for mm=1:12
                    R(mm,yy)=im(mm,yy);
                end
            else
                % for i=1:size(year,2)
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
    
    
end
%
Thp_all=reshape(Thp,12,10*k);monthly_all3=reshape(monthly_all2,12,10*k);
Thp_all=mean(Thp_all,1);monthly_all3=mean(monthly_all3,1);
Thp_all=reshape(Thp_all,10,k);monthly_all3=reshape(monthly_all3,10,k);

% for j=1:size(Thp_all,2)
%     for k=1:2
%     for i=1:5
%         Hp_final_all_percentage(i,j)=[[Thp_all(5*k+1,j)-Thp_all(i,j)]/Thp_all(i,j)]*100;
%         inflow_final_all_percentage(i,j)=[[monthly_all3(5*k+1,j)-monthly_all3(i,j)]/monthly_all3(i,j)]*100;
%     end
%     end
% end
dlmwrite(['/media/iitgn/Elements/data/Results/dams/warming/hp/3degree/mahi_3degree'],Thp_all,' ');
dlmwrite(['/media/iitgn/Elements/data/Results/dams/warming/inflow/3degree/mahi_3degree'],monthly_all3,' ');
