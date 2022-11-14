%%%far period

clc;clear all;
region={'north';'central';'south'};
for r=1:3
    dam_mean_ssp126=[]; dam_mean_ssp585=[];
    dam=importdata(['/media/wcl/Elements/cama_reservoir_data/dams_lalo2/',region{r},'/damname']);
    dam_name=dam;
    model={'BCC-CSM2-MR','INM-CM5-0','MIROC6','NorESM2-MM','TaiESM1'};
    for k=11%:length(dam_name)
        model_mean_ssp126=[]; model_mean_ssp585=[]; historical_temp=[];  ssp126_temp=[]; ssp585_temp=[];
        lalo=importdata(['/media/wcl/Elements/cama_reservoir_data/dams_lalo2/',region{r},'/',dam{k}]);
        disp(k)
        for mm=1:5
            temp_change_ssp126=[]; temp_change_ssp585=[]; model_ssp126=[]; model_ssp585=[]; model_historic=[];
            for i=1:size(lalo,1)
                disp(i)
                forcing_hist=dlmread(['/media/wcl/Pennar_IITGN/Dipesh_BC/VIC_forcing_BC/',model{mm},'/historical/data_',num2str(lalo(i,4)),'_',num2str(lalo(i,3))]);
                t = datetime(1951,01,01):datetime(2014,12,31);
                t=t';
                [y, m, d] =ymd(t);
                time = [y m d];
                hist=[time forcing_hist(:,2:3)];
                temp_hist=hist((16072:end),:);
                temp_hist=[temp_hist(:,1:3) mean(temp_hist(:,4:5),2)];
                
                yearly_temp_hist=mean(temp_hist(:,4),1);
                
                foring_end_ssp126=dlmread(['/media/wcl/Pennar_IITGN/Dipesh_BC/VIC_forcing_BC/',model{mm},'/ssp126/data_',num2str(lalo(i,4)),'_',num2str(lalo(i,3))]);
                foring_end_ssp585=dlmread(['/media/wcl/Pennar_IITGN/Dipesh_BC/VIC_forcing_BC/',model{mm},'/ssp585/data_',num2str(lalo(i,4)),'_',num2str(lalo(i,3))]);
                t = datetime(2015,01,01):datetime(2100,12,31);
                t=t';
                [y, m, d] =ymd(t);
                time = [y m d];
                ssp126=[time foring_end_ssp126(:,2:3)]; ssp585=[time foring_end_ssp585(:,2:3)];
                temp_ssp126=ssp126((24108:end ),:);   temp_ssp585=ssp585((24108:end ),:);
                temp_ssp126=[temp_ssp126(:,1:3) mean(temp_ssp126(:,4:5),2)]; temp_ssp585=[temp_ssp585(:,1:3) mean(temp_ssp585(:,4:5),2)];
                yearly_temp_ssp126=mean(temp_ssp126(:,4),1);    yearly_temp_ssp585=mean(temp_ssp585(:,4),1);
                
                %%%%change in temp
                temp_change_ssp126(i,1)=yearly_temp_ssp126-yearly_temp_hist;
                temp_change_ssp585(i,1)=yearly_temp_ssp585-yearly_temp_hist;
                model_historic(i,1) =   yearly_temp_hist;
            end
            historical_temp(mm,1)= mean(model_historic,1);
            
            model_mean_ssp126(mm,1)=mean(temp_change_ssp126,1);
            model_mean_ssp585(mm,1)=mean(temp_change_ssp585,1);
        end
        
        dam_mean_ssp126(k,1)=(mean(model_mean_ssp126,1));
        dam_mean_ssp585(k,1)=(mean(model_mean_ssp585,1));
        
        dam_mean_ssp126(k,2)=(std(model_mean_ssp126));
        dam_mean_ssp585(k,2)=std((model_mean_ssp585));
    end
    
    % dlmwrite(['/media/wcl/Elements/data/Results/dams/temp_dams/',region{r},'/far_period_ssp126_',char(region{r})],dam_mean_ssp126,' ');
    % dlmwrite(['/media/wcl/Elements/data/Results/dams/temp_dams/',region{r},'/far_period_ssp585_',char(region{r})],dam_mean_ssp585,' ');
end
