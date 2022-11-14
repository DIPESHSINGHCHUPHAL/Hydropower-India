%%%far period

clc;clear all;
region={'north';'central';'south'};
for r=[1 2 3]
dam_mean_ssp126=[]; dam_mean_ssp585=[];
dam=importdata(['/media/wcl/Elements/cama_reservoir_data/dams_lalo2/',region{r},'/damname']);
dam_name=dam;
model={'BCC-CSM2-MR','INM-CM5-0','MIROC6','NorESM2-MM','TaiESM1'};
for k=1:length(dam_name)
    lalo=importdata(['/media/wcl/Elements/cama_reservoir_data/dams_lalo2/',region{r},'/',dam{k}]);
    model_mean_ssp126=[]; model_mean_ssp585=[]; historical_pcp=[];  ssp126_pcp=[]; ssp585_pcp=[];
    disp(k)
    for mm=1:5
            pcp_change_ssp126=[]; pcp_change_ssp585=[]; model_ssp126=[]; model_ssp585=[]; model_historic=[];
       parfor i=1:size(lalo,1)
           disp(i)
            forcing_hist=dlmread(['/media/wcl/Pennar_IITGN/Dipesh_BC/VIC_forcing_BC/',model{mm},'/historical/data_',num2str(lalo(i,4)),'_',num2str(lalo(i,3))]);
            t = datetime(1951,01,01):datetime(2014,12,31);
            t=t';
            [y, m, d] =ymd(t);
            time = [y m d];
            hist=[time forcing_hist(:,1)];
            pcp_hist=hist((16072:end),:);
            
            
            yearly_pcp_hist=(sum(pcp_hist(:,4))/20);
            
            foring_end_ssp126=dlmread(['/media/wcl/Pennar_IITGN/Dipesh_BC/VIC_forcing_BC/',model{mm},'/ssp126/data_',num2str(lalo(i,4)),'_',num2str(lalo(i,3))]);
            foring_end_ssp585=dlmread(['/media/wcl/Pennar_IITGN/Dipesh_BC/VIC_forcing_BC/',model{mm},'/ssp585/data_',num2str(lalo(i,4)),'_',num2str(lalo(i,3))]);
            t = datetime(2015,01,01):datetime(2100,12,31);
            t=t';
            [y, m, d] =ymd(t);
            time = [y m d];
            ssp126=[time foring_end_ssp126(:,1)]; ssp585=[time foring_end_ssp585(:,1)];
            pcp_ssp126=ssp126((24108:end),:);   pcp_ssp585=ssp585((24108:end),:);
            
            yearly_pcp_ssp126=(sum(pcp_ssp126(:,4))/20);     yearly_pcp_ssp585=(sum(pcp_ssp585(:,4))/20);
            %%%%change in pcp
            pcp_change_ssp126(i,1)=yearly_pcp_ssp126-yearly_pcp_hist;
            pcp_change_ssp585(i,1)=yearly_pcp_ssp585-yearly_pcp_hist;
            model_ssp126(i,1)=yearly_pcp_ssp126;       model_ssp585(i,1) = yearly_pcp_ssp585
            model_historic(i,1) =   yearly_pcp_hist;
        end
        historical_pcp(mm,1)= mean(model_historic,1);
        ssp126_pcp(mm,1)= mean(model_ssp126,1);
        ssp585_pcp(mm,1)= mean(model_ssp585,1);

        model_mean_ssp126(mm,1)=mean(pcp_change_ssp126,1);
        model_mean_ssp585(mm,1)=mean(pcp_change_ssp585,1);
    end
    
    dam_mean_ssp126(k,1)=(mean(model_mean_ssp126,1)/mean(historical_pcp,1))*100;
    dam_mean_ssp585(k,1)=(mean(model_mean_ssp585,1)/mean(historical_pcp,1))*100;
    
    dam_mean_ssp126(k,2)=(std(ssp126_pcp)/mean(ssp126_pcp,1))*100;
    dam_mean_ssp585(k,2)=std((ssp585_pcp)/mean(ssp585_pcp,1))*100;
end

dlmwrite(['/media/wcl/Elements/data/Results/dams/precipitation_dams/',region{r},'/far_period_ssp126_',char(region{r})],dam_mean_ssp126,' ');
dlmwrite(['/media/wcl/Elements/data/Results/dams/precipitation_dams/',region{r},'/far_period_ssp585_',char(region{r})],dam_mean_ssp585,' ');
end
