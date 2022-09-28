%setting up environment
clear; clc; close all
addpath(genpath(('C:\Users\Paolo Tasseron\OneDrive - Wageningen University & Research')));

%loading data and storing it to variables
data = readtable('C:\Users\Paolo Tasseron\OneDrive - Wageningen University & Research\Paper projects\Amsterdam Paper\Master_Datasheet_Amsterdam.csv');
load('C:\Users\Paolo Tasseron\OneDrive - Wageningen University & Research\Paper projects\Amsterdam Paper\redtoblue.mat');

%% Computing data

data_means = zeros(1,7);
data_percentages = zeros(7,6);
for i = 1:7
    data_mean = data(data.Day == i, :);
    tot_items = sum(data_mean.TotalCounted);
    data_percentages(i,1:7) = [mean(data_mean.PETf); mean(data_mean.PSf); mean(data_mean.EPSf); mean(data_mean.POHardf);
                               mean(data_mean.POSoftf);mean(data_mean.MLf);mean(data_mean.Otherf)];
    data_means(i) = mean(data_mean.ItemFlux_n_hour_);
end
% ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bridgeFlux = zeros(1,6);
% bridgeFlux_plastic = zeros(1,6);
% data_brige = [];
% data_bridge_month = zeros(6,12);
% %%id_1 = 1, id_2 = 3, id_3 = 2, id_4 = 1, id_5 = 1, id_6 = 3
% 
% %%Correct method (new)
% for j = [1 4 5]
%     data_bridge = data(data.Bridge_ID == j, :);
%     bridgeFlux(j) = sum(data_bridge.ItemFlux_n_hour_)./size(data_bridge,1);
%     data_bridge_plastic = data(data.Bridge_ID == j, :);
%     bridgeFlux_plastic(j) = sum(data_bridge_plastic.itemflux_plastic)./size(data_bridge_plastic,1);
%     for k = 1:13
%         data_bridge_month(j,k) = mean(table2array(data_bridge(data_bridge.month == k,"ItemFlux_n_hour_")));
%     end
% end
% % ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bridgeFlux = zeros(7,6);
bridgeFlux_plastic = zeros(1,6);
data_brige = [];
data_bridge_month = zeros(6,12);

for j = [1 4 5]
    data_bridge = data(data.Bridge_ID == j, :);
    bridgeFlux(1:7,j) = [sum(data_bridge.PETf)./size(data_bridge,1);
                         sum(data_bridge.PSf)./size(data_bridge,1);
                         sum(data_bridge.EPSf)./size(data_bridge,1);
                         sum(data_bridge.POHardf)./size(data_bridge,1);
                         sum(data_bridge.POSoftf)./size(data_bridge,1);
                         sum(data_bridge.MLf)./size(data_bridge,1);
                         sum(data_bridge.Otherf)./size(data_bridge,1)];
    
%     sum(data_bridge.ItemFlux_n_hour_)./size(data_bridge,1);
    data_bridge_plastic = data(data.Bridge_ID == j, :);
    bridgeFlux_plastic(j) = sum(data_bridge_plastic.itemflux_plastic)./size(data_bridge_plastic,1);
    
    for k = 1:13
        data_bridge_month(1:7,j,k) = [mean(table2array(data_bridge(data_bridge.month == k,"PETf")));
                                      mean(table2array(data_bridge(data_bridge.month == k,"PSf")));
                                      mean(table2array(data_bridge(data_bridge.month == k,"EPSf")));
                                      mean(table2array(data_bridge(data_bridge.month == k,"POHardf")));
                                      mean(table2array(data_bridge(data_bridge.month == k,"POSoftf")));
                                      mean(table2array(data_bridge(data_bridge.month == k,"MLf")));
                                      mean(table2array(data_bridge(data_bridge.month == k,"Otherf")))];
    end
end

for j = [2 6]
    data_bridge = data(data.Bridge_ID == j, :);
    bridgeFlux(1:7,j) = [sum(data_bridge.PETf)./(size(data_bridge,1)/3);
                         sum(data_bridge.PSf)./(size(data_bridge,1)/3);
                         sum(data_bridge.EPSf)./(size(data_bridge,1)/3);
                         sum(data_bridge.POHardf)./(size(data_bridge,1)/3);
                         sum(data_bridge.POSoftf)./(size(data_bridge,1)/3);
                         sum(data_bridge.MLf)./(size(data_bridge,1)/3);
                         sum(data_bridge.Otherf)./(size(data_bridge,1)/3)];
    data_bridge_plastic = data(data.Bridge_ID == j, :);
    bridgeFlux_plastic(j) = sum(data_bridge_plastic.itemflux_plastic)./(size(data_bridge_plastic,1)/3);
    for k = 1:13
        selectedMonth = data_bridge(data_bridge.month == k,:);
        data_bridge_month(1:7,j,k) = [mean(sum(reshape(selectedMonth.PETf,3,size(selectedMonth,1)/3)));
                                     mean(sum(reshape(selectedMonth.PSf,3,size(selectedMonth,1)/3)));
                                     mean(sum(reshape(selectedMonth.EPSf,3,size(selectedMonth,1)/3)));
                                     mean(sum(reshape(selectedMonth.POHardf,3,size(selectedMonth,1)/3)));
                                     mean(sum(reshape(selectedMonth.POSoftf,3,size(selectedMonth,1)/3)));
                                     mean(sum(reshape(selectedMonth.MLf,3,size(selectedMonth,1)/3)));
                                     mean(sum(reshape(selectedMonth.Otherf,3,size(selectedMonth,1)/3)));];
    end
end

for j = [3]
    data_bridge = data(data.Bridge_ID == j, :);
    bridgeFlux(1:7,j) = [sum(data_bridge.PETf)./(size(data_bridge,1)/2);
                         sum(data_bridge.PSf)./(size(data_bridge,1)/2);
                         sum(data_bridge.EPSf)./(size(data_bridge,1)/2);
                         sum(data_bridge.POHardf)./(size(data_bridge,1)/2);
                         sum(data_bridge.POSoftf)./(size(data_bridge,1)/2);
                         sum(data_bridge.MLf)./(size(data_bridge,1)/2);
                         sum(data_bridge.Otherf)./(size(data_bridge,1)/2)];

    data_bridge_plastic = data(data.Bridge_ID == j, :);
    bridgeFlux_plastic(j) = sum(data_bridge_plastic.itemflux_plastic)./(size(data_bridge_plastic,1)/2);
    item_flux_corrected = reshape(data_bridge.ItemFlux_n_hour_,2,size(data_bridge.ItemFlux_n_hour_,1)/2);
    for k = 1:13
        selectedMonth = data_bridge(data_bridge.month == k,:);
        data_bridge_month(1:7,j,k) = [mean(sum(reshape(selectedMonth.PETf,2,size(selectedMonth,1)/2)));
                             mean(sum(reshape(selectedMonth.PSf,2,size(selectedMonth,1)/2)));
                             mean(sum(reshape(selectedMonth.EPSf,2,size(selectedMonth,1)/2)));
                             mean(sum(reshape(selectedMonth.POHardf,2,size(selectedMonth,1)/2)));
                             mean(sum(reshape(selectedMonth.POSoftf,2,size(selectedMonth,1)/2)));
                             mean(sum(reshape(selectedMonth.MLf,2,size(selectedMonth,1)/2)));
                             mean(sum(reshape(selectedMonth.Otherf,2,size(selectedMonth,1)/2)))];
    end
end

% ORIGINAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for j = [2 6]
%     data_bridge = data(data.Bridge_ID == j, :);
%     bridgeFlux(j) = sum(data_bridge.ItemFlux_n_hour_)./(size(data_bridge,1)/3);
%     data_bridge_plastic = data(data.Bridge_ID == j, :);
%     bridgeFlux_plastic(j) = sum(data_bridge_plastic.itemflux_plastic)./(size(data_bridge_plastic,1)/3);
%     for k = 1:13
%         selectedMonth = data_bridge(data_bridge.month == k,:);
%         data_bridge_month(j,k) = mean(sum(reshape(selectedMonth.ItemFlux_n_hour_,3,size(selectedMonth,1)/3)));
%     end
% end
% 
% for j = [3]
%     data_bridge = data(data.Bridge_ID == j, :);
%     bridgeFlux(j) = sum(data_bridge.ItemFlux_n_hour_)./(size(data_bridge,1)/2);
%     data_bridge_plastic = data(data.Bridge_ID == j, :);
%     bridgeFlux_plastic(j) = sum(data_bridge_plastic.itemflux_plastic)./(size(data_bridge_plastic,1)/2);
%     item_flux_corrected = reshape(data_bridge.ItemFlux_n_hour_,2,size(data_bridge.ItemFlux_n_hour_,1)/2);
%     for k = 1:13
%         selectedMonth = data_bridge(data_bridge.month == k,:);
%         data_bridge_month(j,k) = mean(sum(reshape(selectedMonth.ItemFlux_n_hour_,2,size(selectedMonth,1)/2)));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INCORRECT METHOD DO NOT USE
% for j = 1:6
%     data_bridge = data(data.Bridge_ID == j, :);
%     bridgeFlux(j) = sum(data_bridge.TotalCounted)./ sum(data_bridge.Duration_min_) * 60;
% end

timeFlux = zeros(1,12);

for k = 7:18
    data_time = data(data.Timestamp == k, :);
    timeFlux(k,1:7) = [mean(data_time.PETf);
                       mean(data_time.PSf);
                       mean(data_time.EPSf);
                       mean(data_time.POHardf);
                       mean(data_time.POSoftf);
                       mean(data_time.MLf);
                       mean(data_time.Otherf)];
end

%%
close all
fg = figure('DefaultAxesFontSize',13); set(gcf,'Position',[0 0 1500 800]);  

subplot(3,2,1);
ba = bar(data_percentages,'stacked','FaceColor','flat'); xticks([1:7]);
ba(1).CData = [140/255 81/255 10/255]; ba(2).CData = [216/255 179/255 101/255]; ba(3).CData = [246/255 232/255 195/255]; 
ba(4).CData = [245/255 245/255 245/255]; ba(5).CData = [199/255 234/255 229/255]; ba(6).CData = [90/255 180/255 172/255];
ba(7).CData = [1/255 102/255 94/255];
labels = {'SUN n=11','MON n=22','TUE n=77','WED n=142','THU n=72','FRI n=22','SAT n=0'}; 
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
xticklabels(labels);
ylim([0 80]); ylabel('Litter flux (items/hr)'); grid on; %title('Mean item fluxes per weekday'); 

subplot(3,2,2)
ba = bar(bridgeFlux','stacked','FaceColor','flat'); xticks([1:6]);
ba(1).CData = [140/255 81/255 10/255]; ba(2).CData = [216/255 179/255 101/255]; ba(3).CData = [246/255 232/255 195/255]; 
ba(4).CData = [245/255 245/255 245/255]; ba(5).CData = [199/255 234/255 229/255]; ba(6).CData = [90/255 180/255 172/255];
ba(7).CData = [1/255 102/255 94/255];
labels = {'WK n=32','WD n=32','WT n=32','GK n=32','OK n=32','PH n=31'};
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
legend('PET','PS','EPS','PO Hard','PO Soft','ML','Other','location','eastoutside')
xticklabels(labels);
ylim([0 80]); ylabel('Litter flux (items/hr)'); grid on;
%title('Mean item fluxes per measurement location'); grid on
%xticklabels({'Westerkanaal','Westerdok (GBB)','Westertoegang','Geldersekade','Oudezijdskolkbrug','Piet Heinkade'})

subplot(3,1,2)
ba = bar(timeFlux,'stacked','FaceColor','flat'); xticks([7:18]); 
ba(1).CData = [140/255 81/255 10/255]; ba(2).CData = [216/255 179/255 101/255]; ba(3).CData = [246/255 232/255 195/255]; 
ba(4).CData = [245/255 245/255 245/255]; ba(5).CData = [199/255 234/255 229/255]; ba(6).CData = [90/255 180/255 172/255];
ba(7).CData = [1/255 102/255 94/255];
labels = {'7:00 n=2','8:00 n=14','9:00 n=30','10:00 n=33','11:00 n=30','12:00 n=28','13:00 n=16','14:00 n=7','15:00 n=11','16:00 n=19','17:00 n=22','18:00 n=13'};
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
xticklabels(labels); grid on
ylim([0 60]); xlim([6 19]); ylabel('Litter flux (items/hr)'); %title('Mean item flux per time of day'); 


subplot(3,1,3)
ba = bar(squeeze(mean(data_bridge_month))','stacked','FaceColor','flat'); xticks([1:13]);
ba(1).CData = [140/255 81/255 10/255]; ba(2).CData = [216/255 179/255 101/255]; ba(3).CData = [246/255 232/255 195/255]; 
ba(4).CData = [245/255 245/255 245/255]; ba(5).CData = [199/255 234/255 229/255]; ba(6).CData = [90/255 180/255 172/255];
ba(7).CData = [1/255 102/255 94/255];
labels = {'Feb-21 n=11','Mar-21 n=55','Apr-21 n=44','May-21 n=11','Jun-21 n=22','Jul-21 n=32','Aug-21 n=11','Sep-21 n=54','Okt-21 n=22','Nov-21 n=33','Dec-21 n=17','Jan-22 n=22','Feb-22 n=11'};
labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);
xticklabels(labels); grid on;
ylim([0 120]); ylabel('Litter flux (items/hr)'); %title('Mean item flux per month'); 
AddLetters2Plots(fg,{'(a)','(b)','(c)','(d)'},'HShift',0,'VShift',0,'Direction','LeftRight','FontSize',15)

%% Correlation analyses
item_flux = data.ItemFlux_n_min_(1:346);

%sunlight
sunlight = data.Sunlight_hrs_(1:346);
corr(item_flux,sunlight);

%rainfall
rainfall_mm = data.Rainfall_mm_(1:346);
corr(item_flux,rainfall_mm);

rainfall_dur = data.RainfallLength_hr_(1:346);
corr(item_flux,rainfall_dur);

%wind
windspd_avg = data.WindSpeedAverage_m_s_(1:346);
winddir = data.WindDirection___(1:346);
windspd_max_hr = data.WindspeedMaxPerHour_MaxWindGusts_m_s_(1:346);

%tide
tide = data.tide_bool(122:346);
% x = [sunlight,rainfall_mm,windspd_avg,winddir];
% mdl = fitlm(x,item_flux);

%%
wind_data_bridge = data(data.Bridge_ID == 6, :); 
%wind_data_bridge = data(data.Bridge_ID == 4 | data.Bridge_ID == 5, :);

Q = wind_data_bridge.ItemFlux_n_hour_;
W = wind_data_bridge.WindDirection;

wind_directions = {'ENE','ESE','N','NE','NNE','NNW','NW','S','SE','SSE','SSW','SW','W','WNW','NW'};
wind_values = [67.5,112.5,0,45,22.5,33.5,315,180,135,157.5,202.5,225,270,202.5,225,270,292.5,247.5];
wind = [];

for i = 1:length(wind_directions)
    wind_direction = wind_directions{i};
    if isnan(mean(Q(strcmp(W,wind_direction))))
        new_input = [];
    else
        new_input = ones(floor(mean(Q(strcmp(W,wind_direction)))),1)*wind_values(i);
    end
    wind = [wind;new_input];
end

figure; set(gcf,'Position',[100 100 400 400]);  ax = polaraxes; handle = polarhistogram(ax, deg2rad(wind),16,'FaceColor','K'); hold on; 
handle.BinEdges = handle.BinEdges + handle.BinWidth/2;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise'; % 90 degrees at the right]
pax = gca; pax.FontSize = 13; pax.ThetaTick = [0,22.5,45,67.5,90,112.5,135,157.5,180,202.5,225,247.5,270,292.5,315,337.5];
labels = {'N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW'};
pax.ThetaTickLabel = labels;
%pax.RTick = [0:25:150];
%pax.RTickLabel = {'','25','50','75','100','125','150'}; 
legend('Outflow (items/hr)','Measurements in wind direction (n)','Location','Southoutside')

% labels = {'N n=15','NNE n=14','NE n=8','ENE n=11','E n=1','ESE n=11','SE n=11','SSE n=43','S n=23','SSW n=2','SW n=33','WSW n=60','W n=51','WNW n=32','NW n=21','NNW n=10'}
% labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);

%% ebb & flood
ebb_outflow = data(data.tide_bool == 1,:); ebb_outflow = ebb_outflow.ItemFlux_n_hour_;
flood_inflow = data(data.tide_bool == 0,:); flood_inflow = flood_inflow.ItemFlux_n_hour_;
flood_inflow = [flood_inflow;NaN(99,1)];
group2 = [
        ones(size(ebb_outflow));
        2*ones(size(ebb_outflow));
];

fg = figure('DefaultAxesFontSize',13); boxplot([ebb_outflow,flood_inflow],group2); ax = gca;
set(ax,'Xtick',[1 2]);
set(ax,'XTickLabel',{'Ebb: Outflow (n = 162)','Flood: Inflow (n = 63)'});
ylabel('Outflow (Items/hour)'); ylim([-10 100]); grid on

 colors = ones(2,3)*0.8;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end


%% Correlation matrices {this NEEDS TO BE REDONE}
pet = data.PETf(1:346); ps = data.PSf(1:346); eps = data.EPSf(1:346); pohard = data.POHardf(1:346); 
posoft = data.POSoftf(1:346); ml = data.MLf(1:346); other = data.Otherf(1:346);
R_test = [pet,ps,eps,pohard,posoft,ml,other,item_flux];

%[Rho,pval] corr

R_rainfall = corr(R_test,rainfall_mm);
R_sunlight = corr(R_test,sunlight);
R_windspeed = corr(R_test,windspd_avg);
R_windspeedmax = corr(R_test,windspd_max_hr);
R_tide = corr(R_test(122:346,:),tide);
R_matrix = [R_rainfall,R_sunlight,R_windspeed,R_windspeedmax,R_tide]';

%generating labels in image
%generating labels in image
t = num2cell(round(R_matrix,3)); 
t = cellfun(@num2str, t,'UniformOutput',false);

t = {'0.048','-0.065','-0.003','0.118*','0.112*','0.105','0.108*','0.129*';
 '0.019','0.053','0.005','0.055','0.098','-0.091','-0.107*','-0.024';
 '-0.058','0.022','0.024','0.048','0.139**','0.076','0.085','0.11*';
 '-0.051','0.022','0.051','0.102*','0.158**','0.085','0.087','0.131*';
 '0.06','-0.194**','0.022','0.001','-0.011','0.057','0.066','0.041'};

y = repmat(1:8,5,1); x = repmat(1:5,8,1)';
figure('DefaultAxesFontSize',13); set(gcf,'Position',[100 100 950 400]);
imagesc(R_matrix); text(y(:), x(:), t, 'HorizontalAlignment', 'Center');
set(gca, 'XTick', 1:8); % center x-axis ticks on bins
set(gca, 'YTick', 1:5); % center y-axis ticks on bins
set(gca, 'XTickLabel', {'PET','PS','EPS','PO Hard','PO Soft','Multilayer','Other','Combined'}); % set x-axis labels
set(gca, 'YTickLabel', {'Rainfall (mm)','Sunlight (hrs)','Windspeed avg (m/s)','Windspeed max (m/s)','Tide'}); % set y-axis labels
colormap(RedToBlue); a = colorbar; caxis([-.2 .2]); ylabel(a, 'r^2')
%% Plotting boxplots
% id_1 = 1, id_2 = 3, id_3 = 2, id_4 = 1, id_5 = 1, id_6 = 3


data_bridge1 = data(data.Bridge_ID == 1,:); A = data_bridge1.ItemFlux_n_hour_;
data_bridge2 = data(data.Bridge_ID == 2,:); data_bridge2b = reshape(data_bridge2.ItemFlux_n_hour_,3,32); B = sum(data_bridge2b)';
data_bridge3 = data(data.Bridge_ID == 3,:); data_bridge3b = reshape(data_bridge3.ItemFlux_n_hour_,2,32); C = sum(data_bridge3b)';
data_bridge4 = data(data.Bridge_ID == 4,:); D = data_bridge4.ItemFlux_n_hour_;
data_bridge5 = data(data.Bridge_ID == 5,:); E = data_bridge5.ItemFlux_n_hour_;
data_bridge6 = data(data.Bridge_ID == 6,:); data_bridge6b = reshape(data_bridge6.ItemFlux_n_hour_,3,31); F = sum(data_bridge6b)';

D = [D; NaN; NaN];E = [E; NaN]; F = [F; NaN];
group = [   ones(size(A))'; 
            2*ones(size(A))';
            3*ones(size(A))';
            4*ones(size(A))';
            5*ones(size(A))';
            6*ones(size(A))'];
        
 figure('DefaultAxesFontSize',13);  boxplot([A,B,C,D,E,F],group); grid on; ax = gca;
 set(ax,'Xtick',[1 2 3 4 5 6]);
 set(ax,'XTickLabel',{'WK','WD','WT','GK','OK','PH'});
 ylabel('Outflow (Items/hour)');

 colors = ones(6,3)*0.8;
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.5);
end        

%%
