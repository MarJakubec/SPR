close all
clear all 

% variables
name = '13pG';
mol_weight = 948.5574;
conc = [4 8 12 16 24 32 48 64 80 96 112 128];

mol_volume = 0.6564;
lipid = 'DMPC';
lipid_weight = 655.7;

%samples to be deleted
CH1_del = [12]; %samples for chanel1;
CH2_del = [];
CH3_del = [];
%----------------------------------------------------------------

% 
% mol_volume = 0.6564*0.90+27.0990*0.10;
% lipid = 'DMPC w 10% LPS';
% lipid_weight = 655.7*0.90+15000*0.10;
CH_del{1} = CH1_del;
CH_del{2} = CH2_del;
CH_del{3} = CH3_del;
% Names of file
files = dir('*Chanel*.txt');
% Changes comma into dots don forget to delete values in the end of file -
% need to bu runned only once 
for i= 1:(size(files,1))
    Data = fileread(files(i).name);
    Data = strrep(Data, ',', '.');
    FID = fopen(files(i).name, 'w');
    fwrite(FID, Data, 'char');
    fclose(FID);
end

% Import files in results grid
for i= 1:(size(files,1))
results{i} = dlmread(files(i).name,'\t',1,0);
end

%move X and Y into separate grid and 
for i= 1:(size(files,1))
   for j = 1:(size(results{i},2))/2
       res_parsed{i}(:,j) = results{i}(:,j*2);
       if results{i}((size(results{i},1)-1),j*2-1) ~= 0
         timeX(:,i) = results{i}(:,j*2-1);
       end
   end
end

% Change 0 in results into NaN - useless 

timeX(timeX == 0) = NaN;

%Substract baseline 

for i = 1:size(res_parsed,2)-1
    res_blank{i} = res_parsed{i+1}; %- res_parsed{1};  
    
end

% RUL calculation
for i = 1:size(res_parsed,2)-1
    res_rul{i} = res_parsed{i+1};
    RUL(i) = mean(res_rul{i}(timeX(:,1)>4000 & timeX(:,1)<5200,1))- ...
         mean(res_rul{i}(timeX(:,1)>0 & timeX(:,1)<30,1));
end



%Filter of only specific time window
filt_ind{1} = [1 1000];
filt_ind{2} = [1155 2154];
filt_ind{3} = [2298 3297];
to_plotX1 = [filt_ind{1}(1):filt_ind{1}(2)]';
to_plotY1 = res_blank{1}(filt_ind{1}(1):filt_ind{1}(2),:);

to_plotX2 = [filt_ind{2}(1):filt_ind{2}(2)]';
to_plotY2 = res_blank{2}(filt_ind{2}(1):filt_ind{2}(2),:);

to_plotX3 = [filt_ind{3}(1):filt_ind{3}(2)]';
to_plotY3 = res_blank{3}(filt_ind{3}(1):filt_ind{3}(2),:);
%Choose injection time
inj1 = 61; 
inj2 = 1215;
inj3 = 2369;

to_plotX1 = to_plotX1 - inj1; 
to_plotX2 = to_plotX2 - inj2; 
to_plotX3 = to_plotX3 - inj3; 

% Move to plot  to zero

    for j = 1:size(to_plotY1,2)
        to_plotY1(:,j)= to_plotY1(:,j) - mean(to_plotY1(to_plotX1 < 0,j));
    end
    
    for j = 1:size(to_plotY2,2)
        to_plotY2(:,j)= to_plotY2(:,j) - mean(to_plotY2(to_plotX2 < 0,j));
    end

    for j = 1:size(to_plotY3,2)
        to_plotY3(:,j)= to_plotY3(:,j) - mean(to_plotY3(to_plotX3 < 0,j));
    end
to_plotX(:,1) = to_plotX1;
to_plotX(:,2) = to_plotX2;
to_plotX(:,3) = to_plotX3;

to_plotY{1} = to_plotY1;
to_plotY{2} = to_plotY2;
to_plotY{3} = to_plotY3;

for o = 1:size(to_plotX,2)
names =sprintf('%s Chanel %d %s',name,o+1,lipid); % add names of samples 

figure('units','normalized','outerposition',[0 0 1 1],Visible='off')
subplot(2,2,[1 3]) 
hold on
if isempty(CH_del{o})
    for j = [4:15]
        plot(to_plotX(:,o),to_plotY{o}(:,j),'LineWidth', 2);      
    end
else
A = CH_del{o} + 3;
j_del = setdiff([4:15],A,'stable');
    for j = j_del
        plot(to_plotX(:,o),to_plotY{o}(:,j),'LineWidth', 2);      
    end
end

  title(names(1,:))
xlabel('time[s]') % x-axis label
ylabel('Relative response unit') % y-axis label
read_point(1) = 190;
if isempty(CH_del{o})
Ext_Values = to_plotY{o}(to_plotX(:,o)==read_point,4:15);
Scatt_X = repmat(read_point,1,size(Ext_Values,2))';
else
Ext_Values = to_plotY{o}(to_plotX(:,o)==read_point,j_del);
Scatt_X = repmat(read_point,1,size(Ext_Values,2))';
end
hline = line([read_point read_point], [-10 max(Ext_Values*1.2)],'Color','red');

ylim([-30 max(Ext_Values)*1.05]) 
xlim([-50 300])


scatter(Scatt_X,Ext_Values, 'LineWidth',2,'MarkerEdgeColor','red')

%Fitting of sample
if isempty(CH_del{o})
X_Conc = conc;
else
X_Conc = conc(setdiff(1:12,CH_del{o},'stable'));
end


y_to_fit = Ext_Values(1:end); 
% y_to_fit1 = y_to_fit/max(y_to_fit);
X_Conc = X_Conc(1:end);

[xData, yData] = prepareCurveData( X_Conc, y_to_fit);
% new custom equaiton - (1099.4*KP*(902/677.5)*x)/(1+sigma*1099.4*KP*x)
% new euation - (0.6564*KP*(0.6775/0.902)*x)/(1+(sigma*0.6564*KP*x))
% Set up fittype and options.
ft = fittype( 'Rmax/(1+Kd/x)+off', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [4000 40 0.157613081677548];
opts.Upper = [Inf Inf Inf];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

subplot(2,2,2)
scatter(X_Conc,y_to_fit,'LineWidth',2,'MarkerEdgeColor','red')
title('K_d')
hold on 
plot(fitresult)

xlabel 'Peptide concentration [\muM]'
ylabel 'Relative response unit'

b = gca; legend(b,'off');
% grid on
% Annotation 
dim = [.9 .4 .3 .3];
str = {sprintf('Kd= %0.3f uM',fitresult.Kd),sprintf('Rmax = %0.3f',fitresult.Rmax),...
    sprintf('off= %0.3f',fitresult.off),sprintf('sse = %0.3f',gof.sse),...
    sprintf('rsquare = %0.3f', gof.rsquare)};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',8)

%Fit Kp

y_to_fit = Ext_Values(1:end); 
% y_to_fit1 = y_to_fit/max(y_to_fit);

y_to_fit = y_to_fit/RUL(1);
[xData_1, yData_1] = prepareCurveData(X_Conc, y_to_fit);
%(0.6564*KP*(900/655.7)*(x))+off
% Set up fittype and options.
ft = fittype( sprintf('(%.4f*KP*(%.2f/%.2f)*(x*1E-6))/(1+sigma*%.4f*KP*x*1E-6)+off',mol_volume,mol_weight,lipid_weight,mol_volume),...
    'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxFunEvals = 500;
opts.StartPoint = [10 0.765516788149002 200];

% Fit model to data.
[fitresult2, gof] = fit( xData_1, yData_1, ft, opts );

subplot(2,2,4)
scatter(X_Conc,y_to_fit,'LineWidth',2,'MarkerEdgeColor','red')
title('K_p')
hold on 
plot(fitresult2)

xlabel 'Peptide concentration [\muM]'
ylabel 'RU_S/RU_L'

b = gca; legend(b,'off');


dim = [.9 .05 .3 .3];
str = {sprintf('Kp= %0.2f',fitresult2.KP),sprintf('sse = %0.3f',gof.sse),...
    sprintf('sigma= %0.2f',fitresult2.sigma),sprintf('rsquare = %0.3f', gof.rsquare)};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',8)


hFig    = gcf;
hAxes   = findobj(allchild(hFig), 'flat', 'Type', 'axes');
axesPos = get(hAxes, 'Position');
set(hAxes(1), 'Position', [0.65 0.1100 0.3347 0.3575]);
set(hAxes(2), 'Position', [0.65 0.5675 0.3347 0.3575]);
set(hAxes(3), 'Position', [0.0500 0.1100 0.5347 0.8150]);

print(sprintf('%s Channel %d %s.png',name,o+1,lipid),'-dpng')
close all 

Kd(o) =  fitresult.Kd;
Kp(o) = fitresult2.KP;
sigma(o) = fitresult2.sigma;
RUSL_exp{o} = y_to_fit;
RUSL_fit_exp{o} = fitresult2(0:max(conc));
end


% %------------------------------
save('Results.mat','to_plotX1','to_plotY1','to_plotX2','to_plotY2',...
    'to_plotX3','to_plotY3','RUL','Kd','Kp','sigma','RUSL_exp',...
    'RUSL_fit_exp','X_Conc','to_plotX','to_plotY')

