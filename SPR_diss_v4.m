clear all 
close all
load('Results.mat')

% variables
name = 'Lulworthinone';
lipid = 'DMPC';


%Conversion to Ru s diss
diss_time = 200;
to_plotX1 = to_plotX1-diss_time;

%extract intensity for each time
to_anal = to_plotY1(:,4:15);



RU_max = (to_anal(to_plotX1==0,1:end));
%figure plotting
figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
subplot(2,2,[1 3]) 
hold on
    for j = [1:12]
        plot(to_plotX1,to_anal(:,j),'LineWidth', 2);
       
    end
  title(sprintf('%s %s Disociation Fit',name,lipid))
xlabel('time[s]') % x-axis label
ylabel('Relative response unit') % y-axis label
xlim([-100 350])
ylim([-Inf Inf])

subplot(2,2,2)
hold on 
for i = 1:to_plotX1(end)
trans_y =  (to_anal(to_plotX1==i,1:end));
trans_x = (RU_max-trans_y);
scatter(trans_x,trans_y,'k')

[xData, yData] = prepareCurveData( trans_x, trans_y );
% Set up fittype and options.
ft = fittype( 'Sl/(1-Sl)*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.8002804688888;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

plot(0:max(trans_x),fitresult(0:max(trans_x)),'k')

Sl(i) = fitresult.Sl;
end

xlabel('RU_{S,max} - RU_{S,diss}') % x-axis label
ylabel('RU_{S,diss}') % y-axis label

subplot(2,2,4)

sl_y = Sl(1:200)+(1-Sl(1));
sl_x = 0:size(sl_y,2)-1;
test_x = 0:size(Sl,2)-1

[xData, yData] = prepareCurveData( sl_x, sl_y );

% Set up fittype and options.
ft = fittype( '(alpha*(2.71828^(-koffa*x)))+(beta*(2.71828^(-koffb*x)))+Sret', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxIter = 4000000;
opts.StartPoint = [0.792207329559554 0.959492426392903 0.655740699156587 0.0357116785741896 0.849129305868777];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

scatter(sl_x,sl_y,'LineWidth',0.8,'MarkerEdgeColor','red')
% title('K_p')
hold on 
plot(fitresult)
% hold off
xlabel 'time[s]'
ylabel 'S_{L}'
xlim([0 200])
b = gca; legend(b,'off');

dim = [.9 .05 .3 .3];
str = {sprintf('S_{ret}= %0.2f',fitresult.Sret),sprintf('alpha = %0.2f',fitresult.alpha),...
    sprintf('beta= %0.2f',fitresult.beta),sprintf('K_{off alpha} = %0.3f s^{-1}', fitresult.koffa),...
    sprintf('K_{off beta} = %0.3f s^{-1}',fitresult.koffb),sprintf('rsquare = %0.3f', gof.rsquare)};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',8)

m = 1;
alpha(m) =  fitresult.alpha;
beta(m) = fitresult.beta;
koffa(m) = fitresult.koffa;
koffb(m) = fitresult.koffb;
Sret(m) = fitresult.Sret;
%calculation mediam koff
% koff = (alpha*koffa + beta*koffb)/(alpha+beta)
koff(m) = (alpha(m)*koffa(m)+beta(m)*koffb(m))/(alpha(m)+beta(m));

dim2 = [.75 .05 .3 .3];
str2 = {sprintf('K_{off} = %0.4f s^{-1}',koff(m))};
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','FontSize',10)

hFig    = gcf;
hAxes   = findobj(allchild(hFig), 'flat', 'Type', 'axes');
axesPos = get(hAxes, 'Position');
set(hAxes(1), 'Position', [0.55 0.1100 0.4347 0.3575]);
set(hAxes(2), 'Position', [0.55 0.5675 0.4347 0.3575]);
set(hAxes(3), 'Position', [0.0500 0.1100 0.4347 0.8150]);
% set(hAxes(1), 'Position', [0.65 0.1100 0.3347 0.3575]);

print(sprintf('%s %s Channel 2 Diss fit.png',name,lipid),'-dpng')
close all 
Sl_exp{m} = Sl;
Sl_fit_exp{m} = fitresult(0:200); 

% slt = ((1-Sret)*2.71828^(-koff*x))+Sret
% slt = (alpha*(2.71828^(-koffa*x)))+(beta*(2.71828^(-koffb*x)))+Sret 

%------------------------------------------------------------------
%---------------------------------------------------------------
to_plotX2 = to_plotX2-diss_time;

%extract intensity for each time
to_anal = to_plotY2(:,4:15);

% % moving to zero
% for i = 1:size(to_anal,2)
%    to_anal(:,i) = to_anal(:,1) - mean(to_anal(to_plotX1>=650,i));  
% end


RU_max = (to_anal(to_plotX2==0,1:end));
%figure plotting
figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
subplot(2,2,[1 3]) 
hold on
    for j = [1:12]
        plot(to_plotX2,to_anal(:,j),'LineWidth', 2);
       
    end
  title(sprintf('%s %s Disociation Fit',name,lipid))
xlabel('time[s]') % x-axis label
ylabel('Relative response unit') % y-axis label
xlim([-100 700])
ylim([-Inf Inf])

subplot(2,2,2)
hold on 
for i = 1:to_plotX2(end)
trans_y =  (to_anal(to_plotX2==i,1:end));
trans_x = (RU_max-trans_y);
scatter(trans_x,trans_y,'k')

[xData, yData] = prepareCurveData( trans_x, trans_y );
% Set up fittype and options.
ft = fittype( 'Sl/(1-Sl)*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.8002804688888;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

plot(0:max(trans_x),fitresult(0:max(trans_x)),'k')

Sl(i) = fitresult.Sl;
end

xlabel('RU_{S,max} - RU_{S,diss}') % x-axis label
ylabel('RU_{S,diss}') % y-axis label

subplot(2,2,4)

sl_y = Sl(1:200)+(1-Sl(1));
sl_x = 0:size(sl_y,2)-1;

[xData, yData] = prepareCurveData( sl_x, sl_y );

% Set up fittype and options.
ft = fittype( '(alpha*(2.71828^(-koffa*x)))+(beta*(2.71828^(-koffb*x)))+Sret', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxIter = 4000000;
opts.StartPoint = [0.792207329559554 0.959492426392903 0.655740699156587 0.0357116785741896 0.849129305868777];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

scatter(sl_x,sl_y,'LineWidth',0.8,'MarkerEdgeColor','red')
% title('K_p')
hold on 
plot(fitresult)
% hold off
xlabel 'time[s]'
ylabel 'S_{L}'
xlim([0 200])
b = gca; legend(b,'off');

dim = [.9 .05 .3 .3];
str = {sprintf('S_ret= %0.3f',fitresult.Sret),sprintf('alpha = %0.3f',fitresult.alpha),...
    sprintf('beta= %0.3f',fitresult.beta),sprintf('K_{off alpha} = %0.3f s^{-1}', fitresult.koffa),...
    sprintf('K_{off beta} = %0.3f s^{-1}',fitresult.koffb),sprintf('rsquare = %0.3f', gof.rsquare)};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',8)

m = 2;
alpha(m) =  fitresult.alpha;
beta(m) = fitresult.beta;
koffa(m) = fitresult.koffa;
koffb(m) = fitresult.koffb;
Sret(m) = fitresult.Sret;
%calculation mediam koff
% koff = (alpha*koffa + beta*koffb)/(alpha+beta)
koff(m) = (alpha(m)*koffa(m)+beta(m)*koffb(m))/(alpha(m)+beta(m));

dim2 = [.75 .05 .3 .3];
str2 = {sprintf('K_{off} = %0.4f s^{-1}',koff(m))};
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','FontSize',10)

hFig    = gcf;
hAxes   = findobj(allchild(hFig), 'flat', 'Type', 'axes');
axesPos = get(hAxes, 'Position');
set(hAxes(1), 'Position', [0.55 0.1100 0.4347 0.3575]);
set(hAxes(2), 'Position', [0.55 0.5675 0.4347 0.3575]);
set(hAxes(3), 'Position', [0.0500 0.1100 0.4347 0.8150]);
% set(hAxes(1), 'Position', [0.65 0.1100 0.3347 0.3575]);

print(sprintf('%s %s Channel 3 Diss fit.png',name,lipid),'-dpng')
close all 
Sl_exp{m} = Sl;
Sl_fit_exp{m} = fitresult(0:200); 

%------------------------------------------------------------------
%---------------------------------------------------------------
to_plotX3 = to_plotX3-diss_time;

%extract intensity for each time
to_anal = to_plotY3(:,4:15);

% % moving to zero
% for i = 1:size(to_anal,2)
%    to_anal(:,i) = to_anal(:,1) - mean(to_anal(to_plotX1>=650,i));  
% end


RU_max = (to_anal(to_plotX3==0,1:end));
%figure plotting
figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
subplot(2,2,[1 3]) 
hold on
    for j = [1:12]
        plot(to_plotX3,to_anal(:,j),'LineWidth', 2);
       
    end
  title(sprintf('%s %s Disociation Fit',name,lipid))
xlabel('time[s]') % x-axis label
ylabel('Relative response unit') % y-axis label
xlim([-100 700])
ylim([-Inf Inf])

subplot(2,2,2)
hold on 
for i = 1:to_plotX3(end)
trans_y =  (to_anal(to_plotX3==i,1:end));
trans_x = (RU_max-trans_y);
scatter(trans_x,trans_y,'k')

[xData, yData] = prepareCurveData( trans_x, trans_y );
% Set up fittype and options.
ft = fittype( 'Sl/(1-Sl)*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.8002804688888;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

plot(0:max(trans_x),fitresult(0:max(trans_x)),'k')

Sl(i) = fitresult.Sl;
end

xlabel('RU_{S,max} - RU_{S,diss}') % x-axis label
ylabel('RU_{S,diss}') % y-axis label

subplot(2,2,4)

sl_y = Sl(1:200)+(1-Sl(1));
sl_x = 0:size(sl_y,2)-1;

[xData, yData] = prepareCurveData( sl_x, sl_y );

% Set up fittype and options.
ft = fittype( '(alpha*(2.71828^(-koffa*x)))+(beta*(2.71828^(-koffb*x)))+Sret', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.MaxIter = 4000000;
opts.StartPoint = [0.792207329559554 0.959492426392903 0.655740699156587 0.0357116785741896 0.849129305868777];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

scatter(sl_x,sl_y,'LineWidth',0.8,'MarkerEdgeColor','red')
% title('K_p')
hold on 
plot(fitresult)
% hold off
xlabel 'time[s]'
ylabel 'S_{L}'
xlim([0 200])
b = gca; legend(b,'off');

dim = [.9 .05 .3 .3];
str = {sprintf('S_ret= %0.3f',fitresult.Sret),sprintf('alpha = %0.3f',fitresult.alpha),...
    sprintf('beta= %0.3f',fitresult.beta),sprintf('K_{off alpha} = %0.3f s^{-1}', fitresult.koffa),...
    sprintf('K_{off beta} = %0.3f s^{-1}',fitresult.koffb),sprintf('rsquare = %0.3f', gof.rsquare)};
annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',8)

m = 3;
alpha(m) =  fitresult.alpha;
beta(m) = fitresult.beta;
koffa(m) = fitresult.koffa;
koffb(m) = fitresult.koffb;
Sret(m) = fitresult.Sret;
%calculation mediam koff
% koff = (alpha*koffa + beta*koffb)/(alpha+beta)
koff(m) = (alpha(m)*koffa(m)+beta(m)*koffb(m))/(alpha(m)+beta(m));

dim2 = [.75 .05 .3 .3];
str2 = {sprintf('K_{off} = %0.4f s^{-1}',koff(m))};
annotation('textbox',dim2,'String',str2,'FitBoxToText','on','FontSize',10)

hFig    = gcf;
hAxes   = findobj(allchild(hFig), 'flat', 'Type', 'axes');
axesPos = get(hAxes, 'Position');
set(hAxes(1), 'Position', [0.55 0.1100 0.4347 0.3575]);
set(hAxes(2), 'Position', [0.55 0.5675 0.4347 0.3575]);
set(hAxes(3), 'Position', [0.0500 0.1100 0.4347 0.8150]);
% set(hAxes(1), 'Position', [0.65 0.1100 0.3347 0.3575]);

print(sprintf('%s %s Channel 4 Diss fit.png',name,lipid),'-dpng')
close all 
Sl_exp{m} = Sl;
Sl_fit_exp{m} = fitresult(0:200); 

save('Diss_ress.mat','koff','Sret','Sl_exp','Sl_fit_exp','alpha','beta',...
    'koffa','koffb','Sret')

