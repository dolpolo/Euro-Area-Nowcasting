clear
clc

% Set Base_Dir
Base_Dir  = cd;
addpath(genpath(Base_Dir));


% start of the estimation sample
P.StartEst = [2000 4];
% first forecast produced in
P.StartEv  = [2017 1];
% last forecast produced in
P.EndEv    = [2024 12];
% forecast horizon
P.fcstH_Q = 1;

% data set
P.Model = 'Medium'; 

% covid period 
P.covid_start = [2020, 3, 1];
P.covid_end = [2020, 12, 1];


P.DataFile = 'Data\IT\Processed\Data_IT';
P.Legend = 'Data\IT\Original\Legend_IT';

P.SerFcst = 'GDP_IT';
funPRT_BenchQ(P)


%--------------------------------------------------------------------------
% Graphical rapresentation 
%--------------------------------------------------------------------------

% Converti DateQQ da numero seriale MATLAB a datetime
DateQQ_dt = datetime(DateQQ, 'ConvertFrom', 'datenum');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forecast Vs Real Values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on;

% Disegna il valore reale del GDP
plot(DateQQ_dt, TrueQQ(:, 1), 'k', 'LineWidth', 2); 

plot(DateQQ_dt, FcstQQ(:, 7, 1), '--', 'LineWidth', 1.5);

plot(DateQQ_dt, FcstARQ(: , 1), '--', 'LineWidth', 1.5);

plot(DateQQ_dt, FcstRW1Q(: , 1), '--', 'LineWidth', 1.5);

plot(DateQQ_dt, FcstRW2Q(: , 1), '--', 'LineWidth', 1.5);

% Legenda e formattazione
xlabel('Time');
ylabel('GDP Level');
title('GDP Nowcasting: Forecast vs Actual');
grid on;
hold off;