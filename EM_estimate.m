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


P.DataFile = 'Data\ES\Processed\Data_ES';
P.Legend = 'Data\ES\Original\Legend_ES';


% Best ex-post parameterisations
% with AR(1) idiosyncratic component
P.idio = '';
% P.p = 1;P.r = [1];   funPRT_ML(P);
P.p = 1;P.r = [2];   funPRT_ML(P);

%--------------------------------------------------------------------------
% Graphical rapresentation 
%--------------------------------------------------------------------------

% Converti DateQQ da numero seriale MATLAB a datetime
DateQQ_dt = datetime(DateQQ, 'ConvertFrom', 'datenum');

% Converti covid_start e covid_end da array [Anno, Mese, Giorno] a datetime
covid_start_dt = datetime(P.covid_start(1), P.covid_start(2), P.covid_start(3));
covid_end_dt = datetime(P.covid_end(1), P.covid_end(2), P.covid_end(3));

% Trova gli indici per il periodo Pre-COVID e Post-COVID
preCovidIdx = DateQQ_dt < datetime(P.covid_start);
postCovidIdx = DateQQ_dt > datetime(P.covid_end);

iCovS = find(DateQQ_V(:,1) == P.covid_start(1) & DateQQ_V(:,2) == P.covid_start(2));
iCovE = find(DateQQ_V(:,1) == P.covid_end(1) & DateQQ_V(:,2) == P.covid_end(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Forecast Vs Real Values %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Funzione anonima per generare le etichette della legenda
generateLabel = @(h) sprintf('Q(%+d)M%d', floor((h-1)/3) - 1, mod(h-1,3) + 1);

% ---------------- Whole Sample ----------------

figure;
hold on;

% Disegna il valore reale del GDP
plot(DateQQ_dt, TrueQQ(:, iVar), 'k', 'LineWidth', 2); 

% Disegna le previsioni per ogni orizzonte temporale con colori diversi
colors = lines(9); % Creiamo una mappa di colori per i 9 forecast horizons
selectedH = 5:7;

for h = selectedH
    plot(DateQQ_dt, FcstQQ(:, h, iVar), '--', 'Color', colors(h, :), 'LineWidth', 1.5);
end

% Evidenzia il periodo COVID con una zona grigia
covid_x = [DateQQ_dt(iCovS), DateQQ_dt(iCovE), DateQQ_dt(iCovE), DateQQ_dt(iCovS)];
covid_y = ylim;
fill(covid_x, [covid_y(1) covid_y(1) covid_y(2) covid_y(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% Legenda e formattazione
xlabel('Time');
ylabel('GDP Level');
title('GDP Nowcasting: Forecast vs Actual');
legend(['True GDP', arrayfun(generateLabel, selectedH, 'UniformOutput', false)], ...
       'Location', 'Best');
grid on;
hold off;


% ---------------- PRE-COVID ----------------

figure;

subplot(2,1,1); % Primo grafico
hold on;
plot(DateQQ_dt(preCovidIdx), TrueQQ(preCovidIdx, 1), 'k', 'LineWidth', 2);

for h = selectedH
    plot(DateQQ_dt(preCovidIdx), FcstQQ(preCovidIdx, h, iVar), '--', 'Color', colors(h, :), 'LineWidth', 1.5);
end

xlabel('Time');
ylabel('GDP Level');
title('GDP Forecast vs Actual (Pre-COVID)');
legend(['True GDP', arrayfun(generateLabel, selectedH, 'UniformOutput', false)], ...
       'Location', 'Best');
grid on;
hold off;

% ---------------- POST-COVID ----------------
subplot(2,1,2); % Secondo grafico
hold on;
plot(DateQQ_dt(postCovidIdx), TrueQQ(postCovidIdx, 1), 'k', 'LineWidth', 2);

for h = selectedH
    plot(DateQQ_dt(postCovidIdx), FcstQQ(postCovidIdx, h, iVar), '--', 'Color', colors(h, :), 'LineWidth', 1.5);
end

xlabel('Time');
ylabel('GDP Level');
title('GDP Forecast vs Actual (Post-COVID)');
legend(['True GDP', arrayfun(generateLabel, selectedH, 'UniformOutput', false)], ...
       'Location', 'Best');
grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-COVID vs Post-COVID Error Boxplot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
boxplot([rmsfeGDP_preCov, rmsfeGDP_postCov], {'Pre-COVID', 'Post-COVID'});
ylabel('RMSFE');
title('Forecast Error Distribution: Pre-COVID vs Post-COVID');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scatter Plot of Errors  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Converti FcstQQ in un vettore bidimensionale corretto
forecast_errors = squeeze(FcstQQ(:, :, iVar)) - TrueQQ(:, iVar);

% Grafico degli errori di previsione
figure;
scatter(DateQQ_dt, forecast_errors, 'filled');
hold on;

% Linee verticali per l'inizio e la fine del COVID
xline(datetime(covid_start(1), covid_start(2), covid_start(3)), '--k', 'COVID Start', 'LineWidth', 1.5);
xline(datetime(covid_end(1), covid_end(2), covid_end(3)), '--k', 'COVID End', 'LineWidth', 1.5);

% Linea orizzontale per lo zero
yline(0, '--r', 'Zero Error', 'LineWidth', 1.5);

% Formattazione asse X
xtickformat('yyyy-Q');
grid on;
xlabel('Time'); ylabel('Forecast Error');
title('Forecast Errors Over Time');
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Istogramm of Errors  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
histogram(FcstQQ(:, iVar, 1) - TrueQQ(:, iVar), 100, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('Forecast Error'); ylabel('Frequency');
title('Histogram of Forecast Errors');
grid on;

