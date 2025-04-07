clear
clc

% Set Base_Dir
Base_Dir  = cd;
addpath(genpath(Base_Dir));


P.DataFile = 'Data\FR\Processed\Data_FR';
P.Legend = 'Data\FR\Original\Legend_FR';

% forecast series
P.SerNews  = 'GDP_FR';
% forecast reference period
P.Qnews    = [2024 10];
% first forecast done in
P.StartEv  = [2024 6];
% last forecast done in
P.EndEv    = [2024 11];
% start of the estimation sample
P.StartEst = [2000 4];


% covid period 
P.covid_start = [2020, 3, 1];
P.covid_end = [2020, 12, 1];

% with AR(1) idiosyncratic component 
P.idio = ''; 

P.r = 1; P.p = 2;
P.Model = 'Medium';funNews_ML(P);
P.r = 1; P.p = 1;
P.Model = 'Small';funNews_ML(P);

%--------------------------------------------------------------------------
% Graphical representation
%--------------------------------------------------------------------------

% Select the covid period
covid_start = P.covid_start;
covid_end = P.covid_end;

% Converti DateQQ da numero seriale MATLAB a datetime (se necessario)
DateQQ_dt = datetime(DateQQ, 'ConvertFrom', 'datenum');

% Creazione del grafico
figure;
hold on;

% Linea per le previsioni vecchie
plot(DateQQ_dt, y_old, '--k', 'LineWidth', 1.5, 'DisplayName', 'Previsione Vecchia');

% Linea per le previsioni nuove
plot(DateQQ_dt, y_new, '-b', 'LineWidth', 2, 'DisplayName', 'Previsione Nuova');

% Punto per il valore reale (se disponibile)
if exist('trueSer', 'var') && exist('DateQQ_V', 'var')
    scatter(DateQQ_dt(end), trueSer, 30, 'r', 'filled', 'DisplayName', 'Valore Reale');
end

% Personalizzazione del grafico
legend('Location', 'Best');
xlabel('Data'); 
ylabel('Valore PIL');
title('Confronto tra Previsioni Vecchie e Nuove');
grid on;
hold off;


%  Impatto delle nuove informazioni (grafico a barre impilate)

figure;
bar(DateQQ, groupnews, 'stacked');
xlabel('Data'); ylabel('Impatto delle News sulla Previsione');
title('Impatto delle Nuove Informazioni sulle Previsioni');
legend(GroupNames, 'Location', 'BestOutside');
grid on;


%  Impatto delle singole variabili (heatmap)

figure;
imagesc(singlenews'); 
colorbar;
xlabel('Data'); ylabel('Variabili');
title('Impatto delle Singole Variabili sulle Previsioni');
set(gca, 'YTick', 1:length(Series), 'YTickLabel', Series);
colormap(jet); % Blu/rosso per variazioni negative/positive

% Impatto delle nuove informazioni (grafico a barre impilate)
figure;
hold on;

% Definizione di colori personalizzati (modificabili)
customColors = [
    0 0.4470 0.7410;   % Blu
    0.8500 0.3250 0.0980; % Arancione
    0.9290 0.6940 0.1250; % Giallo
    0.4940 0.1840 0.5560; % Viola
    0.4660 0.6740 0.1880; % Verde
    0.3010 0.7450 0.9330; % Azzurro
    0.6350 0.0780 0.1840; % Rosso scuro
    0.7350 0.0080 0.0840; % Rosso scuro
    0.0350 0.0780 0.1840; % Rosso scuro
    0.8350 0.1780 0.1840; % Rosso scuro
];

% Adatta la colormap alla dimensione del numero di gruppi
nGroups = size(groupnews, 2);
customColors = customColors(1:nGroups, :);

% Creazione del grafico a barre impilate con colori personalizzati
hb = bar(DateQQ, groupnews, 'stacked');

% Applicazione dei colori definiti sopra
for i = 1:nGroups
    hb(i).FaceColor = customColors(i, :);
end

% Personalizzazione del grafico
xlabel('Data'); 
ylabel('Impatto delle News sulla Previsione');
title('Impatto delle Nuove Informazioni sulle Previsioni');
legend(GroupNames, 'Location', 'BestOutside');
grid on;
hold off;

% Accuratezza delle previsioni (errore di previsione)

if exist('trueSer', 'var')
    figure;
    error_forecast = trueSer - y_new;
    plot(DateQQ, error_forecast, '-r', 'LineWidth', 2);
    yline(0, '--k');
    xlabel('Data'); ylabel('Errore di Previsione');
    title('Accuratezza delle Previsioni');
    grid on;
end


