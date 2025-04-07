%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code performs the following:
%
% (1) Loads data for each Euro Area Country of interest. 
%     The input data are already de-seasonalized and stationarized from Barigozzi and Lissona.
% (2) Process the data and provides new working datasets where
%     - Real Variable for Covid are treated as Missing Values Through the function CovidNaN;
%     - Series are all annualized with respect to each series transformation 
%
%      - TR = 1 ---> c*log(x)                 ---> No Annualization
%	   - TR = 2 ---> \Delta log(c*x)          ---> If M == X12 if Q == X3
%      - TR = 3 ---> \Delta(\Delta log(c*x))  ---> If M == X12 if Q == X3
%      - TR = 4 ---> x (no transformation)    ---> No Annualization
%      - TR = 5 ---> \Delta x                 ---> No Annualization
%      - TR = 6 ---> \Delta(\Delta x)         ---> No Annualization
%
%	whenever TR = 1.5, 2.5, 4.5, the transformation depends on the
%	methodology: if light transformations, take floor(TR), if heavy
%	transformations, take ceil(TR). I Consider just light transformations.
% --------------------------------------------------------------------------
% PARAMETERS TO BE CHOSEN BY THE USER:
%  Preferences for the Covid Period
% ----------------------------------------------------------------------------------------------------------
% OUTPUT:
%		      Excel file (.xlsx) with transformed data annualized exluding
%		      real variables for Covid period for the following cases:
%             For every i in {DE,FR,IT,SP}
%             - Monthly Variables Only Treated and Annualized (TA_Mi)
%             - Quarterly Variables Only (TA_Qi)
%             - Monthly and Quarterly Variables Jointly (TA_MQi)
% ----------------------------------------------------------------------------------------------------------
% SEE ALSO: EAtransform, remove_outliers, EMimputation
% TO NOTE: Data during the Covid perod (2020Q1-2020Q4) for real variables aimed to be 
%     imputed with the EM algorithm exploting the information included in nominal and financial series.
% ----------------------------------------------------------------------------------------------------------
% REFERENCES:
% 
% -----------------------------------------------------------------------------------------------------------
% Author: Davide Delfino (davide.delfino@unibo.it)
% Last update: 30/Mar/2025
% Version: MATLAB 2023b
% Required Toolboxes: /
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; close all; clc;


%% Set Directory and Parameters

% Set Base Directory and Add Paths
Base_Dir = cd;
addpath(genpath(Base_Dir));

% Load Datasets
Legend = 'Data/EA/Original/Legend_EA.xlsx';
data_M = readtable('Data/EA/Original/EAdataM_LT.xlsx');
data_Q = readtable('Data/EA/Original/EAdataQ_LT.xlsx');

% NationCode 
PathParts = strsplit(Legend, filesep);
PathParts = regexp(Legend, '[\/]', 'split');  
NationCode = PathParts{2}; 

% Output directory personalizzata per la nazione
OutputFile = ['Data/', NationCode, '/Processed/'];

% Creare la cartella di output se non esiste
if ~exist(OutputFile, 'dir')
    mkdir(OutputFile);
end

% COVID Period Window
covid_start = datenum(2020,3,1);
covid_end = datenum(2020,12,1);

%% Merge Monthly and Quarterly Data Exluding Covid and Annualizing

[a, b] = xlsread(Legend, 'Description');
Class = b(2:end,11);

time_col = 'Time';

% Back to the ordering 
all_var = ['Time'; b(2:end, 1)];  % Correct variable order with Time
all_var = regexprep(all_var, '\.', '_');  % Sanity Check 

% Find Variables to Add from Quarterly Data
vars_to_add = setdiff(data_Q.Properties.VariableNames, data_M.Properties.VariableNames);

merged_data = outerjoin(data_M, data_Q(:, [time_col, vars_to_add]), ...
                        'Keys', time_col , 'MergeKeys', true, 'Type', 'left');

% Order Variables According to all_var
ordered_vars = intersect(all_var, merged_data.Properties.VariableNames, 'stable');
merged_data = merged_data(:, ordered_vars);

% Excel File with non-treated Quarterly and Monthly Variables Jointly
output_path_MQ = [OutputFile, 'MQ_', NationCode, '.xlsx'];
writetable(merged_data, output_path_MQ);
disp(['Non-treated M and Q dataset saved at: ', output_path_MQ ]);

% Extract the 'Time' column from the table and convert it into detenum
Time = merged_data(:, 1);
TimeDate = merged_data.Time;
TimeNum = datenum(TimeDate);

% Remove the "time" column
merged_data = merged_data(:, 2:end);

% Ensure Light_TR and Frequence are properly assigned
Light_TR = floor(a(:, 7));  % Assuming this is the correct assignment for Light_TR
Frequence = b(2:end, 3);    % Assuming this is the correct assignment for Frequence

% Get the number of variables in the merged dataset
num_vars = size(merged_data, 2);

% Loop through all variables
for i = 1:num_vars
    % Apply the condition based on Light_TR and Frequence
    if Light_TR(i) == 2 || Light_TR(i) == 3 || Light_TR(i) == 4 || Light_TR(i) == 5
        if strcmp(Frequence{i}, 'Q')  % Frequency "Q"
            merged_data{:, i} = merged_data{:, i} * 4;  % Multiply by 4
        elseif strcmp(Frequence{i}, 'M')  % Frequency "M"
            merged_data{:, i} = merged_data{:, i} * 12;  % Multiply by 12
        end
    end
end

% Set Real Variables == NaN during COVID (convert data in double)
merged_data_d = table2array(merged_data);
merged_data = CovidNaN(merged_data_d, TimeNum, Class, covid_start, covid_end);
merged_data = array2table(merged_data, 'VariableNames', ordered_vars(2:end));

merged_data = [Time, merged_data];

% Save Merged Dataset con il nome della nazione
output_path_merged = [OutputFile, 'TA_MQ_', NationCode, '.xlsx'];
writetable(merged_data, output_path_merged);
disp(['Merged dataset saved at: ', output_path_merged]);

%% Disaggregate all quarterly Time seires and Annualize them 

[a, b] = xlsread(Legend, 'QDescriptionFull');
Class = b(2:end,11);

Q_all_var = ['Time'; b(2:end, 1)];  % Correct quarterly variable order with Time
Q_all_var = regexprep(Q_all_var, '\.', '_');  % Sanity Check

data_Q_filtered = data_Q(:, [time_col, vars_to_add]);

Time = data_Q_filtered(:,1);
TimeDate = data_Q_filtered.Time;
TimeNum = datenum(TimeDate);

% Order Variables According to the original pattern
ordered_vars = intersect(Q_all_var, data_Q_filtered.Properties.VariableNames, 'stable');
data_Q_filtered = data_Q_filtered(:, ordered_vars);

% Save Quarterly Variables Only
output_path_MandQ = [OutputFile, 'Data_', NationCode, '.xlsx'];
writetable(data_M, output_path_MandQ, sheet = 'MonthlyLong');
writetable(data_Q_filtered, output_path_MandQ, sheet = 'QuarterlyLong');
disp(['Q and M seprated and NonTreated dataset saved at: ', output_path_MandQ]);

% Remove the "time" column
data_Q_filtered = data_Q_filtered(:, 2:end);

% Ensure Light_TR are properly assigned
Light_TR = floor(a(:, 7));  % Assuming this is the correct assignment for Light_TR

% Get the number of variables in the merged dataset
num_vars = size(data_Q_filtered, 2);

% Loop through all variables
for i = 1:num_vars
    % Apply the condition based on Light_TR
    if Light_TR(i) == 2 || Light_TR(i) == 3 || Light_TR(i) == 4 || Light_TR(i) == 5
        data_Q_filtered{:, i} = data_Q_filtered{:, i} * 4;  % Multiply by 4
    end
end

% Set Real Variables == NaN during COVID (convert data in double)
if istable(data_Q_filtered)
    data_Q_filtered_d = table2array(data_Q_filtered);  % Converte in array numerico
else
    data_Q_filtered_d = data_Q_filtered; 
end

% CovidNaN
data_Q_filtered_d = CovidNaN(data_Q_filtered_d, TimeNum, Class, covid_start, covid_end);

% Se il dataset originale era una tabella, lo riconvertiamo in tabella
if istable(data_Q_filtered)
    data_Q_filtered = array2table(data_Q_filtered_d, 'VariableNames', ordered_vars(2:end));
else
    data_Q_filtered = data_Q_filtered_d;  
end

data_Q_filtered = [Time, data_Q_filtered];

% Save Quarterly Variables Only
output_path_Q = [OutputFile, 'TA_Data_', NationCode, '.xlsx'];
writetable(data_Q_filtered, output_path_Q, sheet = 'QuarterlyLong');
disp(['Quarterly-only dataset saved at: ', output_path_Q]); 

%% Annualize Monthly Time Series

[a, b] = xlsread(Legend, 'MDescriptionFull');
Class = b(2:end,11);
SeriesM = b(2:end,1);

Time = data_M(:,1);
TimeDate = data_M.Time;
TimeNum = datenum(TimeDate);

% Remove the "time" column
data_M = data_M(:, 2:end);

% Ensure Light_TR are properly assigned
Light_TR = floor(a(:, 7));  % Assuming this is the correct assignment for Light_TR

% Get the number of variables in the merged dataset
num_vars = size(data_M, 2);

% Loop through all variables
for i = 1:num_vars
    % Apply the condition based on Light_TR
    if Light_TR(i) == 2 || Light_TR(i) == 3 || Light_TR(i) == 4 || Light_TR(i) == 5
        data_M{:, i} = data_M{:, i} * 12;  % Multiply by 4
    end
end

% Set Real Variables == NaN during COVID (convert data in double)
data_M_d = table2array(data_M);
data_M = CovidNaN(data_M_d, TimeNum, Class, covid_start, covid_end);
data_M = array2table(data_M, 'VariableNames', SeriesM);

data_M = [Time, data_M];

% Save Quarterly Variables Only
output_path_M = [OutputFile, 'TA_Data_', NationCode, '.xlsx'];
writetable(data_M, output_path_M, sheet = 'MonthlyLong');
disp(['Monthly-only dataset saved at: ', output_path_M]);
