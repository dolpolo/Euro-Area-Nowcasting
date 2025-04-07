function [DataCovid, OriginalVar] = CovidNaN(X, T, C, C19_S, C19_E) 
    % Function to substitute the real variables with NaN during Covid-19
    % and preserve the original data for all variables during the COVID period.
    %
    % Input:
    %   X     - Data
    %   T     - Dates in vector form (datenum)
    %   C     - Cell array with Class specification
    %   C19_S - Starting period of COVID (datenum)
    %   C19_E - Ending period of COVID (datenum)
    %
    % Output:
    %   DataCovid    - Data after Covid manipulation
    %   OriginalVar - Original data (all variables) during the COVID period

    % Rows index for COVID period
    covid_idx = (T >= C19_S) & (T <= C19_E);

    % Salva tutte le variabili originali nel periodo COVID
    OriginalVar = X(covid_idx, :);

    % Identifica le variabili reali e sostituisce con NaN solo quelle nel periodo COVID
    real_vars_idx = strcmp(C, 'R'); 
    X(covid_idx, real_vars_idx) = NaN;

    % Output della matrice trasformata
    DataCovid = X;
end
