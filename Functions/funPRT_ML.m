function funPRT_ML(P)

P
r = P.r;
p = P.p;

% loading the details for the forecast evaluation
StartEst = P.StartEst;
StartEv  = P.StartEv;
EndEv    = P.EndEv ;
fcstH_Q = P.fcstH_Q;

% Select the covid period
covid_start = P.covid_start;
covid_end = P.covid_end;

% Vectorize the covid period dates
covid_startV = datenum(covid_start);
covid_endV = datenum(covid_end);

% maximum number of EM iterations
P.max_iter = 300;

%--------------------------------------------------------------------------
% Loading monthly data
%--------------------------------------------------------------------------
DataFile = P.DataFile;
LegendFile = P.Legend;

[a,b] = xlsread(DataFile,'MonthlyLong');

DataM = a(1:end-1,:);
if strcmp(version('-release'),'2006b')
DatesM = datenum(b(2:end-1,1),'dd/m/yy');
else
     DatesM = datenum(b(2:end-1,1));
end
DatesMV = datevec(DatesM);
DatesMV = DatesMV(:,1:2);

T = length(DatesM);

[a,b] = xlsread(LegendFile,'MDescriptionFull');
b = b(2:end,:);
GroupM = b(:,16);
SeriesM = b(:,1);
%Transformation
TransfM = a(:,8);
TransfM = floor(TransfM);
% unbalancedeness patterns
UnbM = a(:,11:13);
% Type of Variable
ClassM = b(:,11);

% Transforming the data Annualizing and Substituting Real variables during
% covid with NaN
DataMTrf = DataM;
for i = 2:5
    DataMTrf(:,TransfM(:,1) == i) = 12*(DataMTrf(:,TransfM(:,1) == i));
end
[DataMTrf, Covid_obsM] = CovidNaN(DataMTrf, DatesM, ClassM, covid_startV, covid_endV);

% first differences
% DataMTrf(2:end,TransfM(:,2) == 1) = (DataMTrf(2:end,TransfM(:,2) == 1) ...
%    - DataMTrf(1:end-1,TransfM(:,2) == 1));
% DataMTrf(1,TransfM(:,2) == 1) = nan;

[tM,nM] = size(DataMTrf);
DataMTrf = [DataMTrf;nan(T-tM,nM)];

%--------------------------------------------------------------------------
% Loading quarterly data
%--------------------------------------------------------------------------
[a,b] = xlsread(DataFile,'QuarterlyLong');

DatesQ = datenum(b(2:end,1));
DataQ = a(:,:);

[a,b] = xlsread(LegendFile,'QDescriptionFull');

ListIdx = find(ismember(b(1,:), P.Model));
List = find(a(:,ListIdx-1)==1);

DataQ = DataQ(1:end,List); 

b = b(2:end,:);

GroupQ = b(:,16);
GroupQ = GroupQ(List);

SeriesQ = b(:,1);
SeriesQ = SeriesQ(List);

%Transformation
TransfQ = a(:,8);
TransfQ = TransfQ(List);
TransfQ = floor(TransfQ);

% unbalancedeness patterns
UnbQ = a(:,11:13);
UnbQ = UnbQ(List,:);

% Type of Variable
ClassQ = b(:,11);
ClassQ = ClassQ(List);

% Transforming the data Annualizing and Substituting Real variables during
% covid with NaN
DataQTrf = DataQ;
for i = 2:5
    DataQTrf(:,TransfQ(:,1) == i) = 4*(DataQTrf(:,TransfQ(:,1) == i));
end
[DataQTrf, Covid_obsQ] = CovidNaN(DataQTrf, DatesQ, ClassQ, covid_startV, covid_endV);

% quarterly at monthly frequency
DataQMTrf = kron(DataQTrf,[nan;nan;1]);

[tQ,nQ] = size(DataQMTrf);
DataQMTrf = [DataQMTrf;nan(T-tQ,nQ)];

% Building the matrix of restrictions on the loadings for the quarterly
% variables
P.Rconstr = [2 -1 0 0 0;...
    3 0 -1 0 0;...
    2 0 0 -1 0;...
    1 0 0 0 -1];
P.q = zeros(4,1);
P.restr = '_restrMQ';

%--------------------------------------------------------------------------
% complete dataset
%--------------------------------------------------------------------------
Data = [DataMTrf DataQMTrf];
Series = [SeriesM;SeriesQ];
Group = [GroupM;GroupQ];
UnbPatt = [UnbM;UnbQ];

% Starting evaluation Date
iEst = find(DatesMV(:,1) == StartEst(1) & DatesMV(:,2) == StartEst(2));

% Final Data Info
Data = Data(iEst:end,:); 
Dates = DatesM(iEst:end,:);
DatesV = DatesMV(iEst:end,:);

nVar = size(Data,2);
nQ = nVar-nM;

%--------------------------------------------------------------------------
% unbalancedness patterns
%--------------------------------------------------------------------------
nn = min(min(UnbPatt));
nn = min(nn,0);

UnbPattM1 = zeros(12-nn,nVar);
UnbPattM2 = zeros(12-nn,nVar);
UnbPattM3 = zeros(12-nn,nVar);

nUnb = 12-nn;

for i = 1:nVar
    UnbPattM1(end-UnbPatt(i,1)+1+nn:end,i) = nan;
    UnbPattM2(end-UnbPatt(i,2)+1+nn:end,i) = nan;
    UnbPattM3(end-UnbPatt(i,3)+1+nn:end,i) = nan;
end

P.nQ = nQ;

%--------------------------------------------------------------------------
% out-of-sample evaluation
%--------------------------------------------------------------------------

iS = find(DatesV(:,1) == StartEv(1) & DatesV(:,2) == StartEv(2));
iE = find(DatesV(:,1) == EndEv(1) & DatesV(:,2) == EndEv(2));

FcstQ = nan(iE,(fcstH_Q+2),nQ);

Month = mod(DatesV(:,2),3);
Month(Month == 0) = 3;

for i = iS:iE

    Date_i = DatesV(i,:);
    Month_i = Month(i);
    disp(['Computing the predictions for the vintages: y', num2str(DatesV(i,1)),' m',num2str(DatesV(i,2))])


    eval(['UnbP = UnbPattM',int2str(Month_i),';']);

    X = Data(1:i-nn,:);
    temp = X(end-nUnb+1:end,:);
    temp(isnan(UnbP)) = nan;
    X(end-nUnb+1:end,:) = temp;
    n_nan = max([0,(fcstH_Q+1)*3-Month_i+nn]);
    X = [X;nan(n_nan,nVar)];
    % estimation and forecast
    eval(['Res = EM_DFM_SS',P.idio,P.restr,'(X,P);'])
    % collecting the forecast for the quarterly variables
    FcstQ(i,:,:) = Res.X_sm([i-Month_i:3:i-Month_i+3*(fcstH_Q+1)],nM+1:end);

end
% For i = 202, FcstQ contiene le previsioni del trimestre passato, presente e futuro nelle righe
% 201,204,207 ( rispettivamente Q(-1)M3,Q(0)M3,Q(+1)M3)). In questo modo per il set di dati 
% della riga 201, un certo vintage, avrò diverse previsioni per i trimestri
% di interesse. quando passerò alla riga 203 avrò sempre le previsioni
% Q(-1)M3,Q(0)M3,Q(+1)M3) ma con un vintage aggiornato. Di tre righe in tre
% righ cambierò l'orizzonte temporale di porevisioni con i trimestri di
% interesse, sulle colonne avrò i le previsioni per i trimestri di
% interesse e su ogni pagina di FcstQ una certa variabile.

FcstQQ = [];
TrueQQ = [];

for k = iS:3:iE-((fcstH_Q+2)*3-1) % si ferma a 8 orizzonti temporali precedenti di modo da prevedere l'ultima osservazione
    f_t = [];
    for j = 1:fcstH_Q+2                                    %     j = 1     |     j = 2     |      j = 3 
        f_t = cat(2,f_t,FcstQ(k+3*(j-1),fcstH_Q+3-j,:),... % 202,3,: = 205 | 205,2,: = 207 |  208,1,: = 209
            FcstQ(k+3*(j-1)+1,fcstH_Q+3-j,:),...           % 203,3,: = 206 | 206,2,: = 208 |  209,1,: = 210
            FcstQ(k+3*(j-1)+2,fcstH_Q+3-j,:));             % 204,3,: = 207 | 207,2,: = 209 | 2010,1,: = 211
    end
    FcstQQ = cat(1,FcstQQ,f_t);
    % F_t = [
    % 202,3,: = 205 | 203,3,: = 206 | 204,3,: = 207 | 205,2,: = 207 | 206,2,: = 208 | 207,2,: = 209 | 208,1,: = 209 | 209,1,: = 210 | 210,1,: = 211
    % 205,3,: = 208 | 206,3,: = 209 | 207,3,: = 210 | 208,2,: = 210 | 209,2,: = 211 | 210,2,: = 212 | 211,1,: = 212 | 212,1,: = 213 | 213,1,: = 214 
    % 208,3,: = 211 | 209,3,: = 212 | 210,3,: = 213 | 211,2,: = 213 | 212,2,: = 214 | 213,2,: = 215 | 214,1,: = 215 | 215,1,: = 216 | 216,1,: = 217 
    % 211,3,: = 214 | 212,3,: = 215 | 213,3,: = 216 | 214,2,: = 216 | 215,2,: = 217 | 216,2,: = 218 | 217,1,: = 218 | 218,1,: = 219 | 219,1,: = 220
    % ... ]
    TrueQQ = [TrueQQ;Data(k+(fcstH_Q )*3 ,nM+1:end)]; % TrueQQ = [TrueQQ;Data(k+(fcstH_Q + 2)*3,nM+1:end)]
end
% Significato delle colonne di FcstQQ
% Colonne 1-3 → Previsioni del trimestre futuro usando dati dei 3 mesi del trimestre corrente.
% Colonne 4-6 → Previsioni del trimestre corrente usando dati dei 3 mesi del trimestre corrente.
% Colonne 7-9 → Previsioni del trimestre passato usando dati dei 3 mesi del trimestre corrente.

count = 1; % Contatore per Covid_obsQ
for i = 1:size(TrueQQ, 1)
    if any(isnan(TrueQQ(i, :))) % Se la riga ha NaN
        TrueQQ(i, :) = Covid_obsQ(count, :); % Sostituisci con la riga corrispondente di Covid_obsQ
        count = count + 1; % Passa alla riga successiva di Covid_obsQ
    end
end

DateQQ = Dates(iS+(fcstH_Q+1)*3-1:3:iE-3);
DateQQ_V = DatesV(iS+(fcstH_Q+1)*3-1:3:iE-3,:);

% From file P.DataFile Extract just the name of the nation
PathParts = strsplit(P.DataFile, filesep); 
NationName = PathParts{2};  

% RMSFE for GDP
nF = size(FcstQQ,2);
iVar = find(ismember(SeriesQ, ['GDP_', NationName]));

iCovS = find(DateQQ_V(:,1) == covid_start(1) & DateQQ_V(:,2) == covid_start(2));
iCovE = find(DateQQ_V(:,1) == covid_end(1) & DateQQ_V(:,2) == covid_end(2));

% Compute the Rmsfe just for quarters before and after covid, ignoring
% missing values during covid (eventually just for real variables)
rmsfeGDP  = sqrt(mean((FcstQQ(1:end,:,iVar) - repmat(TrueQQ(1:end,iVar),1,nF)).^2)');
rmsfeGDP_preCov   = sqrt(mean((FcstQQ(1:iCovS,:,iVar) - repmat(TrueQQ(1:iCovS,iVar),1,nF)).^2)');
rmsfeGDP_postCov  = sqrt(mean((FcstQQ(iCovE+1:end,:,iVar) - repmat(TrueQQ(iCovE+1:end,iVar),1,nF)).^2)');

%--------------------------------------------------------------------------
% saving the results
%--------------------------------------------------------------------------

% File Name
datafile = ['fcst_', NationName, '_', P.idio, strrep(int2str(P.r), ' ', ''), int2str(P.p)];

% Save as a .mat file
eval(['save ', datafile, ...
    ' FcstQQ TrueQQ DateQQ DateQQ_V Data Series SeriesQ Group iVar FcstQ' ...
    ' Dates DatesV P iS iE rmsfeGDP rmsfeGDP_preCov rmsfeGDP_postCov'])

