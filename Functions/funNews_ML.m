function funPRT_ML(P)

P
r = P.r;
p = P.p;

% loading the details for the forecast evaluation
StartEst = P.StartEst;
StartEv  = P.StartEv;
EndEv    = P.EndEv ;
Qnews = P.Qnews;
SerNews = P.SerNews;

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
% computing the forecasts and news
% forecast updates are performed each month
%--------------------------------------------------------------------------

iS = find(DatesV(:,1) == StartEv(1) & DatesV(:,2) == StartEv(2));
iE = find(DatesV(:,1) == EndEv(1) & DatesV(:,2) == EndEv(2));
iQ = find(DatesV(:,1) == Qnews(1) & DatesV(:,2) == Qnews(2));
iSer = find(ismember(Series,SerNews));


Month = mod(DatesV(:,2),3);
Month(Month == 0) = 3;

% data at the start of the forecast sequence
Month_i = Month(iS-1);
% first unbalancedeness pattern
eval(['UnbP = UnbPattM',int2str(Month_i),';']);
X = Data(1:iS-1-nn,:);
temp = X(end-nUnb+1:end,:);
temp(isnan(UnbP)) = nan;
X(end-nUnb+1:end,:) = temp;
X_old = [X;nan(max(0,iQ-(iS-1-nn)),nM+nQ)];


DatesNews = Dates(iS:iE);
y_old = zeros(iE-iS+1,1);
y_new = zeros(iE-iS+1,1);
groupnews = zeros(iE-iS+1,length(unique(Group)));
singlenews = zeros(iE-iS+1,nM+nQ);

for i = iS:iE
       
    Date_i = DatesV(i,:);
    Month_i = Month(i);
    disp(['Computing the news for the vintages: y', num2str(DatesV(i,1)),' m',num2str(DatesV(i,2))])
    
    
    % first unbalancedeness pattern
    eval(['UnbP = UnbPattM',int2str(Month_i),';']);
    
    % updated data set
    X = Data(1:i-nn,:);
    temp = X(end-nUnb+1:end,:);
    temp(isnan(UnbP)) = nan;
    X(end-nUnb+1:end,:) = temp;
    X_new = [X;nan(max(0,iQ-(i-nn)),nM+nQ)];
    
    T_o = size(X_old,1);
    T_n = size(X_new,1);
    X_old = [X_old;nan(T_n-T_o,nM+nQ)];
 
    % new estimates and forecasts
    eval(['R_new = EM_DFM_SS',P.idio,P.restr,'(X_new,P);'])
    R_new.Groups = Group;
    R_new.Series = Series;
    
    % news
    [y_old(i-iS+1,1),y_new(i-iS+1,1),groupnews(i-iS+1,:),singlenews(i-iS+1,:),...
        gainT,serGainT,actual(:,i-iS+1),fcst(:,i-iS+1)] = ...
        News_DFM_ML(X_old,X_new,R_new,iQ,iSer);
    
    X_old = X_new;

  
end

GroupNames = unique(Group)';
trueSer = Data(iQ,iSer);
DateQQ = Dates(iS:iE);
DateQQ_V = DatesV(iS:iE,:);
check = y_new-y_old-sum(groupnews,2)

%--------------------------------------------------------------------------
% saving the results
%--------------------------------------------------------------------------
% From file P.DataFile Extract just the name of the nation
PathParts = strsplit(P.DataFile, filesep); 
NationName = PathParts{2};  

datafile = ['news', NationName, '_', P.idio,strrep(int2str(P.r),' ',''),int2str(P.p)];
eval(['save ',datafile,...
    ' y_old y_new trueSer DateQQ DateQQ_V Data Series GroupNames P groupnews singlenews fcst actual'])




