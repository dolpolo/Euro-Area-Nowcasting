function funPRT_Bench(P)

P
% loading the details for the forecast evaluation
StartEst = P.StartEst;
StartEv  = P.StartEv;
EndEv    = P.EndEv ;
fcstH_Q = P.fcstH_Q;
SerFcst = P.SerFcst;

% Select the covid period
covid_start = P.covid_start;
covid_end = P.covid_end;

% Vectorize the covid period dates
covid_startV = datenum(covid_start);
covid_endV = datenum(covid_end);

P.IC = 'AIC';
P.cflag = 1;
P.lags = 5;

%--------------------------------------------------------------------------
% Loading monthly data
%--------------------------------------------------------------------------
DataFile = P.DataFile;
LegendFile = P.Legend;

[a,b] = xlsread(DataFile,'MonthlyLong');

DataM = a(:,:);
if strcmp(version('-release'),'2006b')
DatesM = datenum(b(4:end,1),'dd/m/yy');
else
    DatesM = datenum(b(4:end,1));
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


%--------------------------------------------------------------------------
% complete dataset
%--------------------------------------------------------------------------

Data = [DataQMTrf];
Series = [SeriesQ];
Group = [GroupQ];
UnbPatt = [UnbQ];

iEst = find(DatesMV(:,1) == StartEst(1) & DatesMV(:,2) == StartEst(2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = Data(iEst:end,:);         %
Dates = DatesM(iEst:end,:);      %
DatesV = DatesMV(iEst:end,:);     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nVar = size(Data,2);
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



%--------------------------------------------------------------------------
% out-of-sample evaluation
%--------------------------------------------------------------------------
iSer = find(ismember(Series,SerFcst));

iS = find(DatesV(:,1) == StartEv(1) & DatesV(:,2) == StartEv(2));
iE = find(DatesV(:,1) == EndEv(1) & DatesV(:,2) == EndEv(2));

Month = mod(DatesV(:,2),3);
Month(Month == 0) = 3;

FcstRW1 = nan(iE,(fcstH_Q+2));
FcstRW2 = nan(iE,(fcstH_Q+2));
FcstAR = nan(iE,(fcstH_Q+2));

for i = iS:iE
       
    Date_i = DatesV(i,:);
    Month_i = Month(i);
    disp(['Computing the predictions for the vintages: y', num2str(DatesV(i,1)),' m',num2str(DatesV(i,2))])
    
    
    % first unbalancedeness pattern
    eval(['UnbP = UnbPattM',int2str(Month_i),';']);
    
    X = Data(1:i-nn,:);
    temp = X(end-nUnb+1:end,:);
    temp(isnan(UnbP)) = nan;
    X(end-nUnb+1:end,:) = temp;
    n_nan = max([0,(fcstH_Q+1)*3-Month_i+nn]);
    X = [X;nan(n_nan,nQ)];
    
    % random walk with drift in levels
    x = X(3:3:end,iSer);
    FcstRW1(i,:) = nanmean(x);

    % random walk without drift for growth rates
    T = length(x);
    iNaN = isnan(x);
    leadNaN = max(find(cumsum(iNaN) == (1:T)'));
    if isempty(leadNaN) leadNaN = 0; end
    h = max(find(cumsum(iNaN(end:-1:1)) == (1:T)'));
    FcstRW2(i,:) = x(end-h);
    
    % AR model
    xFcst = V_AR(x(leadNaN+1:end-h),P,h);
    xFcst = [nan(leadNaN,1);xFcst];
    xFcst = kron(xFcst,[nan;nan;1]);
    FcstAR(i,:) = xFcst([i-Month_i:3:i-Month_i+3*(fcstH_Q+1)]);

   
end

FcstARQ = [];
FcstRW1Q = [];
FcstRW2Q = [];
TrueQQ = [];
for k = iS:3:iE-((fcstH_Q+2)*3-1)
    f_t = [];
    for j = 1:fcstH_Q+2
        f_t = cat(2,f_t,FcstAR(k+3*(j-1),fcstH_Q+3-j),...
            FcstAR(k+3*(j-1)+1,fcstH_Q+3-j),...
            FcstAR(k+3*(j-1)+2,fcstH_Q+3-j));
    end
    FcstARQ = cat(1,FcstARQ,f_t);
    f_t = [];
    for j = 1:fcstH_Q+2
        f_t = cat(2,f_t,FcstRW1(k+3*(j-1),fcstH_Q+3-j),...
            FcstRW1(k+3*(j-1)+1,fcstH_Q+3-j),...
            FcstRW1(k+3*(j-1)+2,fcstH_Q+3-j));
    end
    FcstRW1Q = cat(1,FcstRW1Q,f_t);
    f_t = [];
    for j = 1:fcstH_Q+2
        f_t = cat(2,f_t,FcstRW2(k+3*(j-1),fcstH_Q+3-j),...
            FcstRW2(k+3*(j-1)+1,fcstH_Q+3-j),...
            FcstRW2(k+3*(j-1)+2,fcstH_Q+3-j));
    end
    FcstRW2Q = cat(1,FcstRW2Q,f_t);
    TrueQQ = [TrueQQ;Data(k+(fcstH_Q+1)*3-1,iSer)];
end
DateQQ = Dates(iS+(fcstH_Q+1)*3-1:3:iE-3);
DateQQ_V = DatesV(iS+(fcstH_Q+1)*3-1:3:iE-3,:);

% RMSFE
nF = size(FcstARQ,2);
rmsfeRW1  = sqrt(mean((FcstRW1Q - repmat(TrueQQ(1:end,1),1,nF)).^2)');
rmsfeRW2  = sqrt(mean((FcstRW2Q - repmat(TrueQQ(1:end,1),1,nF)).^2)');
rmsfeAR  = sqrt(mean((FcstARQ - repmat(TrueQQ(1:end,1),1,nF)).^2)');

%--------------------------------------------------------------------------
% saving the results
%--------------------------------------------------------------------------
eval(['save resBench_',strrep(P.SerFcst,' ','_')])
    



