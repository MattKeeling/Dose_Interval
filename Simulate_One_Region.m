function [nDC, nDCH, nDCD, nHospital, inHospital, nICU, inICU, nDeaths, PD_Lockdown, INF, FinalState]=Simulate_One_Region(Region, TAU, ALPHA, INC_P, S_Scale, Factor, h_factor, i_factor, d_factor, h_stretch, i_stretch, Lag, ...
    Start_Date,  WALES_FLAG, ComplianceT, ComplianceO, Run_stop, nV_Beta, nV_Speed, Import_Day, Import_Level, V1, V2, V3, Transmission_Reduction, VEffI, VEffS, VEffH, VEffD, RatioPf, Detection, Susceptibility, ...
    gamma, WaningSpeed, WaningEfficacy, Seasonality, MaxType, FS, W_start)


GETTING_OLDER=1;

Run_ODEs=@LeakyVacc_ODEs; % this function could be compiled to mex code to make everything run faster.
%ODEs=@LeakyVacc_ODEs_mex; % something like this

ODE_to_Obs=@ODE_to_Observables; % this function could be compiled to mex code to make everything run faster.
%ODE_to_Obs=@ODE_to_Observables_mex;  % something like this

m=3; L=21;

Christmas=0;

if Run_stop(1)>1000
    Run_stop=Run_stop+1-datenum(2020,1,1);
end

if length(Lag)==1
    Lag=Lag*ones(MaxType,1);
end

if length(h_factor)<MaxType,  h_factor(end:MaxType)=h_factor(end); end
if length(i_factor)<MaxType,  i_factor(end:MaxType)=i_factor(end); end
if length(d_factor)<MaxType,  d_factor(end:MaxType)=d_factor(end); end
if length(h_stretch)<MaxType,  h_stretch(end:MaxType)=h_stretch(end); end
if length(i_stretch)<MaxType,  i_stretch(end:MaxType)=i_stretch(end); end

if length(nV_Beta)<MaxType-1, nV_Beta(end:(MaxType-1))=nV_Beta(end); end
if length(nV_Speed)<MaxType-1, nV_Speed(end:(MaxType-1))=nV_Speed(end); end
if length(Import_Day)<MaxType-1, Import_Day(end:(MaxType-1))=10e3; end
if length(Import_Level)<MaxType-1, Import_Level(end:(MaxType-1))=0; end

load UK_SetUp_Data

UK_PP=UK_PP'*sum(Region_PP(Region,:))/sum(UK_PP);

UK_from_toH = UK_from_toH .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toW = UK_from_toW .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toS = UK_from_toS .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));
UK_from_toO = UK_from_toO .* ( ones(21,1) * (Region_PP(Region,:)./UK_PP));


Names={'ENGLAND','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland'};

tau=TAU;

Early_R0=3.8;
Early_gamma=gamma/1.5;
Z=1.3;

Susceptibility2=S_Scale*Susceptibility*Sus_Scale(Region);
Susceptibility3=S_Scale*Susceptibility;

if length(Seasonality)==2
    TT=Seasonality(2);
else
    TT=round(datenum(now)+1-datenum(2020,1,1));
end
T=[1:max(Run_stop+7)];
Extra_Seasonality=(1-0.5*Seasonality(1)-0.5*Seasonality(1)*cos((T-593)*2*pi/365)); % RAW
Extra_Seasonality(1:TT)=1;  % IGNORE EARLY STUFF.

if isempty(FS)
    W_start=1;
    
    [T2,S,E,D,nD2,U , Ratio, Da, FinalState]=Run_ODEs(MaxType, UK_from_toH*0 , gamma*(UK_from_toH + UK_from_toW + UK_from_toS + UK_from_toO), INC_P, Early_gamma, Susceptibility, Detection, TAU, nV_Beta, nV_Speed, 0, Region_PP(Region,:)' , V1(1,:), V2(1,:), V3(1,:), Transmission_Reduction, VEffI, VEffS, RatioPf, WaningSpeed, WaningEfficacy, 2, -1);
    
    
    %% EARLY SET-UP
    %Run without controls to 12th March (day 71)
    FinalState(22:end)=FinalState(22:end)*Factor*(1e-6*sum(Region_PP(Region,:)))*INITIAL_IMPORT_SIZE(Region)/sum(FinalState(22:end));
    
    [T,S,E,D,nD,U , Ratio, Da, FinalState]=Run_ODEs(MaxType, gamma*UK_from_toH*Z , gamma*(UK_from_toW + UK_from_toS + UK_from_toO), INC_P, Early_gamma, Susceptibility2, Detection, TAU, nV_Beta, nV_Speed, 0, Region_PP(Region,:)' , V1(1,:), V2(1,:), V3(1,:), Transmission_Reduction, VEffI, VEffS, RatioPf, WaningSpeed, WaningEfficacy, 71-Start_Date, [FinalState]);
    T=T+Start_Date;
    
    %Self Isolation for 4 days
    [t,s,e,d,nd,u , ratio, Da, FinalState]=Run_ODEs(MaxType, gamma*UK_from_toH*Z, gamma*(UK_from_toW + UK_from_toS + UK_from_toO), INC_P, gamma, Susceptibility3, Detection, TAU, nV_Beta, nV_Speed, 0, Region_PP(Region,:)' , V1(1,:), V2(1,:), V3(1,:), Transmission_Reduction, VEffI, VEffS, RatioPf, WaningSpeed, WaningEfficacy, 4, [FinalState]);
    T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];
    
    %Work from Home + Self Isolation for 4 days
    [new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0, 0.3, 0, ComplianceT(:,1), UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
    [t,s,e,d,nd,u , ratio, Da, FinalState]=Run_ODEs(MaxType, gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), INC_P, gamma, Susceptibility3, Detection, TAU, nV_Beta, nV_Speed, 0, Region_PP(Region,:)' , V1(1,:), V2(1,:), V3(1,:), Transmission_Reduction, VEffI, VEffS, RatioPf, WaningSpeed, WaningEfficacy, 4, [FinalState]);
    T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];
    
    %Some x (Social distancing + HHQ + School Closures + Work from Home + Self Isolatio)n for 3 days
    [new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0.8, 0.5, 0.3, ComplianceT(:,1), UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
    [t,s,e,d,nd,u , ratio, Da, FinalState]=Run_ODEs(MaxType, gamma*new_UK_from_toH*Z, gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), INC_P, gamma, Susceptibility3, Detection, tau, nV_Beta, nV_Speed, 0.3*ComplianceT(1,1), Region_PP(Region,:)', V1(1,:), V2(1,:), V3(1,:), Transmission_Reduction, VEffI, VEffS, RatioPf, WaningSpeed, WaningEfficacy, 3, [FinalState]);
    T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];
    
    ShortFlag=0;
else
    FinalState=FS; T=Run_stop(W_start-1); ShortFlag=1;
end
%% THEN RUN SENARIOS

PD_Lockdown = 0;

Flag=0; fixedaC=0;
UF=1; DF=1;

Told=0;

Multi_Factor=1.0;

for W=W_start:length(Run_stop)
    
    m=3;
    mmC=max(max(abs(ComplianceT)));
    if mmC>0
        if length(ComplianceT(:,W))==1
            aC=ones(21,1)*abs(ComplianceT(:,W));
            bC=ones(21,1)*abs(ComplianceO(:,W));
        end
        if length(ComplianceT(:,W))==4
            aC=abs(ComplianceT([1 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4],W));
            bC=abs(ComplianceO([1 2 2 2 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4],W));
        end
        if length(ComplianceT(:,W))==21
            aC=reshape(abs(ComplianceT(:,W)),21,1);
            bC=reshape(abs(ComplianceO(:,W)),21,1);
        end
        aC=aC+[0.5 0.3 0 0 0 0 0.3 0.5 0.6 0.7 0.8 0.9 1 1 1 1 1 1 1 1 1]'.*bC;
        aC(aC>1)=1;
        PD_Lockdown=PD_Lockdown+(Run_stop(W)-T(end)).*(Region_PP(Region,:)*aC)./mmC;
    end
    
    % Add in ALPHA
    if T(end)<Import_Day(1) & Run_stop(W)>=Import_Day(1)
        SP=(9*L+4*m*L)+L;  sL=(9*L+4*m*L)-L;    %sL excludes R class
        FinalState(SP+[1:sL])=Import_Level(1)*FinalState(L+[1:sL]);
        FinalState(L+[1:sL])=(1-Import_Level(1))*FinalState(L+[1:sL]);
        %fprintf(1,'Added in Alpha\n');
    end
    % Add in DELTA
    if T(end)<Import_Day(2) & Run_stop(W)>=Import_Day(2)
        SP=(9*L+4*m*L)+L; SP2=2*(9*L+4*m*L)+L; sL=(9*L+4*m*L)-L;
        FS=FinalState;
        tmp=FinalState(L+[1:sL])+FinalState(SP+[1:sL]);
        if Import_Level(2)>0
            FinalState(SP2+[1:sL])=Import_Level(2)*tmp;
            FinalState(L+[1:sL])=0*FS(L+[1:sL]);  % Get Rid of Wildtypes
            FinalState(SP+[1:sL])=(1-Import_Level(2))*FS(L+[1:sL])+(1-Import_Level(2))*FS(SP+[1:sL]);
        end
    end
    % Add in OMICRON
    if T(end)<Import_Day(3) & Run_stop(W)>=Import_Day(3)
        SP=(9*L+4*m*L)+L; SP2=2*(9*L+4*m*L)+L; SP3=3*(9*L+4*m*L)+L; sL=(9*L+4*m*L)-L;
        tmp=FinalState(SP+[1:sL])+FinalState(SP2+[1:sL]);
        if Import_Level(3)>0
            FinalState(SP3+[1:sL])=Import_Level(3)*tmp;
            FinalState(SP+[1:sL])=0*FinalState(SP+[1:sL]);  % Get Rid of Alpha
            FinalState(SP2+[1:sL])=(1-Import_Level(3))*FinalState(SP2+[1:sL]) + (1-Import_Level(3))*FinalState(SP+[1:sL]);
        else
            FinalState(SP3+[1:sL])=(-Import_Level(3))*tmp/sum(tmp);             % Add in a fixed number
            FinalState(SP+[1:sL])=0*FinalState(SP+[1:sL]);  % Get Rid of Alpha
            FinalState(SP2+[1:sL])=FinalState(SP2+[1:sL]) + FinalState(SP+[1:sL]) - (-Import_Level(3))*tmp/sum(tmp); % Just to keep the population size fixed.
        end
    end
    % Add in BA2
    if T(end)<Import_Day(4) & Run_stop(W)>=Import_Day(4)
        SP=(9*L+4*m*L)+L; SP2=2*(9*L+4*m*L)+L; SP3=3*(9*L+4*m*L)+L; SP4=4*(9*L+4*m*L)+L; sL=(9*L+4*m*L)-L;
        tmp=FinalState(SP2+[1:sL])+FinalState(SP3+[1:sL]);
        if Import_Level(4)>0
            FinalState(SP4+[1:sL])=Import_Level(4)*tmp;
            FinalState(SP2+[1:sL])=0*FinalState(SP2+[1:sL]);  % Get Rid of Delta
            FinalState(SP3+[1:sL])=(1-Import_Level(4))*FinalState(SP3+[1:sL]) + (1-Import_Level(4))*FinalState(SP2+[1:sL]);
        else
            FinalState(SP4+[1:sL])=(-Import_Level(4))*tmp/sum(tmp);             % Add in a fixed number
            FinalState(SP3+[1:sL])=FinalState(SP3+[1:sL]) + FinalState(SP2+[1:sL]) - (-Import_Level(4))*tmp/sum(tmp); % Just to keep the population size fixed.
            FinalState(SP2+[1:sL])=0*FinalState(SP2+[1:sL]);  % Get Rid of Delta
        end
    end
    
    InSchoolFlag=0;
    if Region==10
        if T(end)>=230  % ie Mid-August put schools back in Scotland
            InSchoolFlag=1;
        end
        if T(end)>=292 && T(end)<=294
            InSchoolFlag=0;  % Unless its Half Term
        end
        if (T(end)>=355 && T(end)<=365)
            InSchoolFlag=0;   % ... or Christmas
        end
    else
        if T(end)>=244 % ie 1st Sept put schools back
            InSchoolFlag=1;
        end
        if T(end)>=299 && T(end)<=301
            InSchoolFlag=0;   % Unless its Half Term
        end
        if (T(end)>=355 && T(end)<=365)
            InSchoolFlag=0;   % ... or Christmas
        end
    end
    
    if Region==9 && T(end)>=347 && T(end)<=365
        InSchoolFlag=0;  % Weeks break in Wales
    end
    
    if Region<9 || Region==11   % In England
        if (T(end)>=355 && T(end)<=433-2)
            InSchoolFlag=0;   % ... EXTENDED CHRISTMAS BREAK & LOCKDOWN
        end
    else
        if (T(end)>=355 && T(end)<=440-2)
            InSchoolFlag=0;   % ... EXTENDED CHRISTMAS BREAK & LOCKDOWN
        end
    end
    
    if T(end)>=456 && T(end)<470
        InSchoolFlag=0; % Easter '21!
    end
    
    if T(end)>=datenum(1,7,23)-2 && T(end)<datenum(1,9,3)-2 && Region<9
        InSchoolFlag=0; % Summer holidays '21 England!
    end
    if T(end)>=datenum(1,6,25)-2 && T(end)<datenum(1,8,5)-2 && Region==10
        InSchoolFlag=0; % Summer holidays '21 England!
    end
    if T(end)>=datenum(1,7,20)-2 && T(end)<datenum(1,8,31)-2 && Region==11
        InSchoolFlag=0; % Summer holidays '21 England!
    end
    if T(end)>=datenum(1,7,1)-2 && T(end)<datenum(1,8,12)-2 && Region==12
        InSchoolFlag=0; % Summer holidays '21 England!
    end
    
    if T(end)>=datenum(1,12,17)-2 && T(end)<datenum(2,1,4)-2
        InSchoolFlag=0; % Xmax holidays '21!
    end
    
    if T(end)>=datenum(2,4,8)-2 && T(end)<datenum(2,4,22)-2
        InSchoolFlag=0; % Easter holidays '22!
    end
    
    if T(end)>=datenum(2,7,25)-2 && T(end)<datenum(2,9,3)-2
        InSchoolFlag=0; % Summer holidays '22!
    end
    
    
    y=floor(T(end)/365);
    if GETTING_OLDER
        if Told<=datenum(y,9,4) & T(end)>datenum(y,9,4) % Move up year !!
            FS=FinalState;
            %fprintf('Moving on up! %g\n',T(end))
            SPV=MaxType*(9*L+4*m*L)+L;
            for A=20:-1:2
                FinalState(A:21:end)=0.8*FS(A:21:end)+0.2*FS((A-1):21:end)*Region_PP(Region,A)./Region_PP(Region,A-1);
            end
            A=21;
            FinalState(A:21:end)=(FS(A:21:end)/1.2)+0.2*FS((A-1):21:end)*Region_PP(Region,A)./(1.2*Region_PP(Region,A-1));
            tmpFS=FinalState;
            
            % and reset to get the right population size (but not infection).
            SPV=MaxType*(9*L+4*m*L)+L;
            FS=FinalState;
            for A=1:21
                PopSize=sum(FS(A:L:end))-sum(FS(SPV+A+[0 1]*L))-sum(FS(end+[-20:0]));  %substract numV1, numV2 and Boosters
                FinalState(A:21:end)=FinalState(A:21:end)*Region_PP(Region,A)/PopSize;
            end
        end
        
        % Add in births
        A=1; FS=FinalState;
        B=(Run_stop(W)-T(end))/(5*365);
        FinalState(A:21:end)=FS(A:21:end)*(1-B);
        FinalState(A)=FinalState(A)+B*Region_PP(Region,1);
        %FinalState((A+21):21:end)=FS((A+21):21:end)*(1-B);
    end
    
    if Told<=datenum(y,8,4) & T(end)>datenum(y,8,4)  % Reset booster status every August
        FinalState(end+[-20:0])=0;  %Boosters should be final state of system.
    end
    if Told<=807 & T(end)>807 % Reset booster for Spring Booster
        FinalState(end+[-20:0])=0;  %Boosters should be final state of system.
    end
    if mean(V3((T(end)+1):Run_stop(W),:),'all')==0 % Reset booster status if number of boosters is zero.
        FinalState(end+[-20:0])=0;  %Boosters should be final state of system.
    end
    
    Extra_Beta = mean(Extra_Seasonality(T(end)+1:Run_stop(W)));
    
    Told=T(end);
    
    %Make sure V's don't over run
    SPV=MaxType*(9*L+4*m*L)+L;
    FS_V1=FinalState(SPV+[1:L]);
    FS_V2=FinalState(SPV+L+[1:L]);
    FS_S=FinalState([1:21]);
    FS_B=FinalState(end+[-20:0]);
    mV1=mean(V1((T(end)+1):Run_stop(W),:),1);  MV1=mV1;
    mV2=mean(V2((T(end)+1):Run_stop(W),:),1);  MV2=mV2;
    mV3=mean(V3((T(end)+1):Run_stop(W),:),1);  MV3=mV3;
    mV1=min(mV1, 0.99*(Region_PP(Region,:)-FS_V1-FS_V2)/(Run_stop(W)-T(end)));  mV1(mV1<0)=0;
    mV2=min(mV2, mV1+FS_V1*0.99/(Run_stop(W)-T(end)));      % (Current_V1 / Time) + new_V1
    mV2(mV2<0)=0;
    mV3=min(mV3, 0.99*(Region_PP(Region,:)-FS_S-FS_B)/(Run_stop(W)-T(end))); mV3(mV3<0)=0;
    
    oldaC=aC;
    
    if InSchoolFlag
        if T(end)>datenum(1,8,1)
            aC(2:4)=(0.4-0.003*(T(end)-datenum(1,8,1)))*aC(2:4);
        else
            aC(2:4)=0.5*aC(2:4);
        end
    end
    
    [new_UK_from_toH, new_UK_from_toW, new_UK_from_toS, new_UK_from_toO] = Return_New_Matrices(0.95, 0.8, 0.95, aC, UK_from_toH, UK_from_toW, UK_from_toS, UK_from_toO);
    
    [t,s,e,d,nd,u , ratio, Da, FinalState]=Run_ODEs(MaxType, Extra_Beta*gamma*new_UK_from_toH*Z, Extra_Beta*gamma*(new_UK_from_toW + new_UK_from_toS + new_UK_from_toO), INC_P, gamma, ...
        Susceptibility3, Detection, tau, nV_Beta, nV_Speed, 0.8*abs(ComplianceT(:,W)), Region_PP(Region,:)', mV1, ...
        mV2, mV3, Transmission_Reduction, VEffI, VEffS, RatioPf, WaningSpeed, WaningEfficacy, Run_stop(W)-T(end), [FinalState]);
    
    
    if W==W_start & ShortFlag
        T=t+Run_stop(W_start-1); S=s; E=e; D=d; nD=nd; U=u; Ratio=ratio;
    else
        T=[T; T(end)+t(2:end)]; S=[S; s(2:end,:)]; E=[E; e(2:end,:)]; D=[D; d(2:end,:)]; nD=[nD; nd(2:end,:)]; U=[U; u(2:end,:)]; Ratio=[Ratio; ratio(2:end,:,:)];
    end
    
end

nD(nD<0)=0;


% Detectable Cases accounting for vaccine efficacy.

load Distrib_n_Probs.mat

RatioAZ=1-RatioPf;

for Type=MaxType:-1:1
    a=1:L;
    clear nDC Ratios;
    nDC(T+Lag(Type),:)=nD; nE=E;
    Ratios(T+Lag(Type),:,:)=Ratio; DRatios=Ratios; DRatios(DRatios==0)=1;
    
    RveffH1=(1-VEffH(1,Type,1)*RatioPf(a)-VEffH(2,Type,1)*RatioAZ(a))./(1-VEffS(1,Type,1)*RatioPf(a)-VEffS(2,Type,1)*RatioAZ(a));
    RveffH2=(1-VEffH(1,Type,2)*RatioPf(a)-VEffH(2,Type,2)*RatioAZ(a))./(1-VEffS(1,Type,2)*RatioPf(a)-VEffS(2,Type,2)*RatioAZ(a));
    RveffH3=(1-VEffH(1,Type,3)*RatioPf(a)-VEffH(2,Type,3)*RatioAZ(a))./(1-VEffS(1,Type,3)*RatioPf(a)-VEffS(2,Type,3)*RatioAZ(a));
    RveffH4=(1-VEffH(1,Type,4));
    
    RveffD1=(1-VEffD(1,Type,1)*RatioPf(a)-VEffD(2,Type,1)*RatioAZ(a))./(1-VEffS(1,Type,1)*RatioPf(a)-VEffS(2,Type,1)*RatioAZ(a));
    RveffD2=(1-VEffD(1,Type,2)*RatioPf(a)-VEffD(2,Type,2)*RatioAZ(a))./(1-VEffS(1,Type,2)*RatioPf(a)-VEffS(2,Type,2)*RatioAZ(a));
    RveffD3=(1-VEffD(1,Type,3)*RatioPf(a)-VEffD(2,Type,3)*RatioAZ(a))./(1-VEffS(1,Type,3)*RatioPf(a)-VEffS(2,Type,3)*RatioAZ(a));
    RveffD4=(1-VEffD(1,Type,4));
    
    for a=1:L
        Scale(a,1:7)=[1 RveffH1(a) RveffH2(a) RveffH3(a) WaningEfficacy(Type,4) WaningEfficacy(Type,3) RveffH4];
        Z=Ratios(:,a,1)*Scale(a,1);
        for i=2:7
            Z=Z+Ratios(:,a,i)*Scale(a,i);
        end
        Z=Z./sum(Ratios(:,a,:),3);
        nDCH(1:size(Z,1),a+L*(Type-1))=nDC(:,a+L*(Type-1)).*Z;
        
        Scale(a,1:7)=[1 RveffD1(a) RveffD2(a) RveffD3(a) WaningEfficacy(Type,6) WaningEfficacy(Type,5) RveffD4];
        Z=Ratios(:,a,1)*Scale(a,1);
        for i=2:7
            Z=Z+Ratios(:,a,i)*Scale(a,i);
        end
        Z=Z./sum(Ratios(:,a,:),3);
        nDCD(1:size(Z,1),a+L*(Type-1))=nDC(:,a+L*(Type-1)).*Z;
    end
end

nDCH(isnan(nDCH))=0; nDCD(isnan(nDCD))=0;

aa=1:L;

[NDC(:,aa), nHospital, inHospital, nICU, inICU, nDeaths, INF]=ODE_to_Obs(T, nDCH(:,aa), nE(:,aa), Region, h_factor(1),i_factor(1),d_factor(1),h_stretch(1),i_stretch(1),INC_P,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);
[~, ~, ~, nICU, inICU, nDeaths, ~]=ODE_to_Obs(T, nDCD(:,aa), nE(:,aa), Region, h_factor(1),i_factor(1),d_factor(1),h_stretch(1),i_stretch(1),INC_P,...
    Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp,  WALES_FLAG);

for TYPE=2:MaxType
    LL=(TYPE-1)*L;
    S2H=Sympt_2_hosp;
    if TYPE>=3  S2H=Sympt_2_hosp_Delta; end
    [~, nHospital(:,aa+LL), inHospital(:,aa+LL), nICU(:,aa+LL), inICU(:,aa+LL), nDeaths(:,aa+LL), INF(:,aa+LL)]=ODE_to_Obs(T, nDCH(:,aa+LL), nE(:,aa+LL), Region, h_factor(TYPE),i_factor(TYPE),d_factor(TYPE),h_stretch(TYPE),i_stretch(TYPE),INC_P*nV_Speed(TYPE-1),...
        Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
        Hosp_2_Death, Sympt_2_critcal, S2H,  WALES_FLAG);
    [~, ~, ~, nICU(:,aa+LL), inICU(:,aa+LL), nDeaths(:,aa+LL), ~]=ODE_to_Obs(T, nDCD(:,aa+LL), nE(:,aa+LL), Region, h_factor(TYPE),i_factor(TYPE),d_factor(TYPE),h_stretch(TYPE),i_stretch(TYPE),INC_P*nV_Speed(TYPE-1),...
        Assumed_Delay_Reporting_Deaths, Distribution_Hosp_Time', Distribution_HospICU_Time', Distribution_ICU_Time', Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death,...
        Hosp_2_Death, Sympt_2_critcal, S2H,  WALES_FLAG);
end


end

%%


