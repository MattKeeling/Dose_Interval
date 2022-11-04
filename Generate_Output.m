function []=Generate_Output(LOOP)

%%

Matlab_Save=1;

Rmax=11;
Rtot=[2:11];

RUN_UNTIL=700;

DATE_STR='02_11_22';

MaxType=5; % maximum number of variants.


load(['Data_File_' DATE_STR '.mat']);
Date=datenum(DATE_STR)+1-datenum(2020,1,1);

load Surrogate_Vaccine_Data.mat

load Regional_PP.mat

% LOOP determines the vaccination strategy


%LOOP=0 - use the default vaccination.
%LOOP=1 - same list but at 3 weeks
%LOOP=2 - same list but at 3 weeks + lower vaccine efficacy
%LOOP=3 - reversed list at 12 weeks.
%LOOP=10 - random order of vaccination
if LOOP>0 && LOOP<=3
    VV1=0*V1; VV2=0*V2;
    for R=2:11
        LIST=zeros(ceil(sum(V1(R,:,:),'all')),1);
        k=0; N=zeros(21,1);
        for t=1:size(V1,2)
            for a=1:21
                q=round(sum(V1(R,1:t,a),'all')-N(a));
                if q>0
                    LIST(k+[1:q])=a;
                    k=k+q;
                    N(a)=N(a)+q;
                end
            end
        end

        if LOOP==3
            LIST=LIST(end:-1:1); LIST(LIST==0)=[]; %reverse the list
            DLAY=12*7; % 12 week delay
        else
            DLAY=3*7;  % 3 week delay
        end
        CumVacc=squeeze(cumsum(squeeze(sum(V1(R,:,:)+V2(R,:,:),3))));
        
        k=0;
        for t=(DLAY+1):size(V1,2)
            VV2(R,t,:)=VV1(R,t-DLAY,:);
            Spare=round(CumVacc(t) - sum( VV1(R,1:t,:) + VV2(R,1:t,:) , 'all'));
            if Spare>0
                VV1(R,t,:)=full(sparse(LIST(k+[1:Spare]),1,1,21,1));
                k=k+Spare;
            end
        end
    end
    V1=VV1; V2=VV2;
end

if LOOP==2 % lower the efficacy against alpha and delta.
    %alpha
    VEffI(1,2,1:2)=[63 77]/100; VEffI(2,2,1:2)=[63 65]/100;
    VEffS(1,2,1:2)=[63 82]/100; VEffS(2,2,1:2)=[63 66]/100;
    %delta
    VEffI(1,3,1:2)=[55 80]/100; VEffI(2,3,1:2)=[45 48]/100;
    VEffS(1,3,1:2)=[55 85]/100; VEffS(2,3,1:2)=[45 49]/100;
end

if LOOP==10
    AGE_ORDER=[21:-1:4];
    AGE_ORDER=AGE_ORDER(randperm(length(AGE_ORDER)));
    Num_Vacc=round(squeeze(sum(V2,2))); Num_Vacc(:,1:3)=0; Num_Vacc(:,4)=round(Num_Vacc(:,4)*2/5);
    Total_Daily_Vacc=round(squeeze(sum(V1+V2+V3,3)));
    V1=0*V1; V2=0*V2;
    
    for Region=2:11
        RNV=Num_Vacc(Region,:);
        for T=(Dlay+1):length(Total_Daily_Vacc)
            Tot=Total_Daily_Vacc(Region,T);
            Scnds=squeeze(sum(V1(Region,1:(T-Dlay),:),2)-sum(V2(Region,1:T,:),2));
            if Tot<sum(Scnds)  % not enough vaccine today
                tmp=zeros(sum(Scnds),1); k=0; for i=1:21, tmp(k+[1:Scnds(i)])=i; k=k+Scnds(i); end;
                ptmp=tmp(randperm(length(tmp)));
                ptmp((Tot+1):end)=[];
                V2(Region,T,:)=hist(ptmp,[1:21]);
            else
                V2(Region,T,:)=Scnds;
                Tot=Tot-sum(Scnds);
                o=1;
                while o<=length(AGE_ORDER) & Tot>0
                    A=AGE_ORDER(o);
                    if Tot>RNV(A)
                        V1(Region,T,A)=RNV(A); Tot=Tot-RNV(A); RNV(A)=0;
                    else
                        V1(Region,T,A)=Tot; RNV(A)=RNV(A)-Tot; Tot=0;
                    end
                    o=o+1;
                end
            end
        end
    end
end

%%
clear nALL nDEATHS nHOSP_AD nHOSP_OCC nICU_OCC nICU_AD

Names={'England','East of England','London','Midlands','North East and Yorkshire','North West','South East','South West','Wales','Scotland','Northern Ireland','United Kingdom'};

%%
% generate run_stops by region
Comps = zeros(4,Num_Comp+1,11);

RUN_STOP(RUN_STOP>RUN_UNTIL)=[];
Num_Comp=length(RUN_STOP);

RUN_STOPs = zeros(11,Num_Comp);
for Region = 2:11
    RUN_STOPs(Region,1:length(RUN_STOP)) = RUN_STOP;
end

Keep_RUN_STOPs=RUN_STOPs(2,:);

if size(V1,2)<max(RUN_STOPs,[],'all')
    mx=max(RUN_STOPs,[],'all')+100;
    for r=1:11  for a=1:21
            V1(r,end:mx,a)=V1(r,end,a);
            V2(r,end:mx,a)=V2(r,end,a);
        end
    end
end

%%

maxtime = max(RUN_STOPs(end,:))+30; % needs to be big enough to take account of any lags
nALL = zeros(11,maxtime,105);
nDEATHS = zeros(11,maxtime,105);
nHOSP_AD = zeros(11,maxtime,105);
nHOSP_OCC = zeros(11,maxtime,105);
nICU_AD = zeros(11,maxtime,105);
nICU_OCC = zeros(11,maxtime,105);

% generate Comp by region
clear Comps

RUN_STARTs=[82*ones(11,1) RUN_STOPs(:,1:(end-1))];

% ADD IN DECLINE IN PRECAUTIONARY BEHAVIOUR IF SIMULATION LENGTH
% EXCEEDS DATA.
NC=Num_Comp;
for Region = Rtot
    Comps(Region,1:NC)=nCOMPLIANCE(1:NC,Region);
    Comps(Region,NC:size(RUN_STOPs,2))=Comps(Region,NC);
    CompsO(Region,1:NC)=nCOMPLIANCEO(1:NC,Region);
    CompsO(Region,NC:size(RUN_STOPs,2))=CompsO(Region,NC);
    m=find(RUN_STOPs(Region,:)>Date);  Decline=0.025;
    if length(m)
        for i=1:length(m)
            Comps(Region,m(i))=max(Comps(Region,m(i)-1)-Decline,0.001);  %Drop compliance slowly each week until zero
        end
        CompsO(Region,m)=CompsO(Region,m(1)-1)*Comps(Region,m)/Comps(Region,m(1)-1);
    end
end

MAX_REGION=Rmax;

TXT_STR=['Warwick_Output_Loop' num2str(LOOP)];

%parfor Region=[2:8]
for Region=2:8
    
    fprintf(1,'Simulating %s, Vaccine Loop %d\n',Names{Region},LOOP);
    
    [nDC, ~, ~, nHospital, inHospital, nICU, inICU, nDeaths, ~, All, ~]=Simulate_One_Region(Region, nTAU, nALPHA, nINC_P, nSCALING(Region), nFACTOR(Region), nH_FACTOR(Region:11:end), nI_FACTOR(Region:11:end), nD_FACTOR(Region:11:end),  ...
        nH_STRETCH(Region:11:end), nI_STRETCH(Region:11:end), nLAG(Region:11:end), START_DATE(Region)+1, 1, Comps(Region,:), CompsO(Region,:), RUN_STOPs(Region,:), nNV_BETA(Region:11:end)', nNV_SPEED(Region:11:end)', nIMPORT_DAY(Region:11:end), nIMPORT_LEVEL(Region:11:end),...
        squeeze(V1(Region,:,:)), squeeze(V2(Region,:,:)), squeeze(V3(Region,:,:)), Transmission_Reduction, VEffI, VEffS, VEffH, VEffD, RatioPf,...
        Detection, Susceptibility, nGAMMA, WaningSpeed, WaningEfficacy, [0.1 0], MaxType, [], 0);

    nDeaths=nDeaths.* ((1+nDH_SCALING(Region)*sum(inHospital,2)./sum(Region_PP(Region,:),2))*ones(1,size(nDeaths,2)));
    
    T=1:size(All,1);
    
    padding = maxtime-length(T);
    nALL(Region,:,:)=[All; zeros(padding,105)];
    nDEATHS(Region,:,:)=[nDeaths; zeros(padding,105)];
    nHOSP_OCC(Region,:,:)=[inHospital; zeros(padding,105)];
    nHOSP_AD(Region,:,:) = [nHospital; zeros(padding,105)];
    nICU_OCC(Region,:,:) = [inICU; zeros(padding,105)];
    nICU_AD(Region,:,:) = [nICU; zeros(padding,105)];
end



%% NEED TO COMBINE THE OLD & NEW VARIANT !

A=1:size(nALL,3); a=1:21; b=21+a; c=42+a; d=63+a; e=84+a;
nALL=nALL(:,:,a)+nALL(:,:,b)+nALL(:,:,c)+nALL(:,:,d)+nALL(:,:,e);
nDEATHS=nDEATHS(:,:,a)+nDEATHS(:,:,b)+nDEATHS(:,:,c)+nDEATHS(:,:,d)+nDEATHS(:,:,e);
nHOSP_OCC=nHOSP_OCC(:,:,a)+nHOSP_OCC(:,:,b)+nHOSP_OCC(:,:,c)+nHOSP_OCC(:,:,d)+nHOSP_OCC(:,:,e);
nHOSP_AD=nHOSP_AD(:,:,a)+nHOSP_AD(:,:,b)+nHOSP_AD(:,:,c)+nHOSP_AD(:,:,d)+nHOSP_AD(:,:,e);
nICU_OCC=nICU_OCC(:,:,a)+nICU_OCC(:,:,b)+nICU_OCC(:,:,c)+nICU_OCC(:,:,d)+nICU_OCC(:,:,e);
nICU_AD=nICU_AD(:,:,a)+nICU_AD(:,:,b)+nICU_AD(:,:,c)+nICU_AD(:,:,d)+nICU_AD(:,:,e);


if Matlab_Save
    save([TXT_STR '_' DATE_STR '.mat'],'nALL','nDEATHS','nHOSP_OCC','nHOSP_AD','nICU_OCC','nICU_AD');
end


end

