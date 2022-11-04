function [NDC, nHospital, inHospital, nICU, inICU, nDeaths, INF]=ODE_to_Observables(T, nDC, E, Region, h_factor,i_factor,d_factor,h_stretch,i_stretch,a,...
    Assumed_Delay_Reporting_Deaths, rc_Distribution_Hosp_Time, rc_Distribution_HospICU_Time, rc_Distribution_ICU_Time, Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death, ...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp, Wales_Flag)


nA=size(nDC,2);

Age_Factor=[2.595 2.595 2.595 2.595 2.595 2.595 1.9752 1.9752 1.2771 ...
    1.2771 1.1319 1.1319 0.9547 0.9547 0.7959 0.7959 0.8854 0.8854 ...
    0.8854 0.8854 0.8854];

%new Hospital
nHospital=zeros(max(T)+100,nA);
A=1:nA;
LL=[1:length(Distribution_Symptoms_to_Hospital)];
for t=1:length(T)
    nHospital(T(t)+LL,A)=nHospital(T(t)+LL,A) + Distribution_Symptoms_to_Hospital'*(nDC(T(t),A).*Prob_Modifier_linear(Sympt_2_hosp(A),h_factor)');
end

%new ICU
nICU=zeros(max(T)+100,nA);
LL=[1:length(Distribution_Symptoms_to_ICU)];
for t=1:length(T)
    nICU(T(t)+LL,A)=nICU(T(t)+LL,A) + Distribution_Symptoms_to_ICU' * (nDC(T(t),A).*Prob_Modifier_linear(Sympt_2_critcal(A),i_factor)');
end

%in Hosptial
inHospital=zeros(max(T)+ceil(128*h_stretch),nA);

for A=1:nA
    
    l=length(rc_Distribution_Hosp_Time)-1;
    tmpDist_Hosp_Time=interp1([0:l],rc_Distribution_Hosp_Time,[0:(Age_Factor(A)/h_stretch):l],'linear')';
    l=length(tmpDist_Hosp_Time);
    for d=1:l
        inHospital(T+d-1,A)=inHospital(T+d-1,A)+tmpDist_Hosp_Time(d)*nHospital(T,A);
    end
    
    
    l=length(rc_Distribution_HospICU_Time)-1;
    tmpDist_HospICU_Time=interp1([0:l],rc_Distribution_HospICU_Time,[0:(Age_Factor(A)/h_stretch):l],'linear')';
    l=length(tmpDist_HospICU_Time);
    for d=1:l
        inHospital(T+d-1,A)=inHospital(T+d-1,A)+tmpDist_HospICU_Time(d)*nICU(T,A);
    end
end
nHospital((max(T)+1):end,:)=[];
inHospital((max(T)+1):end,:)=[];

A=1:nA;
%in ICU
inICU=zeros(max(T)+ceil(100*i_stretch),nA);
l=length(rc_Distribution_ICU_Time)-1;
tmpDist_ICU_Time=interp1([0:l],rc_Distribution_ICU_Time,[0:(1/i_stretch):l],'linear')';
l=length(tmpDist_ICU_Time);
for t=1:length(T)
    inICU(T(t)+[1:l]-1,A)=inICU(T(t)+[1:l]-1,A)+tmpDist_ICU_Time*nICU(T(t),A);
end
nICU((max(T)+1):end,:)=[];
inICU((max(T)+1):end,:)=[];

if Region==9 & Wales_Flag % if Wales_Flag set !!  DO IT ALL AGAINST BUT WITH 14 day cut-off to match data feeds
    inHospital2=zeros(max(T)+ceil(100*h_stretch),nA);
    l=length(rc_Distribution_Hosp_Time)-1;
    tmpDist_Hosp_Time=interp1([0:l],rc_Distribution_Hosp_Time,[0:(1/h_stretch):l],'linear')';
    l=length(rc_Distribution_HospICU_Time)-1;
    tmpDist_HospICU_Time=interp1([0:l],rc_Distribution_HospICU_Time,[0:(1/h_stretch):l],'linear')';
    tmpDist_Hosp_Time(15:end)=0;  tmpDist_HospICU_Time(15:end)=0;
    
    l=length(tmpDist_Hosp_Time);
    for t=1:length(T)
        inHospital2(T(t)+[1:l]-1,A)=inHospital2(T(t)+[1:l]-1,A)+tmpDist_Hosp_Time*nHospital(T(t),A);
    end
    l=length(tmpDist_HospICU_Time);
    for t=1:length(T)
        inHospital2(T(t)+[1:l]-1,A)=inHospital2(T(t)+[1:l]-1,A)+tmpDist_HospICU_Time*nICU(T(t),A);
    end
    inHospital2((max(T)+1):end,:)=[];
    MM=min((max(T)+1),147);
    inHospital(1:MM,:)=inHospital2(1:MM,:);
    
    inICU2=zeros(max(T)+ceil(100*i_stretch),nA);
    l=length(rc_Distribution_ICU_Time)-1;
    tmpDist_ICU_Time=interp1([0:l],rc_Distribution_ICU_Time,[0:(1/i_stretch):l],'linear')';
    tmpDist_ICU_Time(15:end)=0;
    
    l=length(tmpDist_ICU_Time);
    for t=1:length(T)
        inICU2(T(t)+[1:l]-1,A)=inICU2(T(t)+[1:l]-1,A)+tmpDist_ICU_Time*nICU(T(t),A);
    end
    inICU2((max(T)+1):end,:)=[];
    inICU(1:MM,:)=inICU2(1:MM,:);
end

% find Deaths
nDeaths=zeros(max(T)+150,nA);
LL=[1:length(Distribution_Hopital_to_Death)];
for t=1:length(T)
    nDeaths(T(t)+LL+Assumed_Delay_Reporting_Deaths(Region),A)=nDeaths(T(t)+LL+Assumed_Delay_Reporting_Deaths(Region),A) + Distribution_Hopital_to_Death'*(nHospital(T(t),A).*Prob_Modifier_linear(Hosp_2_Death(A),d_factor)');
end
NDC=nDC(1:T(end),:);
nDeaths((max(T)+1):end,:)=[];

INF=zeros(T(end),21);
if ~isempty(E)
    INF(T,:)=E*a;
end
end

function [Y] = Prob_Modifier_linear(y,x)

Y=y.*x;

end
