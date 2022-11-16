% Function to apply observation model to tranmsission model outputs

function [NDC, nHospital, inHospital, nICU, inICU, nDeaths, INF]=ODE_to_Observables(T, nDC, E, Region, h_factor,i_factor,d_factor,h_stretch,i_stretch,a,...
    Assumed_Delay_Reporting_Deaths, rc_Distribution_Hosp_Time, rc_Distribution_HospICU_Time, rc_Distribution_ICU_Time, Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death, ...
    Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp, Wales_Flag)

% Inputs:
% T - Output evaluation times 
% nDC - New symptomatic cases array
% E - Latent infected counts
% Region - Region ID
% h_factor,i_factor - Modification factors for symptomatic case to hospitalisation/ICU admission probabilities
% d_factor - Modification factor for hospitalisation to death probabilities
% h_stretch,i_stretch - Modification factors for hospital & ICU length of stay distributions  
% a - Rate out of latent/exposed disease state
% Assumed_Delay_Reporting_Deaths - Region specific death reporting delay
% rc_Distribution_Hosp_Time, rc_Distribution_ICU_Time - Length of stay distributions for hospital admission & ICU admission
% rc_Distribution_HospICU_Time - Distribution of hospital to ICU time.
% Distribution_Symptoms_to_Hospital, Distribution_Symptoms_to_ICU, Distribution_Hopital_to_Death - Distirubtions for time to progress between specified severe health episode stages
% Hosp_2_Death, Sympt_2_critcal, Sympt_2_hosp - Age-group specific probabilities for specified severe health episode progression
% Wales_Flag - If active, runs observation loop with a 14 day cut-off to match data feed

% Outputs: 
% 2D arrays. Row per timestep, column per age group
%   NDC - Symptomatic cases
%   nHospital - Hospital admissions
%   inHospital - Hospital occupancy
%   nICU - ICU admissions
%   inICU - ICU occpancy
%   nDeaths - Deaths
%   INF - Infections

% Assign number of age groups
nA=size(nDC,2);

% Scaling factor for hospital length of stay stretch parameter (per age group)
Age_Factor=[2.595 2.595 2.595 2.595 2.595 2.595 1.9752 1.9752 1.2771 ...
    1.2771 1.1319 1.1319 0.9547 0.9547 0.7959 0.7959 0.8854 0.8854 ...
    0.8854 0.8854 0.8854];

%new hospital admissions
nHospital=zeros(max(T)+100,nA);
A=1:nA;
LL=[1:length(Distribution_Symptoms_to_Hospital)];
for t=1:length(T)
    nHospital(T(t)+LL,A)=nHospital(T(t)+LL,A) + Distribution_Symptoms_to_Hospital'*(nDC(T(t),A).*Prob_Modifier_linear(Sympt_2_hosp(A),h_factor)');
end

%new ICU admissions
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

%in ICU
A=1:nA;
inICU=zeros(max(T)+ceil(100*i_stretch),nA);
l=length(rc_Distribution_ICU_Time)-1;
tmpDist_ICU_Time=interp1([0:l],rc_Distribution_ICU_Time,[0:(1/i_stretch):l],'linear')';
l=length(tmpDist_ICU_Time);
for t=1:length(T)
    inICU(T(t)+[1:l]-1,A)=inICU(T(t)+[1:l]-1,A)+tmpDist_ICU_Time*nICU(T(t),A);
end
nICU((max(T)+1):end,:)=[];
inICU((max(T)+1):end,:)=[];

% Wales: If Wales_Flag set to 1 (true), run again with a 14 day cut-off to match data feed
if Region==9 & Wales_Flag 
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

% Deaths
nDeaths=zeros(max(T)+150,nA);
LL=[1:length(Distribution_Hopital_to_Death)];
for t=1:length(T)
    nDeaths(T(t)+LL+Assumed_Delay_Reporting_Deaths(Region),A)=nDeaths(T(t)+LL+Assumed_Delay_Reporting_Deaths(Region),A) + Distribution_Hopital_to_Death'*(nHospital(T(t),A).*Prob_Modifier_linear(Hosp_2_Death(A),d_factor)');
end
nDeaths((max(T)+1):end,:)=[];

% Infections
INF=zeros(T(end),21);
if ~isempty(E)
    INF(T,:)=E*a;
end

% Symptomatic cases
NDC=nDC(1:T(end),:);

end

%% Supporting functions
function [Y] = Prob_Modifier_linear(y,x)

Y=y.*x;

end
