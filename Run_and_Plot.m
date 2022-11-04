% Script to generate outputs for the study 
% "The Impact of SARS-CoV-2 Vaccine Dose Separation and Dose Targeting on 
% Hospital Admissions and Deaths from COVID-19 in England"

%% Run the different vaccine prioritisation scenarios
% Input parameter value to Generate_Output corresponds to chosen scenario
% Outputs from function call are saved to MAT files
Generate_Output(0); % Default vaccination priorisation order, 12 week dose interval
Generate_Output(1); % Default vaccination priorisation order, 3 week dose interval
Generate_Output(2); % Default vaccination priorisation order, 3 week dose interval & lower vaccine efficacy
Generate_Output(3); % Reversed vaccination prioritisation order, 12 week dose interval

%% Load required output data
clear

% Specify filename suffix
DATE = '02_11_22';

% Data load: Default vaccination priorisation order, 12 week dose interval
load(['Warwick_Output_Loop0_' DATE]);
Hosp0=sum(nHOSP_AD(2:8,:,:),[1 3]); Death0=sum(nDEATHS(2:8,:,:),[1 3]);

% Data load: Default vaccination priorisation order, 3 week dose interval
load(['Warwick_Output_Loop1_' DATE]);
Hosp1=sum(nHOSP_AD(2:8,:,:),[1 3]); Death1=sum(nDEATHS(2:8,:,:),[1 3]);

% Data load: Default vaccination priorisation order, 3 week dose interval & lower vaccine efficacy
load(['Warwick_Output_Loop2_' DATE]);
Hosp2=sum(nHOSP_AD(2:8,:,:),[1 3]); Death2=sum(nDEATHS(2:8,:,:),[1 3]);

% Data load: Reversed vaccination prioritisation order, 12 week dose interval
load(['Warwick_Output_Loop3_' DATE]);
Hosp3=sum(nHOSP_AD(2:8,:,:),[1 3]); Death3=sum(nDEATHS(2:8,:,:),[1 3]);


%% Generate the figure

Xrange=[260 640]; % when boosters start.

figure(1); clf;
X=datenum(2020,1:48,1); XT=X+1-datenum(2020,1,1); XTL=datestr(X,'m'); XTL([3 13 25],1:4)=datestr(X([3 13 25]),'m''yy');

% Plot daily deaths
subplot(2,2,1);
h=plot(1:length(Death0),Death0,'-k',1:length(Death1),Death1,'-r',1:length(Death2),Death2,'--m',1:length(Death3),Death3,':g');
h(4).Color=[0 0.8 0]; set(h,'LineWidth',2);
set(gca,'XLim',Xrange,'XTick',XT,'XTickLabel',XTL,'FontSize',12);
ylabel('Daily Deaths','FontSize',14);

% Plot hospital admissions
subplot(2,2,2);
h=plot(1:length(Hosp0),Hosp0,'-k',1:length(Hosp1),Hosp1,'-r',1:length(Hosp2),Hosp2,'--m',1:length(Hosp3),Hosp3,':g');
h(4).Color=[0 0.8 0];  set(h,'LineWidth',2);
set(gca,'XLim',Xrange,'XTick',XT,'XTickLabel',XTL,'FontSize',12);
ylabel('Hospital Admissions','FontSize',14);

% Plot additional cumulative deaths
subplot(2,2,3);
h=plot(1:length(Death1),cumsum(Death1-Death0)/1e3,'-r',1:length(Death2),cumsum(Death2-Death0)/1e3,'--m',1:length(Death3),cumsum(Death3-Death0)/1e3,':g');
h(3).Color=[0 0.8 0];  set(h,'LineWidth',2);
set(gca,'XLim',Xrange,'XTick',XT,'XTickLabel',XTL,'FontSize',12);
ylabel({'Additional Cumulative Deaths','(thousands)'},'FontSize',14);

% Plot additional cumulative hospital admissions
subplot(2,2,4);
h=plot(1:length(Hosp1),cumsum(Hosp1-Hosp0)/1e3,'-r',1:length(Hosp2),cumsum(Hosp2-Hosp0)/1e3,'--m',1:length(Hosp3),cumsum(Hosp3-Hosp0)/1e3,':g');
h(3).Color=[0 0.8 0];  set(h,'LineWidth',2);
set(gca,'XLim',Xrange,'XTick',XT,'XTickLabel',XTL,'FontSize',12);
ylabel({'Additional Cumulative Hospital','Admissions (thousands)'},'FontSize',14);

% Save combined figure to file
print -dpng Standard_Output.png
