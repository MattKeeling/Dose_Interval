function [T,S,E,D,nD,U,RatioS,Age_structure, FinalState] = LeakyVacc_ODEs(MaxType, M_from_to_H, M_from_to_O, alpha, gamma, sigma, d, tau, nV_Beta, nV_Speed, HHQ, N0, ...
    V1, V2, V3, tmpTransmission_Reduction, tmpVEffI, tmpVEffS, RatioPf, WaningSpeed, WaningEfficacy, MaxTime, InitialState) %#codegen
%%
%

number_E_states=3;

L=length(N0);

if length(sigma)==1
    sigma=sigma+0*N0;
else
    sigma=reshape(sigma,L,1);
end
if length(d)==1
    d=d+0*N0;
else
    d=reshape(d,L,1);
end

if length(tau)==1
    tau=tau+0*N0;
else
    tau=reshape(tau,L,1);
end

if length(HHQ)==1
    hhq=HHQ*ones(size(N0,1),size(N0,2));
end

if length(WaningSpeed)<5
    WaningSpeed(5)=1/200;
end

Transmission_Reduction=zeros(21,MaxType,4); VEffI=zeros(21,MaxType,4); VEffS=zeros(21,MaxType,4);
RatioPf(isnan(RatioPf))=0.5;
for a=1:21
    Transmission_Reduction(a,:,:)=RatioPf(a)*tmpTransmission_Reduction(1,:,:)+(1-RatioPf(a))*tmpTransmission_Reduction(2,:,:);
    VEffI(a,:,:)=RatioPf(a)*tmpVEffI(1,:,:)+(1-RatioPf(a))*tmpVEffI(2,:,:);
    VEffS(a,:,:)=RatioPf(a)*tmpVEffS(1,:,:)+(1-RatioPf(a))*tmpVEffS(2,:,:);
end

M_from_to=M_from_to_H+M_from_to_O;

M_from_to_HAT=zeros(L,L);
for f=1:L
    for t=1:L
        M_from_to_HAT(f,t)=M_from_to(f,t)*d(t)*sigma(t)*(1+tau(f)*(1-d(f))/d(f));
    end
end


if length(InitialState)==1 
    [V,D]=eig(M_from_to_HAT);
    [R0, i]=max(abs(diag(D))); R0=R0/gamma;
    Age_structure=abs(V(:,i));
    
    E0=Age_structure; D0=Age_structure; U0=Age_structure; S0=N0;
    InitialState=zeros(1,(9*L+4*number_E_states*L)+L); L=21;
    InitialState(1:L)=S0'; m=L;
    for i=1:number_E_states
        InitialState(m+[1:3*L])=[E0'/number_E_states 0*E0' 0*E0']; m=m+3*L;
    end
    InitialState(m+[1:4*L])=[D0' zeros(1,length(E0)*2) U0'];
    % extra for newVariant, extra for VoC and 5 more for vaccination status,
    % extra S classes and 3 immune class (V1,V2 & R); and another 8 for
    % Waning & fully protected, and +1 for if you've had a booster
    InitialState=[InitialState zeros(1,(MaxType-1)*(length(InitialState)-L)) zeros(1,L*12)];
    
    for a=1:L
        InitialState(a:L:end)=InitialState(a:L:end)*N0(a)/sum(InitialState(a:L:end));
    end
else
    Age_structure=0;
end

options = odeset('RelTol', 1e-8);

m=number_E_states;

SP=MaxType*(9*L+4*m*L) + L;
NVacc1=InitialState(SP+[1:L]);
NVacc2=InitialState(SP+L+[1:L]);
NVacc0=N0'-NVacc1-NVacc2; NVacc0(NVacc0<0)=0;
V1(V1>NVacc0)=NVacc0(V1>NVacc0);
V2(V2>NVacc1)=NVacc1(V2>NVacc1);


[t, pop]=ode45(@Diff_Equations,[0:1:MaxTime],[InitialState],options,[L number_E_states MaxType reshape(M_from_to_H,1,[]) reshape(M_from_to_O,1,[]) ...
    sigma' d' tau' alpha gamma N0' hhq' nV_Beta nV_Speed V1 V2 V3 reshape(Transmission_Reduction,1,[]) ...
    reshape(VEffI,1,[]) reshape(VEffS,1,[]) WaningSpeed reshape(WaningEfficacy,1,[])]);

T=t; 
S=pop(:,1:L);  
EF=zeros(size(pop,1),MaxType*L); ES1=EF; ES2=EF; EQ=EF;
DF=EF; DS1=EF; DS2=EF;
UF=EF; US=EF;
DQF=EF; DQS=EF; UQ=EF;
a=[1:L];


% DO THE variants in reverse order
for Vtype=MaxType:-1:1
    SP=(Vtype-1)*(9*L+4*m*L);   % offset for VoC
    LL=(Vtype-1)*L;
    for i=1:m
        EF(:,a+LL)=EF(:,a+LL)+pop(:,SP+L+a+3*(i-1)*L);
        ES1(:,a+LL)=ES1(:,a+LL)+pop(:,SP+2*L+a+3*(i-1)*L);
        ES2(:,a+LL)=ES2(:,a+LL)+pop(:,SP+3*L+a+3*(i-1)*L);
        EQ(:,a+LL)=EQ(:,a+LL)+pop(:,SP+6*L+a+3*m*L+(i-1)*L);
    end
    DF(:,a+LL)=pop(:,SP+1*L+a+3*m*L); DS1(:,a+LL)=pop(:,SP+2*L+a+3*m*L); DS2(:,a+L)=pop(:,SP+3*L+a+3*m*L);
    UF(:,a+LL)=pop(:,SP+4*L+a+3*m*L); US(:,a+LL)=pop(:,SP+5*L+a+3*m*L);
    DQF(:,a+LL)=pop(:,SP+6*L+a+4*m*L); DQS(:,a+LL)=pop(:,SP+7*L+a+4*m*L); UQ(:,a+L)=pop(:,SP+8*L+a+4*m*L);
end

E=EF+ES1+ES2+EQ;  D=DF+DS1+DS2+DQF+DQS; U=UF+US+UQ;

%Ratios is a list of all those that can potentially be infected:  S, nV1, nV2+wV, nV3, R+wR, waneV, waneR.
SPV=MaxType*(9*L+4*m*L)+L;
Ratios=zeros(size(pop,1),21,7);  RatioS=Ratios;

Vtype=1; 
mVEffI=squeeze(VEffI(a,Vtype,:)).*sum(sum(D(:,a)+U(:,a))); 
mVEffS=squeeze(VEffS(a,Vtype,:)).*sum(sum(D(:,a)+U(:,a))); 
mWaneEff=squeeze(WaningEfficacy(Vtype,:)).*sum(sum(D(:,a)+U(:,a)));
R1=sum(sum(D(:,a)+U(:,a)));
tmp=sum(sum(D(:,a)+U(:,a)));
for Vtype=2:MaxType
    LL=(Vtype-1)*L;
    mVEffI=mVEffI+squeeze(VEffI(a,Vtype,:)).*sum(sum(D(:,a+LL)+U(:,a+LL))); 
    mVEffS=mVEffS+squeeze(VEffS(a,Vtype,:)).*sum(sum(D(:,a+LL)+U(:,a+LL))); 
    mWaneEff=mWaneEff+squeeze(WaningEfficacy(Vtype,:)).*sum(sum(D(:,a+LL)+U(:,a+LL)));
    tmp=tmp+sum(sum(D(:,a+LL)+U(:,a+LL)));
end
mVEffI=mVEffI/tmp;  mVEffS=mVEffS/tmp; mWaneEff=mWaneEff/tmp; R1=R1/tmp; R2=1-R1;

for A=1:L
      Ratios(:,A,1)=pop(:,A);
    Ratios(:,A,2)=(1-mVEffI(A,1)).*pop(:,SPV+2*L+A);
    Ratios(:,A,3)=(1-mVEffI(A,2)).*pop(:,SPV+3*L+A);
    Ratios(:,A,3)=Ratios(:,A,3) + (1-mVEffI(A,2))+pop(:,SPV+5*L+A);
    Ratios(:,A,3)=Ratios(:,A,3) + R1*(1-mVEffI(A,2))+pop(:,SPV+6*L+A) + R2*(1-mVEffI(A,2))+pop(:,SPV+7*L+A);
    Ratios(:,A,4)=(1-mVEffI(A,3)).*pop(:,SPV+4*L+A);
    Ratios(:,A,5)=(1-mWaneEff(2))*pop(:,SPV+8*L+A);
    Ratios(:,A,5)=Ratios(:,A,5) + R2*(1-mVEffI(A,2))+pop(:,SPV+6*L+A) + R1*(1-mVEffI(A,2))+pop(:,SPV+7*L+A);    
    Ratios(:,A,6)=(1-mWaneEff(1))*pop(:,SPV+8*L+A);
    Ratios(:,A,7)=(1-mVEffI(A,4))*pop(:,SPV+7*L+A);
    for Tp=1:MaxType, SP=Tp*(9*L+4*m*L)+L; Ratios(:,A,7)=Ratios(:,A,7)+(1-mVEffI(A,4))*pop(:,SP-L+A);  end %recovereds 
    
    RatioS(:,A,1:7)=Ratios(:,A,1:7);
    RatioS(:,A,2)=(1-mVEffS(A,1)).*pop(:,SPV+2*L+A);
    Ratios(:,A,3)=(1-mVEffS(A,2)).*pop(:,SPV+3*L+A);
    Ratios(:,A,3)=Ratios(:,A,3) + (1-mVEffS(A,2))+pop(:,SPV+5*L+A);
    Ratios(:,A,3)=Ratios(:,A,3) + R1*(1-mVEffS(A,2))+pop(:,SPV+6*L+A) + R2*(1-mVEffS(A,2))+pop(:,SPV+7*L+A);
    RatioS(:,A,4)=(1-mVEffS(A,3)).*pop(:,SPV+4*L+A);
    RatioS(:,A,7)=(1-mVEffS(A,4))*pop(:,SPV+7*L+A);
    for Tp=1:MaxType, SP=Tp*(9*L+4*m*L)+L; RatioS(:,A,7)=RatioS(:,A,7)+(1-mVEffS(A,4))*pop(:,SP-L+A);  end %recovereds  
end

% Find all new detectable / symptomatic infections
NVS=[1 nV_Speed];
nD=0*EF;
a=[1:L];
for Vtype=MaxType:-1:1
    LL=(Vtype-1)*L;
    SP=(Vtype-1)*(9*L+4*m*L);
    nD(:,a+LL)=(sum(RatioS,3)./sum(Ratios,3)).*(ones(size(pop,1),1)*d(a)').*(alpha*NVS(Vtype)*m*(pop(:,SP+L+a+3*(m-1)*L)+pop(:,SP+2*L+a+3*(m-1)*L)+pop(:,SP+3*L+a+3*(m-1)*L)+pop(:,SP+6*L+a+3*m*L+(m-1)*L)));
end
FinalState=pop(end,:);

end

% Calculates the differential rates used in the integration.
function dPop=Diff_Equations(t, pop, parameter)

L=parameter(1);
m=parameter(2);
MaxType=parameter(3);

s=3; 
l=L*L;  M_from_toH=reshape(parameter(s+[1:l]),L,L);    s=s+l;
l=L*L;  M_from_toO=reshape(parameter(s+[1:l]),L,L);    s=s+l;
l=L;    sigma = parameter(s+[1:l])';           s=s+l;
l=L;    d = parameter(s+[1:l])';               s=s+l;
l=L;    tau = parameter(s+[1:l])';             s=s+l;

l=1;    alpha=parameter(s+[1:l]);                 s=s+l;
l=1;    gamma=parameter(s+[1:l]);                 s=s+l; 
l=L;    N=parameter(s+[1:l])';                     s=s+l;
l=L;    hhq=parameter(s+[1:l])';                   s=s+l;
l=MaxType-1;    nV_Beta=parameter(s+[1:l]);               s=s+l;
l=MaxType-1;    nV_Speed=parameter(s+[1:l]);              s=s+l;

l=L;    V1=parameter(s+[1:l])';                    s=s+l;
l=L;    V2=parameter(s+[1:l])';                    s=s+l;
l=L;    V3=parameter(s+[1:l])';                    s=s+l;

l=21*MaxType*4;    Transmission_Reduction=reshape(parameter(s+[1:l]),21,MaxType,4);               s=s+l;
l=21*MaxType*4;    VEffI=reshape(parameter(s+[1:l]),21,MaxType,4);               s=s+l;
l=21*MaxType*4;    VEffS=reshape(parameter(s+[1:l]),21,MaxType,4);               s=s+l;
l=7;               WaningSpeed=parameter(s+[1:l]);           s=s+l;    % 1 is from REC to 2nd REC, 2 is from VACC_early_var to 2nd Vacc; 3 is VACC_Omicron to 2nd Vacc; 4,5,6 are to waning class, 7 is for boosters.
l=MaxType*8;       WaningEfficacy=reshape(parameter(s+[1:l]),MaxType,8);        s=s+l;    % 1 is waned from Rec, 2 is waned from Vacc. 7 and 8 are reduction in transmission.

% SET UP THE BASICS
S=pop(1:L); a=[1:L];
SPV=MaxType*(9*L+4*m*L)+L;

numV1=pop(SPV+[1:L]);   numV2=pop(SPV+1*L+[1:L]);
nV1=pop(SPV+2*L+[1:L]); nV2=pop(SPV+3*L+[1:L]);   nV3=pop(SPV+4*L+[1:L]);  % these are vaccinated.
WaneV11=pop(SPV+5*L+[1:L]);  WaneV12=pop(SPV+6*L+[1:L]);  % waning vaccination (2x2 to cope with faster waning for Omicron)
WaneV21=pop(SPV+7*L+[1:L]);  WaneV22=pop(SPV+8*L+[1:L]); 
wR=pop(SPV+9*L+[1:L]);  WaneR=pop(SPV+10*L+[1:L]);   % waning after infection.

RRA=zeros(21,1); RR=zeros(21,MaxType);
for Vtype=1:MaxType
    SP=(Vtype-1)*(9*L+4*m*L)+L;
    RR(:,Vtype)=pop(SP+8*L+[1:L]+4*m*L);
    RRA=RRA+RR(:,Vtype);
end

dPop=zeros(length(pop),1);  
nVPos1=SPV+2*L;  nVPos2=SPV+3*L;  nVPos3=SPV+4*L;
WaneVPos11=SPV+5*L; WaneVPos12=SPV+6*L; 
WaneVPos21=SPV+7*L; WaneVPos22=SPV+8*L; 
wRPos=SPV+9*L; WaneRPos=SPV+10*L; 
Boost=SPV+11*L;

% THIS IS ALL FOR THE WT VARIANT
Vtype=1;

SP=(Vtype-1)*(9*L+4*m*L);

%WanedEfficacy for 11, 12, 21, 22 is depended on the first index (ie pre-Omicron)
WEff=[1-VEffI(:,Vtype,2) 1-VEffI(:,Vtype,2) 1-ones(L,1)*WaningEfficacy(Vtype,2) 1-ones(L,1)*WaningEfficacy(Vtype,2)];
TEff=[Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)) Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)) ones(L,1)*(1-WaningEfficacy(Vtype,2))*WaningEfficacy(Vtype,8) ones(L,1)*(1-WaningEfficacy(Vtype,2))*WaningEfficacy(Vtype,8)];
SEff=[(1-VEffS(:,Vtype,2)) (1-VEffS(:,Vtype,2)) (1-ones(L,1)*WaningEfficacy(Vtype,2)) (1-ones(L,1)*WaningEfficacy(Vtype,2))];


SALL=S + (1-VEffI(:,Vtype,1)).*nV1 + (1-VEffI(:,Vtype,2)).*nV2 + (1-VEffI(:,Vtype,3)).*nV3 + (1-VEffI(:,Vtype,4)).*(wR+RRA-RR(:,Vtype)) + WEff(:,1).*WaneV11 + WEff(:,2).*WaneV12 + WEff(:,3).*WaneV21 + WEff(:,4).*WaneV22 + (1-WaningEfficacy(Vtype,1))*WaneR;

EF=zeros(L,m); ES1=zeros(L,m); ES2=zeros(L,m); EQ=zeros(L,m);
for i=1:m
    EF(:,i)=pop(SP+L+[1:L]+3*(i-1)*L);  
    ES1(:,i)=pop(SP+2*L+[1:L]+3*(i-1)*L); 
    ES2(:,i)=pop(SP+3*L+[1:L]+3*(i-1)*L);
    EQ(:,i)=pop(SP+6*L+[1:L]+3*m*L+(i-1)*L);
end
DF=pop(SP+1*L+[1:L]+3*m*L); DS1=pop(SP+2*L+[1:L]+3*m*L); DS2=pop(SP+3*L+[1:L]+3*m*L); 
UF=pop(SP+4*L+[1:L]+3*m*L); US=pop(SP+5*L+[1:L]+3*m*L); 
DQF=pop(SP+6*L+[1:L]+4*m*L); DQS=pop(SP+7*L+[1:L]+4*m*L); UQ=pop(SP+8*L+[1:L]+4*m*L);

IF=DF + tau.*UF;   IS=(DS1+DS2) + tau.*US; 

Trans_Scaling=(1*S + Transmission_Reduction(:,Vtype,1).*(1-VEffI(:,Vtype,1)).*nV1 + ...
    Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)).*nV2 + Transmission_Reduction(:,Vtype,3).*(1-VEffI(:,Vtype,3)).*nV3 + ...
    Transmission_Reduction(:,Vtype,4).*(1-VEffI(:,Vtype,4)).*(wR+RRA-RR(:,Vtype)) + ...
    (1-WaningEfficacy(Vtype,1))*WaningEfficacy(Vtype,7)*WaneR + TEff(:,1).*WaneV11 + TEff(:,2).*WaneV21 + TEff(:,3).*WaneV12 + TEff(:,4).*WaneV22)./SALL;

Symp_Scaling=(1*S + (1-VEffS(:,Vtype,1)).*nV1 + (1-VEffS(:,Vtype,2)).*(nV2) + (1-VEffS(:,Vtype,3)).*nV3 + (1-VEffS(:,Vtype,4)).*(wR+RRA-RR(:,Vtype)) + ...
    (1-WaningEfficacy(Vtype,1))*WaneR + SEff(:,1).*WaneV11 + SEff(:,2).*WaneV12 + SEff(:,3).*WaneV21 + SEff(:,4).*WaneV22)./SALL;


InfF=Trans_Scaling(a).*(sigma(a).*((IF+ IS)'*M_from_toO(:,a))');
InfS1=Trans_Scaling(a).*(sigma(a).*((DF)'*M_from_toH(:,a))');
InfS2=Trans_Scaling(a).*(sigma(a).*((tau.*UF)'*M_from_toH(:,a))');
InfSQ=Trans_Scaling(a).*(sigma(a).*((DQF)'*M_from_toH(:,a))');

s=0;
%dS
dPop(a+s)= - (InfF + InfS1 + InfS2 + InfSQ).*S(a)./N(a);             s=s+L;

dPop(a+nVPos1)= - (1-VEffI(:,Vtype,1)).*(InfF + InfS1 + InfS2 + InfSQ).*nV1(a)./N(a);       
dPop(a+nVPos2)= - (1-VEffI(:,Vtype,2)).*(InfF + InfS1 + InfS2 + InfSQ).*nV2(a)./N(a);       
dPop(a+nVPos3)= - (1-VEffI(:,Vtype,3)).*(InfF + InfS1 + InfS2 + InfSQ).*nV3(a)./N(a);       
dPop(a+wRPos)= - (1-VEffI(:,Vtype,4)).*(InfF + InfS1 + InfS2 + InfSQ).*wR(a)./N(a);       
dPop(a+WaneRPos)= - (1-WaningEfficacy(Vtype,1))*(InfF + InfS1 + InfS2 + InfSQ).*WaneR(a)./N(a);       
dPop(a+WaneVPos11)= - WEff(:,1).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV11(a)./N(a);
dPop(a+WaneVPos12)= - WEff(:,2).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV12(a)./N(a);
dPop(a+WaneVPos21)= - WEff(:,3).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV21(a)./N(a);
dPop(a+WaneVPos22)= - WEff(:,4).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV22(a)./N(a);
for Typ=1:MaxType
    if Typ~=Vtype
        SP=(Typ-1)*(9*L+4*m*L)+L; ps=(SP+8*L+[1:L]+4*m*L);
        dPop(ps) = - (1-VEffI(:,Vtype,4)).*(InfF + InfS1 + InfS2 + InfSQ).*RR(a,Typ)./N(a);  
    end
end
 
%dEF   dES1  and dES2 (first ones)
dPop(a+s)=InfF.*SALL(a)./N(a)  - m*alpha*EF(a,1);                s=s+L;
dPop(a+s)=InfS1.*SALL(a)./N(a) - m*alpha*ES1(a,1);               s=s+L;
dPop(a+s)=InfS2.*SALL(a)./N(a) - m*alpha*ES2(a,1);               s=s+L;

%dEF  dES1  & dES2 (subsequent ones)
for i=2:m
    dPop(a+s)=m*alpha*EF(a,i-1) - m*alpha*EF(a,i);      s=s+L;
    dPop(a+s)=m*alpha*ES1(a,i-1) - m*alpha*ES1(a,i);    s=s+L;
    dPop(a+s)=m*alpha*ES2(a,i-1) - m*alpha*ES2(a,i);    s=s+L;
end

nd=d.*Symp_Scaling;
%dDF  dDS1  and dDS2
dPop(a+s)=nd(a) .* (1-hhq(a)) .* alpha .* EF(a,m) * m - gamma*DF(a);       s=s+L;
dPop(a+s)=nd(a) .* alpha .* ES1(a,m) * m - gamma*DS1(a);                   s=s+L;
dPop(a+s)=nd(a) .* (1-hhq(a)) .* alpha .* ES2(a,m) * m - gamma*DS2(a);     s=s+L;

%dUF and dUS
dPop(a+s)=(1-nd(a)) .* alpha .* EF(a,m) * m - gamma*UF(a);                   s=s+L;
dPop(a+s)=(1-nd(a)) .* alpha .* (ES1(a,m) + ES2(a,m)) * m - gamma*US(a);     s=s+L;

%dEQ
dPop(a+s) = InfSQ.*SALL(a)./N(a) - alpha .* EQ(a,1) * m;       s=s+L;
for i=2:m
    dPop(a+s) = alpha .* EQ(a,i-1) * m - alpha .* EQ(a,i) * m;       s=s+L;
end

%dDQF
dPop(a+s) = nd(a) .* hhq(a) .* alpha .* EF(a,m) * m - gamma .* DQF(a);       s=s+L;

%dDQS
dPop(a+s) = nd(a) .* hhq(a) .* alpha .* ES2(a,m) * m  +  nd(a) .* alpha .* EQ(a,m) * m - gamma*DQS(a);       s=s+L;

%dUQ
dPop(a+s) = (1-nd(a)) .* alpha .* EQ(a,m) * m - gamma*UQ(a);       s=s+L;

%Recovereds
dPop(a+s) = dPop(a+s) + gamma*(DF(a)+DS1(a)+DS2(a)+UF(a)+US(a)+DQF(a)+DQS(a)+UQ(a));  s=s+L;



% REPEAT FOR ALL ADDITIONAL VARIANTS

for Vtype=2:MaxType

if Vtype<=3 % Wildtype, Alpha or Delta
    WEff=[1-VEffI(:,Vtype,2) 1-VEffI(:,Vtype,2) 1-ones(L,1)*WaningEfficacy(Vtype,2) 1-ones(L,1)*WaningEfficacy(Vtype,2)];
    TEff=[Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)) Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)) ones(L,1)*(1-WaningEfficacy(Vtype,2))*WaningEfficacy(Vtype,8) ones(L,1)*(1-WaningEfficacy(Vtype,2))*WaningEfficacy(Vtype,8)];
    SEff=[(1-VEffS(:,Vtype,2)) (1-VEffS(:,Vtype,2)) (1-ones(L,1)*WaningEfficacy(Vtype,2)) (1-ones(L,1)*WaningEfficacy(Vtype,2))];
else
    WEff=[1-VEffI(:,Vtype,2) 1-ones(L,1)*WaningEfficacy(Vtype,2) 1-VEffI(:,Vtype,2) 1-ones(L,1)*WaningEfficacy(Vtype,2)];
TEff=[Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)) ones(L,1)*(1-WaningEfficacy(Vtype,2))*WaningEfficacy(Vtype,8) Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)) ones(L,1)*(1-WaningEfficacy(Vtype,2))*WaningEfficacy(Vtype,8)];
SEff=[(1-VEffS(:,Vtype,2)) (1-ones(L,1)*WaningEfficacy(Vtype,2)) (1-VEffS(:,Vtype,2)) (1-ones(L,1)*WaningEfficacy(Vtype,2))];

end

SALL=S + (1-VEffI(:,Vtype,1)).*nV1 + (1-VEffI(:,Vtype,2)).*nV2 + (1-VEffI(:,Vtype,3)).*nV3 + (1-VEffI(:,Vtype,4)).*(wR+RRA-RR(:,Vtype)) + WEff(:,1).*WaneV11 + WEff(:,2).*WaneV12 + WEff(:,3).*WaneV21 + WEff(:,4).*WaneV22 + (1-WaningEfficacy(Vtype,1))*WaneR;

SP=(Vtype-1)*(9*L+4*m*L);

EF=zeros(L,m); ES1=zeros(L,m); ES2=zeros(L,m); EQ=zeros(L,m);
for i=1:m
    EF(:,i)=pop(SP+L+[1:L]+3*(i-1)*L);  
    ES1(:,i)=pop(SP+2*L+[1:L]+3*(i-1)*L); 
    ES2(:,i)=pop(SP+3*L+[1:L]+3*(i-1)*L);
    EQ(:,i)=pop(SP+6*L+[1:L]+3*m*L+(i-1)*L);
end
DF=pop(SP+1*L+[1:L]+3*m*L); DS1=pop(SP+2*L+[1:L]+3*m*L); DS2=pop(SP+3*L+[1:L]+3*m*L); 
UF=pop(SP+4*L+[1:L]+3*m*L); US=pop(SP+5*L+[1:L]+3*m*L); 
DQF=pop(SP+6*L+[1:L]+4*m*L); DQS=pop(SP+7*L+[1:L]+4*m*L); UQ=pop(SP+8*L+[1:L]+4*m*L); 

IF=DF + tau.*UF;   IS=(DS1+DS2) + tau.*US;

a=[1:L]; nBeta=nV_Beta(Vtype-1)*nV_Speed(Vtype-1);
Trans_Scaling=(1*S + Transmission_Reduction(:,Vtype,1).*(1-VEffI(:,Vtype,1)).*nV1 + ...
    Transmission_Reduction(:,Vtype,2).*(1-VEffI(:,Vtype,2)).*nV2 + Transmission_Reduction(:,Vtype,3).*(1-VEffI(:,Vtype,3)).*nV3 + ...
    Transmission_Reduction(:,Vtype,4).*(1-VEffI(:,Vtype,4)).*(wR+RRA-RR(:,Vtype)) + ...
    (1-WaningEfficacy(Vtype,1))*WaningEfficacy(Vtype,7)*WaneR + TEff(:,1).*WaneV11 + TEff(:,2).*WaneV21 + TEff(:,3).*WaneV12 + TEff(:,4).*WaneV22)./SALL;

Symp_Scaling=(1*S + (1-VEffS(:,Vtype,1)).*nV1 + (1-VEffS(:,Vtype,2)).*(nV2) + (1-VEffS(:,Vtype,3)).*nV3 + (1-VEffS(:,Vtype,4)).*(wR+RRA-RR(:,Vtype)) + ...
     (1-WaningEfficacy(Vtype,1))*WaneR + SEff(:,1).*WaneV11 + SEff(:,2).*WaneV12 + SEff(:,3).*WaneV21 + SEff(:,4).*WaneV22)./SALL;

InfF=nBeta*Trans_Scaling.*(sigma(a).*((IF+ IS)'*M_from_toO(:,a))');
InfS1=nBeta*Trans_Scaling.*(sigma(a).*((DF)'*M_from_toH(:,a))');
InfS2=nBeta*Trans_Scaling.*(sigma(a).*((tau.*UF)'*M_from_toH(:,a))');
InfSQ=nBeta*Trans_Scaling.*(sigma(a).*((DQF)'*M_from_toH(:,a))');

%dS
dPop(a)= dPop(a) - ((InfF + InfS1 +InfS2 + InfSQ).*S(a)./N(a));       %

dPop(a+nVPos1)= dPop(a+nVPos1) - (1-VEffI(:,Vtype,1)).*(InfF + InfS1 + InfS2 + InfSQ).*nV1(a)./N(a);       
dPop(a+nVPos2)= dPop(a+nVPos2) - (1-VEffI(:,Vtype,2)).*(InfF + InfS1 + InfS2 + InfSQ).*nV2(a)./N(a);       
dPop(a+nVPos3)= dPop(a+nVPos3) - (1-VEffI(:,Vtype,3)).*(InfF + InfS1 + InfS2 + InfSQ).*nV3(a)./N(a);       
dPop(a+wRPos)= dPop(a+wRPos) - (1-VEffI(:,Vtype,4)).*(InfF + InfS1 + InfS2 + InfSQ).*wR(a)./N(a);       
dPop(a+WaneRPos)= dPop(a+WaneRPos) - (1-WaningEfficacy(Vtype,1))*(InfF + InfS1 + InfS2 + InfSQ).*WaneR(a)./N(a);       
dPop(a+WaneVPos11)= dPop(a+WaneVPos11) - WEff(:,1).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV11(a)./N(a);
dPop(a+WaneVPos12)= dPop(a+WaneVPos12) - WEff(:,2).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV12(a)./N(a);
dPop(a+WaneVPos21)= dPop(a+WaneVPos21) - WEff(:,3).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV21(a)./N(a);
dPop(a+WaneVPos22)= dPop(a+WaneVPos22) - WEff(:,4).*(InfF + InfS1 + InfS2 + InfSQ).*WaneV22(a)./N(a);

for Typ=1:MaxType
    if Typ~=Vtype
        SP=(Typ-1)*(9*L+4*m*L)+L; ps=(SP+8*L+[1:L]+4*m*L);
        dPop(ps) = dPop(ps) - (1-VEffI(:,Vtype,4)).*(InfF + InfS1 + InfS2 + InfSQ).*RR(a,Typ)./N(a);  
    end
end

nalpha=alpha*nV_Speed(Vtype-1);
ngamma=gamma*nV_Speed(Vtype-1);
nd=d.*Symp_Scaling;

%dEF   dES1  and dES2 (first ones)
dPop(a+s)=InfF.*SALL(a)./N(a)  - m*nalpha*EF(a,1);                s=s+L;
dPop(a+s)=InfS1.*SALL(a)./N(a) - m*nalpha*ES1(a,1);               s=s+L;
dPop(a+s)=InfS2.*SALL(a)./N(a) - m*nalpha*ES2(a,1);               s=s+L;

%dEF  dES1  & dES2 (subsequent ones)
for i=2:m
    dPop(a+s)=m*nalpha*EF(a,i-1) - m*nalpha*EF(a,i);      s=s+L;
    dPop(a+s)=m*nalpha*ES1(a,i-1) - m*nalpha*ES1(a,i);    s=s+L;
    dPop(a+s)=m*nalpha*ES2(a,i-1) - m*nalpha*ES2(a,i);    s=s+L;
end

%dDF  dDS1  and dDS2
dPop(a+s)=nd(a) .* (1-hhq(a)) .* nalpha .* EF(a,m) * m - ngamma*DF(a);       s=s+L;
dPop(a+s)=nd(a) .* nalpha .* ES1(a,m) * m - ngamma*DS1(a);                   s=s+L;
dPop(a+s)=nd(a) .* (1-hhq(a)) .* nalpha .* ES2(a,m) * m - ngamma*DS2(a);     s=s+L;

%dUF and dUS
dPop(a+s)=(1-nd(a)) .* nalpha .* EF(a,m) * m - ngamma*UF(a);                   s=s+L;
dPop(a+s)=(1-nd(a)) .* nalpha .* (ES1(a,m) + ES2(a,m)) * m - ngamma*US(a);     s=s+L;

%dEQ
dPop(a+s) = InfSQ.*SALL(a)./N(a) - nalpha .* EQ(a,1) * m;       s=s+L;
for i=2:m
    dPop(a+s) = nalpha .* EQ(a,i-1) * m - nalpha .* EQ(a,i) * m;       s=s+L;
end

%dDQF
dPop(a+s) = nd(a) .* hhq(a) .* nalpha .* EF(a,m) * m - ngamma .* DQF(a);       s=s+L;

%dDQS
dPop(a+s) = nd(a) .* hhq(a) .* nalpha .* ES2(a,m) * m  +  nd(a) .* nalpha .* EQ(a,m) * m - ngamma*DQS(a);       s=s+L;

%dUQ
dPop(a+s) = (1-nd(a)) .* nalpha .* EQ(a,m) * m - ngamma*UQ(a);       s=s+L;

% Recovereds
dPop(a+s) = dPop(a+s) + ngamma*(DF(a)+DS1(a)+DS2(a)+UF(a)+US(a)+DQF(a)+DQS(a)+UQ(a));  s=s+L;
end

% NOW DO VACCINATION
SPV=MaxType*(9*L+4*m*L)+L; s=SPV;

%mean waning speed vaccine
WS=WaningSpeed(2)*WaningSpeed(5)./(WaningSpeed(2)+WaningSpeed(5));
if WaningSpeed(2)*WaningSpeed(5)==0  WS=0; end

%mean waning speed from recovery
WSR=WaningSpeed(1)*WaningSpeed(4)./(WaningSpeed(1)+WaningSpeed(4));
if WaningSpeed(1)*WaningSpeed(4)==0  WSR=0; end

%mean waning speed from booster
WSB=WaningSpeed(7);

%%% Subtract vaccination from S.
dPop(a) = dPop(a) - (V1.*S./(1+N-numV1-numV2)) ;

% -- everyone that has been vaccinated, numV1
dPop(a+s) = V1 - V2;          s=s+L;   %dPop(a+SPV)
% -- everyone that has been vaccinated, numV2
dPop(a+s) = V2;          s=s+L;        %dPop(a+SPV+L)

% -- vaccinated once and not recovered, nV1
dPop(a+s) = dPop(a+s) +  (V1.*S./(1+N-numV1-numV2))  - (V2.*nV1./(1+numV1));        s=s+L;  %dPop(a+SPV+2L)
% -- vaccinated twice and not recovered, nV2
dPop(a+s) = dPop(a+s) +  (V2.*nV1./(1+numV1));            s=s+L;    %dPop(a+SPV+3L)

% DO WANING.
% first recovereds
 NumRec=zeros(L,MaxType);
for Vtype=1:MaxType
    RPos=(Vtype-1)*(9*L+4*m*L)+L + (8*L+4*m*L);   % should be location of R.
    NumRec(a,Vtype)=pop(a+RPos);
    dPop(a+RPos)=dPop(a+RPos) - pop(a+RPos)*WaningSpeed(1);
    dPop(a+wRPos)=dPop(a+wRPos) + pop(a+RPos)*WaningSpeed(1);
end
NumRec(sum(NumRec(a,:),2)==0,:)=1;
dPop(a+wRPos)=dPop(a+wRPos) - pop(a+wRPos)*WaningSpeed(4);
dPop(a+WaneRPos)=dPop(a+WaneRPos) + pop(a+wRPos)*WaningSpeed(4);

%now 1 dose vaccine
dPop(a+nVPos1)=dPop(a+nVPos1) - pop(a+nVPos1)*WS;
dPop(a+WaneVPos22)=dPop(a+WaneVPos22) + pop(a+nVPos1)*WS;

%and 2 doses vaccine
dPop(a+nVPos2)=dPop(a+nVPos2) - pop(a+nVPos2)*WaningSpeed(2);
dPop(a+WaneVPos11)=dPop(a+WaneVPos11) + pop(a+nVPos2)*WaningSpeed(2) - pop(a+WaneVPos11)*WaningSpeed(5) - pop(a+WaneVPos11)*WaningSpeed(6);
dPop(a+WaneVPos12)=dPop(a+WaneVPos12) + pop(a+WaneVPos11)*WaningSpeed(6) - pop(a+WaneVPos12)*WaningSpeed(5);
dPop(a+WaneVPos21)=dPop(a+WaneVPos21) + pop(a+WaneVPos11)*WaningSpeed(5) - pop(a+WaneVPos21)*WaningSpeed(6);
dPop(a+WaneVPos22)=dPop(a+WaneVPos22) + pop(a+WaneVPos12)*WaningSpeed(5) + pop(a+WaneVPos21)*WaningSpeed(6);

%and 3 doses vaccine (if boosters act like third dose)
dPop(a+nVPos3)=dPop(a+nVPos3) - pop(a+nVPos3)*WSB;
dPop(a+WaneVPos22)=dPop(a+WaneVPos22) + pop(a+nVPos3)*WSB;

% NOW ADD IN BOOSTERS - assume everyone except susceptible and recovered
dPop(a+Boost) = V3;

% Vaccinated and boosting
Frac_UnBoosted=(1+N-S-pop(a+Boost));

% OPTION 1 BOOSTERS TAKE YOU TO V3.
dPop(a+nVPos3) = dPop(a+nVPos3) + V3.*(pop(a+nVPos2)+pop(a+WaneVPos11)+pop(a+WaneVPos12)+pop(a+WaneVPos21)+pop(a+WaneVPos22)+pop(a+WaneRPos))./Frac_UnBoosted;
dPop(a+nVPos2) = dPop(a+nVPos2) - V3.*pop(a+nVPos2)./Frac_UnBoosted;
dPop(a+WaneVPos11) = dPop(a+WaneVPos11) - V3.*pop(a+WaneVPos11)./Frac_UnBoosted;
dPop(a+WaneVPos12) = dPop(a+WaneVPos12) - V3.*pop(a+WaneVPos12)./Frac_UnBoosted;
dPop(a+WaneVPos21) = dPop(a+WaneVPos21) - V3.*pop(a+WaneVPos21)./Frac_UnBoosted;
dPop(a+WaneVPos22) = dPop(a+WaneVPos22) - V3.*pop(a+WaneVPos22)./Frac_UnBoosted;
% booster nV1 takes you to nV2
dPop(a+nVPos2) = dPop(a+nVPos2) + V3.*pop(a+nVPos1)./Frac_UnBoosted;
dPop(a+nVPos1) = dPop(a+nVPos1) - V3.*pop(a+nVPos1)./Frac_UnBoosted;

% Vaccine action on waning recovereds assume action occurs with first and
% booster doses (due to the short interval between 1st and 2nd).
VaccDoses=(V1./(1+N-numV1-numV2)) + V3./Frac_UnBoosted;
dPop(a+WaneRPos) = dPop(a+WaneRPos) - VaccDoses.*pop(a+WaneRPos);
for Vtype=1:MaxType
    RPos=(Vtype-1)*(9*L+4*m*L)+L + (8*L+4*m*L);
    dPop(a+RPos) = dPop(a+RPos) + VaccDoses.*pop(a+wRPos).*(NumRec(a,Vtype)./sum(NumRec(a,:),2));
    dPop(a+RPos) = dPop(a+RPos) + V1.*pop(a+WaneRPos).*(NumRec(a,Vtype)./sum(NumRec(a,:),2))./(1+N-numV1-numV2);
end
dPop(a+wRPos) = dPop(a+wRPos) - VaccDoses.*pop(a+wRPos);

end
