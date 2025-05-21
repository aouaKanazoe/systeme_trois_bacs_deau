% KANAZOE Aoua Asmaa

clear all 
close all

% Actionneur de la pompe 1

Kq = 1.6e5;              %coefficients de conversion du CNA
Offsetq = -9.2592;       % autour du point d'Ã©quilibre

% Capteurs de niveau (A COMPLETER)
% Avril 2001
% ----------
Kh1 = -3.72429e-2;    
OffsetH1 = 0.3586575 ;
Kh2 = -3.7037e-2;
OffsetH2 = 0.3700;
Kh3 = -3.7037e-2;
OffsetH3 = 0.3770;

%____________________________
%DONNEES DU PROCEDE DE 3 BACS
S  = 0.0154;      %surface d'un rÃ©servoir
Sn = 5e-5;        %surface d'un tuyau entre 2 bacs
g  = 9.81;        %acceleration terrestre

%Coefficients de debit :

az13 = 0.4753;  az32 = 0.4833;  az20 = 0.9142; % identification autour du point d'Ã©quilibre
a13 = az13*Sn*sqrt(2*g);    
a32 = az32*Sn*sqrt(2*g);    
a20 = az20*Sn*sqrt(2*g);  

%______________________________________________________
%RESISTANCES A L'ECOULEMENT AUTOUR DU POINT D'EQUILIBRE

%Q10 = 3.5e-5;   %(m3/s) debit de la pompe 1 au point d'Ã©quilibre
%Q20 = 0.5e-5; %(m3/s) debit de la pompe 2 au point d'Ã©quilibre
Q10 = 5e-5;      %(m3/s) debit de la pompe 1 au point d'Ã©quilibre
Q20 = 2e-5
%H10 = 0.247;  H20 = 0.03;  H30 = 0.137; % ancienne valeurs niveaux d'eau Ã  l'Ã©quilibre H0 

%H10 = 0.1923 ; H20 = 0.0299; H30 = 0.1111 % nouvelle valeurs 
%%%%%%%%%%%%% calcul des niveaux d'eau à l'équilibre %%%%%%%%%%%%%%%

H20 = ((Q10 + Q20) / a20)^2
H30 = (Q10 / a13)^2 + H20
H10 = 2*(Q10 / a13)^2+ H20

%Resistances a l'ecoulement
R13 = 2*sqrt(H10-H30)/a13;
R32 = 2*sqrt(H30-H20)/a32;
R20 = 2*sqrt(H20)/a20;

%________________________
%Consigne de niveau d'eau

w1 = 0.03;
Kdebit=Q10/H10 ;
Q1max=12e-5 ;
Te=1;

%DÃ©finition des matrices du systÃ¨me linÃ©arisÃ©
a= -1/(S*R13);
b= 1/(S*R13);
c= 0;
d= 1/(S*R13);
e= (-1/S)*((1/R13)+(1/R32));
f= 1/(S*R32);
g= 0;
h= f;
i= (-1/S)*((1/R32)+(1/R20));

%Definition  de la forme linÃ©arisÃ©e du systÃ¨me Ã  3 bacs d'eau
A=[a b c;d e f;g h i];

B=[1/S 0;0 0;0 1/S];

C=[1 0 0];

C46=[1 0 0;0 1 0;0 0 1]

D=0;

Syst= ss(A,B,C,D)
Syst46=ss(A,B,C46,D)

% Analyse de la stabilité asymptotique du système linéarisation
vp= eig(A) % le système linéarisation est stable car les valeurs propres valent -0.0525,-0.0241,-0.0033
% cela était il prévisible au vu des considérations prises sur le procédé reel?

% Analyse de la commandabilité du système

Com1= ctrb(A,B(:,1))
rangCom1 =rank(Com1) % le rang de la matrice de commandabilitÃ© est de 3 donc le systÃ¨me est commandable par la 1Ã¨re entrÃ©e


Com2= ctrb(A,B(:,2))
rangCom2 =rank(Com2) % le rang de la matrice de commandabilitÃ© est de 3 donc le systÃ¨me est commandable par avec la 2iÃ¨me entrÃ©e

Com12= ctrb(A,B)
rangCom12 =rank(Com12) % le rang de la matrice de commandabilitÃ© est de 3 donc le systÃ¨me est commandable par les 2 entrÃ©es

% Les indices de commandabilite 

w = [64.9351 0 -0.6691; 0 0 0.6691; 0 64.9351 0]
rang_w = rank(w)

% Analyse de l'observabilité

Obsv = obsv(A,C)
rang_obsv = rank(Obsv)
%Analyse de la fonction de transfert%
F=tf(Syst) %pour obtenir la matrice de transfert
poles1=zpk(Syst)
poles = pole(F)
gainSys= dcgain(F) %calcul du gain statique

F1=F(1,1)
zeros1 = zero(F1)
gain1= dcgain(F1)
%%%%%%%%%%%%%%%%%%%%%
F2=F(1,2)
zeros2 = zero(F2)
gain2= dcgain(F2)

%%%%%%%%%%%%%%%%%%Etude des caractéristique temporelles de F pour Q10=5*10^-5 et q20 = 2*10^-5

% figure(5)
% step(F)
% % Obtenir les informations sur la réponse indicielle
% stepInfo = stepinfo(F);
% 
% % Afficher les temps de réponse dans la fenêtre de commande
% fprintf('Temps de montée: %.4f\n', stepInfo.RiseTime);
% fprintf('Temps de réponse: %.4f\n', stepInfo.SettlingTime);
% fprintf('Temps de dépassement: %.4f\n', stepInfo.Overshoot);
figure(6)
step(F1)
% Obtenir les informations sur la réponse indicielle
stepInfo1 = stepinfo(F1);

% Afficher les temps de réponse dans la fenêtre de commande
fprintf('Le temps de montée du transfert 1: %.4f\n', stepInfo1.RiseTime);
fprintf('Le temps de réponse du tansfert 1 : %.4f\n', stepInfo1.SettlingTime);
fprintf('Le temps de dépassement du transfert 1: %.4f\n', stepInfo1.Overshoot);
figure(7)
step(F2)
% Obtenir les informations sur la réponse indicielle
stepInfo2 = stepinfo(F2);

% Afficher les temps de réponse dans la fenêtre de commande
fprintf('Le temps de montée du transfert 2: %.4f\n', stepInfo2.RiseTime);
fprintf('Le temps de réponse du tansfert 2 : %.4f\n', stepInfo2.SettlingTime);
fprintf('Le temps de dépassement du transfert 2: %.4f\n', stepInfo2.Overshoot);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F46 =tf(Syst46)
figure(8)
step(Syst46)

% stepInfo46 = stepinfo(F46);
% 
% % Afficher les temps de réponse dans la fenêtre de commande
% fprintf('Le temps de montée du transfert : %.4f\n', stepInfo2.RiseTime);
% fprintf('Le temps de réponse du tansfert : %.4f\n', stepInfo2.SettlingTime);
% fprintf('Le temps de dépassement du transfert : %.4f\n', stepInfo2.Overshoot);



% Commande par retour d'etat par la methode de Bass-Gura q1

vpdes=[-0.033;-0.033;-0.033]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q1=[1;1] % 
B_tilde1 = B*q1

Co1=ctrb (A, B_tilde1 )
rang_Co1 =rank(Co1)
K1=acker(A,B_tilde1,vpdes)
K11=q1*K1
% N1 le prÃ©compensateur
N1=1/(C*inv(B_tilde1*K1 - A)*B_tilde1)% calcul de N1
systCorrige1 = ss(A-B_tilde1*K1,B_tilde1,C,D)
ComSystcor1 = ctrb(systCorrige1)
rangCor1 = rank(ComSystcor1)

SystCompCor1= ss(A-B*q1*K1,B*q1*N1,C,D)


%figure(1)
%step(SystCompCor1)
% opt = RespConfig;
% opt.Amplitude = 0.05;
% yl = step(SystCompCor1, t, opt)

opt = stepDataOptions;
opt.InputOffset = 0;
opt.StepAmplitude = 0.05;

figure(1)
step(SystCompCor1,opt)

%vpRecompense = eig(A-B*K)
%vpRetourE= eig(A)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q2 = [0.8 ; 0.5] 
B_tilde2 = B*q2

Co2=ctrb (A, B_tilde2 )
rang_Co2 =rank(Co2)
K2=acker(A,B_tilde2,vpdes)
K22=q2*K2
% N2 le prÃ©compensateur
N2=1/(C*inv(B_tilde2*K2 - A)*B_tilde2)% calcul de N2
systCorrige2 = ss(A-B_tilde2*K2,B_tilde2,C,D)
ComSystcor2 = ctrb(systCorrige2)
rangCor2 = rank(ComSystcor2)

SystCompCor2 = ss(A-B*q2*K2,B*q2*N2,C,D)


opt = stepDataOptions;
opt.InputOffset = 0;
opt.StepAmplitude = 0.05;
figure(2)
step(SystCompCor2,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q3=[1;0]
B_tilde3 = B*q3

Co3=ctrb (A, B_tilde3 )
rang_Co2 =rank(Co3)
K3=acker(A,B_tilde3,vpdes)
K33=q2*K2

N3=1/(C*inv(B_tilde3*K3 - A)*B_tilde3)% calcul de N3

systCorrige3 = ss(A-B_tilde3*K3,B_tilde3,C,D)
ComSystcor3 = ctrb(systCorrige3)
rangCor3 = rank(ComSystcor3)
SystCompCor3= ss(A-B*q3*K3,B*q3*N3,C,D)

opt = stepDataOptions;
opt.InputOffset = 0;
opt.StepAmplitude = 0.05;
figure(3)
step(SystCompCor3,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q4 = [0;1]
B_tilde4 = B*q4

Co4=ctrb (A, B_tilde4 )
rang_Co2 =rank(Co4)
K4=acker(A,B_tilde4,vpdes)
K44=q4*K4

N4=1/(C*inv(B_tilde3*K4 - A)*B_tilde4)% calcul de N3

systCorrige4 = ss(A-B_tilde4*K4,B_tilde4,C,D)
ComSystcor4 = ctrb(systCorrige4)
rangCor4 = rank(ComSystcor4)
SystCompCor4= ss(A-B*q4*K4,B*q4*N4,C,D)

opt = stepDataOptions;
opt.InputOffset = 0;
opt.StepAmplitude = 0.05;
figure(4)
step(SystCompCor4,opt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Commande par retour d'Ã©tat par mise sous forme commandable 
Cob= [64.9351 -0.6691 0;0 0.6691 0; 0 0 64.9351]
InvCob= inv(Cob)
P=[Cob(2,:);Cob(2,:)*A;Cob(3,:)]

Ac=P*A*inv(P)
Bc=P*B
Cc=C*inv(P)
Dc=D
%Bcc= diag(Bc)
%Bcc= diag(0,5)
Bcc= [0 0 ;4 9 ; 0 4000] 
R= inv(Bc'* Bc)*Bc'*Bcc
vpdes2=[-0.033;-0.030;-0.020]

Command=ctrb(Ac, Bcc )
rang_Command =rank(Command)
%K=acker(Ac,Bcc,vpdes)
K_tilde=place(Ac,Bcc,vpdes2)
%%%%%% Pour avoir K dans la base initiale %%%%%%%%%
P_inv = inv(P)
Bcc_init = B*R
K_init=place(A,Bcc_init,vpdes2)