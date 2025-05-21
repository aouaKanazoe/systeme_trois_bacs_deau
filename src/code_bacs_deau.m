
% Ce code decris les éléments du système à trois bacs d'eau %

clc;
clear all; 
close all;

s = tf('s');

S = 0.0154;
Sn = 5e-5;
az13 = 0.4753;
az32 = 0.4833;
az20 = 0.9142;
g = 9.81;
Q10 = 3e-5;
Q20 = 0.5e-5;

a13 = az13 * Sn * sqrt(2 * g);
a32 = az32 * Sn * sqrt(2 * g);
a20 = az20 * Sn * sqrt(2 * g);

H20 = ((Q10 + Q20) / a20)^2;
H30 = (Q10 / a13)^2 + H20;
H10 = 2*(Q10 / a13)^2+ H20;

R13 = (2 * sqrt(abs(H10 - H30)) / a13);
R20 = (2 * sqrt(abs(H20)) / a20);
R32 = (2 * sqrt(abs(H30 - H20)) / a32);

A = [-1/(S*R13), 1/(S*R13), 0;
     1/(S*R13), -1/(S*(1/R13 + 1/R32)), 1/(S*R32);
     0, 1/(S*R32), -1/(S*(1/R32 + 1/R20))];

B = [1/S, 0;
     0, 0;
     0, 1/S];

C = [1, 0, 0];

D = 0;

SYS =ss(A,B,C,D);