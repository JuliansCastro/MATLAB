close all; clc; clear;
T = 1;
L = 0.5;
C = 0.25;

index = 1;
for R = [0.1, 1.0, 10.0, 100.0, 1000.0]
    T = tf([R/L 0],[1 R/L (1/(L*C))]);
    figure(index)
    bode(T, 'k-');
    grid on;
    index = index + 1;
end


%{

% Original:
figure(index+1)
asymptoteBode(T); %Genera Bode asintotíco
%}