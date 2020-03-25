% Edited: Julian A. Castro [juacastropa@unal.edu.co]
% The function bode_asymptotic() corresponds to bode(),
% but it also plots asymptotes for the magnitude and phase graphs.
% bode_asymptotic() only accepts monovariable transfer functions.
% Credits: Valerio Scordamaglia (2020). MATLAB Central File Exchange. Retrieved March 24, 2020.

function asymptoticBode(system)

%% SETTINGS:
activated = 'False';
%activated = 'True';
color = 'k.-';
asymptoteLineWidth = 1.5; %Line width for graphs

%%
disp(['Continuous Bode activated: <<', activated, '>>']);

[numerator,denominator] = tfdata(system,'v');
rootsZeros=roots(numerator);
rootsPoles=roots(denominator);
[iz,jz]=find(abs(rootsZeros)==0);
numz=length(iz);
[ip,jp]=find(abs(rootsPoles)==0);
nump=length(ip);
slope = numz-nump;
sz=[1:length(rootsZeros)];
sp=[1:length(rootsPoles)];
rootsZeros=rootsZeros(setdiff(sz,iz),:);
rootsPoles=rootsPoles(setdiff(sp,ip),:);
zroots2=abs(rootsZeros)-j*real(rootsZeros)./abs(rootsZeros);
proots2=abs(rootsPoles)-j*real(rootsPoles)./abs(rootsPoles);
zz = ones(size(zroots2,1),3);
zz(:,2) = real(zroots2);
zz(:,3) = imag(zroots2);
pp = -ones(size(proots2,1),3);
pp(:,2) = real(proots2);
pp(:,3) = imag(proots2);

%magnitude tracking
pp_mag = pp;
zz_mag = zz;
vect = [zz_mag; pp_mag];
vect = sortrows(vect,2);
wPoints = zeros(1,size(vect,1)+2); %Frecuency points
wPoints(:,[2:(size(vect,1)+1)]) = vect(:,2)';

try % Setting limits on x:
    wPoints(1,1)=0.01*vect(1,2);
    wPoints(1,end)=100*vect(end,2);
catch
    wPoints(1,1)=0.01*1;
    wPoints(1,2)=100*1;
end
[decimalGain, phase] = bode(numerator,denominator,wPoints(1));
grid on

yMagnitude(1) = 20*log10(decimalGain);
index = 2;
y_old = yMagnitude(1);
w_old = wPoints(1);
for w = wPoints(:,[2:length(wPoints)-1])
    yMagnitude(index) = 20*slope*log10(w) + y_old - 20*slope*log10(w_old);
    slope = slope + vect(index-1,1);
    w_old = w;
    y_old = yMagnitude(index);
    index = index+1;
end
yMagnitude(index) = 20*slope*log10(wPoints(end)) + y_old - 20*slope*log10(w_old);


%% MAGNITUDE ASYMPTOTIC GRAPH
npoints = length(wPoints);
subplot(2,1,1)
semilogx(wPoints(2:npoints-1), yMagnitude(2:npoints-1),color,'MarkerSize',14,'LineWidth', asymptoteLineWidth);
hold on; grid on;
semilogx(wPoints(1:2),yMagnitude(1:2),color,'LineWidth', asymptoteLineWidth); %increase
semilogx(wPoints(npoints-1:npoints),yMagnitude(npoints-1:npoints),color,'LineWidth', asymptoteLineWidth); %decrease
hold off;

% Setting limits of magnitude in the graph:
magnitudeMax = max([yMagnitude, decimalGain']);
magnitudeMin = min([yMagnitude, decimalGain']);
ySpan = magnitudeMax - magnitudeMin;
for yDelta = [2 4 10 20 40]
    yhlp = ySpan/yDelta;
    if yhlp < 8
        break;
    end
end
magnitudeMin = yDelta*floor(magnitudeMin/yDelta);
magnitudeMax = yDelta*ceil(magnitudeMax/yDelta);
if magnitudeMin == magnitudeMax
    magnitudeMin = magnitudeMin - 2; magnitudeMax = magnitudeMax + 2;
end
magLimitMargin =(magnitudeMax - magnitudeMin)*0.1;
set(get(gcf, 'CurrentAxes'),'YTick',magnitudeMin:yDelta:magnitudeMax);
set(get(gcf, 'CurrentAxes'),'YLim',[magnitudeMin magnitudeMax]);
yScaleMag = axis;
yScaleMag(3) = magnitudeMin - magLimitMargin;
yScaleMag(4) = magnitudeMax + magLimitMargin;
axis(yScaleMag);
ylabel('Magnitude (dB)')

% MAGNITUDE CONTINUOUS GRAPH
if strcmp(activated, 'True')
    [AAA,BBB,CCC,DDD] = tf2ss(numerator,denominator);
    sys = ss(AAA,BBB,CCC,DDD);
    [m1, m2, m3] = bode(numerator,denominator,{wPoints(1,1),wPoints(1,end)});
    hold on;
    plot(m3,20*log10(m1), 'r'); % continuous Bode graph
    hold off;
end

clear interval;

%% PHASE ASYMPTOTIC GRAPH
[ii,jj] = find(vect(:,3)<0);
vect(ii,1) = -(vect(ii,1));

index=1;
for w = 1:size(vect,1)
    wPoints(1,index) = 0.1*vect(w,2);
    wPoints(2,index) = vect(w,1);
    index = index + 1;
    wPoints(1,index) = vect(w,2);
    wPoints(2,index) = 0;
    index = index + 1;
    wPoints(1,index) = 10*vect(w,2);
    wPoints(2,index) = -vect(w,1);
    index = index + 1;
end
try
    wPoints = wPoints'; wPoints = sortrows(wPoints,1); wPoints = wPoints';
catch
    wPoints(:,1) = [0.01*1;0]; wPoints(:,2) = [100*1;0];
    wPoints = wPoints'; wPoints = sortrows(wPoints,1); wPoints = wPoints';
end
wPoints = [[0.1*wPoints(1,1);0] wPoints [10*wPoints(1,end);0]];
%[decimalGain, phase, w3]=bode(numerator,denominator,wPoints(1));


pp = phase;
x = round(phase/90);
phase = 90*x;
yPhase(1) = phase;
yP_old = yPhase(1);
w_old = wPoints(1,1);
slope = 0;
index = 2;
for w = wPoints(1,[2:size(wPoints,2)-1])
    yPhase(index) = 45*slope*log10(w) + yP_old - 45*slope*log10(w_old);
    slope = slope + wPoints(2,index);
    w_old = w;
    yP_old = yPhase(index);
    index = index+1;
end
yPhase(index) = yP_old;

subplot(2,1,2)
npoints = length(wPoints);
wX = wPoints(1,:);
semilogx(wX(2:npoints-1),yPhase(2:npoints-1),color,'MarkerSize',14,'LineWidth',asymptoteLineWidth)
hold on;
% No dots at start and end points:
semilogx(wX(1:2),yPhase(1:2),color,'LineWidth',asymptoteLineWidth); % Head
semilogx(wX(npoints-1:npoints),yPhase(npoints - 1:npoints),color,'LineWidth',asymptoteLineWidth); % Tail
grid on
xlabel('Frequency (rad/s)')
ylabel('Phase (deg)')

% Setting limits of phase in the graph:
h = get(findall(get(gcf,'Children'),'String','Phase (deg)'),'Parent');
axes(h);

% Phase asymptotes may exceed current upper and lower limit
% values in phase diagram. Find max and min phase asymptote
% values, and re-scale phase diagram:

phaseMax = max([yPhase,phase']);
phaseMin = min([yPhase,phase']);
ySpan = phaseMax - phaseMin;
for yDelta = [10 15 30 45 90 180]
    yhlp = ySpan/yDelta;
    if yhlp < 8
        break;
    end
end
phaseMin = yDelta*floor(phaseMin/yDelta);
phaseMax = yDelta*ceil(phaseMax/yDelta);
if phaseMin == phaseMax
    phaseMin = phaseMin - 45; phaseMax = phaseMax + 45;
end
phLimitMargin = (phaseMax-phaseMin)*0.1;
set(get(gcf, 'CurrentAxes'),'YLim',[phaseMin phaseMax]);
yScaleph = axis;
yScaleph(3) = phaseMin - phLimitMargin;
yScaleph(4) = phaseMax + phLimitMargin;
axis(yScaleph);

%% PHASE CONTINUOUS GRAPH
if strcmp(activated, 'True')
    [m1,m2,m3]=bode(numerator,denominator,{wPoints(1,1),wPoints(1,end)});
    diff = m2(1) - pp;
    m2 = m2 - diff;
    hold on;
    semilogx(m3,m2);
    hold off;
end
end
