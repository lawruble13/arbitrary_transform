close all;
S_FREQ = 1000;
MAX = 100;
RUN_TIME = 1;
L = RUN_TIME*S_FREQ;
N_DRAW = 4;

PRIMES = primes(MAX);

t = 1/S_FREQ * (0:L-1);
y = ((-pi/2)*sawtooth(2*pi*t));

sawtoothIndices = 1:MAX*MAX;
sawtoothBaseMagnitudes = 1./sawtoothIndices;
sawtoothMagnitudes = sawtoothBaseMagnitudes;

used = NaN(size(PRIMES));
S = y;
active = 1;
td = figure();
hold on;
while ~isequal(used, ones(size(used)))
    if N_DRAW > 1 && sum(used(1+ceil(log2(N_DRAW-1)):end),'omitnan') == 0
        plot(t, ((-1)^(sum(used, 'omitnan')))*(-pi/2)*sawtooth(2*pi*t*active)/active, 'DisplayName', sprintf('Sawtooth %d',active));
    end
    i = 1;
    while used(i) == 1
        used(i) = NaN;
        i = i+1;
    end
    used(i) = 1;
    active = prod(PRIMES.*used,'omitnan');
    if active > 2*S_FREQ
        continue
    end
    S = S + ((-1)^(sum(used, 'omitnan')))*(-pi/2)*sawtooth(2*pi*t*active)/active; 
end
if(N_DRAW > 0)
    plot(t, ((-1)^(sum(used, 'omitnan')))*(-pi/2)*sawtooth(2*pi*t*active)/active, 'DisplayName', sprintf('Sawtooth %d',active));
end
hold on;
plot(t, S);
title("Sine wave approximation using sawtooth waves");
figure;
hold off;
T = fft(S);
P2 = abs(T/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = S_FREQ*(0:(L/2))/L;
plot(f, P1);
title("Single-Sided Amplitude Spectrum of sine wave approximation")


