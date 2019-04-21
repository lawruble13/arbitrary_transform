ITU_R;
close all;

recent = datetime('now');
start = recent;

S_FREQ = 44100;
MAX = 50;
RUN_TIME = 1;
W_FREQ = 100;
L = RUN_TIME*S_FREQ;

prime_arr = primes(MAX);
prime_arr = prime_arr(2:end);

t = 1/S_FREQ * (0:L-1);
squareWaves = containers.Map('KeyType','double','ValueType','double');
squareWaves(W_FREQ) = pi/4;
used = NaN(size(prime_arr));
active = 1;
curmax = length(prime_arr);
i=0;
while true
    count = 0;
    while i == 0 || sum(used,'omitnan') >= curmax
        count = count+1;
        i = 1;
        if isequal(used, ones(size(used)))
            used = zeros(size(used));
            break;
        end
        while used(i) == 1
            used(i) = NaN;
            i = i+1;
        end
        used(i) = 1;
        if mod(count, 100000) == 0 && (datetime('now')-recent > seconds(10))
            break;
        end
    end
    if(datetime('now')-recent > seconds(10)) && ~isequal(used, zeros(size(used)))
        recent = datetime('now');
        done = 0;
        for j=1:length(used)
            if ~isnan(used(j))
                done = done + used(j)*2^(j-1);
            end
        end
        percent_done = 100*done/(2^length(prime_arr));
        trem = (datetime('now')-start)*(100-percent_done)/percent_done;
        fprintf('Processed %d numbers, %.2f%% done. (ETA: %d minutes, %.2f seconds)\n', done, percent_done,idivide(seconds(trem),int32(60)),mod(seconds(trem),60));
    end
    if isequal(used, zeros(size(used)))
        break;
    end
    active = prod(prime_arr.*used,'omitnan');
    if(active > S_FREQ*2)
        curmax = sum(used, 'omitnan');
        continue;
    end
    if(isKey(squareWaves,active*W_FREQ))
        squareWaves(active*W_FREQ) = squareWaves(active*W_FREQ)+((-1)^(sum(used, 'omitnan')))*(pi/4)/active;
    else
        squareWaves(active*W_FREQ) = ((-1)^(sum(used, 'omitnan')))*(pi/4)/active;
    end
    i=0;
end
S = zeros(size(t));
for k = cell2mat(keys(squareWaves))
    %if(k<30000)
        S = S+squareWaves(k)*square(2*pi*k*t);
    %end
end
td = figure();
movegui(td, 'west');
hold on;
sound(S, S_FREQ);
plot(t, S,'DisplayName','Approximation');
plot(t, sin(2*pi*W_FREQ*t), 'DisplayName', 'Actual');
title("Sine wave approximation using square waves");
fd = figure();
movegui(fd, 'east');
T = fft(S);
P2 = abs(T/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = S_FREQ*(0:(L/2))/L;
P1d = mag2db(abs(P1));
plot(f, P1d,'DisplayName','Unfiltered');
hold on;
for i = 1:length(f)
    j = 1;
    while ~(ITU_R_468(j,1) < f(i) && ITU_R_468(j+1,1) > f(i)) && j < length(ITU_R_468)
        j=j+1;
    end
    if j == length(ITU_R_468)
        P1d(i) = P1d(i) + ITU_R_468(j,2);
    else
        P1d(i) = P1d(i) + (f(i)-ITU_R_468(j,1)) * (ITU_R_468(j+1,2)-ITU_R_468(j,2)) / (ITU_R_468(j+1,1)-ITU_R_468(j,1));
    end
end
plot(f, P1d,'DisplayName','Filtered by human hearing');
title("Single-Sided Amplitude Spectrum of sine wave approximation");
hold off;


