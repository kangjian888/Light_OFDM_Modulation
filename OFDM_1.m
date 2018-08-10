clear all;
close all; 

IFFT_Length = 1024;
carrier_count = 200;
cp_count = 64; % 1/16 of IFFT Length
symbols_per_carrier = 50;
bit_per_symbol = 2;
SNR = 10;
baseband_out_length = carrier_count * symbols_per_carrier * bit_per_symbol;

carrier = 2:(carrier_count+1);
conjugate_carriers = IFFT_Length - carrier + 2;

% generate baseband signal
baseband_out = round(rand(1,baseband_out_length));
convert_matrix = reshape(baseband_out,bit_per_symbol,length(baseband_out)/bit_per_symbol);
% transform binary to decimal
for k = 1:(length(baseband_out)/bit_per_symbol)
	modulo_baseband(k) = 0;
	for i = 1:bit_per_symbol
		modulo_baseband(k) = modulo_baseband(k)+ convert_matrix(i,k)*2^(bit_per_symbol-i);
	end
end

% Serial to Parallel
carrier_matrix_1 = reshape(modulo_baseband,carrier_count,symbols_per_carrier)'; %"'" means transposed matrix
% Differential Coding for each carrier symbol
carrier_matrix = [zeros(1,carrier_count);carrier_matrix_1]; %add zero on the first row
for i = 2:(symbols_per_carrier + 1)
	carrier_matrix(i,:) = rem(carrier_matrix(i,:)+carrier_matrix(i-1,:),2^bit_per_symbol);
end
%I think this is PSK modulation scheme
carrier_matrix = carrier_matrix * (2*pi)/(2^bit_per_symbol);
[X,Y] = pol2cart(carrier_matrix,ones(size(carrier_matrix,1),size(carrier_matrix,2)));
complex_carrier_matrix = complex(X,Y);

IFFT_modulation = zeros(symbols_per_carrier+1,IFFT_Length);
IFFT_modulation(:,carrier) = complex_carrier_matrix;
IFFT_modulation(:,conjugate_carriers) = conj(complex_carrier_matrix);


figure(1)
stem(0:IFFT_Length-1,abs(IFFT_modulation(2,1:IFFT_Length)),'b*-')
grid on
axis ([0 IFFT_Length -0.5 1.5])
ylabel('Magnitude')
xlabel('IFFT Bin')
title('OFDM Carrier Frequency Magnitude')

figure(2)
plot(0:IFFT_Length-1,(180/pi)*angle(IFFT_modulation(2,1:IFFT_Length)),'go')
hold on
stem(carrier-1,(180/pi)*angle(IFFT_modulation(2,carrier)),'b*-')
stem(conjugate_carriers-1,(180/pi)*angle(IFFT_modulation(2,conjugate_carriers)),'b*-')
axis([0 IFFT_Length -200 +200])
grid on
ylabel('Phase(degree)')
xlabel('IFFT Bin')
title('OFDM Carrier Phase')
 
%from frequency domain to time domain
time_wave_matrix = ifft(IFFT_modulation');
time_wave_matrix = time_wave_matrix';

figure(3)
plot(0:IFFT_Length-1,time_wave_matrix(2,:))
grid on
ylabel('Amplitude')
xlabel('Time')
title('OFDM Time Signal, One Symbol Period')

for f = 1 : carrier_count
	temp(1:IFFT_Length) = 0 + 0j;
	temp(carrier(f)) =  IFFT_modulation(2,carrier(f));
	temp(conjugate_carriers(f)) = IFFT_modulation(2,conjugate_carriers(f));
	temp_time = ifft(temp');
	figure(4)
	plot(0:IFFT_Length-1,temp_time)
	hold on
end
grid on
ylabel('Amplitude')
xlabel('Time')
title('Separated Time Waveforms Carrier')
%We do not add windows here
for i = 1:symbols_per_carrier+1
	%windowed_time_wave_matrix(i,:) = real(time_wave_matrix(i,:)).*hamming(IFFT_Length)';
	windowed_time_wave_matrix(i,:) = time_wave_matrix(i,:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The Below Code has no CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%paralell to serial
ofdm_modulation = reshape(windowed_time_wave_matrix',1,IFFT_Length*(symbols_per_carrier+1));
%figure of all period OFDM
temp_time = IFFT_Length*(symbols_per_carrier+1);
figure(5)
plot(0:temp_time-1,ofdm_modulation)
grid on
ylabel('Amplitude(volts)')
xlabel('Time(samples)')
title('OFDM Time Signal')
%frequency domain of OFDM signal
symbols_per_average = ceil(symbols_per_carrier/5);
avg_temp_time = IFFT_Length*symbols_per_average;
averages = floor(temp_time/avg_temp_time);
average_fft(1:avg_temp_time) = 0;
for a = 0:(averages-1)
	subset_ofdm = ofdm_modulation((a*avg_temp_time)+1:(a+1)*avg_temp_time);
	subset_ofdm_f = abs(fft(subset_ofdm));
	average_fft = average_fft + (subset_ofdm_f/averages);
end
average_fft_log = 20*log10(average_fft);
figure(6)
plot((0:(avg_temp_time - 1))/avg_temp_time,average_fft_log);
hold on
plot(0:1/IFFT_Length:1,-35,'rd')
grid on
axis([0 0.5 -40 max(average_fft_log)])
ylabel('Magnitude(dB)')
xlabel('Normalized Frequency(0.5 = fs/2)')
title('OFDM Signal Spectrum')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add CP in windowed_time_wave_matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
windowed_time_wave_matrix_add_cp = zeros((symbols_per_carrier+1),(IFFT_Length + cp_count));
cp_prefix = (IFFT_Length - cp_count + 1) : IFFT_Length;
windowed_time_wave_matrix_add_cp = [windowed_time_wave_matrix(:,cp_prefix), windowed_time_wave_matrix];

%paralell to serial
ofdm_modulation_cp = reshape(windowed_time_wave_matrix_add_cp',1,(IFFT_Length + cp_count) *(symbols_per_carrier+1));

%figure of all period OFDM
temp_time_cp = (IFFT_Length + cp_count)*(symbols_per_carrier+1);
figure(7)
plot(0:temp_time_cp-1,ofdm_modulation_cp)
grid on
ylabel('Amplitude(volts)')
xlabel('Time(samples)')
title('OFDM Time Signal Add CP')

%frequency domain of OFDM signal which is added CP
symbols_per_average_cp = ceil(symbols_per_carrier/5);
avg_temp_time_cp = (IFFT_Length+cp_count)*symbols_per_average_cp;
averages_cp = floor(temp_time_cp/avg_temp_time_cp);
average_fft_cp(1:avg_temp_time_cp) = 0;
for a = 0:(averages_cp - 1)
	subset_ofdm_cp = ofdm_modulation_cp((a*avg_temp_time_cp)+1:(a+1)*avg_temp_time_cp);
	subset_ofdm_f_cp = abs(fft(subset_ofdm_cp));
	average_fft_cp = average_fft_cp + (subset_ofdm_f_cp/averages_cp);
end
average_fft_log_cp = 20*log10(average_fft_cp);
figure(8)
plot((0:(avg_temp_time_cp - 1))/avg_temp_time_cp,average_fft_log_cp);
hold on
plot(0:1/IFFT_Length:1,-35,'rd')
grid on
axis([0 0.5 -40 max(average_fft_log_cp)])
ylabel('Magnitude(dB)')
xlabel('Normalized Frequency(0.5 = fs/2)')
title('OFDM Signal Spectrum Added CP')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Put Data to the channel
%Tx_data = ofdm_modulation;
%Tx_signal_power = var(Tx_data);
%linear_SNR = 10^(SNR/10);
%noise_sigma = Tx_signal_power/linear_SNR;
%noise_sigma_factor = sqrt(noise_sigma);
%noise = randn(1,length(Tx_data))*noise_sigma_factor;
%Rx_data = Tx_data + noise;

% Put Data added CP to channel
Tx_data_cp = ofdm_modulation_cp;
%Simulate the multi-path signal channel
d1 = 4;
a1 = 0.2;
d2 = 5;
a2 = 0.3;
d3 = 6;
a3 = 0.4;
d4 = 7;
a4 = 0.5;

copy1 = zeros(size(Tx_data_cp));
for i = (1 + d1) : length(Tx_data_cp)
	copy1(i) = a1*Tx_data_cp(i - d1);
end

copy2 = zeros(size(Tx_data_cp));
for i = (1 + d2) : length(Tx_data_cp)
	copy1(i) = a2*Tx_data_cp(i - d2);
end

copy3 = zeros(size(Tx_data_cp));
for i = (1 + d3) : length(Tx_data_cp)
	copy1(i) = a3*Tx_data_cp(i - d3);
end

copy4 = zeros(size(Tx_data_cp));
for i = (1 + d4) : length(Tx_data_cp)
	copy1(i) = a1*Tx_data_cp(i - d4);
end
Tx_data_cp = copy1 + copy2 + copy3 + copy4;

Tx_signal_power_cp = var(Tx_data_cp);
linear_SNR = 10^(SNR/10);
noise_sigma_cp = Tx_signal_power_cp/linear_SNR;
noise_sigma_factor_cp = sqrt(noise_sigma_cp);
noise_cp = randn(1,length(Tx_data_cp))*noise_sigma_factor_cp;
Rx_data_cp = Tx_data_cp + noise_cp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Below is Recieve without CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve
%Rx_data_matrix = reshape(Rx_data,IFFT_Length,symbols_per_carrier+1);
%Rx_spectrum = fft(Rx_data_matrix);% FFT to the colomn
%
%figure(9)
%stem(0:IFFT_Length-1,abs(Rx_spectrum(1:IFFT_Length,2)),'b*-')
%grid on
%axis([0 IFFT_Length -0.5 1.5])
%ylabel('Magnitude')
%xlabel('FFT Bin')
%title('OFDM Recieve Spectrum, Magnitude')
%
%figure(10)
%plot(0:IFFT_Length-1,(180/pi)*angle(Rx_spectrum(1:IFFT_Length,2)),'go')
%hold on
%stem(carrier-1,(180/pi)*angle(Rx_spectrum(carrier,2)),'b*-')
%stem(conjugate_carriers-1,(180/pi)*angle(Rx_spectrum(conjugate_carriers,2)),'b*-')
%axis([0 IFFT_Length -200 +200])
%grid on
%ylabel('Phase(degree)')
%xlabel('FFT Bin')
%title('OFDM Recieve Spectrum, Phase')
%
%%Distribution of Recieved Signal
%figure(11)
%Rx_carriers = Rx_spectrum(carrier,:)';
%Rx_phase_P = angle(Rx_carriers);
%Rx_mag_P = abs(Rx_carriers);
%polar(Rx_phase_P,Rx_mag_P,'bd');
%
%Rx_phase = angle(Rx_carriers)*(180/pi);
%Rx_phase = rem(Rx_phase+360,360);
%
%Rx_decoded_phase = diff(Rx_phase);
%Rx_decoded_phase = rem(Rx_decoded_phase+360,360);
%
%% Judgement
%base_phase = 360/2^bit_per_symbol;
%delta_phase = base_phase/2;
%Rx_decoded_symbols = zeros(size(Rx_decoded_phase,1),size(Rx_decoded_phase,2));
%
%for i = 1:(2^bit_per_symbol - 1)
%	cenner_phase = base_phase * i;
%	plus_delta = cenner_phase + delta_phase;
%	minus_delta = cenner_phase - delta_phase;
%	decoded = find((Rx_decoded_phase <= plus_delta) & (Rx_decoded_phase > minus_delta));
%	Rx_decoded_symbols(decoded) = i;
%end
%
%Rx_serial_symbols = reshape(Rx_decoded_symbols',1,size(Rx_decoded_symbols,1)*size(Rx_decoded_symbols,2));
%
%%Convert symbol to Binary Stream
%for i = bit_per_symbol:-1:1
%	if i~=1 %i is not equal to 1
%		Rx_binary_matrix(i,:) = rem(Rx_serial_symbols,2);
%		Rx_serial_symbols = floor(Rx_serial_symbols/2);
%	else
%		Rx_binary_matrix(i,:) = Rx_serial_symbols;
%	end
%end
%baseband_in = reshape(Rx_binary_matrix,1,size(Rx_binary_matrix,1)*size(Rx_binary_matrix,2));
%%Calculate the error bit
%bit_errors = find(baseband_in ~= baseband_out);
%bit_error_count = size(bit_errors,2);
%bit_error_rate = bit_error_count/(size(baseband_in,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Below is Recieve witH CP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recieve
Rx_data_matrix_cp = reshape(Rx_data_cp,(IFFT_Length+cp_count),symbols_per_carrier+1);
Rx_data_matrix_decp = zeros(IFFT_Length,symbols_per_carrier+1);
Rx_data_matrix_decp = Rx_data_matrix_cp((1+cp_count):(IFFT_Length+cp_count),:);

Rx_spectrum_decp = fft(Rx_data_matrix_decp);% FFT to the colomn

figure(12)
stem(0:IFFT_Length-1,abs(Rx_spectrum_decp(1:IFFT_Length,2)),'b*-')
grid on
axis([0 IFFT_Length -0.5 1.5])
ylabel('Magnitude')
xlabel('FFT Bin')
title('OFDM Recieve Spectrum, Magnitude')

figure(13)
plot(0:IFFT_Length-1,(180/pi)*angle(Rx_spectrum_decp(1:IFFT_Length,2)),'go')
hold on
stem(carrier-1,(180/pi)*angle(Rx_spectrum_decp(carrier,2)),'b*-')
stem(conjugate_carriers-1,(180/pi)*angle(Rx_spectrum_decp(conjugate_carriers,2)),'b*-')
axis([0 IFFT_Length -200 +200])
grid on
ylabel('Phase(degree)')
xlabel('FFT Bin')
title('OFDM Recieve Spectrum, Phase')

%Distribution of Recieved Signal
figure(14)
Rx_carriers_decp = Rx_spectrum_decp(carrier,:)';
Rx_phase_P_decp = angle(Rx_carriers_decp);
Rx_mag_P_decp = abs(Rx_carriers_decp);
polar(Rx_phase_P_decp,Rx_mag_P_decp,'bd');

Rx_phase_decp = angle(Rx_carriers_decp)*(180/pi);
Rx_phase_decp = rem(Rx_phase_decp+360,360);

Rx_decoded_phase_decp = diff(Rx_phase_decp);
Rx_decoded_phase_decp = rem(Rx_decoded_phase_decp+360,360);

% Judgement
base_phase = 360/2^bit_per_symbol;
delta_phase = base_phase/2;
Rx_decoded_symbols_decp = zeros(size(Rx_decoded_phase_decp,1),size(Rx_decoded_phase_decp,2));

for i = 1:(2^bit_per_symbol - 1)
	cenner_phase = base_phase * i;
	plus_delta = cenner_phase + delta_phase;
	minus_delta = cenner_phase - delta_phase;
	decoded_decp = find((Rx_decoded_phase_decp <= plus_delta) & (Rx_decoded_phase_decp > minus_delta));
	Rx_decoded_symbols_decp(decoded_decp) = i;
end

Rx_serial_symbols_decp = reshape(Rx_decoded_symbols_decp',1,size(Rx_decoded_symbols_decp,1)*size(Rx_decoded_symbols_decp,2));

%Convert symbol to Binary Stream
for i = bit_per_symbol:-1:1
	if i~=1 %i is not equal to 1
		Rx_binary_matrix_decp(i,:) = rem(Rx_serial_symbols_decp,2);
		Rx_serial_symbols_decp = floor(Rx_serial_symbols_decp/2);
	else
		Rx_binary_matrix_decp(i,:) = Rx_serial_symbols_decp;
	end
end
baseband_in_decp = reshape(Rx_binary_matrix_decp,1,size(Rx_binary_matrix_decp,1)*size(Rx_binary_matrix_decp,2));
%Calculate the error bit
bit_errors_decp = find(baseband_in_decp ~= baseband_out);
bit_error_count_decp = size(bit_errors_decp,2);
bit_error_rate_decp = bit_error_count_decp/(size(baseband_in_decp,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%