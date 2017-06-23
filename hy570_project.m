close all;
clear;
clc;

% CHANNEL ESTIMATION USING ADAPTIVE FILTERING AND OFDM MODULATION         %

% Gkiolias A. 940, Karamolegkos N. 943                                    %
% HY 570 Project                                                          %
% Advisor: Mouchtaris A.
% June 2017                                                               %

% OFDM approach based on books of 
% 1) D. Tse, and P. Viswanath, Fundamentals of wireless communications,   %
%    Cambridge university press, August 13, 2004                          %
% 2) G. Stuber and Y. Li, Orthogonal frequency division multiplexing      %
%    for wireless communications, Springer, 2006                          %

% Adaptive filtering base on book of
% S. Haykin, Adaptive Filter Theory, 1996                                 %


% %                   initializations                                   % %

qam=4;
% number of SYMBOLS to be transmitted, should be multiple of log2(qam)
N=512;
% number of training seq SYMBOLS, should be multiple of log2(qam) 
M=10; 
% taps of the channel and length of cyclic prefix
CP_len=8;
L=4;
channel_var=1;
% white gaussian noise var (on receiver)
noise_var=0.0;

% for LMS
mu=0.01;

% for RLS
delta=10^-4;
lambda=0.98;

% %                   end of initializations                             %%
                                                                                   
% create the packet and the trainning sequence
[data_bitsIn ,s_tilda, train_seq_tilda] = create_symbol_packet(N,M,qam);

% inverse fft to train sequence symbols 
train_seq=sqrt(M)*ifft(train_seq_tilda);

% inverse fft to data symbols 
s=sqrt(N)*ifft(s_tilda);

% add to CP of length CP_len > L (taps of the channel) to trainnig sequence
% and to data symbols
train_seq_cp=[train_seq(M-CP_len+1:M); train_seq];
s_cp=[s(N-CP_len+1:N); s];

% create OFDM frame
x=[train_seq_cp; s_cp];

% transmit
[y, h] = transmit_packet(x, noise_var, channel_var, L, N, M, CP_len);
h


% the desired responce is the received signal for the trainning symbols
d=y(CP_len+1:M+CP_len);
% input to the filter
u=train_seq;

% LS
                                                                                                                                              
h_LS = LS(u, d, L)

H_est = fft(h_LS,N);

H_est=reshape(H_est,N,1);

% using OFDM, we convert a wideband channel into a set of N parallel
% narrowband channels. As a result, no complex equalization is required
Y_data= fft(y(M+2*CP_len+1:end),N)*(1/sqrt(N));

% cancel the effect of the channel divide with the conjugate and the abs^2
% symbol-by-symbol decision for each information symbol (no complex
% equalization as mentioned)
r= Y_data.*conj(H_est)./(abs(H_est).^2);

% scatterplot(r);

num_bits_wrong = demodulate(r, qam, N, data_bitsIn)

% LMS
[e_lms, h_lms] = LMS(d, u, mu, L);
h_lms

% figure
% plot(abs(e_lms));
% title('LMS error');

% RLS


[xi_rls, h_rls] = RLS(d, u, delta, lambda, L);
h_rls
% 
% figure
% plot(abs(xi_rls));
% title('RLS error');

% % BER

loops=1000;
SNRdb_min=-10;
SNRdb_max=20;
step_db=2;
SNR_vector = SNRdb_min:step_db:SNRdb_max;

% % flat fading channel
% disp('Flat fading channel');
% L=1;
% [BER_LS, BER_LMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector,BER_RLS,'g');
% grid on;
% title ('BER for flat fading channel'); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));
% 
% % frequency selective channel with 2 taps
% L=2;
% disp('Frequency selective channel L=2');
% [BER_LS, BER_LMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));

% frequency selective channel with 4 taps
L=4;
disp('Frequency selective channel L=4');
[BER_LS, BER_LMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda);

figure;
semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector,BER_RLS,'g');
grid on; 
title (sprintf('BER for frequency selective channel with %g taps',L)); 
xlabel ('SNR in dB');
ylabel('Bit error probability');
legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));

% % frequency selective channel with 4 taps but different mu
% mu=0.001;
% L=4;
% disp('Frequency selective channel L=4 for different mu');
% [BER_LS, BER_LMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));
% 
% % frequency selective channel with 4 taps but different delta
% mu=0.01;
% delta=10^-2;
% L=4;
% disp('Frequency selective channel L=4 for different delta');
% [BER_LS, BER_LMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));
% 
% % frequency selective channel with 4 taps but different lambda
% mu=0.01;
% delta=10^-4;
% lambda=0.3;
% L=4;
% disp('Frequency selective channel L=4 for different lambda');
% [BER_LS, BER_LMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));




