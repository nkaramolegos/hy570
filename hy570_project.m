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
M=256; 
% taps of the channel and length of cyclic prefix
CP_len=8;
L=1;
channel_var=1;
% white gaussian noise var (on receiver)
noise_var=0.0;

% for LMS
% we saw experimentally that mu should be in space 0<mu<0.05 (Conservative
% approach)
mu=0.01; 

% for NLMS
epsilon=10^-4;
beta=0.8;

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


% NLMS
[e_nlms, h_nlms] = NLMS(d, u, beta, L, epsilon);
h_nlms


% RLS
[xi_rls, h_rls] = RLS(d, u, delta, lambda, L);
h_rls

figure
plot(1:M, abs(e_lms),1:M, abs(e_nlms),'r', 1:M, abs(xi_rls),'g');
grid on;
title('Error abs');
xlabel('Iteration');
ylabel('Error');
legend('LMS', 'NLMS', 'RLS');

figure
t(1)=subplot(3,2,1);
stem(1:L,abs(h),'c','LineWidth',1.5);
title('Theoretical');
xlim([0 5]);
grid on;
t(2)=subplot(3,2,2);
stem(1:L,abs(h_LS),'b','LineWidth',1.5);
title('LS');
xlim([0 5]);
grid on;
t(3)=subplot(3,2,3);
stem(1:L,abs(h_lms),'r','LineWidth',1.5);
title('LMS');
xlim([0 5]);
grid on;
t(4)=subplot(3,2,4);
stem(1:L,abs(h_nlms),'k','LineWidth',1.5);
title('NLMS');
xlim([0 5]);
grid on;
t(5)=subplot(3,2,5);
stem(1:L,abs(h_rls),'g','LineWidth',1.5);
title('RLS');
xlim([0 5]);
grid on;

pos = get(t,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(t(5),'Position',[new,pos{end}(2:end)])

% return;

% % BER

loops=1000;
SNRdb_min=-10;
SNRdb_max=20;
step_db=2;
SNR_vector = SNRdb_min:step_db:SNRdb_max;

% flat fading channel
% disp('Flat fading channel');
% L=1;
% [BER_LS, BER_LMS, BER_NLMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda, epsilon, beta);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector, BER_NLMS,'k', SNR_vector,BER_RLS,'g');
% grid on; 
% title ('BER for Flat fading channel'); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('NLMS \\beta=%g',beta), sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));


% % frequency selective channel with 2 taps
% L=2;
% disp('Frequency selective channel L=2');
% [BER_LS, BER_LMS, BER_NLMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda, epsilon, beta);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector, BER_NLMS,'k', SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('NLMS \\beta=%g',beta), sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));
% 
% frequency selective channel with 4 taps
% L=4;
% disp('Frequency selective channel L=4');
% [BER_LS, BER_LMS, BER_NLMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda, epsilon, beta);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector, BER_NLMS,'k', SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('NLMS \\beta=%g',beta), sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));

% 
% 
% % frequency selective channel with 4 taps but different mu
% mu=0.05;
% L=4;
% disp('Frequency selective channel L=4 for different mu');
% [BER_LS, BER_LMS, BER_NLMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda, epsilon, beta);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector, BER_NLMS,'k', SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('NLMS \\beta=%g',beta), sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));
% 
% 
% % frequency selective channel with 4 taps but different delta
% mu=0.01;
% delta=10^-2;
% L=4;
% disp('Frequency selective channel L=4 for different delta');
% [BER_LS, BER_LMS, BER_NLMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda, epsilon, beta);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector, BER_NLMS,'k', SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('NLMS \\beta=%g',beta), sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));
% 
% % frequency selective channel with 4 taps but different lambda
% mu=0.01;
% delta=10^-4;
% lambda=0.3;
% L=4;
% disp('Frequency selective channel L=4 for different lambda');
% [BER_LS, BER_LMS, BER_NLMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda, epsilon, beta);
% 
% figure;
% semilogy(SNR_vector,BER_LS,'b',SNR_vector,BER_LMS,'r',SNR_vector, BER_NLMS,'k', SNR_vector,BER_RLS,'g');
% grid on; 
% title (sprintf('BER for frequency selective channel with %g taps',L)); 
% xlabel ('SNR in dB');
% ylabel('Bit error probability');
% legend('LS',sprintf('LMS \\mu=%g',mu),sprintf('NLMS \\beta=%g',beta), sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));

for SNR=-10:10:10

    [MSE_LMS, MSE_NLMS, MSE_RLS] = calculate_MSE(loops, qam, L, N, M, CP_len, channel_var, SNR, mu, delta, lambda, epsilon, beta);

    figure;
    semilogy(1:M, MSE_LMS,'r',1:M, MSE_NLMS,'k',1:M, MSE_RLS,'g');
    grid on; 
    title (sprintf('MSE for SNR=%d',SNR));
    ylabel ('MSE');
    xlabel('Iteration')
    legend(sprintf('LMS \\mu=%g',mu),sprintf('NLMS \\beta=%g',beta), sprintf('RLS \\delta=%g, \\lambda=%g',delta, lambda));
end
