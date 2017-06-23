function [y, h] = transmit_packet(x, noise_var, channel_var, L, N, M, CP_len)
%transmit_packet
%  Transmits the packet and return the output on receiver

% create channel                                                         
h=sqrt(channel_var/2)*randn(L,1) + 1i*randn(L,1)*sqrt(channel_var/2);

% convolution with the channel
y=conv(h,x);

% create additive white gaussian noise on receiver  
w=sqrt(noise_var/2)*randn(length(y),1) + 1i*sqrt(noise_var/2)*randn(length(y),1); 

% add noise to receiver
y=y+w;
% keep the useful information of the convolution, discard the last ones
y=y(1:N+M+2*CP_len);
end

