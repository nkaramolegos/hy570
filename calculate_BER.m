function [BER_LS, BER_LMS, BER_RLS] = calculate_BER(loops, qam, L, N, M, CP_len, channel_var, SNR_vector, mu, delta, lambda)
%calculate_BER Calculate the Bit Error Rate (BER) using the channel
%estimation made by each one of the adaptive algorithm

P=length(SNR_vector);
num_bits_wrong_LS=zeros(P,1);
num_bits_wrong_LMS=zeros(P,1);
num_bits_wrong_RLS=zeros(P,1);

for j=1:P
        noise_var= (2/3*(qam-1))./  (10^ (SNR_vector(j)/10)); % for QAM modulation, note that real and image part of the noise have noise_var/2. So the final noise is noise_var

        for NumPackets=1:loops

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


            % the desired responce is the received signal for the trainning symbols
            d=y(CP_len+1:M+CP_len);
            % input to the filter
            u=train_seq;

            %% LS
                                                                                                                                 
            h_LS = LS(u, d, L);

            H_est = fft(h_LS,N);

            H_est=reshape(H_est,N,1);

            % using OFDM, we convert a wideband channel into a set of N parallel
            % narrowband channels. As a result, no complex equalization is required
            Y_data= fft(y(M+2*CP_len+1:end),N)*(1/sqrt(N));

            % cancel the effect of the channel divide with the conjugate and the abs^2
            % symbol-by-symbol decision for each information symbol (no complex
            % equalization as mentioned)
            r= Y_data.*conj(H_est)./(abs(H_est).^2);


            num_bits_wrong_LS(j) = num_bits_wrong_LS(j) + demodulate(r, qam, N, data_bitsIn);

            %% LMS
            if  ( mu < 2/(norm(u).^2) )
                [e_lms, h_lms] = LMS(d, u, mu, L);

                H_est = fft(h_lms,N);

                H_est=reshape(H_est,N,1);

                % using OFDM, we convert a wideband channel into a set of N parallel
                % narrowband channels. As a result, no complex equalization is required
                Y_data= fft(y(M+2*CP_len+1:end),N)*(1/sqrt(N));

                % cancel the effect of the channel divide with the conjugate and the abs^2
                % symbol-by-symbol decision for each information symbol (no complex
                % equalization as mentioned)
                r= Y_data.*conj(H_est)./(abs(H_est).^2);


                num_bits_wrong_LMS(j) = num_bits_wrong_LMS(j) + demodulate(r, qam, N, data_bitsIn);
            else
                disp('!!!......................Voitheia......................!!!');
            end

             %% RLS

            [xi_rls, h_rls] = RLS(d, u, delta, lambda, L);

            H_est = fft(h_rls,N);

            H_est=reshape(H_est,N,1);

            % using OFDM, we convert a wideband channel into a set of N parallel
            % narrowband channels. As a result, no complex equalization is required
            Y_data= fft(y(M+2*CP_len+1:end),N)*(1/sqrt(N));

            % cancel the effect of the channel divide with the conjugate and the abs^2
            % symbol-by-symbol decision for each information symbol (no complex
            % equalization as mentioned)
            r= Y_data.*conj(H_est)./(abs(H_est).^2);


            num_bits_wrong_RLS(j) = num_bits_wrong_RLS(j)+ demodulate(r, qam, N, data_bitsIn);


        end
        
BER_LS(j)=num_bits_wrong_LS(j)/(log2(qam)*N*loops);
BER_LMS(j)=num_bits_wrong_LMS(j)/(log2(qam)*N*loops);
BER_RLS(j)=num_bits_wrong_RLS(j)/(log2(qam)*N*loops);       
        
end

end

