function [MSE_LMS, MSE_NLMS, MSE_RLS] = calculate_MSE(loops, qam, L, N, M, CP_len, channel_var, SNR, mu, delta, lambda, epsilon, beta)
%calculate_MSE Calculate the Mean Squared Error (MSE) using the channel
%estimation made by each one of the adaptive algorithm

mse_mat_LMS = zeros(loops,M);
mse_mat_NLMS = zeros(loops,M);
mse_mat_RLS = zeros(loops,M);

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
            [y, h] = transmit_packet(x, SNR, channel_var, L, N, M, CP_len);


            % the desired responce is the received signal for the trainning symbols
            d=y(CP_len+1:M+CP_len);
            % input to the filter
            u=train_seq;
            
            %% LMS

            [e_lms, h_lms] = LMS(d, u, mu, L);
            
            %% NLMS

            [e_nlms, h_nlms] = NLMS(d, u, beta, L, epsilon);
            
            %% RLS

            [xi_rls, h_rls] = RLS(d, u, delta, lambda, L);
            
            %% gather MSE data
               
            mse_mat_LMS(NumPackets,:) = e_lms;
            mse_mat_NLMS(NumPackets,:) = e_nlms;
            mse_mat_RLS(NumPackets,:) = xi_rls;
        end

MSE_LMS = mean( abs(mse_mat_LMS).^2 );
MSE_NLMS = mean( abs(mse_mat_NLMS).^2 );
MSE_RLS = mean( abs(mse_mat_RLS).^2 );


