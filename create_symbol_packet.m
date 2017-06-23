function [data_bitsIn ,s, train_seq] = create_symbol_packet(N,M,qam)
%create_symbol_packet
%  Create bit stream and convert it to symbols for a specific modulation

% Generate vector of trainning bits to be transmitted
train_bitsIN = randi([0 1],log2(qam)*M,1);         

% Generate vector of data bits to be transmitted
data_bitsIn = randi([0 1],log2(qam)*N,1);                                          

% match the train sequence bits to QAM symbols

train_bitsInMatrix = reshape(train_bitsIN,M,log2(qam));

x=bi2de(train_bitsInMatrix);                                                              

train_seq=qammod(x,qam); % the trainning sequence symbols to be transmmitted

% match the data sequence bits to QAM symbols

data_bitsInMatrix = reshape(data_bitsIn,N,log2(qam));

x=bi2de(data_bitsInMatrix);                                                              

s=qammod(x,qam); % the trainning sequence symbols to be transmmitted

end

