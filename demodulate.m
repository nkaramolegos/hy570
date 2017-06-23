function num_bits_wrong = demodulate(r, qam, N, data_bitsIn)
%demodulate returns the number of bits that are wrong

dataSymbolsOut = qamdemod(r,qam);
dataOutMatrix = de2bi(dataSymbolsOut,log2(qam));

train_bitsInMatrix = reshape(dataOutMatrix,log2(qam)*N,1);

num_bits_wrong=nnz(data_bitsIn -train_bitsInMatrix);

end

