function x = qpsk_generation(x0, beta, sym_num, samp_per_sym)
% qpsk generation: row vector in, row vector out

rct_filt = comm.RaisedCosineTransmitFilter('Shape','Normal',...
    'RolloffFactor',beta,...
    'FilterSpanInSymbols',sym_num,...
    'OutputSamplesPerSymbol',samp_per_sym);
x = step(rct_filt,x0')';

end