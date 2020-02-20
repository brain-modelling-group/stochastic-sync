function [sig_symbol, sig_text] = get_sig_star(pval)

if pval > 0.05
    sig_symbol = '';
    sig_text = 'ns';
elseif pval > 0.01 && pval <= 0.05
    sig_symbol = '*';
    sig_text = '\leq0.05';
elseif pval > 0.001 && pval <= 0.01
    sig_symbol = '**';
    sig_text = '\leq0.01';
elseif pval > 0.0001 && pval <= 0.001
    sig_symbol = '***';
    sig_text = '\leq0.001';
elseif pval <= 0.0001
    sig_symbol = '****';
    sig_text = '\leq0.0001';
end