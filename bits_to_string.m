%% EGB342 Assignment 2B
%% str = bits_to_string(b)
% Converts an array of bits into a 7-bit ASCII character string.

function str = bits_to_string(b)

if rem(numel(b), 7)
    error('Input length must be a factor of 7');
end

b = reshape(b > 0, 7, numel(b)/7)';
str = char(bi2de(b, 'left-msb')');
