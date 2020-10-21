% Filter MA Right
% nama: Jihan
% prodi: S1 TT 
 
function y = MovingAverageright(data, win_size)
 
% untuk menambahkan elemen di kanan data
wpad = win_size-1
data2 = padarray(data, [0, wpad], 'replicate', 'post');
% buffer
y = zeros(1, length(data));
% Proses MA
for i = 1 : length(data2) - win_size + 1
    y(i) = mean(data2(i : (i + win_size - 1)), 'omitnan');
end