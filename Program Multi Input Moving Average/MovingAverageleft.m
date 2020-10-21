% Filter MA Left
% nama: Jihan
% prodi: S1 TT 
 
function y =  MovingAverageleft(data, win_size)
 
% untuk menambahkan elemen di kiri
wpad = win_size-1;
data2 = padarray(data, [0, wpad], 'replicate', 'pre');
% buffer
y = zeros(1, length(data));
% proses MA 
for i = win_size : length(data2)
    y(i - win_size + 1) = mean(data2((i - win_size + 1) : i), 'omitnan');
end