% Program Filter MA Symmetry
% nama: Jihan
% prodi: S1 TT

function y = MovingAveragesym(data, win_size)

% untuk menambahkan elemen di kiri-kanan data
wpad = (win_size - 1)/2;
data2 = padarray(data, [0, wpad], 'replicate', 'both');

% buffer
y = zeros(1, length(data));
for i = (win_size - 1) / 2 + 1 : length(data2) - (win_size - 1) / 2;
    y(i-wpad) = mean(data2((i-(win_size-1)/2) : (i+((win_size-1)/2))), 'omitnan');
end