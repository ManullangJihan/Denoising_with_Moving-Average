% Program utama untuk filter MA terhadap sinyal PCG - single input
% nama: Tuah Jihan
% prodi: S1 TT

% Environment
warning off;
clear all;
close all;
clc;

%% Memilih folder untuk menyimpan

[fname, pname] = uigetfile('*.wav', 'Pilih sebuah data PCG');

if ~isequal(fname, 0) || ~isequal(pname, 0)
    
    %% Import data
    pcgfile = fullfile(pname, fname);
    [x, fs] = audioread(pcgfile);
    fprintf('Processing: %s\n', fname);
    
    % karena datanya stereo, pilih 1 data saja
    M = 5;
        
    %% Index data selection
    % Gunakan data ntuk t detik saja
    t1 = 1 / fs;
    t2 = length(x) / fs;
    N1 = round(t1 * fs);
    N2 = round(t2 * fs);
    x = x(N1 : N2-1);
    
    %% Index Data dengan menggunakan Panjang Data
    % x = x(1:50000);
    
    %% Preprocessing
    % normalisasi data mentah agar berada pada -1 hingga +1 volt
    x = x ./ max(abs(x));
    
    % centering
    x = x - mean(x);

    % Noise added to signal with certain SNR value
    % tambahkan noise acak N(0,1)
    snrawgn = 5;
    datan = awgn(x, snrawgn, 'measured'); %Input Signal+Noise
    xnoise = x+datan;

    %% MA Processing
    % % MA left
    y  = MovingAverageleft(xnoise, M);

    % MA Right
    % y  = MovingAverageright(xnoise, M);

    % % MA Symmetri
    % y = MovingAveragesym(xnoise, M);
    
      
    %% Post processing - needed for performance measurement
    % transpose data
    yt = y';
    
    % centering output
    yt = yt - mean(yt);
    
    % normalisasi output
    y = yt ./ max(abs(yt));
       
    
    %% Analisis Parameter
    % Menghitung error dengan MSE, SNR dan RMSE
    
    % Hitung MSE
    err1 = (norm(x(:)-y(:),2).^2)/numel(x);
    fprintf('>> The Mean-squared Error is %0.4f\n', err1);

    % Hitung SNR
    noiseampestimation = x-xnoise;
    snr1 = 20*log10(rms(x)/rms(noiseampestimation));
    fprintf('>> The Signal Noise to ratio is %0.4f\n', snr1);
    
    % Hitung RMSE
    RMSE = sqrt(err1);
    fprintf('>> The RMSE is %0.4f\n', RMSE);

    
    %% Menampilkan hasil setiap langkah
    addpath('./plots');
    
    outfolder = 'Output Plots';
    if ~exist(outfolder, 'dir')
        mkdir(outfolder);
    end
    sname = fname(1:length(fname)-4);
    
    %Plot Sinyal Asli
    foname = sprintf('%s_pcgAsli.jpg', sname);
    onam1 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(x))/fs, x,'LineWidth', .5);
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal PCG - \it{Ground Truth}');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam1, '-jpg', '-r200', '-a4', '-painters', '-transparent');
     
    % Sinyal Ditambah Noise AWGN
    foname = sprintf('%s_pcgNoisy_orde%d_%s.jpg', sname, M);
    onam2 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(xnoise))/fs, xnoise,'LineWidth', .5);
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal PCG - Tercampur Noise Acak');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
     
    % Denoised Signal
    foname = sprintf('%s_pcgDenoised_orde%d_%s.jpg', sname, M);
    onam3 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(y))/fs, y,'LineWidth', .5);
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal PCG - Hasil Denoising');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam3, '-jpg', '-r200', '-a4', '-painters', '-transparent');
     
    % Noisy Signal vs Denoised Signal
    foname = sprintf('%s_pcgAslivsDenoised_orde%d_%s.jpg', sname, M);
    onam4 = fullfile(outfolder, foname);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    plot((1:length(xnoise))/fs, xnoise,'LineWidth', .5); hold on;
    plot((1:length(y))/fs, y,'r','LineWidth', 2); hold off;    
    xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Sinyal Noise vs Hasil \it{Denoising}');
    legend('Noisy Signal','Denoised Signal')
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam4, '-jpg', '-r200', '-a4', '-painters', '-transparent');
        
    % Frequency Response
    foname = sprintf('%s_pcgFreqResp_orde%d_%s.jpg', sname, M);
    onam5 = fullfile(outfolder, foname);
    h = ones(1, M) / M;
    [H, W] = freqz(h, 1, 1024);
    figure;
    plot(W/pi, 20*log10(abs(H)), 'LineWidth', .5);
    xlabel('Frekuensi Ternormalisasi', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel('Magnitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title('Tanggapan Frekuensi Filter');
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    export_fig (onam5, '-jpg', '-r200', '-a4', '-painters', '-transparent');
end

% % Menghitung tanggapan frekuensi
% Frequency Respone
% k = 0 : pi / 400 : pi;
% H5 = (1/M) * (1-exp(-1i*k*M))./(1-exp(-1i*k)); 
% H10 = (1/(2*M)) * (1-exp(-1i*k*2*M))./(1-exp(-1i*k)); 
% H15 = (1/(3*M)) * (1-exp(-1i*k*3*M))./(1-exp(-1i*k));
  
%% Menampilkan hasil setiap langkah
  
% figure;
% ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
% ax2.ActivePositionProperty = 'position';
% plot(k,[abs(H5);abs(H10);abs(H15)]);
% xlabel('Frekuensi Ternormalisasi', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
% title('Sinyal PCG - Hasil Denoising');
% set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
% set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
% export_fig (onam3, '-jpg', '-r200', '-a4', '-painters', '-transparent');
% axis([0, pi, 0, 1]);
% legend('M=5','M=10','M=15');
% hold off;

%% END