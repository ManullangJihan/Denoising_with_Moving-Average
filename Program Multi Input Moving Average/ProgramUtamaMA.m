% Program utama Moving Average Multi Input terhadap sinyal PCG
% nama: Tuah Jihan
% prodi: S1 TT 

% Environment 
warning off;
clear all;
close all;
clc;

%% Memilih folder untuk menyimpan
direk =  uigetdir('Choose a folder where you store the data');

if ~isequal(direk, 0)
    
    Nfiles = dir(fullfile(direk, '*.wav'));
    M = 5;
    for id = 5
        fprintf('Moving Average symmetric 5 point');
        for ix = 1: numel(Nfiles);
            
            % Import data ke Matlab
            fname = Nfiles(ix).name;
            dname = fullfile(direk, fname);
            [x, fs] = audioread(dname);
            
            fprintf('%d) File: %s', ix, fname);
            
            % karena datanya stereo, pilih 1 data saja
            data = x(:, 1)';
            
            % Gunakan data ntuk t detik saja
            % t1 = 1 / fs;
            % t2 = length(x) / fs;
            % N1 = round(t1 * fs);
            % N2 = round(t2 * fs);
            % x = x(N1 : N2-1);
            x = x(1:100);
            
            % Normalisasi data mentah agar berada pada -1 hingga +1 volt
            x = x ./ max(abs(x));
            
            % centering
            x = x - mean(x);
            
            % tambahkan noise acak N(0,1)
            % datan = wgn(length(data), 1, 0)';
            snrawgn = 5;
            datan = awgn(x, snrawgn, 'measured'); %Input Signal+Noise
            xnoise = x+datan;
            
            %% Proses Moving Average
            % MA left
            y  = MovingAverageleft(xnoise, M);

            % % MA Right
            % y  = MovingAverageright(xnoise, M);

            % % MA Symmetri
            % y = MovingAveragesym(xnoise, M);

            %% Post Processing
            % transpose data
            yt = y';

            % Centering output
            yt = yt - mean(yt);

            % Normalisasi output
            y = yt ./ max(abs(yt));
            
    %% Analisis Parameter
    % Menghitung error dengan MSE, SNR dan PSNR
    
    % Hitung MSE
    err1 = (norm(x(:)-y(:),2).^2)/numel(x);
    fprintf('>> The Mean-squared Error is %0.4f\n', err1);

    % Hitung SNR
    noiseampestimation = x-y;
    snr1 = 20*log10(rms(x)/rms(noiseampestimation));
    fprintf('>> The Signal Noise to ratio is %0.4f\n', snr1);
    
    % Hitung RMSE
    RMSE = sqrt(err1);
    fprintf('>> The RMSE is %0.4f\n', RMSE);
    
%     %% Menampilkan hasil setiap langkah
%     addpath('./plots');
%     
%     outfolder = 'Output Plots';
%     if ~exist(outfolder, 'dir')
%         mkdir(outfolder);
%     end
%     sname = fname(1:length(fname)-4);
%     
%     % Plot Sinyal Asli
%     foname = sprintf('%s_pcgAsli.jpg', sname);
%     onam1 = fullfile(outfolder, foname);
%     figure;
%     ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
%     ax2.ActivePositionProperty = 'position';
%     plot((1:length(x))/fs, x,'LineWidth', .5);
%     xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     title('Sinyal PCG - \it{Ground Truth}');
%     set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
%     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
%     export_fig (onam1, '-jpg', '-r200', '-a4', '-painters', '-transparent');
%    
%     % Sinyal Kena Noise
%     foname = sprintf('%s_pcgNoisy_%d.jpg', sname);
%     onam2 = fullfile(outfolder, foname);
%     figure;
%     ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
%     ax2.ActivePositionProperty = 'position';
%     plot((1:length(xnoise))/fs, xnoise,'LineWidth', .5);
%     xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     title('Sinyal PCG - Tercampur Noise Acak');
%     set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
%     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
%     export_fig (onam2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
%     
%     % Denoised Signal
%     foname = sprintf('%s_pcgDenoised_%d.jpg', sname);
%     onam3 = fullfile(outfolder, foname);
%     figure;
%     ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
%     ax2.ActivePositionProperty = 'position';
%     plot((1:length(y))/fs,y,'LineWidth', .5);
%     xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     title('Sinyal PCG - Hasil Denoising');
%     set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
%     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
%     export_fig (onam3, '-jpg', '-r200', '-a4', '-painters', '-transparent');
%     
%     % Sinyal Noise vs Sinyal Denoised
%     foname = sprintf('%s_pcgAslivsDenoised_%d.jpg', sname);
%     onam4 = fullfile(outfolder, foname);
%     figure;
%     ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
%     ax2.ActivePositionProperty = 'position';
%     plot((1:length(xnoise))/fs,xnoise,'LineWidth', .5); hold on;
%     plot((1:length(y))/fs,y,'LineWidth', 2); hold off;    
%     xlabel('Waktu (sec)', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     ylabel('Amplitudo', 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
%     title('Sinyal PCG vs Hasil \it{Denoising}');
%     set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 2, 'GridAlpha', 0.1);
%     set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
%     export_fig (onam4, '-jpg', '-r200', '-a4', '-painters', '-transparent');
%     legend('Noisy Signal','Denoised Signal')
    
%% Menampilkan Grafik MSE, SNR dan RMSE Sebelum Denoising
    % MSE
    %% Menampilkan Grafik Analisis Parameter
    addpath('./plots');
    bgdir = pwd;
    out_folder = 'PARAMETER';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('%s_mse all', numel(Nfiles));
    oname2 = fullfile(out_full, ot2);
    xlbl = 'Data number'; 
    ylbl = 'MSE'; 
    til1=('MSE Reduksi Noise AWGN')
    y = err1; 
    x = 1 : numel(Nfiles);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    for iix = 1 : length(err1)
        plot(x, y(iix), 'LineWidth', 1.5);
        hold on;
    end
    hold off;
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    
    % xy label
    xlabel(xlbl, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel(ylbl, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title(til1, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal')
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
  
    % SNR 
    bgdir = pwd;
    out_folder = 'PARAMETER';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('%s_snr after',numel(Nfiles));
    oname2 = fullfile(out_full, ot2);
    xlb1= 'Data number'; 
    ylb1 = 'SNR'; 
    til1=('SNR Sesudah Reduksi Noise AWGN');
    y = snr1; x = 1 : numel(Nfiles);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    for iix = 1 : length(err1)
        plot(x, y(iix), 'LineWidth', 1.5);
        hold on;
    end
    hold off;
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    % xy label
    xlabel(xlb1, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel(ylb1, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title(til1, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal')
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
    % RMSE
    bgdir = pwd;
    out_folder = 'PARAMETER';
    if ~exist(out_folder, 'dir');
        mkdir(out_folder);
    end
    out_full = fullfile(bgdir, out_folder);
    ot2 = sprintf('%s_snr after',numel(Nfiles));
    oname2 = fullfile(out_full, ot2);
    xlb1= 'Data number'; 
    ylb1 = 'RMSE'; 
    til1=('SNR Sesudah Reduksi Noise AWGN');
    y = RMSE; 
    x = 1 : numel(Nfiles);
    figure;
    ax2 = axes('Position',[0.14 0.17 0.78 0.74]);
    ax2.ActivePositionProperty = 'position';
    for iix = 1 : length(err1)
        plot(x, y(iix), 'LineWidth', 1.5);
        hold on;
    end
    hold off;
    pbaspect([3.6 2 1]);
    set(gca, 'XColor', 'black', 'YColor', 'black', 'LineWidth', 1.5, 'GridAlpha', 0);
    set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
    % xy label
    xlabel(xlb1, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    ylabel(ylb1, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal');
    title(til1, 'FontName', 'Rockwell', 'FontSize', 12, 'FontWeight', 'normal')
    export_fig (oname2, '-jpg', '-r200', '-a4', '-painters', '-transparent');
    close;
    
        end
    end
end
%% END