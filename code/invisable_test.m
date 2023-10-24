% Invisibility test

% Clear workspace variables and command window
clc
clear all
close all
% Image paths
carrierImages = {'carrier_1.jpg', 'carrier_2.jpg', 'carrier_3.jpg', 'carrier_4.jpg'};
embedImages = {'carrier_embed_1.jpg', 'carrier_embed_2.jpg', 'carrier_embed_3.jpg', 'carrier_embed_4.jpg'};

% Iterate over each pair of images
for i = 1:length(carrierImages)
    % Read images
    carrier = imread(carrierImages{i});
    embed = imread(embedImages{i});
    
    % Compute PSNR
    psnrValue = psnr(carrier, embed);
    fprintf('PSNR between %s and %s: %.2f dB\n', carrierImages{i}, embedImages{i}, psnrValue);
    
    % Compute SSIM
    ssimValue = ssim(carrier, embed);
    fprintf('SSIM between %s and %s: %.4f\n', carrierImages{i}, embedImages{i}, ssimValue);
end
