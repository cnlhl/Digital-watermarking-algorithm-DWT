clc
clear all
close all

% read the carrier image

carrier_embedded = imread('carrier_embedded.jpg');
t = 350;    % threshold of the embedded watermark

% extract the watermark

[secret_extract] = extract_secret( carrier_embedded,t);
figure;
subplot(1,2,1);
imshow(carrier_embedded);   
title('carrier embedded');
subplot(1,2,2);
imshow(secret_extract);      
title('extract watermark');

function [water_extract] = extract_secret(imagemark, t)
    % Convert the input image to YCbCr color space and extract the Cb channel
    YCbCrmark = rgb2ycbcr(imagemark);
    Cbmark = double(YCbCrmark(:,:,2));
    
    % Apply 2D discrete wavelet transform (DWT) to the Cb channel
    cAw = dwt2(Cbmark, 'haar');
    [M, N, ~] = size(cAw);
    
    % Prepare the partition vectors for dividing cAw into cells
    x = uint8(4 * ones(1, M / 4));
    y = uint8(4 * ones(1, N / 4));
    cAw_ = mat2cell(cAw, x, y);
    
    % Initialize the water extraction matrix
    water_extract = ones(64, 64);
    
    % Iterate over each cell in cAw_
    for i = 1:M/4
        for j = 1:N/4
            % Apply the Schur decomposition to the current cell
            [~, T] = schur(cAw_{i, j});
            
            % Find the maximum value in T
            Tmax = max(T(:));
            
            % Determine the extracted water bit based on Tmax and threshold t
            if mod(Tmax, t) < t/2
                water_extract(i, j) = 0;
            else
                water_extract(i, j) = 1;
            end
        end
    end
end
