clc
clear all
close all

% read the carrier image and secret image

image = imread('lena.jpg');       % read carrier image
secret_image = imread('carrier_1.bmp');    % reand the secret
[M, N, C] = size(image);
t = 350;  % threshold of the embedded watermark

if ~(M == 512 && N == 512 && C == 3)
    error('input size required: 512x512x3');
end

% embed the watermark

[carrier_embedded] = embed_secret(image,secret_image,t);
subplot(1,3,1);
imshow(image);
title('original image');
subplot(1,3,2);
imshow(secret_image);
title('secret image');
subplot(1,3,3);
imshow(carrier_embedded);
title('after embedding');

imwrite(carrier_embedded, 'carrier_embedded.jpg', 'jpg');

function [imagemark] = embed_secret(carrier_image, secret_image, t)
    % Convert the carrier image to YCbCr color space
    YCbCr = rgb2ycbcr(carrier_image);
    
    % Extract the Cb channel and convert it to double
    Cb = double(YCbCr(:,:,2));
    
    % Apply 2D discrete wavelet transform (DWT) to the Cb channel
    [cA, cB, cC, cD] = dwt2(Cb, 'haar');
    
    % Get the size of cA matrix
    [M, N, ~] = size(cA);
    
    % Get the size of the secret image
    [wM, wN] = size(secret_image);

    % Prepare the partition vectors for dividing cA into cells
    x = uint8(4 * ones(1, M / 4));
    y = uint8(4 * ones(1, N / 4));
    
    % Divide cA into a cell array cA_
    cA_ = mat2cell(cA, x, y);

    % Iterate over each cell in cA_
    for i = 1:M/4
        for j = 1:N/4
            % Check if i and j are within the range of the secret image
            if i <= wM && j <= wN
                % Apply the Schur decomposition to the current cell
                [U, T] = schur(cA_{i, j});
                
                % Find the index of the maximum value in T
                [m, n] = find(T == max(max(T)));
                
                % Get the maximum value from T
                Tmax = T(m, n);
                
                % Adjust the maximum value based on the corresponding pixel in the secret image
                if secret_image(i, j) == 1
                    Tmax = Tmax - mod(Tmax, t) + 0.75 * t;
                else
                    Tmax = Tmax - mod(Tmax, t) + 0.25 * t;
                end
                
                % Update the maximum value in T
                T(m, n) = Tmax;
                
                % Reconstruct the cell using the updated T matrix
                cA_{i, j} = U * T * U';
            end
        end
    end

    % Convert the cell array cA_ back to a matrix cA
    cA = cell2mat(cA_);
    
    % Perform inverse DWT to obtain the modified Cb channel
    Cbmark = uint8(idwt2(cA, cB, cC, cD, 'haar'));
    
    % Replace the original Cb channel in YCbCr with the modified Cb channel
    YCbCr(:,:,2) = Cbmark;
    
    % Convert YCbCr image back to RGB color space
    imagemark = ycbcr2rgb(YCbCr);
end


