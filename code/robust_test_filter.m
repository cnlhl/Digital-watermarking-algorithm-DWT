clc
clear all
close all

carriers = {'carrier_1.jpg', 'carrier_2.jpg', 'carrier_3.jpg', 'carrier_4.jpg'};
t = 100;

water = imread('secret_fudan.bmp');

figure;
for carrierIndex = 1:numel(carriers)
    carrier = imread(carriers{carrierIndex});
    imagemark = embed_secret(carrier, water, t);
    water_extract = extract_water(imagemark, t);
    
    H = fspecial('gaussian', [5,5], 0.1);
    im_filter = imfilter(imagemark, H);
    
    subplot(2, numel(carriers), carrierIndex);
    imshow(im_filter);
    title('image-filter');
    
    water_extract_filter = extract_water(im_filter, t);
    subplot(2, numel(carriers), carrierIndex+4);
    imshow(water_extract_filter);
    title('water-extract-filter');
    [CDR] = cdr(water, water_extract_filter);
    disp('filter CDR:');
    disp(CDR);
end



function [ CDR ] = cdr( X,Y)
[MX,NX] = size(X);
sum =0.0;
for i = 1:1:MX
    for j = 1:1:NX
        if X(i,j) == Y(i,j)
            sum = sum +1;
        end
    end
end
CDR = sum/MX/NX;
end

function [water_extract] = extract_water(imagemark, t)
    YCbCrmark = rgb2ycbcr(imagemark);
    Cbmark = double(YCbCrmark(:,:,2));
    cAw = dwt2(Cbmark, 'haar');
    [M, N, ~] = size(cAw);
    x = uint8(4 * ones(1, M / 4));
    y = uint8(4 * ones(1, N / 4));
    cAw_ = mat2cell(cAw, x, y);
    water_extract = ones(64, 64);
    for i = 1:M/4
        for j = 1:N/4
            [~, T] = schur(cAw_{i, j});
            Tmax = max(T(:));
            if mod(Tmax, t) < t/2
                water_extract(i, j) = 0;
            else
                water_extract(i, j) = 1;
            end
        end
    end
end

function [imagemark] = embed_secret(carrier_image, secret_image, t)
    YCbCr = rgb2ycbcr(carrier_image);
    Cb = double(YCbCr(:,:,2));
    [cA, cB, cC, cD] = dwt2(Cb, 'haar');
    [M, N, ~] = size(cA);
    [wM, wN] = size(secret_image);
    x = uint8(4 * ones(1, M / 4));
    y = uint8(4 * ones(1, N / 4));
    cA_ = mat2cell(cA, x, y);
    for i = 1:M/4
        for j = 1:N/4
            if i <= wM && j <= wN
                [U, T] = schur(cA_{i, j});
                [m, n] = find(T == max(max(T)));
                Tmax = T(m, n);
                if secret_image(i, j) == 1
                    Tmax = Tmax - mod(Tmax, t) + 0.75 * t;
                else
                    Tmax = Tmax - mod(Tmax, t) + 0.25 * t;
                end
                T(m, n) = Tmax;
                cA_{i, j} = U * T * U';
            end
        end
    end
    cA = cell2mat(cA_);
    Cbmark = uint8(idwt2(cA, cB, cC, cD, 'haar'));
    YCbCr(:,:,2) = Cbmark;
    imagemark = ycbcr2rgb(YCbCr);
end
