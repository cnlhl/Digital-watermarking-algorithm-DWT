clc
clear all
close all

image = imread('carrier_1.jpg');
water = imread('secret_fudan.bmp');
[M, N, C] = size(image);
t = 100;

if (M ~= 512 || N ~= 512 || C ~= 3)
    error('请输入512X512x3大小的图像');
end

[imagemark] = embed_secret(image, water, t);

[water_extract] = extract_water(imagemark, t);

compress = [90 85 75 65 45];
figure;

for i = 1:1:5
    str = strcat('image_compressed_', num2str(compress(i)), '.jpg');
    imwrite(imagemark, str, 'Quality', compress(i));
end

for i = 1:1:5
    str = strcat('image_compressed_', num2str(compress(i)), '.jpg');
    image_compress = imread(str);
    [image_compress_water_extract] = extract_water(image_compress, t);
    [CDR] = cdr(water, image_compress_water_extract);
    
    subplot(5, 2, 2*i-1);
    imshow(image_compress_water_extract);
    str = strcat(num2str(compress(i)), '-cmpress-extract-watermark');
    title(str);
    
    subplot(5, 2, 2*i);
    text(0.5, 0.5, sprintf('%s\n%.2f', strcat(num2str(compress(i)), '倍压缩率下的正确解码率'), CDR), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12);
    axis off;
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
