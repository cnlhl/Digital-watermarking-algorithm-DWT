clc
clear all
close all

carriers = {'carrier_1.jpg', 'carrier_2.jpg', 'carrier_3.jpg', 'carrier_4.jpg'};
t = 100;

water = imread('secret_fudan.bmp');
water = arnold(water,1,1,10);

figure;
for carrierIndex = 1:numel(carriers)
    carrier = imread(carriers{carrierIndex});
    [M, N, C] = size(carrier);
    imagemark = embed_secret(carrier, water, t);
    water_extract = extract_water(imagemark, t);
    im_cut = ones(M,N,C);
    for i=1:1:3
      im_cut(:,:,i) = uint8(imagemark(:,:,i));
      im_cut(150:250,100:250,i) = 0; 
    end
    im_cut = uint8(im_cut);

    subplot(2, numel(carriers), (carrierIndex));
    imshow(im_cut);
    title('image-cut');
    water_extract_cut = extract_water( im_cut,t ); 
    
    subplot(2, numel(carriers), (carrierIndex)+4);
    imshow(water_extract_cut);
    title('water-extract-cut');
    [CDR] = cdr(water, water_extract_cut);
    disp('CDR under cut:');
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
    water_extract = rearnold(water_extract,1,1,10);
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


% Function to apply Arnold transform
function arnoldImg = arnold(img,a,b,n)
    [h,w] = size(img);
    N = h;
    arnoldImg = zeros(h,w);
    for i=1:n
        for y=1:h
            for x=1:w
                xx = mod((x-1) + b*(y-1), N) + 1;
                yy = mod(a*(x-1) + (a*b+1)*(y-1), N) + 1;  
                arnoldImg(yy,xx) = img(y,x);              
            end
        end
        img = arnoldImg;
    end
end

% Function to reverse Arnold transform
function img = rearnold(arnoldImg,a,b,n)
    [h,w] = size(arnoldImg);
    img = zeros(h,w);
    N = h;
    for i=1:n
        for y=1:h
            for x=1:w           
                xx = mod((a*b+1)*(x-1) - b*(y-1), N) + 1;
                yy = mod(-a*(x-1) + (y-1), N) + 1;      
                img(yy,xx) = arnoldImg(y,x);              
            end
        end
        arnoldImg = img;
    end
end

