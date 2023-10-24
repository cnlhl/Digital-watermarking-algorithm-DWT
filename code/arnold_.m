secret_image = imread('secret_fudan.bmp');

% Applying Arnold transform
arnoldim = arnold(secret_image,1,1,10);

% Recovering the original image
recover = rearnold(arnoldim,1,1,10);

% Displaying the transformed image
subplot(1,2,1);
imshow(arnoldim);   
title('After Disturbance');

% Displaying the recovered image
subplot(1,2,2);
imshow(recover);   
title('After Recovery');

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
