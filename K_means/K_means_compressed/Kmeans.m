clear 
clc

train_img = double(imread('timg2.png'));   %导入数据（图片）

[h,w,~] = size(train_img);

kn = 30;                        %聚类数为16
key = zeros(kn, 3);
%选取中心点
for i = 1 : kn
    x = round(random('unif', 1, h));
    y = round(random('unif', 1, w));
    for j = 1 : 3
        key(i, j) = train_img(x, y ,j);
    end
end

max_iterate = 50;      %最大迭代次数

for t = 1 : max_iterate
    sumloc = zeros(kn, 3);
    count = zeros(kn, 1);

    for i = 1 : h
        for j = 1 : w
            min_d = dis(key(1, :), train_img(i, j, :));
            min_k = 1;
            for k = 2 : kn
                d = dis(key(k, :), train_img(i, j, :));
                if d < min_d
                    min_d = d;
                    min_k = k;
                end
            end
            count(min_k) = count(min_k) + 1;
            for c = 1 : 3
                sumloc(min_k, c) = sumloc(min_k, c) + train_img(i, j, c);
            end
        end
    end
%更新聚类中心
old_key = key;
    for k = 1 : kn
        key(k, :) = round(sumloc(k, :) / count(k));
    end
end


large = double(imread('timg2.png'));
compressed = large;

[h,w,~] = size(train_img);

for i = 1 : h
    for j = 1 : w
        min_d = dis(key(1, :), large(i, j, :));
        min_k = 1;
        for k = 2 : kn
            d = dis(key(k, :), large(i, j, :));
            if d < min_d
                min_d = d;
                min_k = k;
            end
        end
        for c = 1 : 3
            compressed(i, j, c) = key(min_k, c);
        end
    end
end

figure
subplot(1,2,1)
imshow(uint8(round(large)))
title('original')

subplot(1,2,2)
imshow(uint8(round(compressed)))
title('compressed')

imwrite(uint8(round(compressed)), 'compressed.bmp');
