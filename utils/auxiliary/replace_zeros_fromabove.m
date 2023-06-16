%replace zeros with nearest value above it

function[image2] = replace_zeros_fromabove(image)

image2 = image;

zero_indices = find(image==0); %index of zeros

maxr = size(image,1);
maxc = size(image,2);
maxdepth = 10;
minr = 1;
minc = 1;

for i = 1:length(zero_indices)
    ind = zero_indices(i);
    [r,c] = ind2sub(size(image),ind);

    %find (if any) nearest non-zero above
    r2 = r;
    notfound = 1;
    while r2 > max(0,r2-maxdepth) && notfound
        if image(r2,c) ~= 0
            image2(r,c) = image(r2,c);
            notfound = 0;
        end
        r2 = r2-1;
    end
end
