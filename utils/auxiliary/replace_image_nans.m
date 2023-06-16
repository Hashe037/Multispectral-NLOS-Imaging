%replace nans with weighted nearest non-nan neighbors. If all neighbors are
%also nan, leave as nan and rerun after

function[image2] = replace_image_nans(image)

image2 = image;

nan_indices = find(isnan(image)); %index of nans

maxr = size(image,1);
maxc = size(image,2);
minr = 1;
minc = 1;

iteration = 1;
while (iteration < 10) && (sum(isnan(image(:)))>0)
    for i = 1:length(nan_indices)
        ind = nan_indices(i);
        [r,c] = ind2sub(size(image),ind);
       
        a1 = [min(maxr,r+1),max(minc,c-1)];
        a2 = [max(minr,r-1),max(minc,c-1)];
        a3 = [min(maxr,r+1),min(maxc,c+1)];
        a4 = [max(minr,r-1),min(maxc,c+1)];
        a5 = [r,min(maxc,c+1)];
        a6 = [r,max(minc,c-1)];
        a7 = [min(maxr,r+1),c];
        a8 = [max(minr,r-1),c];
        
        nvals = [image(a1(1),a1(2)),image(a2(1),a2(2)),image(a3(1),a3(2)),image(a4(1),a4(2)),...
            image(a5(1),a5(2)),image(a6(1),a6(2)),image(a7(1),a7(2)),image(a8(1),a8(2))];
        
        image2(r,c) = mean(nvals,'omitnan');
    end
    nan_indices = find(isnan(image)); %index of nans
    iteration = iteration + 1;
    
    image = image2;
end
