%turn multispectral data into RGB colors with given color list

function[image_colored] = transform_multitocolor(image_multi,color_list)

image_colored = zeros(size(image_multi,1),size(image_multi,2),3);

for spec=1:size(image_multi,3)
    image_colored(:,:,1) = image_colored(:,:,1)+image_multi(:,:,spec)*color_list(spec,1)/norm(color_list(spec,:));
    image_colored(:,:,2) = image_colored(:,:,2)+image_multi(:,:,spec)*color_list(spec,2)/norm(color_list(spec,:));
    image_colored(:,:,3) = image_colored(:,:,3)+image_multi(:,:,spec)*color_list(spec,3)/norm(color_list(spec,:));
end