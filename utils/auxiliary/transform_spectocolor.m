%turn spectrum array into RGB colors with given color list

function[spec_rgb] = transform_spectocolor(spec_list,color_list)

spec_rgb = [0,0,0];

for spec=1:length(spec_list)
    spec_rgb = spec_rgb + spec_list(spec)*color_list(spec,:);%/norm(color_list(spec,:));
end
spec_rgb = rescale(spec_rgb);