%cut the vantage positions to limit noise/bad vantages
%
function[vant_pos,lfield,lfield_ground,lfield_clutter] ...
    = cutVantPos(vant_pos,vant_cut,lfield,lfield_ground,lfield_clutter)


num_specs = length(lfield); %number of spectra measurements
num_ground = length(lfield_ground); %number of ground measurements
num_clutter = length(lfield_clutter); %number of clutter measurements

vant_pos = vant_pos(vant_cut(1):vant_cut(2)); %update vantages
for ii=1:num_specs
    lfield{ii} = lfield{ii}(vant_cut(1):vant_cut(2),:);
end
for ground_i=1:num_ground
    lfield_ground{ground_i} = lfield_ground{ground_i}(vant_cut(1):vant_cut(2),:);
end
for clutter_i=1:num_clutter
    lfield_clutter{clutter_i} = lfield_clutter{clutter_i}(vant_cut(1):vant_cut(2),:);
end
