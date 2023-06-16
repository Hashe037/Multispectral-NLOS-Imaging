%Normalize the Lincs and return new structure for better comparisons
%Assume depth is 2

function[Sinc_all_norm] = normalize_Sinc_struct(Sinc_all)

Sinc_all_norm = struct();

top_fields = fieldnames(Sinc_all);
for top_ind = 1:length(top_fields)
    top_field = convertCharsToStrings(top_fields{top_ind}); 
    mid_fields = fieldnames(Sinc_all.(top_field));
    for mid_ind = 1:length(mid_fields) 
        mid_field = convertCharsToStrings(mid_fields{mid_ind}); 
        arr = Sinc_all.(top_field).(mid_field);
        if length(size(arr))==2 %2d
            for i=1:size(arr,1) %each run
                Sinc_all_norm.(top_field).(mid_field)(i,:) = arr(i,:)/norm(squeeze(arr(i,:)));
            end 
        else %3d
            for i=1:size(arr,1)
                for j=1:size(arr,2)
                    Sinc_all_norm.(top_field).(mid_field)(i,:) = arr(i,j,:)/norm(squeeze(arr(i,j,:)));
                end 
            end
        end
    end
end
