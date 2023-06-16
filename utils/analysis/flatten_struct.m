%Flatten multirun struct based on its order
%
%NOTE: each structure must have the same depth for all fields

function[dict_all] = flatten_struct(results)


%find order
order = 0;
tempresults = results{1};
try
    while 1==1
        names = fieldnames(tempresults);
        tempresults = tempresults.(names{1});
        order = order+1;
    end
catch
    1==1;
end

%perform based on order (max order is 3)
dict_all = struct();
if order == 1
    bot_fields = fieldnames(results{1});
    for bot_ind = 1:length(bot_fields)
        bot_field = convertCharsToStrings(bot_fields{bot_ind});
        dict_all.(bot_field) =  return_results_from_dict(results,[bot_field]);
    end
elseif order==2
    mid_fields = fieldnames(results{1});
    for mid_ind = 1:length(mid_fields) 
        mid_field = convertCharsToStrings(mid_fields{mid_ind}); 
        bot_fields = fieldnames(results{1}.(mid_field));
        for bot_ind = 1:length(bot_fields)
            bot_field = convertCharsToStrings(bot_fields{bot_ind});
            dict_all.(mid_field).(bot_field) = return_results_from_dict(results,[mid_field,bot_field]);
        end
    end
elseif order==3
    top_fields = fieldnames(results{1});
    for top_ind = 1:length(top_fields)
        top_field = convertCharsToStrings(top_fields{top_ind}); 
        mid_fields = fieldnames(results{1}.(top_field));
        for mid_ind = 1:length(mid_fields) 
            mid_field = convertCharsToStrings(mid_fields{mid_ind}); 
            bot_fields = fieldnames(results{1}.(top_field).(mid_field));
            for bot_ind = 1:length(bot_fields)
                bot_field = convertCharsToStrings(bot_fields{bot_ind});
                dict_all.(top_field).(mid_field).(bot_field) = return_results_from_dict(results,[top_field,mid_field,bot_field]);
            end
        end
    end
else
    error('Error: order not supported')
end