


function[recmet_all,dist_all,predmet_all,Sinc_all,separated_all]...
    = load_across_runs(recmet_results,dist_results,predmet_results,Sinc_results,separated_results)



%do reconstruction metrics
if iscell(recmet_results)
    recmet_all = flatten_struct(recmet_results);
else
    recmet_all = nan;
end
if iscell(dist_results)
    dist_all = flatten_struct(dist_results);
else
    dist_all = nan;
end
if iscell(predmet_results)
    predmet_all = flatten_struct(predmet_results);
else
    predmet_all = nan;
end
% if iscell(spec_results)
%     spec_all = flatten_struct(spec_results);
% else
%     spec_all = nan;
% end
% if iscell(ldif_results)
%     ldif_all = flatten_struct(ldif_results);
% else
%     ldif_all = nan;
% end
if iscell(Sinc_results)
    Sinc_all = flatten_struct(Sinc_results);
else
    Sinc_all = nan;
end
if iscell(separated_results)
    separated_all = flatten_struct(separated_results);
else
    separated_all = nan;
end


% 
% num_exps = length(remat_results); %number of runs
% 
% 
% 
% recmet_all = struct(); %reconstruction metrics
% top_fields = fieldnames(recmet_results{1});
% for top_ind = 1:length(top_fields) %each top-level field (total, each)
%     top_field = convertCharsToStrings(top_fields{top_ind}); 
%     mid_fields = fieldnames(recmet_results{1}.(top_field));
%     for mid_ind = 1:length(mid_fields) %each mid-level field (metric)
%         mid_field = convertCharsToStrings(mid_fields{mid_ind}); 
%         bot_fields = fieldnames(recmet_results{1}.(top_field).(mid_field));
%         for bot_ind = 1:length(bot_fields) %each bot-level field (method)
%             bot_field = convertCharsToStrings(bot_fields{bot_ind}); %each top-level field (total, each)
%             if iscell(recmet_results{1}.(top_field).(mid_field).(bot_field))
%                 recmet_all.(top_field).(mid_field).(bot_field) = return_results_from_dict_cellversion(recmet_results,[top_field,mid_field,bot_field]);
%             else
%                 recmet_all.(top_field).(mid_field).(bot_field) = return_results_from_dict(recmet_results,[top_field,mid_field,bot_field]);
%             end
%         end
%     end
% end
%     
% %do distance metrics
% distance_all = struct();
% bot_fields = fieldnames(distance_results{1});
% for bot_ind = 1:length(bot_fields)
%     bot_field = convertCharsToStrings(bot_fields{bot_ind});
%     if iscell(distance_results{1}.(bot_field))
%         distance_all.(bot_field) =  return_results_from_dict_cellversion(distance_results,[bot_field]);
%     else
%         distance_all.(bot_field) =  return_results_from_dict(distance_results,[bot_field]);
%     end
% end
% 
% %do predictive metrics
% predmet_all = struct();
% bot_fields = fieldnames(predmet_results{1});
% for bot_ind = 1:length(bot_fields)
%     bot_field = convertCharsToStrings(bot_fields{bot_ind});
%     if iscell(predmet_results{1}.(bot_field))
%         predmet_all.(bot_field) =  return_results_from_dict_cellversion(predmet_results,[bot_field]);
%     else
%         predmet_all.(bot_field) =  return_results_from_dict(predmet_results,[bot_field]);
%     end
% end
% 
% %do spec metrics
% spec_all = struct();
% mid_fields = fieldnames(spec_results{1});
% for mid_ind = 1:length(bot_fields)
%     bot_field = convertCharsToStrings(bot_fields{bot_ind});
%     if iscell(predmet_results{1}.(bot_field))
%         predmet_all.(bot_field) =  return_results_from_dict_cellversion(predmet_results,[bot_field]);
%     else
%         predmet_all.(bot_field) =  return_results_from_dict(predmet_results,[bot_field]);
%     end
% end
% 
% 
% recmet_all.each = struct(); %errors between each object and its best reconstruction
% recmet_all.each.mse = struct(); %non-normalized MSE
% recmet_all.each.shape = struct(); %normalized MSE
% recmet_all.each.powdif = struct(); %difference in "power" where power is defined by us
% 
% 
% 
% recmet_all.total = struct(); %errors between all objects and all reconstructions
% recmet_all.total.mse = struct(); %non-normalized MSE
% recmet_all.total.shape = struct(); %normalized MSE
% recmet_all.total.powdif = struct(); %difference in "power" where power is defined by us
% 
% 
% 
% 
% distance_all = struct();
% predmet_all = struct();
% spec_all = struct();
% lfielddif_all = struct();
% L_inc_all = struct();
% 
