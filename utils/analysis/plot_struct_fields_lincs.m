%Plot the Lincs of structures

function[] = plot_struct_fields_lincs(lincs_all, inc_angs, inc_angs_cut, xvals, xinds, title_name, xname, yname)

%make figure
figure('units','normalized','outerposition',[.1 .1 .9 .9]), hold on;
num_row = 3; num_col = 3;
lwidth = 2;
sgtitle(title_name)
xlabel(xname),ylabel(yname)
    
%plot fields
for ii = 1:length(xinds) %each xvalue
    ind = xinds(ii); %index of run
    subplot(floor(length(xinds)/num_row),num_row,ii), hold on; %subplot to focus on
    
    exclude_fields = []; %fields to exclude
    names = fieldnames(lincs_all); %names of field
    %plot lfields
    for i=1:length(names)
        name = names{i};
        if isnan(lincs_all.(name))
            exclude_fields(end+1) = i; %add to fields to exclude
        else
            if isnan(inc_angs) %use normal xvalues
                plot(lincs_all.(name)(ind,:),'LineWidth',lwidth) %plot results otherwise
            else %use inc_angs
                if length(lincs_all.(name)(ind,:)) ~= length(inc_angs) %not same size
                    plot(linspace(inc_angs(1),inc_angs(end),length(lincs_all.(name)(ind,:))) ...
                        , lincs_all.(name)(ind,:),'LineWidth',lwidth)
                else
                    plot(inc_angs,lincs_all.(name)(ind,:),'LineWidth',lwidth)
                end
                if ~isnan(inc_angs_cut) %cut inc angles to limit
                   xlim([inc_angs_cut(1),inc_angs_cut(2)]);
                end
            end
        end
    end
    names(exclude_fields) = []; %remove fields to exclude
    legend(names);
    title(sprintf('Plot at run value %.2f',xvals(ind)))
end




% figure('units','normalized','outerposition',[.1 -.6 1.4 1.5]); per_row = 3;
% for exp_ind = 1:length(lfield_locations_list)
%     Lincs = Linc_results{exp_ind};
%     inc_angles = Lincs.inc_angles;
%     subplot(floor(length(lfield_locations_list)/per_row),per_row,exp_ind)
%     title(strcat("Object at ",num2str(exp_xvals(exp_ind))))
%     hold on, plot(inc_angles,-Lincs.ground), plot(inc_angles,-Lincs.spatial), plot(inc_angles,-Lincs.red), plot(inc_angles,-Lincs.cp),
%     plot(inc_angles,-Lincs.cpr2), plot(inc_angles,-Lincs.bss)
%     legend('Ground','Spatial','Just Red','CP Alg','Relaxed CP','BSS')
% end