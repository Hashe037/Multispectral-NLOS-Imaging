%Plot the fields of a structure (one-deep)

function[] = plot_struct_fields(results, xvals, do_log, title_name,xname,yname)

lwidth = 3;
msize = 6;
tag = '-o';

%make figure
figure,hold on
title(title_name)
xlabel(xname),ylabel(yname)
    
%plot fields
exclude_fields = []; %fields to exclude
names = fieldnames(results); %names of field
for i=1:length(names)
    name = names{i};
    if isnan(results.(name))
        exclude_fields(end+1) = i; %add to fields to exclude
    else
        if isnan(xvals) %use normal xvalues
            if do_log
                plot(log10(results.(name)),tag,'Linewidth',lwidth,'MarkerSize',msize) %plot results otherwise
            else
                plot(results.(name),tag,'Linewidth',lwidth,'MarkerSize',msize) %plot results otherwise
            end
        else
            if do_log
                plot(xvals,log10(results.(name)),tag,'Linewidth',lwidth,'MarkerSize',msize)
            else
                plot(xvals,results.(name),tag,'Linewidth',lwidth,'MarkerSize',msize)
            end
        end
    end
end

names(exclude_fields) = []; %remove fields to exclude
legend(names);