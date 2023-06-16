%Find the mixture to use for the BSS method.
%
%--------------------------------------------------------------------------
% Inputs
% -lfield,alpha_meas,dalpha_meas,Sinc_vant -- defined in previous
% functions
% -precon: which preconditionig on mixture to use (string with only several supported
% values)
%
%--------------------------------------------------------------------------
% Outputs
% -X linear mixture measurements (numSpecs x numElements)


function[X] = find_bss_mixture(lfield,alpha_meas,dalpha_meas,Sinc_vant, ...
    precon)

num_specs = length(lfield);

if strcmp(precon,"alpha") %on alpha values
    X = zeros(num_specs,numel(alpha_meas{1})); %hold all light fields
    for spec = 1:num_specs
        temp = alpha_meas{spec};
        X(spec,:) = temp(:);
    end
elseif strcmp(precon,"diff") %differential light field
    X = zeros(num_specs,numel(lfield{1}(2:end,:))); %hold all light fields
    for spec = 1:num_specs
        temp = smoothdata(diff(lfield{spec},1),1);
        X(spec,:) = temp(:);
    end
elseif strcmp(precon,"dalpha") %on dalpha values
    X = zeros(num_specs,numel(dalpha_meas{1})); %hold all light fields
    for spec = 1:num_specs
        temp = dalpha_meas{spec};
        X(spec,:) = temp(:);
    end
% elseif strcmp(precon,"Sinc_vant") %on Sinc_vant reconstructions
%     X = zeros(num_specs,numel(Sinc_vant{1})); %hold all light fields
%     for spec = 1:num_specs
%         X(spec,:) = Sinc_vant{spec};
%     end
% elseif strcmp(precon,"Sinc_vant_inverted") %invert Sinc_vant and add to end
%     X = zeros(num_specs,numel(Sinc_vant{1})*2);
%     for spec = 1:num_specs
%         X(spec,:) = [Sinc_vant{spec},-1*Sinc_vant{spec}];
%     end 
elseif strcmp(precon,"dalpha_inverted") %invert dalpha and add to end
    X = zeros(num_specs,numel(dalpha_meas{1})*2); %hold all light fields
    for spec = 1:num_specs
        temp = dalpha_meas{spec}';
        X(spec,:) = [temp(:);-1*temp(:)];
    end
elseif strcmp(precon,"opt_precon") %run with opt precon as preconditioner
    X = zeros(num_specs,numel(lfield{1}(2:end,:))); %hold all light fields
    for spec = 1:num_specs
        X(spec,:) = run_params.l_field_precon{spec}(:);
    end
else %on light fields
    X = zeros(num_specs,numel(lfield{1})); %hold all light fields
    for spec = 1:num_specs
        temp = lfield{spec};
        X(spec,:) = temp(:);
    end
end