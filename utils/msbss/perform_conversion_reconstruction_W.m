% MAKE COMMENTS
%
%convert W matrix (BSS components) to light fields and convert to
%reconstructions

function[lfields,alpha_meas,dalpha_meas,Sinc_vant,Sinc,dist] ...
    = perform_conversion_reconstruction_W(W,lfields_normal,alpha_normal,dalpha_normal,...
    msbss_params,scene_params,recon_params,brdf_params)

U_vant_t = brdf_params.U_vant_t;
vector_num = size(U_vant_t{1},2); %how many SVD vectors
precon = msbss_params.precon;

lfields = {}; %light fields of components
alpha_meas = {};
dalpha_meast = {};
dalpha_meas = {};

for ii = 1:size(W,2) %each unmixed component
    
    %extract dalpha for reconstructions
    if strcmp(precon,"alpha") %BSS was performed on alpha coef
        alpha_meas{ii} = reshape(W(:,ii),size(alpha_normal{1}));
        dalpha_meast = ([alpha_meas{ii}',zeros(vector_num,1)]-[zeros(vector_num,1),alpha_meas{ii}'])';
        dalpha_meas{ii} = dalpha_meast(2:end-1,:);
        lfields{ii} = zeros(size(alpha_meas{ii},1),size(U_vant_t{1},1));
        for vant=1:(size(alpha_meas{ii},1))
            lfields{ii}(vant,:) = lfields{ii}(vant,:)-(U_vant_t{vant}*alpha_meas{ii}(vant,:)')';
        end
        lfields{ii}(end,:) = lfields{ii}(end-1,:); %extrapolate last value
        
    elseif strcmp(precon,"diff") || strcmp(precon,"opt_precon")
        lfields_diff{ii} = reshape(W(:,ii),size(lfields_normal{ii}(2:end,:)));
        lfields{ii} = cumsum(lfields_diff{ii},'reverse');
        lfields{ii} = vertcat(lfields{ii}(1,:),lfields{ii}); %add row to be same size as lfields
        for i=1:length(r_params.vant_pos)
            alpha_meas{ii}(i,:) = project_lfield_singular2(r_params.U_vant_t{i},r_params.D_tot,lfields{ii},i);
        end
        dalpha_meast = ([alpha_meas{ii}',zeros(vector_num,1)]-[zeros(vector_num,1),alpha_meas{ii}'])';
        dalpha_meas{ii} = dalpha_meast(2:end-1,:);
        
        
    elseif strcmp(precon,"dalpha") %BSS was performed on dalpha coef
        dalpha_meas{ii} = reshape(W(:,ii),size(dalpha_normal{1}'))';
        alpha_meas{ii} = cumsum(dalpha_meas{ii},'reverse');
        lfields{ii} = zeros(size(alpha_meas{ii},1)+1,size(U_vant_t{1},1));
        for vant=1:(size(alpha_meas{ii},1)-1)
            lfields{ii}(vant,:) = lfields{ii}(vant,:)-(U_vant_t{vant}*alpha_meas{ii}(vant,:)')';
        end
        lfields{ii}(end,:) = lfields{ii}(end-1,:); %extrapolate last value

    elseif strcmp(precon,"dalpha_inverted") %last portion is inverted so discard
        half_size = size(W,1)/2;
        dalpha_meas{ii} = reshape(W(1:half_size,ii),size(dalpha_normal{1}'))';
        alpha_meas{ii} = cumsum(dalpha_meas{ii},'reverse');
        lfields{ii} = zeros(size(alpha_meas{ii},1)+1,size(U_vant_t{1},1));
        for vant=1:(size(alpha_meas{ii},1)-1)
            lfields{ii}(vant,:) = lfields{ii}(vant,:)-(U_vant_t{vant}*alpha_meas{ii}(vant,:)')';
        end
        lfields{ii}(end,:) = lfields{ii}(end-1,:); %extrapolate last value

    % elseif strcmp(precon,"Sinc_vant") %BSS performed on Sinc_vant
    %     dalpha_meas{ii} = repmat(W(:,ii),1,size(r_params.dalpha_occ,1))'.*r_params.dalpha_occ;
    %     alpha_meas{ii} = cumsum(dalpha_meas{ii}','reverse');
    %     lfields{ii} = zeros(size(alpha_meas{ii},1)+1,size(U_vant_t{1},1));
    %     for vant=1:(size(alpha_meas{ii},1)-1)
    %         lfields{ii}(vant,:) = lfields{ii}(vant,:)-(U_vant_t{vant}*alpha_meas{ii}(vant,:)')';
    %     end
    %     lfields{ii}(end,:) = lfields{ii}(end-1,:); %extrapolate last value
    % 
    % elseif strcmp(precon,"Sinc_vant_inverted") %BSS performed on Sinc_vant inverted
    %     half_size = size(W,1)/2;
    %     dalpha_meas{ii} = repmat(W(1:half_size,ii),1,size(r_params.dalpha_occ,1))'.*r_params.dalpha_occ;
    %     alpha_meas{ii} = cumsum(dalpha_meas{ii}','reverse');
    %     lfields{ii} = zeros(size(alpha_meas{ii},1)+1,size(U_vant_t{1},1));
    %     for vant=1:(size(alpha_meas{ii},1)-1)
    %         lfields{ii}(vant,:) = lfields{ii}(vant,:)-(U_vant_t{vant}*alpha_meas{ii}(vant,:)')';
    %     end
    %     lfields{ii}(end,:) = lfields{ii}(end-1,:); %extrapolate last value

    else %BSS performed on light field
        
        lfields{ii} = reshape(W(:,ii),size(lfields_normal{ii}));
        for i=1:length(r_params.vant_pos)
            alpha_meas{ii}(i,:) = project_lfield_singular2(r_params.U_vant_t{i},r_params.D_tot,lfields{ii},i);
        end
        dalpha_meast = ([alpha_meas{ii}',zeros(vector_num,1)]-[zeros(vector_num,1),alpha_meas{ii}'])';
        dalpha_meas{ii} = dalpha_meast(2:end-1,:);
    end

    %smooth dalpha
    %normally want to smooth alpha but cannot in these methods since we use
    %dalpha already
    for i=1:vector_num
        f = fit((1:size(dalpha_meas{ii},1))',dalpha_meas{ii}(:,i),'smoothingspline','SmoothingParam',recon_params.alpha_smooth);
        dalpha_meas{ii}(:,i) = f(1:length(dalpha_meas{ii}));
    end
end

%perform reconstruction
[Sinc_vant,Sinc,dist] = performLfieldRecon(dalpha_meas,scene_params,brdf_params,recon_params);
