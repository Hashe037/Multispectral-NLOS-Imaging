%IMPROVE COMMENT
%
%Perform incident radiance with CP algorithm, which is the closest point
%between two sets. Takes an estimated spectrum of the object
%
%--------------------------------------------------------------------------
% Inputs
% gamma_cfov -- estimated spectrum of revealed CFOV slice
% msdfov_params -- holds information needed for MS-DFOV reconstruction
% r_params -- holds information about scene geometry needed for DFOV
% reconstruction
% components
%--------------------------------------------------------------------------
% Outputs
% L_inc,L_inc_dw -- reconstructed incident radiance
% dist -- distance for each dw (high means high residual = bad recon)
% projs -- holds information about the estimated strengths of the clutter
% objects
% extra_info -- structure that holds information that might be useful
%--------------------------------------------------------------------------
%
function[Sinc,Sinc_vant,dist,projs,extra_info] = ...
    perform_cpa_reconstruction(gamma_cfov,scene_params,brdf_params,meas_params,recon_params,mscpa_params,weighted_spec,dis_method,do_diff)


%instantiate some params
dalpha = meas_params.dalpha_meas_mscpa;
dalpha_occ = recon_params.dalpha_occ_mscpa;
diff_occ = recon_params.diff_occ;
diff_meas = meas_params.diff_meas;
% null_size_max = mscpa_params.nullmax;
% nullsize_method = mscpa_params.nullsize_method;
nullcut = mscpa_params.nullcut;
brdf_inc_angs_up = brdf_params.brdf_inc_angs_up;
vant_pos = scene_params.vant_pos;
occ_x = scene_params.occ_x;
occ_y = scene_params.occ_y;
num_specs = length(dalpha);
num_vants = size(dalpha{1},1);

extra_info = struct(); %holds extra parameters that could be useful

vecs = {}; %list of vectors from difference in points
a_list = {}; %keeps track of null vectors
A_list = {}; %keeps track of clutter space vectors
%vecs will be length M, where M is the vector number.

%make the object spectrally white and subtract each measurement from the first
%measurement to create num_specs-1 vectors
for i=2:num_specs
    vecs{end+1} = (dalpha{i}./gamma_cfov(:,i)-dalpha{1}./gamma_cfov(:,1))'; 
end

%go through each vantage point and run algorithm
dist = []; %distances
projs = []; %projection coefficients into A

rez = 1; %how many spatial vantages to look at a time (larger means more averaging)
for i=1:rez:num_vants

    %define vantages to look over
    if (i+rez) <= num_vants
        vants = i:(i+rez-1);
    else
        vants = i:num_vants;
    end

    %make the difference matrix. If the object is spectrally white, this will
    %remove the influence of the CFOV and only keep clutter contributions
    A = [];
    if do_diff %diff of light field
        for j=2:num_specs
            temp = diff_meas{j}(vants,:)./repmat(gamma_cfov(vants,j),1,size(diff_meas{j},2)) ...
                - diff_meas{1}(vants,:)./repmat(gamma_cfov(vants,1),1,size(diff_meas{1},2)) ;
            A(:,j-1) = sum(temp,1);
        end
    else %dalpha
        for j=1:length(vecs)
            temp = vecs{j}(:,vants); A(:,j) = temp(:);
        end
    end
    %construct the hyperplane equivalent by taking the nullspace of the
    %matrix A plus offset which is from the weighted_spec (dalpha{weighted_spec}).
    %Nullspace is found by taking SVD of composition and finding the right
    %singular vectors that are greater than the rank.
    [Utemp,Dtemp,Vtemp] = svd(A');

    %find null_size by the number of dimensions that contain 1-nullcut
    %amount of the power in the matrix, where nullcut is set by the user
    %(typically .05)
    if ~do_diff 
        null_size(i) = find_nullsize(diag(Dtemp),nullcut); %not all vantages have same amount of clutter
        if size(A,1) > size(A,2) %numM > (N-1) so add difference
            null_size(i) = null_size(i) + size(A,1)-size(A,2);
        end
    
        %make sure nullsize is not too small
        null_size(i) = max(null_size(i),size(Vtemp,2)-2);%?
    else %if doing differential, set nullsize to 4 (autmoatic doesn't work well)
        null_size(i) = 4;%null_size_max;
    end

    % null_size(i) = 2;
    nulla2 = Vtemp(:,(end-null_size(i)+1):end)';

    %finding dalpha_occ_j which is the dalpha_occ for this specific set of
    %vantages
    if ~do_diff %dalpha
        temp = dalpha{weighted_spec}(vants,:)/gamma_cfov(vants(1),weighted_spec); x0 = temp(:);
        temp = dalpha_occ(vants,:); occ_j = temp(:);
    else %differential
        x0 = sum(diff_meas{weighted_spec}(:,vants),1)'/gamma_cfov(vants(1),weighted_spec);
        occ_j = sum(diff_occ(vants,:),2);
    end

    %perform actual CPA algorithm
    [dist(vants(1)),Sinc_vant(vants(1)),temp_projs,otherpoint(vants(1),:)] = cp_alg(occ_j,A',nulla2,x0,dis_method);
    dist(vants) = dist(vants(1));
    Sinc_vant(vants) = Sinc_vant(vants(1));
    projs(vants) = norm(temp_projs(:));

    %some other metrics
    ortho(i) = norm(nulla2/norm(nulla2(:))*occ_j/norm(occ_j(:))); %orthogonality of nullspace to beta
    met1(i) = norm(nulla2'*nulla2);
    a_list{i} = nulla2; %save nullspace 
    A_list{i} = A; %save clutter space  
    U_list{i} = Utemp; %save utemp
    V_list{i} = Vtemp; %save utemp
end
%convert to incident radiance
Sinc = incvant_to_incang(Sinc_vant,brdf_inc_angs_up,vant_pos,occ_x,occ_y); %convert to angles

%save extra parameters
extra_info.ortho = A_list; %list of clutter spaces
extra_info.ortho = a_list; %list of null spaces
extra_info.ortho = ortho; %orthogonality of beta vector to nullspace
extra_info.met1 = met1; %
extra_info.null_size = null_size;





