%Perform incident radiance with CP algorithm, which is the closest point
%between two sets. Takes an estimated spectrum of the object
%
%--------------------------------------------------------------------------
% Inputs
% estspec -- estimated spectrum of revealed CFOV slice
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
function[L_inc,L_inc_dw,dist,projs,extra_info] = ...
    perform_cp_reconstruction3(estspec,msdfov_params,r_params,weighted_spec)


%instantiate some params
diff = msdfov_params.diff_meas;
dalpha = msdfov_params.dalpha_meas;
beta_inc = msdfov_params.beta_inc;
beta_diff_inc = msdfov_params.beta_diff_inc;
null_size_max = msdfov_params.nullmax;
dis_method = msdfov_params.dis_method;
nullsize_method = msdfov_params.nullsize_method;
nullcut = msdfov_params.nullcut;
inc_angles_up = r_params.inc_angles_up;
vant_pos = r_params.vant_pos;
d = r_params.d;
y_l = r_params.y_l;
numSpecs = length(dalpha);
numVants = size(dalpha{1},2);

extra_info = struct(); %holds extra parameters that could be useful

vecs = {}; %list of vectors from difference in points
a_list = {}; %keeps track of null vectors
A_list = {}; %keeps track of clutter space vectors
%vecs will be length M, where M is the vector number.

%go through each vantage point and run algorithm
dist = []; %distances
projs = []; %projection coefficients into A
max_svd = numSpecs;%size(dalpha{1},2);

%make the difference matrix
% for i=2:numSpecs
%     vecs{end+1} = diff{i}./estspec(:,i)'-diff{1}./estspec(:,1)'; %MAY NEED TO DOUBLE CHECK THIS
% end

rez = 1; %how many vantages to look at a time
for i=1:rez:numVants

    %define vantages to look over
    if (i+rez) <= numVants
        vants = i:(i+rez);
    else
        vants = i:numVants;
    end

    %construct the difference matrix A
    A = [];
    for j=2:numSpecs
        temp = diff{j}(vants,:)./repmat(estspec(vants,j),1,size(diff{j},2)) ...
            - diff{1}(vants,:)./repmat(estspec(vants,1),1,size(diff{1},2)) ;
        A(:,j-1) = sum(temp,1);
    end
    %construct the hyperplane equivalent by taking the nullspace of the
    %matrix A plus offset which is from the first BRDF (dalpha_meas{1}).
    %Nullspace is found by taking SVD of composition and finding the right
    %singular vectors that are greater than the rank.
%     if size(A,2) >= size(A,1) %no nullspace, so shorten using PCA? HAVE
%     TO DO THIS FOR BETA AND EVERYTHING ELSE
%         des_size = size(A,1)-1;
%         [Utemp,Dtemp,Vtemp] = svd(A');
%         A = Utemp(:,1:des_size)*Dtemp(1:des_size,1:des_size);%*Vtemp(:,1:des_size)';
%         A = A';
%     end
    [Utemp,Dtemp,Vtemp] = svd(A','econ');

    %nullsize ideally = numM-numClutter
    if nullsize_method == 1 %OLED, two clutter one object
        if i<90
            null_size(i) = max_svd-2;
        else
            null_size(i) = max_svd-1;
        end
    elseif nullsize_method == 2 %OLED, two clutter one object v2
        if i<10
            null_size(i) = max_svd-2;
        else
            null_size(i) = max_svd-1;
        end
    elseif nullsize_method == 3 %OLED, 2 clut, multiple object
        null_size(i) = max_svd-2;%null_size_max;
    elseif nullsize_method == 4 %OLED, 1 clut, multiple object
        null_size(i) = max_svd-1;%null_size_max;
    elseif nullsize_method == 5 %OLED, 3 clut, multiple object
        null_size(i) = max_svd-3;%null_size_max;
    else %use automatic method
        null_size(i) = find_nullsize(diag(Dtemp),nullcut); %not all vantages have same amount of clutter
%         if size(A,1) > size(A,2) %numM > (N-1) so add difference
%             null_size(i) = null_size(i) + size(A,1)-size(A,2);
%         end
%         null_size(i) = max(null_size(i),null_size_max);
    end

%     null_size(i) = 4;
%     null_size(i) = max(null_size(i),null_size_max);
%     null_size(i) = min(null_size(i),4);
%     null_size(i) = 1;

    %hardcoded two clutter
%     if i<90
%         null_size(i) = 2;
%     else
%         null_size(i) = 3;
%     end

    %hardcoded two clutter v2
%     if i<100
%         null_size(i) = max_svd-2;
%     else
%         null_size(i) = max_svd-1;
%     end

%     %hardcoded multiobject (flat run3)
    null_size(i) = 4;%null_size_max;

    %hardcoded multiobject (flat 10ints2clut)
%     null_size(i) = 2;%null_size_max;
% 
%     null_size(231) = 2;

    nulla2 = Vtemp(:,(end-null_size(i)+1):end)';
%     nulla2 = Vtemp(:,numSpecs:end)'; %hyperplane of nullspace
    
    %find closest points
%     x0 = sum(diff{1}(vants,:),1)';%/estspec(vants(1),1);
    x0 = sum(diff{weighted_spec}(vants,:),1)'/estspec(vants(1),weighted_spec);
%     x0 = sum(diff{4}(vants,:),1)'/estspec(vants(1),4);
    beta_j = sum(beta_diff_inc(:,vants),2);
    [dist(vants(1)),L_inc_dw(vants(1)),temp_projs,otherpoint(vants(1),:)] = runCPalg(beta_j,A',nulla2,x0,dis_method);
    dist(vants) = dist(vants(1));
    L_inc_dw(vants) = L_inc_dw(vants(1));
%     otherpoint(vants,:) = otherpoint(vants(1),:)
    projs(vants) = norm(temp_projs(:));


    %some other metrics
    ortho(i) = norm(nulla2/norm(nulla2(:))*beta_j/norm(beta_j(:))); %orthogonality of nullspace to beta
    met1(i) = norm(nulla2'*nulla2);
    a_list{i} = nulla2; %save nullspace 
    A_list{i} = A; %save clutter space  
    U_list{i} = Utemp; %save utemp
    V_list{i} = Vtemp; %save utemp
end

%convert to incident radiance
L_inc = dw_to_inc2(L_inc_dw,inc_angles_up,vant_pos,d,y_l);
% figure,plot(L_inc_dw),title('lincdw')
% figure,plot(dist),title('Dist')
% figure,plot(null_size),title('Nullsize')

%save extra parameters
extra_info.ortho = A_list; %list of clutter spaces
extra_info.ortho = a_list; %list of null spaces
extra_info.ortho = ortho; %orthogonality of beta vector to nullspace
extra_info.met1 = met1; %
extra_info.null_size = null_size;




