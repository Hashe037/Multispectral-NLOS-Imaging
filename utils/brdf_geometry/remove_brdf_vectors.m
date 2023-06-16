%remove SVD vectors from the BRDF components

function[U_tot2,V_tot2,D_tot2,U_vant2,V_vant2,vector_num] ...
    = remove_brdf_vectors(vectors_remove,U_tot,V_tot,D_tot,U_vant,V_vant)

%remove certain svd vectors
U_tot2 = U_tot;
V_tot2 = V_tot;
D_tot2 = D_tot;
U_vant2 = U_vant;
V_vant2 = V_vant;
%account for removal
vectors_remove = vectors_remove-(0:(length(vectors_remove)-1));
for ind = length(vectors_remove)
    U_tot2(:,vectors_remove(ind)) = [];
    D_tot2(vectors_remove(ind)) = [];
    V_tot2(:,vectors_remove(ind)) = [];
    for i = 1:length(U_vant)
        U_vant2{i}(:,vectors_remove(ind)) = [];
        % V_vant2{i}(:,vectors_remove(ind)) = [];
    end
end
vector_num = size(U_tot2,2);