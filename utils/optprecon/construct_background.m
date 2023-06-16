%construct background measurements for preconditioning

function[B_mat] = construct_background(l_field_clutter,w,do_diff)

%B is going to be the clutter light fields
B = {};
for clut_ind = 1:length(l_field_clutter)
    B{clut_ind} = l_field_clutter{clut_ind};
end

%consider differential
if do_diff
    Btemp = {};
    for clut_ind = 1:size(B,2)
        Btemp{clut_ind} = smoothdata(diff(B{clut_ind},1),1);
        Btemp{clut_ind} = Btemp{clut_ind}/norm(Btemp{clut_ind}(:));
    end
    B = Btemp;
end

%define the preconditioner
B_mat = zeros(numel(B{1}),length(B));
for clut_ind = 1:length(l_field_clutter)
    B_mat(:,clut_ind) = w*B{clut_ind}(:)/norm(B{clut_ind}(:));
end