% project lfield onto singular vectors of BRDF. Projection will be
% alpha_i, where i is the singular vector. This is done just for one
% vantage point.
%
% Input:
% V -- left singular vectors of VANTAGE POINT (outgoing angles)
% D -- diagonal values
% l_field -- measured light field
% vant -- vantage position

function[alpha] = project_lfield_onto_vector(U,D,l_field,vant)

alpha = zeros(size(U,2),1);
for i=1:size(U,2) %each singular vector
    alpha(i) = 1/D(i)*dot(U(:,i),l_field(vant,:)); %projection
%     alpha(i) = dot(U(:,i),l_field(vant,:))/length(V(:,i)); %projection
end
