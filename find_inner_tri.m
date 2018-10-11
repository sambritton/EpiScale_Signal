function [v12,v13,v23] = find_inner_tri(e12,e13,e23)
D12_13 = pdist2(e12,e13);
[min_vec,ind_12] = min(D12_13,[],1);
[min_val,ind_13] = min(min_vec);
v12 = e12(ind_12(ind_13),:);
v13 = e13(ind_13,:);

D23_13 = pdist2(e23,e13);
[min_vec,ind_23] = min(D23_13,[],1);
[min_val,ind_13] = min(min_vec);
v23 = e23(ind_23(ind_13),:);
