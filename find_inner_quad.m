function [v12,v23,v34,v41] = find_inner_quad(e12,e23,e34,e41)
D12_23 = pdist2(e12,e23);
[min_vec,ind_12] = min(D12_23,[],1);
[min_val,ind_23] = min(min_vec);
v12 = e12(ind_12(ind_23),:);
v23 = e23(ind_23,:);

D34_41 = pdist2(e34,e41);
[min_vec,ind_34] = min(D34_41,[],1);
[min_val,ind_41] = min(min_vec);
v34 = e34(ind_34(ind_41),:);
v41 = e41(ind_41,:);

