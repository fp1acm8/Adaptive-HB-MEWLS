function [ QI_coeff, condition] = getcoeff_weighted_least_squares_pen(hspace,hmsh,data,f,lambda,weight)
% coefficients computed by global least squares (with smoothing) in the B-spline space

mass=op_gradgradu_gradgradv_hier (hspace, hspace, hmsh);
mass_matrix=mass;
col_matrix=sp_eval_alt(speye(hspace.ndof), hspace, data);
bvector=f;
bbbvector=col_matrix*diag(weight)*bvector(:);
S=col_matrix*diag(weight)*col_matrix'+ lambda*mass_matrix; %lambda scalar
QI_coeff=S\bbbvector; %pinv(S)*bbbvector;   %A/b
QI_coeff=QI_coeff';
condition=condest(S);

end