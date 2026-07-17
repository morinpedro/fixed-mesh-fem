% ---- data of the manufactured problem -------------------------------------

% diffusion coefficient (a) of the equation
coeff_a = @(x) 1.0;
% reaction coefficient (c) of the equation
coeff_c = @(x) 0.0;

u_ex = @(x) exp(-10*norm(x)^2);
grad_u_ex = @(x) -20*u_ex(x).*x;
fc_f = @(x) 20*exp(-10*norm(x)^2)*(2-20*norm(x)^2);
fc_gD = @(x) u_ex(x);
fc_gN = @(x) 0;

%% ---- resolution loop ----------

prevH = -1;
prevL = -1;

fprintf("-----------------------------------------\n");
fprintf("  N     Np    H1err   rate   L2err   rate\n");
fprintf("-----------------------------------------\n");

for N = 2.^[2:6]
    gen_mesh_rectangle(N);
    fem
    l2 = L2_err(elem_vertices,vertex_coordinates,uh,u_ex);
    h1 = H1_err(elem_vertices,vertex_coordinates,uh,grad_u_ex);
    if (prevH < 0)
        str = " %3d  %5d  %4.2e ----  %4.2e ----\n";
        fprintf(str,N,size(vertex_coordinates,1),h1,l2);
    else
        rateH1 = log2(prevH / h1);
        rateL2 = log2(prevL / l2);
        str = " %3d  %5d  %4.2e %4.2f  %4.2e %4.2f\n";
        fprintf(str,N,size(vertex_coordinates,1),h1,rateH1,l2,rateL2);
    end
    prevH = h1;
    prevL = l2;
end
fprintf("-----------------------------------------\n");
