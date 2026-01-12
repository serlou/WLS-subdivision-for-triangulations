function theta = get_noise_reduction_coefficient(w,L,L_stencil,basis)

    A = compute_mask(w,L,L_stencil,basis,false);
    theta = 0;

    for i=1:2
        for j=1:2
            a = A(i:2:end,j:2:end);
            theta = max(theta,sum(a(:).^2));
        end
    end
end