function f = apply_mask(f, A, N)
    % N iterations of subdivision in uniform regular grids, with mask A
    % f : initial data (2D array)
    
    for i=1:N
        % 2D upscaling with zeros in the even positions
        f1 = zeros(2*size(f));
        f1(1:2:end, 1:2:end) = f;
        
        % 2D valid convolution
        f = conv2(f1, A, 'valid');
    end
end