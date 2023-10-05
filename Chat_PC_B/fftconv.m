function res = fftconv(f,g)
% fftconv  Convolve two vectors. Adds negligeble time cost (~1ms) for small
% arrays (both lengths < 10 000) but improves significantly for larger arrays. 

    assert(ismatrix(f))
    assert(min(size(f))==1)
    if size(f,1) ~=1
        f = f';
    end
    
    assert(ismatrix(g))
    assert(min(size(g))==1)
    if size(g,1) ~=1
        g = g';
    end
    
    f_len = length(f);
    g_len = length(g);
    
    f = [f zeros(1,g_len-1)];
    g = [g zeros(1,f_len-1)];

    res = ifft(fft(g) .* fft(f));
end