function Lphi = Lap(phi,k2)
    Lphi = real(ifft2(k2 .* fft2_filtered(phi)));
end