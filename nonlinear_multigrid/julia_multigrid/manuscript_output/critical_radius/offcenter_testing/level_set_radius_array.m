function [] = level_set_radius_array(R0, m, Nx,indir, suffix)
    folder = "plots";
    total_time=10; 
    everyR=10;
    eps = m * (1 / Nx) / (2 * sqrt(2) * atanh(0.9))
    epsilon_name = sprintf("%.5g", eps)
    R0_name = sprintf('%.5g', R0)

    [rr,tt] = level_set_plot(2.5e-5, indir, total_time, everyR, epsilon_name, R0_name, folder, Nx, suffix);
    R0_vector = repmat(R0, 1,length(tt));
    length(rr)
    length(R0_vector)
    length(tt)
    tt=tt(1:length(rr));
    length(tt)
    R0_vector=R0_vector(1:length(rr));
    length(R0_vector)

    T = table(transpose(rr), transpose(tt), transpose(R0_vector),'VariableNames',{'radius', 'time', 'R0'});
    % writetable(T,sprintf('%s/radius_0.5_level_set_epsilon_%s_%d_%s.txt',indir, epsilon_name, Nx, suffix),'WriteMode','append')

    % M = [[rr; nan], [tt; nan], R0_vector]
end
