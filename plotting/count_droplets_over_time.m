% Script to calculate the number of droplets at each timestep for
% Ostwald ripening simulations. This uses the trackedDroplets.mat file that
% is generated for a simulation when level_set_radius_multiple_droplets.m
% is run. 
% Sarah Groves July 25, 2024
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/nonlinear_multigrid/plotting/radii_over_time_level_set_plots/domain_0_2_e_0.0075_noisy_cohesin";
total_time = 0.03;
%name = "phi_512_19660_1.0e-5__CPC_20_cohesin_8_eps_0.007504684956431058";
filePattern = fullfile(indir, 'phi*'); % Change to whatever pattern you need.
% filePattern = fullfile(myFolder, '**/*.exe'); % To look in subfolders too, add **/.
folders = dir(filePattern);
for k = 1 : length(folders)
    name = folders(k).name
    name_split  = regexp(name,"_",'split');
    % phi_512_19661_1.0e-5__CPC_0.12_cohesin_0.2_eps_0.0075_alpha_0_domain_0_2
    tmp = name_split(8);
    CPC = tmp{1};
    tmp = name_split(10);
    cohesin = tmp{1};
    eps = "0.0075";
    alpha = "0";
    % CPC = "0.15";
    % cohesin = "0.1";
    % eps = "0.0096";
    % alpha = "0";
    % name = sprintf("phi_512_19661_1.0e-5__CPC_%s_cohesin_%s_eps_%s_alpha_%s", CPC, cohesin, eps, alpha);
    load(sprintf("%s/%s/trackedDroplets.mat", indir, name));
    dt = .1*(1/256)^2;
    num_timepoints = round(total_time/(10*dt))+1;
    num_droplets = zeros(num_timepoints,1);

    for t = 1:num_timepoints
        time = (t-1)*dt*10;
        for d = 1:max([trackedDroplets.id])
            if ismember(d, [1 2]) %ignore IC droplets
            elseif ismembertol(time,trackedDroplets(d).time, 1e-3)
                % only look within 1.6um of IC
                if trackedDroplets(d).center(2) <= 384 && trackedDroplets(d).center(2) >= 128 
                    num_droplets(t) = num_droplets(t)+ 1;
                end
            end
        end
    end

    % just calculating 1/2 of the cohesin axis, for 1.6 away from IC
    num_droplets = num_droplets/2;

    % cumulative_count_data = cumsum(num_droplets);
    % Create a bar chart
    [GC,GR] = groupcounts(num_droplets) ;
    % bar(GR(2:end), GC(2:end)/sum(GC(2:end)));
    norm_gc = GC(2:end)/sum(GC(2:end));
    fid = fopen('simulated_droplet_distributions_e_0.0075_noisy_cohesin.csv', 'a+');
    %assumes there are never more than 4 droplets
    % fprintf(fid, '%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n', CPC, cohesin, eps, alpha, ...
    %     norm_gc(1),norm_gc(2),norm_gc(3),norm_gc(4), ...
    %     GC(1),GC(2),GC(3),GC(4),GC(5),trackedDroplets(1).radius(end));
    % fclose(fid);
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%d\n', CPC, cohesin, eps, alpha, ...
        mat2str(norm_gc), ...
        mat2str(GC),6.4*trackedDroplets(1).radius(end));
    fclose(fid);
end
% title(sprintf('Number of non-IC droplets within 1.6um from IC \n %s', name));
% ylabel("Frequency")
% xlabel("Number of non-IC droplets in one direction from IC")
% print(gcf,sprintf('%s/%s/num_droplets_1600nm.png', indir, name),"-dpng")
% 
% bar(GR(2:end), GC(2:end)*dt*10);
% title(sprintf('Number of non-IC droplets within 1.6um from IC \n %s', name));
% ylabel("Seconds")
% xlabel("Number of non-IC droplets in one direction from IC")
% print(gcf,sprintf('%s/%s/num_droplets_1600nm_seconds.png', indir, name),"-dpng")
