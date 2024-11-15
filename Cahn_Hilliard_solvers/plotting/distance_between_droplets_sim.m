% Script to calculate the distance between droplets at each timestep for
% Ostwald ripening simulations. This uses the trackedDroplets.mat file that
% is generated for a simulation when level_set_radius_multiple_droplets.m
% is run. 
% Sarah Groves August 28, 2024
% indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/domain_0_2_noisy_cohesin_sd_0.25/";
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/domain_0_1_eps_0.0075/";

total_time = 0.03;
spacing = 0.000001525878906;
x = 0:spacing:0.03;
%name = "phi_512_19660_1.0e-5__CPC_20_cohesin_8_eps_0.007504684956431058";
cohesin = "0.08";
eps = "0.0075";
chr_lengths = load("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/plotting/image_analysis/chromosome_lengths.csv");
num_chr = length(chr_lengths);
for CPC = ["0.1","0.12","0.125", "0.15", "0.22","0.25","0.3", "0.35"]
    sampled_times = datasample(x, num_chr, 'Replace', false);

    Nx = 512;
    domain = 6.4;
    name = sprintf("phi_256_19661_1.0e-5__CPC_%s_cohesin_%s_eps_%s_domain_0_1", CPC, cohesin, eps);
    load(sprintf("%s/%s/trackedDroplets.mat", indir, name));
    dt = .1*(1/256)^2;
    num_timepoints = round(total_time/(10*dt))+1;
    all_distances = cell(1, num_timepoints);
    all_dist_list = [];
    for t = 1:num_timepoints
        time = (t-1)*dt*10;
        droplet_centers = [];
        for d = 1:max([trackedDroplets.id])
            if ismember(d, [2,3,4]) %ignore second IC droplet
            elseif ismembertol(time,trackedDroplets(d).time, 1e-3)
                if ismembertol(time, .001, 1e-3)
                    d
                end
                droplet_centers(end+1) = trackedDroplets(d).center(2);
            end
        end
        droplet_centers = sort(droplet_centers);
        distances = domain*diff(droplet_centers)/Nx;
        if ismembertol(time, [sampled_times], 1e-3)
            % name
            % disp(distances);
            % histogram(distances);
            % saveas(gcf,sprintf('%s/%s/distances_between_droplets_t_0.001.png', indir, name))

            fid = fopen('simulated_droplet_distances_e_0.0075_domain_0_1.csv', 'a+');
            fprintf(fid, '%s,%s,%s,%s,%s \r\n', CPC, cohesin, eps, time, mat2str(distances));
            fclose(fid);

            
        end
        % all_distances{t} = distances;
        % all_dist_list = [all_dist_list distances];
    end

    % histogram(all_dist_list);
    % saveas(gcf,sprintf('%s/%s/distances_between_droplets_all_time.png', indir, name))

end