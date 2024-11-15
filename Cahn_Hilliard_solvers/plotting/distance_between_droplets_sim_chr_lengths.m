% Script to calculate the distance between droplets at each timestep for
% Ostwald ripening simulations. This uses the trackedDroplets.mat file that
% is generated for a simulation when level_set_radius_multiple_droplets.m
% is run. 
% Sarah Groves August 28, 2024
% indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/domain_0_2_noisy_cohesin_sd_0.11/";
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_model/Cahn_Hilliard_solvers/plotting/radii_over_time_level_set_plots/domain_0_2_from_rivanna_kymographs_e_0.0075/";

total_time = 0.03;
spacing = 0.000001525878906;
x = 0:spacing:0.03;
%name = "phi_512_19660_1.0e-5__CPC_20_cohesin_8_eps_0.007504684956431058";
cohesin = "0.08";
eps = "0.0075";
chr_lengths = load("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_model/Cahn_Hilliard_solvers/plotting/image_analysis/chromosome_lengths.csv");
num_chr = length(chr_lengths);
% seed = "1111";
% for CPC = "0.12"
for CPC = ["0.1","0.12","0.125", "0.15", "0.22","0.25","0.3", "0.35"]
    sampled_times = datasample(x, num_chr, 'Replace', false);

    Nx = 512;
    domain = 6.4;
    grid_spacing = domain/Nx;
    midpoint = Nx/2;
    name = sprintf("phi_512_19661_1.0e-5__CPC_%s_cohesin_%s_eps_%s_alpha_0_domain_0_2",CPC, cohesin, eps);
    load(sprintf("%s/%s/trackedDroplets.mat", indir, name));
    dt = .1*(1/256)^2;
    for i = 1:num_chr
        t = sampled_times(i);
        arm_length = chr_lengths(i);

        droplet_centers = [];
        %check if each droplet exists at time t and if so, add to
        %droplet_centers
        for d = 1:max([trackedDroplets.id])
            if ismember(d, [2]) %ignore second IC droplet because it is in the same y location as droplet 1
            elseif ismembertol(t,trackedDroplets(d).time, 1e-3)
                if trackedDroplets(d).center(2) <= midpoint+(arm_length/grid_spacing) && trackedDroplets(d).center(2) >= midpoint-(arm_length/grid_spacing)
                    droplet_centers(end+1) = trackedDroplets(d).center(2);
                end
            end
        end
        droplet_centers = sort(droplet_centers);
        distances = domain*diff(droplet_centers)/Nx;
        fid = fopen('simulated_droplet_distances_e_0.0075_domain_0_2_chr_lengths.csv', 'a+');
        fprintf(fid, '%s,%s,%s,%f,%f,%s \r\n', CPC, cohesin, eps, t,arm_length, mat2str(distances));
        fclose(fid);
         

    end


end

