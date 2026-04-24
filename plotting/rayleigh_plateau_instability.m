% Calculate the characteristic length between droplets compared to width of cohesion stripe
% and the droplet radius compared to width of cohesion stripe to evaluate whether Rayleigh-Plateau instability occurs.
% Sarah Groves April 17, 2026
% Load in the trackedDroplets.mat file that is generated for a simulation when level_set_radius_multiple_droplets.m
% is run.
indir = "/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/plotting/radii_lineplots_kymographs/domain_0_2_e_0.0067/";
CPC = "0.35";
for cohesin = ["0.1", "0.09", "0.08"]
    eps = "0.0067";
    Nx = 512;
    domain_to_nm = 6.4;
    grid_spacing = domain_to_nm/Nx;
    midpoint = Nx/2;
    name = sprintf("phi_%d_19661_1.0e-5__CPC_%s_cohesin_%s_eps_%s_domain_0_2",Nx,CPC, cohesin, eps);
    load(sprintf("%s/%s/trackedDroplets.mat", indir, name));
    droplet_centers = [];
    droplet_radii = [];

    % droplet_centers(end+1) = trackedDroplets(1).center(2);

    for d = [7,8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28] %only keep 7,8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28 because those are the ones that are truly circular droplets in cohesin = 0.09
        % grab first nonzero radius value for each droplet, which is the radius at the first timepoint when the droplet is detected. This is because the radius can fluctuate over time due to noise, but we want to capture the characteristic radius of the droplet when it first forms.
        droplet_centers(end+1) = trackedDroplets(d).center(2);
        droplet_radii(end+1) = trackedDroplets(d).radius(find(trackedDroplets(d).radius>0,1));

    end
    droplet_centers = sort(droplet_centers)
    distances = domain_to_nm*diff(droplet_centers)/Nx;
    % get rid of max of distances because that is the distance between the two droplets on either side of the central droplet 
    distances = distances(distances<max(distances))
    cohesion_radius_stripe = str2double(cohesin)/2; %width of cohesion stripe is the fraction of the domain that has cohesin multiplied by the total domain size
    characteristic_length = mean(distances);
    characteristic_length_over_width_cohesion_stripe = characteristic_length/cohesion_radius_stripe;
    droplet_radii*domain_to_nm
    mean_radius = mean(droplet_radii)*domain_to_nm;
    radius_over_width_cohesion_stripe = mean_radius/cohesion_radius_stripe;
    fprintf('Simulation parameters: CPC = %s, cohesin = %s, eps = %s\n', CPC, cohesin, eps);
    fprintf('Characteristic length between droplets: %f \n', characteristic_length);
    fprintf('Characteristic length between droplets / cohesion radius stripe: %f \n', characteristic_length_over_width_cohesion_stripe);
    fprintf('Mean droplet radius: %f \n', mean_radius);
    fprintf('Mean droplet radius over cohesion radius stripe: %f \n', radius_over_width_cohesion_stripe); 

    % append to a file
    fid = fopen('rayleigh_plateau_instability_results.csv', 'a+');
    fprintf(fid, '%s,%s,%s,%f,%f,%f,%f,%s \r\n', CPC, cohesin, eps, characteristic_length, characteristic_length_over_width_cohesion_stripe, mean_radius, radius_over_width_cohesion_stripe, name);
    fclose(fid);
    % droplets that are truly circular droplets in cohesin = 0.09: data 7,8, 11, 12, 15, 16, 19, 20, 23, 24, 27, 28
end