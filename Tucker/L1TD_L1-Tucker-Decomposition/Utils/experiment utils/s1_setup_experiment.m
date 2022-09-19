function [] = s1_setup_experiment(D,d,G_std,parent_folder)
%% Select base experiment folder, for storing the simulation results
    
    if ~exist(parent_folder,'dir')
        mkdir(parent_folder)
    end
    
    G_mu = 0;
    I = length(D);
    G_true = tensor(normrnd(G_mu,G_std,d));
    Un_true = generate_orth_basis(I,D,d);
    file_path = fullfile(parent_folder,'init.mat');
    save(file_path,'D','d','G_true','G_std','Un_true');
end