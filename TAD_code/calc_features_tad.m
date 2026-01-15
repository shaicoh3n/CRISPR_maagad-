%resolution changes from cell_type to cell_type
% T-100K res
% one can find out res by dividing length of chromosome by number of bins

function calc_features_tad(sites_name)

    switch sites_name
        case {'T','K562','U937','ICS','Tomato'}
            reso = 100000;
    end
    sites = readtable(sprintf("C://Users//shaic//Desktop//PhD//CRISPR_maagad//code//features//TAD_code//sites//%s_w_isana.csv",sites_name)); 
    if ~strcmp(sites_name,'Tomato')
        % Split each string by '_'
        split_data = cellfun(@(x) strsplit(x, '_'), sites.g_rna_info, 'UniformOutput', false);
        
        % Convert to structure
        sites.chrm = cellfun(@(x) strcat('chr',x{1}), split_data, 'UniformOutput', false); % Chromosome
        sites.start_site = cellfun(@(x) str2double(x{2}), split_data);       % Start site
        sites.end_site = cellfun(@(x) str2double(x{3}), split_data);         % End site
        sites.strand = cellfun(@(x) x{4}, split_data, 'UniformOutput', false); % Strand
    else
        split_data = cellfun(@(x) strsplit(x, ':'), sites.g_rna_id, 'UniformOutput', false);
        sites.chrm = cellfun(@(x) x{7}, split_data, 'UniformOutput', false);      
        sites.start_site = cellfun(@(x) str2double(x{8}), split_data);       
        sites.end_site = sites.start_site+20;
        %sites.strand exists already
    end
    if strcmp(sites_name,'U937')
        load(sprintf('tad_%s.mat','THP'),'tad');
    elseif strcmp(sites_name,'ICS')
        load(sprintf('tad_%s.mat','iPSC'),'tad');
    elseif strcmp(sites_name,'Tomato')
        load(sprintf('tad_%s.mat','leaf_tomato'),'tad');
    else
        load(sprintf('tad_%s.mat',sites_name),'tad');
    end
    num_particle = [ceil(sites.start_site / reso), ceil(sites.end_site / reso)];
    for i = 1:height(sites)
        if mod(i,50)==0
            disp(i)
        end
        curr_particle = num_particle(i,:);
        if strcmp(sites_name,'Tomato')
            inds = [find(strcmp(sites.chrm{i}, tad.chr) & tad.part == curr_particle(1)),...
                    find(strcmp(sites.chrm{i}, tad.chr) & tad.part == curr_particle(2))];
        else
            inds = [find(strcmp(sites.chrm{i}, strcat('chr',tad.chr)) & tad.part == curr_particle(1)),...
                    find(strcmp(sites.chrm{i}, strcat('chr',tad.chr)) & tad.part == curr_particle(2))];
        end
        if isempty(inds)
            sites.tad_density(i) = 0;
            sites.tad_angle(i) =  0;
            sites.tad_interactions(i) =  0;
        else
            sites.tad_density(i) = mean([tad.density(inds(1)),tad.density(inds(2))]);
            sites.tad_angle(i) =  mean([tad.angle(inds(1)),tad.angle(inds(2))]);
            sites.tad_interactions(i) =  mean([tad.interactions(inds(1)),tad.interactions(inds(2))]);
        end
    end
    writetable(sites,sprintf("C://Users//shaic//Desktop//PhD//CRISPR_maagad//code//features//TAD_code//sites//%s_w_tad.csv",sites_name));
    disp('done')
end