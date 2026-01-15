files = dir(); % Get all files and folders in the current directory

for i = 1:length(files)
    filename = files(i).name;
    if contains(filename, '.csv') && ~files(i).isdir
        disp(filename); % Print the file name
        sites = readtable(files(i).name);
        foldings = cellfun(@oligoprop, strcat(sites.up_seq, sites.target_seq, sites.down_seq));
        foldings = {foldings(:).Thermo};
        foldings = cellfun(@(x) x(:,3),foldings,'UniformOutput',false);
        sites.DNA_fold = transpose(cellfun(@mean,foldings));
        writetable(sites,strcat(filename(1:end-4),'_w_DNAfold.csv'));
    end
end

for i=1:166
    try
        if ~isequal(T_post{:,i},T_pre{:,i})
            disp(T_post.Properties.VariableNames{i})
        end
    catch
        sprintf('bad %s',T_post.Properties.VariableNames{i})
    end
end
disp('end')
% foldings = cellfun(@oligoprop, strcat(sites.up_seq, sites.target_seq, sites.down_seq));
% foldings = {foldings(:).Thermo};
% foldings = cellfun(@(x) x(:,3),foldings,'UniformOutput',false);
% sites.DNA_fold = transpose(cellfun(@mean,foldings));





% for i=10:20
% %     seqs = cellfun(@(x,y,z), strcat(x(end-i,end), y, z(1:i),);
%     up_seqs_i = cellfun(@(x) x(end-i:end), up_seqs,'uni',0);
%     down_seqs_i = cellfun(@(x) x(1:j), up_seqs,'uni',0);
%     foldings = cellfun(@oligoprop, strcat(up_seqs_i, target_seqs, down_seqs_i));
%     foldings = {foldings(:).Thermo};
%     foldings = cellfun(@(x) x(:,3),foldings,'UniformOutput',false);
%     res = transpose(cellfun(@mean,foldings));
%     res=res(1);
%     sprintf('%f %f %d %d',res, a_sites.DNA_fold(1) ,i,j)
% end

i=17
j=17
up_seqs_i = cellfun(@(x) x(end-i:end), up_seqs,'uni',0);
down_seqs_i = cellfun(@(x) x(1:j), up_seqs,'uni',0);
foldings = cellfun(@oligoprop, strcat(up_seqs_i, target_seqs, down_seqs_i));
foldings = {foldings(:).Thermo};
foldings = cellfun(@(x) x(:,3),foldings,'UniformOutput',false);
res = transpose(cellfun(@mean,foldings));
res=res(1);
sprintf('%f %f %d %d',res, a_sites.DNA_fold(1) ,i,j)