% 
% Edited so that eigenvalues can be extracted to separate aligned(+ve eigenvalues) and liberal(-ve eigenvalues) eigenvectors
% Then we calculate the number of eigenvectors (k) that can give us 80%, 90% and 95% variance explained
% Do this for all subjects then find group-level average number of k for 80%, 90% and 95% variance explained
% For longitudinal subjects, k is first averaged within subject --> this average is used for averaging at group level

clear;clc

p_dim = 126; % nROIs in parcellation


for i = 1:length(subjects)
    subj_idx = find(unique_idx == i);
    for j = 1:length(subj_idx) %suject has multiple visits (longituidnal data); so here calculate for each visit of each subject
        phase = all_subjects{subj_idx(j)}(end);

        load('/xx/xx/p_mat.mat');%load SC: p_dim * p_dim

        A = p_mat;

        [V, D] = eig(A); % eigendecomposition on A %% we only want this step from this portion

        % This part makes sure the eigenvalues are associated with ascending
        % eigenvalues. %% modifications start from here

        [lambda{i}(:,j), order] = sort(diag(D), 'ascend');
        V = V(:, order);

        % use eigenvalues to define align/liberal eigenvectors (there should not be lambda = 0)
        align_idx = find(lambda{i}(:,j) > 0);
        liberal_idx = find(lambda{i}(:,j) < 0);

        % get total variance explained for aligned and liberal eigenvectors (separately) --> sum to zero?
        aligned_var = sum(lambda{i}(align_idx,j));
        liberal_var = sum(abs(lambda{i}(liberal_idx,j)));

        % cumulative sum of variance explained for aligned and liberal eigenvectors (separately) from largest to smallest variance explained
        aligned_sum_var{i}{j} = cumsum(sort(lambda{i}(align_idx,j),'descend')); % need to reverse order for aligned signals to get largest to smallest
        liberal_sum_var{i}{j} = cumsum(abs(lambda{i}(liberal_idx,j)));

        k_aligned_subj(j,1) = min(find(aligned_sum_var{i}{j} >= 0.8*aligned_var));
        k_aligned_subj(j,2) = min(find(aligned_sum_var{i}{j} >= 0.9*aligned_var));
        k_aligned_subj(j,3) = min(find(aligned_sum_var{i}{j} >= 0.95*aligned_var));
        k_liberal_subj(j,1) = min(find(liberal_sum_var{i}{j} >= 0.8*liberal_var));
        k_liberal_subj(j,2) = min(find(liberal_sum_var{i}{j} >= 0.9*liberal_var));
        k_liberal_subj(j,3) = min(find(liberal_sum_var{i}{j} >= 0.95*liberal_var));

    end

    % for each subject, sum across timepoints (because CDE subjects have longitudinal data)
    k_aligned(i,:) = mean(k_aligned_subj,1);
    k_liberal(i,:) = mean(k_liberal_subj,1);

    clear k_aligned_subj k_liberal_subj

end


% generate some figures for visualisation
figure; boxplot(k_aligned); title('k aligned'); set(gca,'XTickLabel',{'80% var','90% var','95% var'});
figure; boxplot(k_liberal); title('k liberal'); set(gca,'XTickLabel',{'80% var','90% var','95% var'});

figure; hold on; title('aligned cumulative variance normalised')
for i = 1:length(aligned_sum_var) % subject
    for j = 1:length(aligned_sum_var{i}) % phase
        plot(aligned_sum_var{i}{j}/max(aligned_sum_var{i}{j}))
    end
end
figure; hold on; title('liberal cumulative variance normalised')
for i = 1:length(liberal_sum_var) % subject
    for j = 1:length(liberal_sum_var{i}) % phase
        plot(liberal_sum_var{i}{j}/max(liberal_sum_var{i}{j}))
    end
end

