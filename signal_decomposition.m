
% this script performs the signal vector decomposition using graph signal
% processing approach (Graph Fourier Transform) given Adjacency matrix

%%%change this:
p_dim = 126; % dimension of probability matrix (nROIs)

%%%%aligned means coupled; liberal means decoupled
%%%change this:
choose_aligned_range = 12; % number of most aligned (coupled) values to extract 
choose_liberal_range = 67; % number of least aligned (most liberal) values to extract 

frequency_range_A_align = (p_dim-choose_aligned_range+1):(p_dim);
frequency_range_A_liberal = 1:choose_liberal_range;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%signal decomposition:
for subj_num = 1:length(subjects)

    load('/xx/xx/myTS.mat');%rsfMRI timeseries
    load('/xx/xx/p_mat.mat');%SC
    A=p_mat; % nROIs * nROIs
    x=myTS; % nROIs * timepoints(length)

    %%%%%%
    frequency_range=frequency_range_A_align;

    [V, D] = eig(A); % eigendecomposition on A
    % This part makes sure the eigenvalues are associated with ascending
    % eigenvalues.
    %
    [~, order] = sort(diag(D), 'ascend');
    V = V(:, order);

    % This part defines the graph frequency filter, which is a diagonal matrix
    % with entires 1 on frequency_range, and 0 elsewhere.
    %
    n = size(A, 1);
    filter = zeros(n, 1);
    filter(frequency_range) = 1;
    filter = diag(filter); % this is the filter in graph frequency domain
    filter_nodespace = V * filter * V'; % this is the filter in graph node domain

    % This part performs the decomposition using the graph filter applied onto
    % x, resulting in the output x_decomposed
    %
    x_decomposed = filter_nodespace * x;
    aligned_values(:,:,subj_num) = x_decomposed;
    %%%aligned_values: Roi*timepoints*subjects

    %%%%%%%%%
    frequency_range=frequency_range_A_liberal;

    [V, D] = eig(A); % eigendecomposition on A
    % This part makes sure the eigenvalues are associated with ascending
    % eigenvalues.
    %
    [~, order] = sort(diag(D), 'ascend');
    V = V(:, order);

    % This part defines the graph frequency filter, which is a diagonal matrix
    % with entires 1 on frequency_range, and 0 elsewhere.
    %
    n = size(A, 1);
    filter = zeros(n, 1);
    filter(frequency_range) = 1;
    filter = diag(filter); % this is the filter in graph frequency domain
    filter_nodespace = V * filter * V'; % this is the filter in graph node domain

    % This part performs the decomposition using the graph filter applied onto
    % x, resulting in the output x_decomposed
    %
    x_decomposed = filter_nodespace * x;
    liberal_values(:,:,subj_num) = x_decomposed;
    %%%liberal_values: Roi*timepoints*subjects
end

% average aligned and liberal values across timepoints/TRs
aligned_avg = squeeze(mean(aligned_values,2))'; %%for each subject, each ROI has a aligned value
liberal_avg = squeeze(mean(liberal_values,2))'; %%for each subject, each ROI has a liberal value

%%%you can average aligned and liberal values across whole brain or within
%%%each network, for the following statistical analysis