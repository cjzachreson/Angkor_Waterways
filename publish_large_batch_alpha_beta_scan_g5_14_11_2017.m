% main script for network dynamics simulation

% Code and model by Cameron Zachreson, University of
% Sydney, civil engineering, complex systems unit.

% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

% This code simulates the erosion and sedimentation of earthen water
% channel networks. A phenomenological approach approximates fluid velocity
% as the channel flow divided by the channel capacity. Sedimentation or
% erosion occur when fluid velocity falls below or above the sedimentation
% or erosion thresholds, respectively. Erosion increases the capacity of
% the eroded edge, sedimentation decreases it. The functions describing
% this update system are scripted in the Edge class.

% the model runs as follows:

% 1) the adjacency matrix describing the network is used to create lists of
% edges and nodes.

% 2) the topology and edge weights are used to determine equilibrium flow
% distribution. These flows are then assigned as initial edge capacities.
% This can also be interpreted as setting an initial fluid velocity of 1
% for all edges.

% 3) gamma-distributed random flows are added to each edge. The scale
% parameter is the flood magnitude multiplied by the initial capacity
% (equilibrium flow)

% 4) edge flows are re-calculated with the flood perturbations added.
% In addition to the random additions to each edge, source flow is increased
% by a factor of gamma multiplied by the equilibrium input flow.

% 5) Once the flow in the network converges (sink input equals total flow),
% the fluid velocities and erosion/sedimentation updates are computed.
% Using these new edge weights, the flow distribution is
% re-calculated from source to sink (flood perturbation ignored), and
% compared to the initial distribution to calculate the normalized change
% in flow distribution on each edge (dF).

% 6) the process of flooding the network and updating capacities as
% appropriate continues until dF stops changing (up to a pre-set
% convergence threshold). During these iterations, the same random flood
% configuration is applied in each round, to ensure convergence.

%%%%%%%%

%initialization. set path to adj. matrix and set up alpha beta sweep

path(path,'C:\Users\Cameron\Documents\MATLAB\flood_model_final_14_08_2017\Adjacency_matrices')

clearvars

numreps = 100 % number of times to sweep through alpha and beta values (different random perturbation each time)

flood_factor = 5; %determines mean of fluctuation distribution on each edge (this is the 'gamma' parameter)

% fluctuations are gamma distributed:
for k_gam = [1]%2, 1, 0.5, 0.25]; % shape parameter controls whether the distribution is long tailed or , scale parameter is defined later to give the desired mean fluctuation
% (for gamma distribution, k*theta = mean = flood_factor * f_initial)

d_alpha = 0.05; % alpha interval
d_beta = 0.05; % beta interval

alpha_max = 1 - d_alpha;
alpha_min = 0 + d_alpha;
beta_max = 1 - d_alpha;
beta_min = 0 + d_alpha;

alphas = alpha_min : d_alpha : alpha_max; %erosion threshold = 1/alpha
betas =  beta_min : d_beta : beta_max; %sedimentation threshold = beta

Q_average = zeros(numel(alphas), numel(betas)); %output matrices for the current repetition
tf_average = zeros(numel(alphas), numel(betas));

%set up output figure axis
num_x_ticks = 10;
num_y_ticks = 10;
x_delta = ceil(numel(betas) / num_x_ticks);
y_delta = ceil(numel(alphas)/num_y_ticks);
x_tick_indices = [x_delta:x_delta:numel(betas)];
y_tick_indices = [y_delta:y_delta:numel(alphas)];
x_axis_labels = betas(x_tick_indices);
y_axis_labels = alphas(y_tick_indices);
%*****


%initialize network: create edges and nodes by going through
%adjacency matrix element by element and making a new edge for each entry.


%pick an adjacency matrix
adj_mat = Adjacency_Matrix_Angkor_directed_ST_CZ('Angkor_Adj_mat_ST_CZ_20_09_2017.txt');
%adj_mat = Adjacency_Matrix_Cdi_21;
%adj_mat = Adjacency_Matrix_simple_loop;
%adj_mat = Adjacency_Matrix_py_N_val(10);
%adj_mat = Adjacency_Matrix_simple_4_val;
%adj_mat = Adjacency_Matrix_di_21_val;

% if a real network is used with virtual super source and super sink,
% exclude the source and sink edges from the damage calculation:
ST_exclusion_flag = 1;

% set output root
root = fullfile('output',  '14_11_2017', 'flood_damage', adj_mat.name,...
    ['gamma_', num2str(flood_factor)], ['k_' num2str(k_gam)], ['n_reps_' num2str(numreps)]);


convergence_threshold = 5.0 * 10^(-3);
% set convergence threshold, this is the percentage change in dF
% that will signal convergence of network damage


tic

parfor rep = 1:numreps
    
    output_dir = fullfile(root, ['rep_' num2str(rep)]);
    
    if ~isdir(output_dir)
        mkdir(output_dir)
    end
    
    Q_mat = zeros(numel(alphas), numel(betas)); % holds damage values for each alpha and beta set
    
    tf_mat = zeros(numel(alphas), numel(betas));
    
    source_node = 1; %for initialization, source node is super source
    sink_node = adj_mat.N;
    
    flood_flows = []; % this vector will contain the random perturbations associated with a flood configuration
    
    for alpha_i = 1:numel(alphas)
        for beta_i = 1:numel(betas)
            
            alpha = alphas(alpha_i);
            beta = betas(beta_i);
            
            % organize output directories
            % *******
            dirname_edge_data = fullfile(output_dir, 'edge_data');
            if ~isdir(dirname_edge_data)
                mkdir(dirname_edge_data)
            end
            
            %********
            
            % creates lists of edges and nodes from the selected adjacency
            % matrix.
            [edges, nodes] = initialize_network_valued(adj_mat, alpha, beta);
            
            % set an arbitrary source flow input value (this will be
            % re-scaled later so that all flows are > 1, a step that is
            % only important if distributions used for floods require
            % integer inputs or inputs > 1, for gamma-distributed floods
            % this is not necessary).
            
            f_eq = 1;
            
            nodes{source_node}.inFlow = 0;
            nodes{source_node}.outFlow = f_eq;
            nodes{sink_node}.inFlow = 0;
            
            % starting from super-source, divides flow based on relative edge weights
            [edges, nodes] = initialize_flows_loops_ok...
                (adj_mat, edges, nodes, source_node, sink_node, f_eq);
            
            % tabulate initial flows for comparison to damaged topology
            initial_flows = zeros(1, numel(edges));
            for eij = 1:numel(edges)
                initial_flows(eij) = edges{eij}.tot_flow_t;
            end
            
            sink_edges = nodes{sink_node}.inEdges;
            source_edges = nodes{source_node}.outEdges;
            
            % if the edges of the super source and super sink do not
            % represent physical channels, exclude them from the damage
            % calculation here.
            if ST_exclusion_flag == 1
                excluded_edges = [source_edges, sink_edges];
                initial_flows(excluded_edges) = NaN;
                %initial_flows = initial_flows(~isnan(initial_flows));
            else
                excluded_edges = [];
            end
            
            %***********
            % rescale all flows to the minimum flow, this is only necessary
            % if the flood distribution requires flow values greater than
            % 1, this is the case if the flood is poisson distributed for
            % example.
            
            min_flow = f_eq; % placeholder, larger than any flow
            
            for eij = 1:numel(edges)
                if edges{eij}.tot_flow_t < min_flow
                    min_flow = edges{eij}.tot_flow_t;
                end
            end
            
            initial_flows = initial_flows / min_flow;
            
            f_eq = f_eq / min_flow;
            
            for eij = 1:numel(edges)
                edges{eij}.tot_flow_t = edges{eij}.tot_flow_t / min_flow;
                edges{eij}.initial_flow = edges{eij}.initial_flow / min_flow;
                edges{eij}.flow = edges{eij}.flow / min_flow;
            end
            
            for n_i = 1:numel(nodes)
                nodes{n_i}.inFlow = nodes{n_i}.inFlow / min_flow;
                nodes{n_i}.outFlow = nodes{n_i}.outFlow / min_flow;
            end
            %***********
            
            
            % compute capacities, velocities, and erosion/sedimentation thresholds
            for eij = 1:numel(edges)
                edges{eij}.initial_flow = edges{eij}.tot_flow_t;
                edges{eij}.initial_capacity = edges{eij}.tot_flow_t;
                edges{eij}.capacity = edges{eij}.tot_flow_t;
                edges{eij} = edges{eij}.set_velocity(edges{eij});
                edges{eij} = edges{eij}.set_thresholds(edges{eij});
            end
            
            % output initial equilibrium state as text file
            filename = fullfile(dirname_edge_data, 'eq_graph.txt');
            if rep == 1
                save_edge_list(edges, f_eq, filename)
            end
            
            
            % with equilibrium flows and thresholds set, dynamics can be initialized
            
            % reset flows
            [edges, nodes] = reset_flows(edges, nodes);
            
            %*****
            % this is where random flood flows are computed. Random
            % perturbations (delta > 0) are added to each edge of the network
            % to simulate a flood scenario.
            
            % apply flood, if it's the first run, find the random values,
            % otherwise keep them constant for all values of alpha and beta
            % tested
            
            flood_edges = 1:numel(edges);

            %flow from super-source increases by the flood factor, this
            %reflects the tendency for all rainfall in the watershed region
            %to increase
            f_tot = f_eq + f_eq * flood_factor;
            nodes{source_node}.outFlow = f_eq + f_eq * flood_factor;
            
            if isempty(flood_flows)
                for fe_i = flood_edges
                    
                    % lambda is average fluctuation
                    lambda = edges{fe_i}.initial_flow * flood_factor;
                    
                    %flood_flow = rand() * lambda;%uniform random fluctuations
                    %flood_flow = poissrnd(lambda);%poisson random fluctuation from equilibrium
                    
                    % gamma-distributed
                    theta = lambda / k_gam;
                    flood_flow = gamrnd(k_gam, theta);
                    
                    % f_tot is used for flow convergence detection so it
                    % must increase if flow is added as flood fluctuations
                    f_tot = f_tot + flood_flow;
                    
                    flood_flows = [flood_flows, flood_flow];
                    edges{fe_i}.tot_flow_t = flood_flow;
                    nodes{edges{fe_i}.outNode}.inFlow = ...
                        nodes{edges{fe_i}.outNode}.inFlow + flood_flow;
                    nodes{edges{fe_i}.outNode}.outFlow = ...
                        nodes{edges{fe_i}.outNode}.outFlow + flood_flow;
                end
            else
                index = 0;
                for fe_i = flood_edges
                    index = index + 1;
                    flood_flow = flood_flows(index);
                    f_tot = f_tot + flood_flow;
                    edges{fe_i}.tot_flow_t = flood_flow; % the equilibrium flow will be added during the flow distribution loop
                    nodes{edges{fe_i}.outNode}.inFlow = ...
                        nodes{edges{fe_i}.outNode}.inFlow + flood_flow;
                    nodes{edges{fe_i}.outNode}.outFlow = ...
                        nodes{edges{fe_i}.outNode}.outFlow + flood_flow;
                end
            end
            
            if alpha == alphas(1) && beta == betas(1)
                filename = fullfile(dirname_edge_data, 'flood_graph.txt');
                save_edge_list_minimal(edges, f_tot, filename, source_node, sink_node)
            end
            
            
            convergence_flag = 0; % switch for converged dF (damage) values
            
            in_nodes = source_node; %nodes providing input (sending flow)
            
            out_nodes = [];
            
            Q_this = 0;
            
            iterations = 0;
            
            tic
            while ~ convergence_flag
                
                Q_last = Q_this;
                
                % loop through nodes sending flow (current_nodes)
                
                [edges, nodes, out_nodes] = distribute_flow(adj_mat, edges, nodes, in_nodes);
                
                in_nodes = unique(out_nodes);
                
                test_val =  abs(nodes{sink_node}.inFlow - f_tot) / f_tot;
                
                % if flows have converged, calculate damage
                if test_val < 1 * 10^(-12)
                    
                    iterations = iterations + 1;
                    
                    % recalculate capacities
                    for eij = 1:numel(edges)
                        % velocity changes with flow
                        edges{eij} = edges{eij}.set_velocity(edges{eij});
                    end
                    
                    for eij = 1:numel(edges)
                        %update capacities based on new velocities
                        edges{eij} = edges{eij}.update_capacity(edges{eij});
                    end
                    
                    % reset capacities on super source edges, as their
                    % topology should remain fixed (super source does not
                    % represent a physical node)
                    if ST_exclusion_flag
                        for e_i = source_edges
                            edges{e_i}.capacity = edges{e_i}.initial_capacity;
                        end
                    end
                    
                    % calculate average flow differential between current
                    % state and initial state
                    %**********
                    
                    % duplicate the edges and nodes with zero-ed flow, to
                    % test equilibrium behavior of damaged network.
                    [edges_tst, nodes_tst] = reset_flows(edges, nodes);
                    
                    nodes_tst{source_node}.inFlow = 0;
                    nodes_tst{source_node}.outFlow = f_eq;
                    nodes_tst{sink_node}.inFlow = 0;
                    
                    [edges_tst, nodes_tst] = test_flows_loops_ok(adj_mat, edges_tst, nodes_tst, source_node, sink_node, f_eq);
                    
                    test_flows = zeros(1, numel(edges));
                    for eij = 1:numel(edges)
                        test_flows(eij) = edges_tst{eij}.tot_flow_t;
                    end
                    
                    test_flows(excluded_edges) = NaN;
                    
                    Q_i = zeros(1, numel(test_flows)); Q_i(excluded_edges) = NaN;
                    
                    depleted_flows = test_flows(test_flows <= initial_flows);
                    init_depleted = initial_flows(test_flows <= initial_flows);
                    
                    flooded_flows = test_flows(test_flows > initial_flows);
                    init_flooded = initial_flows(test_flows > initial_flows);
                    
                    Q_i(test_flows > initial_flows) = (flooded_flows - init_flooded)./(f_eq - init_flooded);
                    Q_i(test_flows <= initial_flows) = (init_depleted - depleted_flows)./init_depleted;
                    
                    Q_this = sum(Q_i(~isnan(Q_i)) .* initial_flows(~isnan(initial_flows)))...
                        / sum(initial_flows(~isnan(initial_flows)));
                    
                    if abs(Q_last - Q_this)/Q_this < convergence_threshold || Q_this == 0
                        
                        convergence_flag = 1; % stops iteration loop
                        
                        % save graph of converged flows for visualization
                        % and post-processing of edge data
                        Q_label = round(Q_this * 10) / 10;
                        Q_label = num2str(Q_label);
                        
                        filename = fullfile(dirname_edge_data, ...
                            ['state_graph_a_ ' num2str(alpha) '_b_' num2str(beta)...
                            '_tf_' num2str(iterations) '_Q_' Q_label '.txt']);
                        save_edge_list_minimal(edges_tst, f_eq, filename, source_node, sink_node)
                        
                        Q_mat(alpha_i, beta_i) = Q_this;
                        
                        tf_mat(alpha_i, beta_i) = iterations;
                        
                        disp(['tf = ' num2str(iterations), ' Q = ' num2str(Q_this), ' alpha = ' num2str(alpha), ' beta = ' num2str(beta)])
                        
                        toc
                        
                    end
                    
                    % reset network for next iteration
                    [edges, nodes] = reset_flows(edges, nodes);
                    
                    % continue applying flood
                    index = 0;
                    for fe_i = flood_edges
                        index = index + 1;
                        flood_flow = flood_flows(index);
                        edges{fe_i}.tot_flow_t = flood_flow;
                        nodes{edges{fe_i}.outNode}.inFlow = ...
                            nodes{edges{fe_i}.outNode}.inFlow + flood_flow;
                        nodes{edges{fe_i}.outNode}.outFlow = ...
                            nodes{edges{fe_i}.outNode}.outFlow + flood_flow;
                    end
                    
                    nodes{source_node}.outFlow = f_eq + f_eq * flood_factor;
                    
                    in_nodes = source_node;
                    
                end
            end
        end
    end
    
    a = figure;
    subplot(1, 2, 1)
    imagesc(Q_mat)
    ax = gca;
    ax.XTick = x_tick_indices;
    ax.YTick = y_tick_indices;
    ax.XTickLabel = x_axis_labels;
    ax.YTickLabel = y_axis_labels;
    ax.DataAspectRatio = [d_alpha/d_beta, 1, 1];
    title('average relative deviation Q')
    xlabel('\beta')
    ylabel('\alpha')
    colormap gray
    
    subplot(1, 2, 2)
    imagesc(tf_mat)
    ax = gca;
    ax.XTick = x_tick_indices;
    ax.YTick = y_tick_indices;
    ax.XTickLabel = x_axis_labels;
    ax.YTickLabel = y_axis_labels;
    ax.DataAspectRatio = [d_alpha/d_beta, 1, 1];
    
    ct_pow = floor(log10(convergence_threshold))
    ct_base = convergence_threshold / (1 * 10^ct_pow);
    title(['convergence time $(dQ <' num2str(ct_base) '\times 10^{' num2str(ct_pow)...
        '})$'], 'Interpreter', 'Latex')
    
    xlabel('\beta')
    ylabel('\alpha')
    colormap gray
    
    
    saveas(a, fullfile(output_dir, ['rep_' num2str(rep) '_Q.fig']))
    dlmwrite(fullfile(output_dir, ['rep_' num2str(rep) '_Q.txt']), Q_mat)
    dlmwrite(fullfile(output_dir, ['rep_' num2str(rep) '_conv_t.txt']), tf_mat)
    
    
end


toc


for rep = 1:numreps
    
    data_dir = fullfile(root, ['rep_' num2str(rep)]);
    
    filename = fullfile(data_dir, ['rep_' num2str(rep) '_Q.txt']);
    Q_rep = dlmread(filename);
    
    filename = fullfile(data_dir, ['rep_' num2str(rep) '_conv_t.txt']);
    tf_rep = dlmread(filename);
    
    Q_average = Q_average + Q_rep;
    tf_average = tf_average + tf_rep;
    
end

Q_average = Q_average / numreps;

b = figure;

imagesc(Q_average)
ax = gca;
ax.XTick = x_tick_indices;
ax.YTick = y_tick_indices;
ax.XTickLabel = x_axis_labels;
ax.YTickLabel = y_axis_labels;
ax.DataAspectRatio = [d_alpha/d_beta, 1, 1];
title('average relative deviation Q')
xlabel('\beta')
ylabel('\alpha')
colormap gray

saveas(b, fullfile(root, ['Q_avg_N_' num2str(numreps) '.fig']))
dlmwrite(fullfile(root, ['Q_avg_N_' num2str(numreps) '.txt']), Q_average)

tf_average = tf_average / numreps;

c = figure;
imagesc(tf_average)
ax = gca;
ax.XTick = x_tick_indices;
ax.YTick = y_tick_indices;
ax.XTickLabel = x_axis_labels;
ax.YTickLabel = y_axis_labels;
ax.DataAspectRatio = [d_alpha/d_beta, 1, 1];
ct_pow = floor(log10(convergence_threshold));
ct_base = convergence_threshold / (1 * 10^ct_pow);
title(['convergence time $(df <' num2str(ct_base) '\times 10^{' num2str(ct_pow)...
    '})$'], 'Interpreter', 'Latex')
xlabel('\beta')
ylabel('\alpha')
colormap gray

saveas(c, fullfile(root, ['tf_avg_N_' num2str(numreps) '.fig']))
dlmwrite(fullfile(root, ['tf_avg_N_' num2str(numreps) '.txt']), tf_average)


end





