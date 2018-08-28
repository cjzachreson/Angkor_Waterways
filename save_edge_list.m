% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)


% makes an edge list for import to Cytoscape for visualization

function save_edge_list(edges, f_tot, filename)

edge_list = fopen(filename, 'w');

fprintf(edge_list, ['source \t target \t directed \t initial_flow \t change_in_flow \t transparency_flow \t'...
    'rel_flow \t norm_damage \t norm_abs_damage \t transparency_damage \t norm_rel_damage'...
    '\t transp_rel_damage \t rel_capacity \t ID \n']);


rel_C = [];
abs_dam = [];
rel_abs_dam = [];

for i = 1:numel(edges)
    rel_C = [rel_C edges{i}.capacity / edges{i}.initial_capacity];
    abs_dam = [abs_dam, abs(edges{i}.tot_flow_t - edges{i}.initial_flow)];
    rel_abs_dam = [rel_abs_dam, abs(edges{i}.tot_flow_t - edges{i}.initial_flow) / edges{i}.initial_flow];
end

dam_max = max(abs_dam);
rel_dam_max = max(rel_abs_dam);
rel_C_max = max(rel_C);


for i = 1:numel(edges)
    source = num2str(edges{i}.inNode);
    target = num2str(edges{i}.outNode);
    rel_flow = num2str(edges{i}.tot_flow_t / f_tot); % 0 to 1
    transparency_flow = num2str(edges{i}.tot_flow_t / f_tot * 255); % 0 to 255
    
    damage = edges{i}.tot_flow_t - edges{i}.initial_flow; %unbounded
    init_flow = edges{i}.initial_flow;
    norm_abs_damage = abs_dam(i)/dam_max; %0 to 1
    transparency_damage = norm_abs_damage * 255; %0 to 255
    
    norm_rel_damage = (damage / edges{i}.initial_flow) / rel_dam_max; % -1 to 1
    transp_rel_damage = (rel_abs_dam(i) / rel_dam_max) * 255; % 0 to 255
    
    rel_capacity = (edges{i}.capacity ./ edges{i}.initial_capacity) / rel_C_max; % 0 to 1
    
    norm_damage = num2str(damage / dam_max);
    
    norm_abs_damage = num2str(norm_abs_damage);
    
    transparency_damage = num2str(transparency_damage);
    
    norm_rel_damage = num2str(norm_rel_damage);
    
    transp_rel_damage = num2str(transp_rel_damage);
    
    rel_capacity = num2str(rel_capacity);
    
    ID = num2str(i);
    
    change_in_flow = num2str(damage);
    
    init_flow = num2str(init_flow);
    
    fprintf(edge_list, [source, '\t', target, '\t', 'TRUE', '\t', init_flow, '\t' change_in_flow '\t' transparency_flow,...
        '\t'  rel_flow, '\t' norm_damage, '\t' norm_abs_damage '\t' transparency_damage,...
        '\t' norm_rel_damage, '\t' transp_rel_damage, '\t', rel_capacity '\t' ID '\n']);
end

fclose(edge_list);
    
    
    