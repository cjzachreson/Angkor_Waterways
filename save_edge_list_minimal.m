% Copyright (C) August 28th, 2018. Cameron Zachreson
% This work is licensed under the Creative Commons Attribution 4.0
% International (CC BY 4.0)
% [https://creativecommons.org/licenses/by/4.0/legalcode]
% Any use of this code, or the associated data must be accompanied by
% citation of the work entitled:
% "The demise of Angkor: systemic vulnerability of urban infrastructure to 
% climatic variations" by Dan Penny et al. (2018)

% makes an edge list for import to Cytoscape for visualization

function save_edge_list_minimal(edges, f_tot, filename, Super_source, Super_sink)

edge_list = fopen(filename, 'w');

fprintf(edge_list, ['source \t target \t directed \t initial_flow \t'...
                    'flow_t \t rel_flow \t ID \t source_edge \t sink_edge \n']);


for i = 1:numel(edges)
    
    source = num2str(edges{i}.inNode);
    
    target = num2str(edges{i}.outNode);
    
    init_flow = num2str(edges{i}.initial_flow);
    
    flow_t = num2str(edges{i}.tot_flow_t);
    
    rel_flow = num2str(edges{i}.tot_flow_t / f_tot);
    
    ID = num2str(i);
    
    if edges{i}.inNode ~= Super_source && edges{i}.outNode ~= Super_sink
        source_edge = 'FALSE';
        sink_edge = 'FALSE';
    else if edges{i}.inNode == Super_source
            source_edge = 'TRUE';
            sink_edge = 'FALSE';
        else
            source_edge = 'FALSE';
            sink_edge = 'TRUE';
        end
    end
    
    
    fprintf(edge_list, [source '\t' target '\t' 'TRUE' '\t' init_flow...
                        '\t' flow_t '\t' rel_flow '\t' ID '\t' source_edge...
                        '\t' sink_edge '\n']);
end

fclose(edge_list);
    
    
    