% reg_list
% keeped 115 regions
load reg_keep.mat
% single units in these regions
reg_list=deblank(h5read('transient_6.hdf5','/reg'));
% allen CCF structure tree
structure_tree_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\structure_tree_safe_2017.csv';
st = loadStructureTree(structure_tree_location);

keep = reg_set(1:115); % discard white matter, and unlabeled
neuron_num = cell(length(keep),1);
name = cell(length(keep),1);
for i = 1:length(keep)
    acronym = keep{i};
    neuron_num{i} = nnz(cellfun(@(x) strcmp(x,acronym),reg_list));
    name{i} = st.name{cellfun(@(x) strcmp(x,acronym),st.acronym)};
end
% construct table
T = table(name,keep,neuron_num,'VariableNames',{'name','acronym','neuron number'});
writetable(T,'reg_tabel.csv');


