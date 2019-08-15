function [nuc_mask,bad_mask] = excludelargeandwarped(nuc_mask,boulderarea)
% antimask=bwareaopen(nuc_mask,boulderarea);
% nuc_mask=nuc_mask-antimask;
nuc_label=bwlabel(nuc_mask);
bad_label = nuc_label;
nuc_props=regionprops(nuc_label>0,'Solidity','Area');
nuc_solidity = [nuc_props.Solidity];
nuc_area = [nuc_props.Area];
cellnum = numel(nuc_area);
warpedobjects= find(nuc_solidity<0.9 | nuc_area>=boulderarea); %default 10xbin1:0.9
% goodobjects=find(nuc_solidity>=0.9 | nuc_area<boulderarea); 
goodobjects = setdiff(1:cellnum,warpedobjects);
bad_label(ismember(bad_label,goodobjects)) = 0;
bad_mask = bad_label>0;
nuc_label(ismember(nuc_label,warpedobjects))=0;
nuc_mask=nuc_label>0;
