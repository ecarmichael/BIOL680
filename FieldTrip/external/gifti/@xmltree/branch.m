function subtree = branch(tree,uid)
% XMLTREE/BRANCH Branch Method
% FORMAT uid = parent(tree,uid)
% 
% tree    - XMLTree object
% uid     - UID of the root element of the subtree
% subtree - XMLTree object (a subtree from tree)
%_______________________________________________________________________
%
% Return a subtree from a tree.
%_______________________________________________________________________
% Copyright (C) 2002-2008  http://www.artefact.tk/

% Guillaume Flandin <guillaume@artefact.tk>
% $Id: branch.m 3261 2011-03-31 15:06:56Z roboos $

error(nargchk(2,2,nargin));

if uid > length(tree) || ...
   numel(uid)~=1 || ...
   ~strcmp(tree.tree{uid}.type,'element')
    error('[XMLTree] Invalid UID.');
end

subtree = xmltree;
subtree = set(subtree,root(subtree),'name',tree.tree{uid}.name);
%- fix by Piotr Dollar to copy attributes for the root node:
subtree = set(subtree,root(subtree),'attributes',tree.tree{uid}.attributes); 

child = children(tree,uid);

for i=1:length(child)
    l = length(subtree);
    subtree = sub_branch(tree,subtree,child(i),root(subtree));
    subtree.tree{root(subtree)}.contents = [subtree.tree{root(subtree)}.contents l+1];
end

%=======================================================================
function tree = sub_branch(t,tree,uid,p)

    l = length(tree);
    tree.tree{l+1} = t.tree{uid};
    tree.tree{l+1}.uid = l + 1;
    tree.tree{l+1}.parent = p;
    tree.tree{l+1}.contents = [];
    if isfield(t.tree{uid},'contents')
        contents = get(t,uid,'contents');
        m = length(tree);
        for i=1:length(contents)
            tree.tree{l+1}.contents = [tree.tree{l+1}.contents m+1];
            tree = sub_branch(t,tree,contents(i),l+1);
            m = length(tree);
        end
    end
