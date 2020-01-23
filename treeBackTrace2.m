function es= treeBackTrace2(G,A,leavesSeq,p0)
%This function computes the expected time spent in each state for each
%lineage, conditional on the entire final population and tree structure.
%
%Inputs:
%G= the (directed) graph structure of the tree, in adjacency matrix form (nxn),
%   where node 1 is the root node, n is the total number of nodes in the tree.
%A= the markov chain (4x4)
%leavesSeq = the states of all the leaf nodes, formatted as a (1xn) vector, with
%             0s in the indices of non-leaf nodes.
%p0 = the probability distribution of the initial state (a 4x1 vector)
%
%Output:
%es = an (nx4) matrix, where 
    n=length(G);
    es=zeros(n,4);
    p=0;
    [ess,ps]=treeBackTraceHelper(G,A,leavesSeq);
    for i=1:4
        es=es+p0(i)*ess(:,:,i);
        p=p+p0(i)*ps(i);
    end
    es=es/p;
end

function [es,ps,usedInds]=treeBackTraceHelper(G,A,leavesSeq)
    
    n=size(G,1);
    leavesinds=leavesSeq>0;
    es=zeros(n,4,4);
    ps=zeros(1,4);
    usedInds=[];
    if leavesSeq(1)
        mystate=leavesSeq(1);
        ps(mystate)=1;
        es(1,mystate,mystate)=1;
        usedInds=1;
    else    
        children=find(G(1,:));
        nchild=length(children);
        childes=zeros(n,4,4,nchild);
        childps=zeros(4,nchild);
        for j=1:nchild
            child=children(j); 
            subinds=(child:n);
            subG=G(subinds,subinds);
            subSeq=leavesSeq(subinds);
            [myes,myps,myUsed]=treeBackTraceHelper(subG,A,subSeq);
            childes(subinds,:,:,j)=myes;
            childps(:,j)=myps;
            usedInds=[usedInds,subinds(myUsed)];
        end
        startSeqs=getstartSeqs(children,leavesSeq);
        usedLeaves=usedInds(leavesSeq(usedInds)>0);
        for state0=1:4
            for i=1:size(startSeqs,1)
                childseq=startSeqs(i,:);
                myps0=zeros(1,nchild);
                myps1=zeros(1,nchild);
                for j=1:nchild
                    myps0(j)=A(childseq(j),state0);
                    myps1(j)=childps(childseq(j),j);
                end
                mytotalp=prod(myps0)*prod(myps1);
                if mytotalp>0
                    ps(state0)=ps(state0)+mytotalp;
                    for j=1:nchild
                        myes=childes(:,:,childseq(j),j);
                        backgroundP=mytotalp/myps1(j);
                        myes=myes*backgroundP;                
                        es(:,:,state0)=es(:,:,state0)+myes;
                    end
                    es(usedLeaves,state0,state0)=es(usedLeaves,state0,state0)+mytotalp;
                end
            end
        end
    end
end

function startSeqs = getstartSeqs(children,leavesSeq)
    startSeqs=[];
    for child=children
        if leavesSeq(child)
            states=leavesSeq(child);
        else
            states=(1:4);
        end
        if isempty(startSeqs)
            startSeqs=states';
        else
            l=size(startSeqs,1);
            newSeqs=[];
            for state=states
                newSeqs=[newSeqs;startSeqs,state*ones(l,1)];
            end
            startSeqs=newSeqs;
        end
    end
end