 %% Determines time spent in each state, conditional on the tree structure. Also generates Figure S6

%final data table is storMat which is time in each state (columns 1:4 in
%state order G1/2H G2H LES and P1H), column 5 is number of generations
%along the trajectory, and column 6 is the state classification of
%the endpoint cell (the leaf of the trajectory).

clear all
load('KCA_datamatrix_withStateAssignments.mat')
btStrp = 500;
ctr = 0
load('colonyIDs.mat')

%colonyIDs = colonyIDs;   %ColIDsP, ColIDsG, colonyIDs,ColsLowG2
for mdl= 1:2        %compare the inferred model versus irreversible model
    for aIdx =1:btStrp    %repeat the backpropagation 1000 times to get error for the time spent in each state
        
        storMat = zeros(1,6);
        ctr = ctr+1
        for k = 1:length(colonyIDs)
            
            
            colI = colonyIDs(k);
            
            
            load(['Colony_' num2str(colI) '.mat']);
            
            clear nodes leavesSeq tmp
            
            
            nodes = unique(colony.lineage); %number of nodes in the tree)
            nodes(nodes == 0) =[]; %remove 0 (used for padding colony.lineage)
            
            p0 = [0 0 1 0]'; %initial distribution
            leaves = zeros(1,length(nodes)); %storage for states of the leaves
            clrs = zeros(length(nodes),3);
            if mdl == 1
                muA =[0.88    0.02    0.0000    0.0000;
                    0.12    0.86    0.025    0.0000;
                    0.0000    0.13    0.95    0.28;
                    0.0000    0.0000    0.022    0.73];
                stdA =[0.04    0.0008    0.0000    0.0000;
                    0.05    0.02    0.007    0.0000;
                    0.0000    0.02    0.008    0.03;
                    0.0000    0.0000    0.005    0.04];
            else
                muA =[0.88    0.02    0    0;
                    0    0.86    0.025    0;
                    0    0    0.95    0;
                    0    0    0.022    0.73];
                stdA =[0.04    0.0008    0.0000    0.0000;
                    0.05    0.02    0.007    0.0000;
                    0.0000    0.02    0.008    0.03;
                    0.0000    0.0000    0.005    0.04];
            end
            
            
            A = normrnd(muA,stdA);
            
            for jj =1:4
                A(:,jj) = A(:,jj) / sum(A(:,jj));
            end
            
            G = zeros(size(nodes));  %storage for directed graph of the pedigree. Rows are outputs from each node, columns are inputs
            tmp = KCA(KCA.colonyID ==colonyIDs(k),:);
            for i = 1:length(nodes)
                %first assign the state to each leaf in leavesSeq
                if ismember(nodes(i),colony.CellID)
                    if tmp.State(tmp.CellID == nodes(i)) == {'Low'}
                        leaves(i) = 3;
                        clrs(i,:) = [0.3 0.3 0.3];
                    elseif tmp.State(tmp.CellID == nodes(i)) == {'G1/2H'}
                        leaves(i) = 1;
                        clrs(i,:) = [1 0 0];
                    elseif tmp.State(tmp.CellID == nodes(i)) == {'G2H'}
                        leaves(i) = 2;
                        clrs(i,:) = [0.8 0.6 0.8];
                    elseif tmp.State(tmp.CellID == nodes(i)) == {'P1H'}
                        leaves(i) = 4;
                        clrs(i,:) = [0 0 1];
                    else
                        leaves(i) = 0;
                        clrs(i,:) = [1 1 0.00];
                    end
                else
                    leaves(i) = 0;
                    clrs(i,:) = [0.75 0.75 0.75];
                end
                
                
                %then do the outputs from each node
                outIds = [nodes(i)*2, nodes(i)*2+1];
                isANode = double(ismember(outIds,nodes));
                for j =1:2
                    if isANode(j) == 1
                        [r,~] = find(nodes == outIds(j));
                        G(i,r) = 1;
                    else
                    end
                end
                
                %finally do the inputs
                if mod(nodes(i),2) == 0
                    inId = (nodes(i))/2;
                    [ri,~] = find(nodes == inId);
                    G(ri,i) = 1;
                elseif nodes(i) ~= 1
                    inId = (nodes(i)-1)/2;
                    [ri,~] = find(nodes == inId);
                    G(ri,i) = 1;
                else
                end
            end
            es = treeBackTrace2(G,A,leaves,p0);
            
            if sum(isnan(es) == 0)
                es(:,5) = sum(es,2);
                nu = horzcat(es,leaves');
                storMat= vertcat(storMat,nu);
                
            else
            end
            
            
            
            
            
%             if aIdx == 1
%                 figure(1)
%                 subplot(15,8,k)
%                 P = graph(G,'upper');
%                 plot(P,'NodeColor',clrs,'MarkerSize',3,'NodeLabel',[],'EdgeColor','k')
%                 
%                 set(gca,'visible','off')
%                 set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96],'color','w');
%             else
%             end
            
            
        end
        xs = find(storMat(:,5)>4 & storMat(:,6) ==1);
        g12h = storMat(xs,:);
        G12Hcells(aIdx,:) = mean(g12h);
        
        xs = find(storMat(:,5)>4 & storMat(:,6) ==2);
        g2h = storMat(xs,:);
        G2Hcells(aIdx,:) = mean(g2h);
        
        xs = find(storMat(:,5)>4 & storMat(:,6) ==3);
        les = storMat(xs,:);
        LEScells(aIdx,:) = mean(les);
        
        xs = find(storMat(:,5)>4 & storMat(:,6) ==4);
        p1h = storMat(xs,:);
        P1Hcells(aIdx,:) = mean(p1h);
    end
    
    
    %% This will generate the figure comparing time spent in each state for reversible and irreversible chains
    
    G12Hcells_time= ((G12Hcells(:,1:4)./G12Hcells(:,5)).*96);
    G2Hcells_time= ((G2Hcells(:,1:4)./G2Hcells(:,5)).*96);
    LEScells_time= ((LEScells(:,1:4)./LEScells(:,5)).*96);
    P1Hcells_time= ((P1Hcells(:,1:4)./P1Hcells(:,5)).*96);
    
   
    x = 1:19; %4 per state with 5, 11
    mu = [mean(G12Hcells_time), 0, mean(G2Hcells_time), 0, mean(LEScells_time), 0, mean(P1Hcells_time)];
    sdevs = [std(G12Hcells_time), 0, std(G2Hcells_time), 0, std(LEScells_time), 0, std(P1Hcells_time)];

    figure(1)
    subplot(1,2,mdl)
    bar(x,mu)
    hold on
    er = errorbar(x,mu, 2*sdevs);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    hold off
    
    
    
end

