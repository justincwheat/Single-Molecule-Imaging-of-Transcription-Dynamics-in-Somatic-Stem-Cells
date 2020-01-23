%% Will produce all simulation panels of figure 2 of the manuscript
%Parameters for each state need to be pasted into the appropriate area
%below

clear all
load('KLdatamatrix.mat');


[r,c] = find(fn2.State=={'Low','PU1 High', 'Gata2 High','Gata1 High'});
fnCMP = fn2(r,:);
Tend = 10000;    %total simulation time
npath =1000;      %number of path trajectories to simulate
MasterStateSeq = zeros(Tend,npath,4);  %Storage for state seq
mRNAcounts = zeros(Tend,3,npath);
%% Cycle through and determine the fraction of time spent in each state
for st_ctr=1:4
    StateIndex =st_ctr;  %Low G2H G1H P1H
    
    
    %% DiffusionMap (don't run everytime)
    %     fnCMP = sortrows(fnCMP,'State','ascend');
    %     data = (table2array(fnCMP(:,[11,19,27])));
    %     sigma=50;
    %     n = 40;
    %     [T,phi0] = transition_matrix(data,'nn',[sigma, n]);
    %
    %
    %     root=[1788:1850]; %index of root cell(s)
    %     branching=1; %detect branching? 0/1
    %     num_tips=2;
    %
    %     [M, tips] = dpt_input(T, phi0, branching, 'maxdptdist', root);
    %     %  [M, tips] = dpt_input(T, phi0, branching, 'manual',num_tips);% labels;
    %     [Branch,DPT]=dpt_analyse(M,branching,tips);
    %
    %     [phi, lambda] = diffusionmap.eig_decompose_normalized(T,3);
    %
    %     fnCMP.DPT= DPT;
    %     fnCMP.phiX = phi(:,3)+abs(min(phi(:,3)));
    %     fnCMP.phiY = phi(:,2)+abs(min(phi(:,2)));
    %
    %     figure(1)
    %     idL = find(fnCMP.State == 'Low');
    %     idG2 = find(fnCMP.State == 'Gata2 High');
    %     idG1 = find(fnCMP.State == 'Gata1 High');
    %     idP1 = find(fnCMP.State == 'PU1 High');
    %
    %     scatter(fnCMP.phiX(idL),fnCMP.phiY(idL),'MarkerFaceAlpha',0.2,'Marker','o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
    %     hold on
    %     scatter(fnCMP.phiX(idG2),fnCMP.phiY(idG2),'MarkerFaceAlpha',0.2,'Marker','o','MarkerFaceColor',[1 0.6 0.78],'MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
    %     hold on
    %     scatter(fnCMP.phiX(idG1),fnCMP.phiY(idG1),'MarkerFaceAlpha',0.2,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
    %     hold on
    %     scatter(fnCMP.phiX(idP1),fnCMP.phiY(idP1),'MarkerFaceAlpha',0.2,'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerEdgeAlpha',0.1);
    %     hold off
    %
    %     [gataIds,~] = find(fnCMP.State =={'Low','Gata2 High','Gata1 High','MEP'});
    %
    % fng = fnCMP(gataIds,:);
    % fng = sortrows(fng,'DPT','ascend');
    % gTSbursts = cell(3,1);
    % gTSbursts{1} = fng.p1TS;
    % gTSbursts{2} = fng.g2TS;
    % gTSbursts{3} = fng.g1TS;
    %
    % [pu1Ids,~] = find(fnCMP.State =={'Low','PU1 High','GMP'});
    %
    % fnp = fnCMP(pu1Ids,:);
    % fnp = sortrows(fnp,'DPT','ascend');
    %
    % pTSbursts = cell(3,1);
    % pTSbursts{1} = fnp.p1TS;
    % pTSbursts{2} = fnp.g2TS;
    % pTSbursts{3} = fnp.g1TS;
    %
    % figure(1)
    % subplot(3,1,1)
    % plot(pTSbursts{1},'Color','b')
    % subplot(3,1,2)
    % plot(pTSbursts{2},'Color','k')
    % subplot(3,1,3)
    % plot(pTSbursts{3},'Color','r')
    %
    % figure(2)
    % subplot(3,1,1)
    % plot(gTSbursts{1},'Color','b')
    % subplot(3,1,2)
    % plot(gTSbursts{2},'Color','k')
    % subplot(3,1,3)
    % plot(gTSbursts{3},'Color','r')
    
    
    %% Parameters for the Low State
    clear fn3
    idsim=find(fnCMP.State == 'Low');
    fn3=fnCMP(idsim,:);
    
    ParamsL =   [0.0104	0.0163	0.0125;
        0.0831	0.1027	0.372;
        0.3153	0.0387	1.967;
        0.0039	0.0058	0.023];
    %PU1    %Gata1  %Gata2
    
    activP=find(fn3.p1TS>=1);
    activG1=find(fn3.g1TS>=1);
    activG2=find(fn3.g2TS>=1);
    
    %Compute the average burst size (for those >= 1TS/gene) and TS freq
    bL=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
    bL(isnan(bL)) = 0;
    fL=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
    ICL=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];
    
    %% Parameters for the Gata2 High State
    clear fn3
    idsim=find(fnCMP.State == 'Gata2 High');
    fn3=fnCMP(idsim,:);
    
    ParamsG2 =   [0.009	0.0034	0.1026
        0.0799	0.2188	0.1247
        0.3189	0.8649	0.9534
        0.0039	0.0058	0.0231];
    %PU1    %Gata1  %Gata2
    activP=find(fn3.p1TS>=1);
    activG1=find(fn3.g1TS>=1);
    activG2=find(fn3.g2TS>=1);
    
    %Compute the average burst size (for those >= 1TS/gene) and TS freq
    bG2=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
    bG2(isnan(bG2)) = 0;
    fG2=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
    ICG2=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];
    
    %% Parameters for the Gata1 High State
    clear fn3
    idsim=find(fnCMP.State == 'Gata1 High');
    fn3=fnCMP(idsim,:);
    
    ParamsG1 =   [0.006	0.0172	0.3268;
        0.0704	0.0604	0.1995;
        0.321	1	1.7;
        0.0039	0.0058	0.0231];
    %PU1    %Gata1  %Gata2
    
    activP=find(fn3.p1TS>=1);
    activG1=find(fn3.g1TS>=1);
    activG2=find(fn3.g2TS>=1);
    
    %Compute the average burst size (for those >= 1TS/gene) and TS freq
    bG1=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
    bG1(isnan(bG1)) = 0;
    fG1=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
    ICG1=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];
    
    %% Parameters for the PU1 High State
    clear fn3
    idsim=find(fnCMP.State == 'PU1 High');
    fn3=fnCMP(idsim,:);
    
    ParamsP1 =   [0.0468	0.000782	0.0135;
        0.068	0.0333	0.191;
        0.5	0.2933	1.0015;
        0.0039	0.0058	0.0231];
    %PU1    %Gata1  %Gata2
    
    activP=find(fn3.p1TS>=1);
    activG1=find(fn3.g1TS>=1);
    activG2=find(fn3.g2TS>=1);
    
    %Compute the average burst size (for those >= 1TS/gene) and TS freq
    bP1=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
    bP1(isnan(bP1)) = 0;
    fP1=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
    ICP1=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];
    
    %% Storage for all parameters and average burst and frequencies
    Allparams = zeros(4,3,4);
    
    Allparams(:,:,1) = ParamsL;
    Allparams(:,:,2) = ParamsG2;
    Allparams(:,:,3) = ParamsG1;
    Allparams(:,:,4) = ParamsP1;
    
    Allb = vertcat(bL',bG2',bG1',bP1');
    Allf = vertcat(fL',fG2',fG1',fP1');
    IC = vertcat(ICL,ICG2,ICG1,ICP1);
    
    %% Initialize the storage matrices
    
    Params = Allparams(:,:,StateIndex);
    b = [mean(fnCMP.p1b1(fnCMP.p1b1 >2 )), mean(fnCMP.g1b1(fnCMP.g1b1 >2 )), mean(fnCMP.g2b1(fnCMP.g2b1 >2 ))];
    f = Allf(StateIndex,:);
    
    gamma=Params(4,:); %mRNA half lives
    
    timeEvolutionT=zeros(Tend,npath);   %Stores time for each path
    timeEvolutionS=zeros(Tend,3,npath);      %storage vector for the time
    timeEvolutionG=zeros(Tend,3,npath);     % Stores state of the promoters
    timeEvolutionU=zeros(Tend,3,npath);     % Stores state of the promoters
    
    %% GSSA using exact reaction path
    gmax=2;
    rxnSeq = zeros(Tend,1,npath);
    
    for n=1:npath
        timeEvolutionS(1,:,n) = IC(st_ctr,:);   %Starting conditions for the trajectory
        timeEvolutionG(1,:,n)=[0 0 0];
        
        for t=1:Tend-1
            
            
            %Calculte Propensities
            kon=[Params(1,1)*(gmax-timeEvolutionG(t,1,n)),Params(1,2)*(gmax-timeEvolutionG(t,2,n)),Params(1,3)*(gmax-timeEvolutionG(t,3,n))];
            koff=[Params(2,1)*(timeEvolutionG(t,1,n)),Params(2,2)*(timeEvolutionG(t,2,n)),Params(2,3)*(timeEvolutionG(t,3,n))];
            kini=[Params(3,1),Params(3,2),Params(3,3)];
            kr=[kini(1)*timeEvolutionG(t,1,n),kini(2)*timeEvolutionG(t,2,n),kini(3)*timeEvolutionG(t,3,n)];
            kd=gamma.*timeEvolutionS(t,:,n);
            
            %Sum of the propensities
            asum=sum(kon)+sum(koff)+sum(kr)+sum(kd);
            
            
            timeEvolutionT(t+1,n)=timeEvolutionT(t,n)+log(1/rand(1))/asum;
            %Exact GSSA
            mu=rand(1);
            
            if 0 <= mu && mu < kon(1)/asum        %PU1 on
                rxnSeq(t,1,n) = 1;
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n)+1;
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= (timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= (timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= (timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif kon(1)/asum <= mu && mu  < (kon(1) + kon(2))/asum %Gata1 on
                rxnSeq(t,1,n) = 2;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n)+1;
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= (timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= (timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= (timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (kon(1) + kon(2))/asum <= mu && mu  < (sum(kon))/asum %Gata2 on
                rxnSeq(t,1,n) = 3;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n)+1;
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon))/asum <= mu && mu  < (sum(kon)+koff(1))/asum %PU1 off
                rxnSeq(t,1,n) = 4;
                
                timeEvolutionG(t+1,1,n)= max(timeEvolutionG(t,1,n)-1,0);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+koff(1))/asum <= mu && mu  < (sum(kon)+koff(1)+koff(2))/asum %Gata1 off
                rxnSeq(t,1,n) = 5;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= max(timeEvolutionG(t,2,n)-1,0);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+koff(1)+koff(2))/asum <= mu && mu  < (sum(kon)+sum(koff))/asum %Gata2 off
                rxnSeq(t,1,n) = 6;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= max(timeEvolutionG(t,3,n)-1,0);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
            elseif (sum(kon)+sum(koff))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1))/asum %PU1 burst
                rxnSeq(t,1,n) = 7;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n)+1);
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=1;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+sum(koff)+kr(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1)+kr(2))/asum %Gata1 burst
                rxnSeq(t,1,n) = 8;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n)+1);
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=1;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+sum(koff)+kr(1)+kr(2))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr))/asum %Gata2 burst
                rxnSeq(t,1,n) = 9;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n)+1);
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=1;
                
                
            elseif (sum(kon)+sum(koff)+sum(kr))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum %PU1 decay
                rxnSeq(t,1,n) = 10;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= max(timeEvolutionS(t,1,n)-1,0);
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1)+kd(2))/asum %Gata1 decay
                rxnSeq(t,1,n) = 11;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= max(timeEvolutionS(t,2,n)-1,0);
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            else  %Gata2 decay
                rxnSeq(t,1,n) = 12;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= timeEvolutionS(t,1,n);
                timeEvolutionS(t+1,2,n)= timeEvolutionS(t,2,n);
                timeEvolutionS(t+1,3,n)= max(timeEvolutionS(t,3,n)-1,0);
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            end
            
            
        end
        
        path_num=n
        
    end
    
    
    %%
    rpath = randi(npath);
    
    figure(1)
    subplot(1,4,st_ctr) %mRNA counts
    stairs(timeEvolutionT(:,rpath),timeEvolutionS(:,1,rpath),'LineWidth',1,'Color','b');
    hold on
    stairs(timeEvolutionT(:,rpath),timeEvolutionS(:,2,rpath),'LineWidth',1,'Color','r');
    hold on
    stairs(timeEvolutionT(:,rpath),timeEvolutionS(:,3,rpath),'LineWidth',1,'Color','k');
    hold off
    
    figure(2)   %TS activity
    subplot(3,4,st_ctr)
    plot(timeEvolutionG(:,1,rpath),'Color','b')
    subplot(3,4,st_ctr+4)
    plot(timeEvolutionG(:,2,rpath),'Color','r')
    subplot(3,4,st_ctr+8)
    plot(timeEvolutionG(:,3,rpath),'Color','k')
    
    
    cumMuRP1 = reshape(cumsum(timeEvolutionU(:,1,:),'forward'),[Tend,npath]);
    cumMuRG1 = reshape(cumsum(timeEvolutionU(:,2,:),'forward'),[Tend,npath]);
    cumMuRG2 = reshape(cumsum(timeEvolutionU(:,3,:),'forward'),[Tend,npath]);
    
    muP1(:,st_ctr) = mean(cumMuRP1,2);
    muG1(:,st_ctr) = mean(cumMuRG1,2);
    muG2(:,st_ctr) = mean(cumMuRG2,2);
    
    stdP1(:,st_ctr) = std(cumMuRP1,1,2);
    stdG1(:,st_ctr) = std(cumMuRG1,1,2);
    stdG2(:,st_ctr) = std(cumMuRG2,1,2);
    
    
    %% Assign cells to a state at each time point
    
    clear cState
    
    StateSeq=zeros(Tend,npath);
    histCs = zeros(4,npath);
    for i=1:size(timeEvolutionS,3)
        for j=1:size(timeEvolutionS,1)
            
            if timeEvolutionS(j,2,i) >= 10 && timeEvolutionS(j,1,i) <= 20
                StateSeq(j,i) = 1;
            elseif timeEvolutionS(j,1,i) >= 40 && timeEvolutionS(j,3,i)<50
                StateSeq(j,i) = 4;
            elseif timeEvolutionS(j,3,i)>25
                StateSeq(j,i) = 2;
            else
                StateSeq(j,i) = 3;
            end
        end
        %         h = histcounts(StateSeq(:,i),'Normalization','probability');
        %         histCs(1:numel(h),i) = h;
    end
    x= [1 2 3 4];
    colors = [0.5 0.5 0.5; 1 0.6 0.78; 1 0 0; 0 0 1];
    figure(3)
    subplot(1,4,st_ctr)
    %     bar(x,mean(histCs,2),'FaceColor',colors(st_ctr,:),'FaceAlpha',1,'LineWidth',2);
    %     hold on
    %     er = errorbar(x,mean(histCs,2),std(histCs,1,2));
    %     er.Color = [0 0 0];
    %     er.LineStyle = 'none'
    %     hold off
    histogram(StateSeq,'Normalization','probability','FaceColor',colors(st_ctr,:))
    
    MasterStateSeq(:,:,st_ctr) = StateSeq;
   
end
hold off

%% Cumulative mRNA

figure(4)

subplot(1,4,1)

plot(muG2(:,1),'LineWidth',2,'Color','k');
hold on
area(stdG2(:,1)+muG2(:,1),'FaceColor','k','FaceAlpha',0.4,'LineStyle','none','EdgeColor','k');
hold on;
area(muG2(:,1)-stdG2(:,1),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','k');
hold on

plot(muP1(:,1),'LineWidth',2,'Color','b');
hold on;
area(stdP1(:,1)+muP1(:,1),'FaceColor','b','FaceAlpha',0.4,'LineStyle','none','EdgeColor','b');
hold on;
area(muP1(:,1)-stdP1(:,1),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','b');
hold on

plot(muG1(:,1),'LineWidth',2,'Color','r');
hold on
area(stdG1(:,1)+muG1(:,1),'FaceColor','r','FaceAlpha',0.4,'LineStyle','none','EdgeColor','r');
hold on;
area(muG1(:,1)-stdG1(:,1),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','r');
hold off

subplot(1,4,2)
plot(muG2(:,2),'LineWidth',2,'Color','k');
hold on
area(stdG2(:,2)+muG2(:,2),'FaceColor','k','FaceAlpha',0.4,'LineStyle','none','EdgeColor','k');
hold on;
area(muG2(:,2)-stdG2(:,2),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','k');
hold on

plot(muP1(:,2),'LineWidth',2,'Color','b');
hold on;
area(stdP1(:,2)+muP1(:,2),'FaceColor','b','FaceAlpha',0.4,'LineStyle','none','EdgeColor','b');
hold on;
area(muP1(:,2)-stdP1(:,2),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','b');
hold on

plot(muG1(:,2),'LineWidth',2,'Color','r');
hold on
area(stdG1(:,2)+muG1(:,2),'FaceColor','r','FaceAlpha',0.4,'LineStyle','none','EdgeColor','r');
hold on;
area(muG1(:,2)-stdG1(:,2),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','r');
hold off

subplot(1,4,3)
plot(muG2(:,3),'LineWidth',2,'Color','k');
hold on
area(stdG2(:,3)+muG2(:,3),'FaceColor','k','FaceAlpha',0.4,'LineStyle','none','EdgeColor','k');
hold on;
area(muG2(:,3)-stdG2(:,3),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','k');
hold on

plot(muG1(:,3),'LineWidth',2,'Color','r');
hold on
area(stdG1(:,3)+muG1(:,3),'FaceColor','r','FaceAlpha',0.4,'LineStyle','none','EdgeColor','r');
hold on;
area(muG1(:,3)-stdG1(:,3),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','r');
hold on

plot(muP1(:,3),'LineWidth',2,'Color','b');
hold on;
area(stdP1(:,3)+muP1(:,3),'FaceColor','b','FaceAlpha',0.4,'LineStyle','none','EdgeColor','b');
hold on;
area(muP1(:,3)-stdP1(:,3),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','b');
hold off

subplot(1,4,4)
plot(muP1(:,4),'LineWidth',2,'Color','b');
hold on;
area(stdP1(:,4)+muP1(:,4),'FaceColor','b','FaceAlpha',0.4,'LineStyle','none','EdgeColor','b');
hold on;
area(muP1(:,4)-stdP1(:,4),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','b');
hold on

plot(muG2(:,4),'LineWidth',2,'Color','k');
hold on
area(stdG2(:,4)+muG2(:,4),'FaceColor','k','FaceAlpha',0.4,'LineStyle','none','EdgeColor','k');
hold on;
area(muG2(:,4)-stdG2(:,4),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','k');
hold on

plot(muG1(:,4),'LineWidth',2,'Color','r');
hold on
area(stdG1(:,4)+muG1(:,4),'FaceColor','r','FaceAlpha',0.4,'LineStyle','none','EdgeColor','r');
hold on;
area(muG1(:,4)-stdG1(:,4),'FaceColor','w','FaceAlpha',1,'LineStyle','none','EdgeColor','r');
hold off

%% Initialize in the LES state and transition through fluctuations alone to other states

clear Tend npath MasterStateSeq
Tend = 1000;    %total simulation time
npath =10000;      %number of path trajectories to simulate
MasterStateSeq = zeros(Tend,npath);  %Storage for state seq for each path
mRNAcounts = zeros(Tend,3,npath);

% Parameters for the Low State
clear fn3
idsim=find(fnCMP.State == 'Low');
fn3=fnCMP(idsim,:);

ParamsL =   [0.0104	0.0163	0.0125;
    0.0831	0.1027	0.372;
    0.3153	0.0387	1.967;
    0.0039	0.0058	0.023];
%PU1    %Gata1  %Gata2

activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
% bL=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
% bL(isnan(bL)) = 0;
fL=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICL=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Parameters for the Gata2 High State
clear fn3
idsim=find(fnCMP.State == 'Gata2 High');
fn3=fnCMP(idsim,:);

ParamsG2 =   [0.0129	0.0027	0.129;
    0.097	0.1343	0.1415;
    0.3895	0.667	1.059;
    0.0039	0.0058	0.023];
%PU1    %Gata1  %Gata2
activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
bG2=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
bG2(isnan(bG2)) = 0;
fG2=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICG2=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Parameters for the Gata1 High State
clear fn3
idsim=find(fnCMP.State == 'Gata1 High');
fn3=fnCMP(idsim,:);

ParamsG1 =   [0.006	0.0172	0.3268;
    0.0704	0.0604	0.1995;
    0.321	1	1.7;
    0.0039	0.0058	0.0231];
%PU1    %Gata1  %Gata2

activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
bG1=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
bG1(isnan(bG1)) = 0;
fG1=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICG1=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Parameters for the PU1 High State
clear fn3
idsim=find(fnCMP.State == 'PU1 High');
fn3=fnCMP(idsim,:);

ParamsP1 =   [0.0468	0.000782	0.0135;
    0.068	0.0333	0.191;
    0.5	0.2933	1.0015;
    0.0039	0.0058	0.0231];
%PU1    %Gata1  %Gata2

activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
bP1=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
bP1(isnan(bP1)) = 0;
fP1=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICP1=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Storage for all parameters and average burst and frequencies
Allparams = zeros(4,3,4);

Allparams(:,:,1) = ParamsG1;
Allparams(:,:,2) = ParamsG2;
Allparams(:,:,3) = ParamsL;
Allparams(:,:,4) = ParamsP1;

Allb = vertcat(bG1',bG2',bL',bP1');
Allf = vertcat(fG1',fG2',fL',fP1');
IC = vertcat(ICG1,ICG2,ICL,ICP1);

% Initialize the storage matrices

clear timeEvolutionG timeEvolutionU timeEvolutionS timeEvolutionT

timeEvolutionT=zeros(Tend,npath);   %Stores time for each path
timeEvolutionS=zeros(Tend,3,npath);      %storage vector for the time
timeEvolutionG=zeros(Tend,3,npath);     % Stores state of the promoters
timeEvolutionU=zeros(Tend,3,npath);     % Stores state of the promoters

% GSSA using exact reaction path
gmax=2;
rxnSeq = zeros(Tend,1,npath);

for n=1:npath
    timeEvolutionS(1,:,n) = [0 0 0];   %All start in LES
    timeEvolutionG(1,:,n)=[0 0 0];
    
    for t=1:Tend-1
        %determine the current state of the cell
        
        if timeEvolutionS(t,2,n) >= 10 && timeEvolutionS(t,1,n) <= 20
            MasterStateSeq(t,n) = 1;
        elseif timeEvolutionS(t,1,n) >= 40 && timeEvolutionS(t,3,n)<50
            MasterStateSeq(t,n) = 4;
        elseif timeEvolutionS(t,3,n)>25
            MasterStateSeq(t,n) = 2;
        else
            MasterStateSeq(t,n) = 3;
        end
        StateIndex = MasterStateSeq(t,n);
        %set the parameters for this time step
        
        Params = Allparams(:,:,StateIndex);
%         b = Allb(StateIndex,:);
        f = Allf(StateIndex,:);
        
        gamma=Params(4,:); %mRNA half lives
        
        %Calculte Propensities
        kon=[Params(1,1)*(gmax-timeEvolutionG(t,1,n)),Params(1,2)*(gmax-timeEvolutionG(t,2,n)),Params(1,3)*(gmax-timeEvolutionG(t,3,n))];
        koff=[Params(2,1)*(timeEvolutionG(t,1,n)),Params(2,2)*(timeEvolutionG(t,2,n)),Params(2,3)*(timeEvolutionG(t,3,n))];
        kini=[Params(3,1),Params(3,2),Params(3,3)];
        kr=[kini(1)*timeEvolutionG(t,1,n),kini(2)*timeEvolutionG(t,2,n),kini(3)*timeEvolutionG(t,3,n)];
        kd=gamma.*timeEvolutionS(t,:,n);
        
        %Sum of the propensities
        asum=sum(kon)+sum(koff)+sum(kr)+sum(kd);
        
        
        timeEvolutionT(t+1,n)=timeEvolutionT(t,n)+log(1/rand(1))/asum;
        %Exact GSSA
        mu=rand(1);
        
        if 0 <= mu && mu < kon(1)/asum        %PU1 on
            rxnSeq(t,1,n) = 1;
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n)+1;
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= (timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= (timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= (timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif kon(1)/asum <= mu && mu  < (kon(1) + kon(2))/asum %Gata1 on
            rxnSeq(t,1,n) = 2;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n)+1;
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= (timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= (timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= (timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (kon(1) + kon(2))/asum <= mu && mu  < (sum(kon))/asum %Gata2 on
            rxnSeq(t,1,n) = 3;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n)+1;
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon))/asum <= mu && mu  < (sum(kon)+koff(1))/asum %PU1 off
            rxnSeq(t,1,n) = 4;
            
            timeEvolutionG(t+1,1,n)= max(timeEvolutionG(t,1,n)-1,0);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+koff(1))/asum <= mu && mu  < (sum(kon)+koff(1)+koff(2))/asum %Gata1 off
            rxnSeq(t,1,n) = 5;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= max(timeEvolutionG(t,2,n)-1,0);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+koff(1)+koff(2))/asum <= mu && mu  < (sum(kon)+sum(koff))/asum %Gata2 off
            rxnSeq(t,1,n) = 6;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= max(timeEvolutionG(t,3,n)-1,0);
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
        elseif (sum(kon)+sum(koff))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1))/asum %PU1 burst
            rxnSeq(t,1,n) = 7;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n)+1);
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=1;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+sum(koff)+kr(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1)+kr(2))/asum %Gata1 burst
            rxnSeq(t,1,n) = 8;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n)+1);
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=1;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+sum(koff)+kr(1)+kr(2))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr))/asum %Gata2 burst
            rxnSeq(t,1,n) = 9;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n)+1);
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=1;
            
            
        elseif (sum(kon)+sum(koff)+sum(kr))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum %PU1 decay
            rxnSeq(t,1,n) = 10;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= max(timeEvolutionS(t,1,n)-1,0);
            timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1)+kd(2))/asum %Gata1 decay
            rxnSeq(t,1,n) = 11;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
            timeEvolutionS(t+1,2,n)= max(timeEvolutionS(t,2,n)-1,0);
            timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        else  %Gata2 decay
            rxnSeq(t,1,n) = 12;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS(t+1,1,n)= timeEvolutionS(t,1,n);
            timeEvolutionS(t+1,2,n)= timeEvolutionS(t,2,n);
            timeEvolutionS(t+1,3,n)= max(timeEvolutionS(t,3,n)-1,0);
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        end
        
        
    end
    
    path_num=n
    
end
%repeat without the state transitions

clear timeEvolutionG timeEvolutionU timeEvolutionS2 timeEvolutionT2

timeEvolutionT2=zeros(Tend,npath);   %Stores time for each path
timeEvolutionS2=zeros(Tend,3,npath);      %storage vector for the time
timeEvolutionG=zeros(Tend,3,npath);     % Stores state of the promoters
timeEvolutionU=zeros(Tend,3,npath);     % Stores state of the promoters

% GSSA using exact reaction path
gmax=2;
rxnSeq = zeros(Tend,1,npath);

for n=1:npath
    timeEvolutionS2(1,:,n) = [0 0 0];   %All start in LES
    timeEvolutionG(1,:,n)=[0 0 0];
    
    for t=1:Tend-1
        %determine the current state of the cell
        
        if timeEvolutionS2(t,2,n) >= 10 && timeEvolutionS(t,1,n) <= 20
            MasterStateSeq2(t,n) = 1;
        elseif timeEvolutionS2(t,1,n) >= 40 && timeEvolutionS(t,3,n)<50
            MasterStateSeq2(t,n) = 4;
        elseif timeEvolutionS2(t,3,n)>25
            MasterStateSeq2(t,n) = 2;
        else
            MasterStateSeq2(t,n) = 3;
        end
        %   StateIndex = MasterStateSeq2(t,n);
        %keep the parameters always in LES
        
        Params = Allparams(:,:,3);
%         b = Allb(3,:);
        f = Allf(3,:);
        
        gamma=Params(4,:); %mRNA half lives
        
        %Calculte Propensities
        kon=[Params(1,1)*(gmax-timeEvolutionG(t,1,n)),Params(1,2)*(gmax-timeEvolutionG(t,2,n)),Params(1,3)*(gmax-timeEvolutionG(t,3,n))];
        koff=[Params(2,1)*(timeEvolutionG(t,1,n)),Params(2,2)*(timeEvolutionG(t,2,n)),Params(2,3)*(timeEvolutionG(t,3,n))];
        kini=[Params(3,1),Params(3,2),Params(3,3)];
        kr=[kini(1)*timeEvolutionG(t,1,n),kini(2)*timeEvolutionG(t,2,n),kini(3)*timeEvolutionG(t,3,n)];
        kd=gamma.*timeEvolutionS2(t,:,n);
        
        %Sum of the propensities
        asum=sum(kon)+sum(koff)+sum(kr)+sum(kd);
        
        
        timeEvolutionT2(t+1,n)=timeEvolutionT2(t,n)+log(1/rand(1))/asum;
        %Exact GSSA
        mu=rand(1);
        
        if 0 <= mu && mu < kon(1)/asum        %PU1 on
            rxnSeq(t,1,n) = 1;
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n)+1;
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= (timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= (timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= (timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif kon(1)/asum <= mu && mu  < (kon(1) + kon(2))/asum %Gata1 on
            rxnSeq(t,1,n) = 2;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n)+1;
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= (timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= (timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= (timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (kon(1) + kon(2))/asum <= mu && mu  < (sum(kon))/asum %Gata2 on
            rxnSeq(t,1,n) = 3;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n)+1;
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon))/asum <= mu && mu  < (sum(kon)+koff(1))/asum %PU1 off
            rxnSeq(t,1,n) = 4;
            
            timeEvolutionG(t+1,1,n)= max(timeEvolutionG(t,1,n)-1,0);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+koff(1))/asum <= mu && mu  < (sum(kon)+koff(1)+koff(2))/asum %Gata1 off
            rxnSeq(t,1,n) = 5;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= max(timeEvolutionG(t,2,n)-1,0);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+koff(1)+koff(2))/asum <= mu && mu  < (sum(kon)+sum(koff))/asum %Gata2 off
            rxnSeq(t,1,n) = 6;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= max(timeEvolutionG(t,3,n)-1,0);
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
        elseif (sum(kon)+sum(koff))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1))/asum %PU1 burst
            rxnSeq(t,1,n) = 7;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n)+1);
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=1;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+sum(koff)+kr(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1)+kr(2))/asum %Gata1 burst
            rxnSeq(t,1,n) = 8;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n)+1);
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=1;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+sum(koff)+kr(1)+kr(2))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr))/asum %Gata2 burst
            rxnSeq(t,1,n) = 9;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n)+1);
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=1;
            
            
        elseif (sum(kon)+sum(koff)+sum(kr))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum %PU1 decay
            rxnSeq(t,1,n) = 10;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= max(timeEvolutionS2(t,1,n)-1,0);
            timeEvolutionS2(t+1,2,n)= ceil(timeEvolutionS2(t,2,n));
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        elseif (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1)+kd(2))/asum %Gata1 decay
            rxnSeq(t,1,n) = 11;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= ceil(timeEvolutionS2(t,1,n));
            timeEvolutionS2(t+1,2,n)= max(timeEvolutionS2(t,2,n)-1,0);
            timeEvolutionS2(t+1,3,n)= ceil(timeEvolutionS2(t,3,n));
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        else  %Gata2 decay
            rxnSeq(t,1,n) = 12;
            
            timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
            timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
            timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
            
            timeEvolutionS2(t+1,1,n)= timeEvolutionS2(t,1,n);
            timeEvolutionS2(t+1,2,n)= timeEvolutionS2(t,2,n);
            timeEvolutionS2(t+1,3,n)= max(timeEvolutionS2(t,3,n)-1,0);
            
            timeEvolutionU(t+1,1,n)=0;
            timeEvolutionU(t+1,2,n)=0;
            timeEvolutionU(t+1,3,n)=0;
            
        end
        
        
    end
    
    path_num=n
    
end

%Plot a histogram of the resulting endpoint distribution
figure(5)

histogram(MasterStateSeq2(end-1,:),'Normalization','probability','FaceColor','k','FaceAlpha',0.8,'LineWidth',2)
hold on
histogram(MasterStateSeq(end-1,:),'Normalization','probability','FaceColor',[0.35 0.35 0.35],'FaceAlpha',0.8,'LineWidth',2)
hold off
%% 
%%Determine the probability of not co-expressing both instantaneously and
%%integreated over time
figure(6)
ProbIN = zeros(npath,1);
TotalT = zeros(npath,1);
for j =1:npath
    for i = 1:Tend-1
        if timeEvolutionS2(i,1,j)>=1 && timeEvolutionS2(i,2,j)>=1 && timeEvolutionS2(i,3,j)>=1
            ProbIN(j,1) = timeEvolutionT2(i,j);
            break
        else
        end
        if i == Tend-1
            ProbIN(j,1) = 1000;
        end
    end
end
ProbIN(ProbIN >720) = 1000;



for j =1:npath
    clear Y dT iDx X
    X = timeEvolutionT2(:,j)>720;
    iDx = find(X,1,'first');
    dT = diff(timeEvolutionT2(1:iDx,j));
    Y = zeros(iDx-1,1);
    
    for i = 1:iDx-1
             if timeEvolutionS2(i,1,j)>=1 && timeEvolutionS2(i,2,j)>=1 && timeEvolutionS2(i,3,j)>=1
                Y(i) = 1;  
             else
             end
    end
        TotalT(j,1) = sum(dT(logical(Y)))/720;
        
end
subplot(2,1,1)
histogram(ProbIN./60,'Normalization','probability','FaceAlpha',1,'FaceColor',[0.5 0.5 0.5], 'LineWidth',1)
subplot(2,1,2)
histogram(TotalT,'Normalization','probability','FaceAlpha',1,'FaceColor',[0.5 0.5 0.5], 'LineWidth',1)
%% 

%Plot a trajectory plot with associated marginal distribution of states
%as a function of the PU1-(Gata1+Gata2)
figure (7)
Ridx = randi(npath,1000,1);

clrs = [1,0,0; 0.8 0.4 0.6; 0.5 0.5 0.5; 0 0 1];
for ctr = 1:1000
    kk = Ridx(ctr);
    deltaM = timeEvolutionS(:,1,kk)-(timeEvolutionS(:,2,kk)+ timeEvolutionS(:,3,kk));
    State = MasterStateSeq(end-1,kk);
    if State ~= 3
         stairs(deltaM,'Color',clrs(State,:),'LineWidth',2)
        hold on
    else
       stairs(deltaM,'Color',clrs(State,:),'LineWidth',0.15)
       hold on
    end       %for a plot versus cumulative time, use
    %stairs(timeEvolutionT(:,kk),deltaM,'Color',clrs(State,:),'LineWidth',1)
    hold on
end

hold off
%% 

%create histograms of state 
     figure (8)
    Alldeltas = timeEvolutionS(end-1,1,:)-(timeEvolutionS(end-1,2,:)+ timeEvolutionS(end-1,3,:));
    histogram(Alldeltas,'Normalization','probability','FaceColor','k', 'LineWidth',1,'NumBins',50)
%     camroll(-90)

     figure (9)
    Alldeltas2 = timeEvolutionS2(end-1,1,:)-(timeEvolutionS2(end-1,2,:)+ timeEvolutionS2(end-1,3,:));
    histogram(Alldeltas2,'Normalization','probability','FaceColor','k', 'LineWidth',1,'NumBins',50)
%     camroll(-90)

figure(10)
subplot(1,3,3)
scatter(timeEvolutionS(end-1,1,:),(timeEvolutionS(end-1,2,:)+ timeEvolutionS(end-1,3,:)),'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerFaceAlpha',0.6,'LineWidth',1)
hold on 
scatter(timeEvolutionS2(end-1,1,:),(timeEvolutionS2(end-1,2,:)+ timeEvolutionS2(end-1,3,:)),'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4,'LineWidth',1)
hold off 
subplot(1,3,1)
scatter(timeEvolutionS(Tend/10,1,:),(timeEvolutionS(Tend/10,2,:)+ timeEvolutionS(Tend/10,3,:)),'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerFaceAlpha',0.6,'LineWidth',1)
hold on 
scatter(timeEvolutionS2(Tend/10,1,:),(timeEvolutionS2(Tend/10,2,:)+ timeEvolutionS2(Tend/10,3,:)),'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4,'LineWidth',1)
hold off 
subplot(1,3,2)
scatter(timeEvolutionS(Tend/2,1,:),(timeEvolutionS(Tend/2,2,:)+ timeEvolutionS(Tend/2,3,:)),'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerFaceAlpha',0.6,'LineWidth',1)
hold on 
scatter(timeEvolutionS2(Tend/2,1,:),(timeEvolutionS2(Tend/2,2,:)+ timeEvolutionS2(Tend/2,3,:)),'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'MarkerFaceAlpha',0.4,'LineWidth',1)
hold off 

%% Repeat simulations with different half lives instead of changes in bursting

%% Initialize in the LES state and transition through fluctuations alone to other states

clear Tend npath MasterStateSeq
Tend = 1000;    %total simulation time
npath =10000;      %number of path trajectories to simulate
MasterStateSeq = zeros(Tend,npath);  %Storage for state seq for each path
mRNAcounts = zeros(Tend,3,npath);

% Parameters for the Low State
clear fn3
idsim=find(fnCMP.State == 'Low');
fn3=fnCMP(idsim,:);

ParamsL =   [0.0104	0.0163	0.0125;
    0.0831	0.1027	0.372;
    0.3153	0.0387	1.967;
    0.0039	0.0058	0.023];
%PU1    %Gata1  %Gata2

activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
% bL=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
% bL(isnan(bL)) = 0;
fL=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICL=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Parameters for the Gata2 High State
clear fn3
idsim=find(fnCMP.State == 'Gata2 High');
fn3=fnCMP(idsim,:);

ParamsG2 =   [0.0129	0.0027	0.129;
    0.097	0.1343	0.1415;
    0.3895	0.667	1.059;
    0.0039	0.0058	0.023];
%PU1    %Gata1  %Gata2
activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
bG2=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
bG2(isnan(bG2)) = 0;
fG2=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICG2=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Parameters for the Gata1 High State
clear fn3
idsim=find(fnCMP.State == 'Gata1 High');
fn3=fnCMP(idsim,:);

ParamsG1 =   [0.006	0.0172	0.3268;
    0.0704	0.0604	0.1995;
    0.321	1	1.7;
    0.0039	0.0058	0.0231];
%PU1    %Gata1  %Gata2

activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
bG1=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
bG1(isnan(bG1)) = 0;
fG1=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICG1=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Parameters for the PU1 High State
clear fn3
idsim=find(fnCMP.State == 'PU1 High');
fn3=fnCMP(idsim,:);

ParamsP1 =   [0.0468	0.000782	0.0135;
    0.068	0.0333	0.191;
    0.5	0.2933	1.0015;
    0.0039	0.0058	0.0231];
%PU1    %Gata1  %Gata2

activP=find(fn3.p1TS>=1);
activG1=find(fn3.g1TS>=1);
activG2=find(fn3.g2TS>=1);

%Compute the average burst size (for those >= 1TS/gene) and TS freq
bP1=[nanmean(fn3.pu1nas(activP)./fn3.p1TS(activP));nanmean(fn3.g1nas(activG1)./fn3.g1TS(activG1));mean(fn3.g2nas(activG2)./fn3.g2TS(activG2))];
bP1(isnan(bP1)) = 0;
fP1=[mean(fn3.p1TS); mean(fn3.g1TS); mean(fn3.g2TS)];
ICP1=[mean(fn3.pu1), mean(fn3.gata1), mean(fn3.gata2)];

% Storage for all parameters and average burst and frequencies
Allparams = zeros(4,3,4);

Allparams(:,:,1) = ParamsG1;
Allparams(:,:,2) = ParamsG2;
Allparams(:,:,3) = ParamsL;
Allparams(:,:,4) = ParamsP1;

Allb = vertcat(bG1',bG2',bL',bP1');
Allf = vertcat(fG1',fG2',fL',fP1');
IC = vertcat(ICG1,ICG2,ICL,ICP1);

% Initialize the storage matrices
for ctr =1:4 %multiplier for the half life
    clear timeEvolutionG timeEvolutionU timeEvolutionS timeEvolutionT
    
    timeEvolutionT=zeros(Tend,npath);   %Stores time for each path
    timeEvolutionS=zeros(Tend,3,npath);      %storage vector for the time
    timeEvolutionG=zeros(Tend,3,npath);     % Stores state of the promoters
    timeEvolutionU=zeros(Tend,3,npath);     % Stores state of the promoters
    
    % GSSA using exact reaction path
    gmax=2;
    rxnSeq = zeros(Tend,1,npath);
    
    for n=1:npath
        timeEvolutionS(1,:,n) = [0 0 0];   %All start in LES
        timeEvolutionG(1,:,n)=[0 0 0];
        
        for t=1:Tend-1
            %determine the current state of the cell
            
            if timeEvolutionS(t,2,n) >= 10 && timeEvolutionS(t,1,n) <= 20
                MasterStateSeq(t,n) = 1;
            elseif timeEvolutionS(t,1,n) >= 40 && timeEvolutionS(t,3,n)<50
                MasterStateSeq(t,n) = 4;
            elseif timeEvolutionS(t,3,n)>25
                MasterStateSeq(t,n) = 2;
            else
                MasterStateSeq(t,n) = 3;
            end
            StateIndex = MasterStateSeq(t,n);
            %set the parameters for this time step
            
            Params = Allparams(:,:,3);
            %         b = Allb(StateIndex,:);
            f = Allf(StateIndex,:);
            
            if StateIndex == 1
                gamma=[log(2)/180 log(2)/(ctr*120) log(2)/20];
            elseif StateIndex ==2
                gamma=[log(2)/180 log(2)/120 log(2)/(ctr*20)]; %mRNA half lives
            elseif StateIndex ==4
                gamma=[log(2)/(ctr*180) log(2)/120 log(2)/20];
            else
                gamma=[log(2)/180 log(2)/120 log(2)/20];
            end
            
            
            %Calculte Propensities
            kon=[Params(1,1)*(gmax-timeEvolutionG(t,1,n)),Params(1,2)*(gmax-timeEvolutionG(t,2,n)),Params(1,3)*(gmax-timeEvolutionG(t,3,n))];
            koff=[Params(2,1)*(timeEvolutionG(t,1,n)),Params(2,2)*(timeEvolutionG(t,2,n)),Params(2,3)*(timeEvolutionG(t,3,n))];
            kini=[Params(3,1),Params(3,2),Params(3,3)];
            kr=[kini(1)*timeEvolutionG(t,1,n),kini(2)*timeEvolutionG(t,2,n),kini(3)*timeEvolutionG(t,3,n)];
            kd=gamma.*timeEvolutionS(t,:,n);
            
            %Sum of the propensities
            asum=sum(kon)+sum(koff)+sum(kr)+sum(kd);
            
            
            timeEvolutionT(t+1,n)=timeEvolutionT(t,n)+log(1/rand(1))/asum;
            %Exact GSSA
            mu=rand(1);
            
            if 0 <= mu && mu < kon(1)/asum        %PU1 on
                rxnSeq(t,1,n) = 1;
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n)+1;
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= (timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= (timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= (timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif kon(1)/asum <= mu && mu  < (kon(1) + kon(2))/asum %Gata1 on
                rxnSeq(t,1,n) = 2;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n)+1;
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= (timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= (timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= (timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (kon(1) + kon(2))/asum <= mu && mu  < (sum(kon))/asum %Gata2 on
                rxnSeq(t,1,n) = 3;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n)+1;
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon))/asum <= mu && mu  < (sum(kon)+koff(1))/asum %PU1 off
                rxnSeq(t,1,n) = 4;
                
                timeEvolutionG(t+1,1,n)= max(timeEvolutionG(t,1,n)-1,0);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+koff(1))/asum <= mu && mu  < (sum(kon)+koff(1)+koff(2))/asum %Gata1 off
                rxnSeq(t,1,n) = 5;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= max(timeEvolutionG(t,2,n)-1,0);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+koff(1)+koff(2))/asum <= mu && mu  < (sum(kon)+sum(koff))/asum %Gata2 off
                rxnSeq(t,1,n) = 6;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= max(timeEvolutionG(t,3,n)-1,0);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
            elseif (sum(kon)+sum(koff))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1))/asum %PU1 burst
                rxnSeq(t,1,n) = 7;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n)+1);
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=1;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+sum(koff)+kr(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+kr(1)+kr(2))/asum %Gata1 burst
                rxnSeq(t,1,n) = 8;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n)+1);
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=1;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+sum(koff)+kr(1)+kr(2))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr))/asum %Gata2 burst
                rxnSeq(t,1,n) = 9;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n)+1);
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=1;
                
                
            elseif (sum(kon)+sum(koff)+sum(kr))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum %PU1 decay
                rxnSeq(t,1,n) = 10;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= max(timeEvolutionS(t,1,n)-1,0);
                timeEvolutionS(t+1,2,n)= ceil(timeEvolutionS(t,2,n));
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            elseif (sum(kon)+sum(koff)+sum(kr)+kd(1))/asum <= mu && mu  < (sum(kon)+sum(koff)+sum(kr)+kd(1)+kd(2))/asum %Gata1 decay
                rxnSeq(t,1,n) = 11;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= ceil(timeEvolutionS(t,1,n));
                timeEvolutionS(t+1,2,n)= max(timeEvolutionS(t,2,n)-1,0);
                timeEvolutionS(t+1,3,n)= ceil(timeEvolutionS(t,3,n));
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            else  %Gata2 decay
                rxnSeq(t,1,n) = 12;
                
                timeEvolutionG(t+1,1,n)= timeEvolutionG(t,1,n);
                timeEvolutionG(t+1,2,n)= timeEvolutionG(t,2,n);
                timeEvolutionG(t+1,3,n)= timeEvolutionG(t,3,n);
                
                timeEvolutionS(t+1,1,n)= timeEvolutionS(t,1,n);
                timeEvolutionS(t+1,2,n)= timeEvolutionS(t,2,n);
                timeEvolutionS(t+1,3,n)= max(timeEvolutionS(t,3,n)-1,0);
                
                timeEvolutionU(t+1,1,n)=0;
                timeEvolutionU(t+1,2,n)=0;
                timeEvolutionU(t+1,3,n)=0;
                
            end
            
            
        end
        
        path_num=n
        
    end
    
    %Plot a histogram of the resulting endpoint distribution
    figure(1)
    subplot(4,1,ctr)
    histogram(MasterStateSeq(end-1,:),'Normalization','probability','FaceColor',[0.35 0.35 0.35],'FaceAlpha',0.8,'LineWidth',2)
    
    %Plot a trajectory plot with associated marginal distribution of states
    %as a function of the PU1-(Gata1+Gata2)
    figure (2)
    subplot(4,1,ctr)
    Ridx = randi(npath,1000,1);
    
    clrs = [1,0,0; 0.8 0.4 0.6; 0.5 0.5 0.5; 0 0 1];
    for ctrx = 1:1000
        kk = Ridx(ctrx);
        deltaM = timeEvolutionS(:,1,kk)-(timeEvolutionS(:,2,kk)+ timeEvolutionS(:,3,kk));
        State = MasterStateSeq(end-1,kk);
        if State ~= 3
            stairs(deltaM,'Color',clrs(State,:),'LineWidth',2)
            hold on
        else
            stairs(deltaM,'Color',clrs(State,:),'LineWidth',0.15)
            hold on
        end       %for a plot versus cumulative time, use
        %stairs(timeEvolutionT(:,kk),deltaM,'Color',clrs(State,:),'LineWidth',1)
        hold on
    end
    
    hold off
    
end
