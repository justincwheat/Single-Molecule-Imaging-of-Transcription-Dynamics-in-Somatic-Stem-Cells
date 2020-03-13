%Figure 4e: testing 3-cell state correlations


clear all
%Positions with State assignments/counts
load(['bdry.mat']);
load(['colonyIDs.mat']);
load(['20180114_analysis.mat']);
load(['transMatMean.mat']); %stored from the KCA run
load(['transMatSTD.mat']);  %stored from the KCA run
load(['SteadyStateDist.mat']);

%% Gata2 fitting and Q matrix calculation
clear fd sd pdf_nbinommixture x mixingPrct

x=KCA.gata2(KCA.gata1<100);

kp1=[2 0.0721 11.239 0.10975];  %Derived from Distribution Fitter
pdf_nbinommixture = @(x,p) ...
    p*nbinpdf(x,kp1(1),kp1(2)) + (1-p)*nbinpdf(x,kp1(3),kp1(4));

firstcomp = @(x,p) ...
    p*nbinpdf(x,kp1(1),kp1(2));
secondcomp= @(x,p,c,d) ...
    (1-p)*nbinpdf(x,kp1(3),kp1(4));

p0=[0.6];
options = statset('MaxIter',10000, 'MaxFunEvals',100000);
lb = [0];
ub = [1];
[mixingPrct, pci2] = mle(x, 'pdf',pdf_nbinommixture, 'start',p0, ...
    'lower',lb, 'upper',ub, 'options',options)
G2params=[mixingPrct kp1];

figure(1)

histogram(x,'Normalization','pdf')
hold on
xgrid = linspace(0,900,900);
pdfgrid = pdf_nbinommixture(1:900,mixingPrct(1));
g2fd= firstcomp(1:900,mixingPrct(1));
g2sd= secondcomp(1:900,mixingPrct(1));
plot(xgrid,pdfgrid,'-')
hold on
plot(xgrid,g2fd,'.')
hold on
plot(xgrid,g2sd,'.')
hold off
xlabel('x')
ylabel('Probability Density')

G2Qmat(1,1) = sum( (g2fd ./ (g2fd + g2sd)) .* (g2fd ./ mixingPrct) );
G2Qmat(2,1) = sum( (g2sd ./ (g2fd + g2sd)) .* (g2fd ./ mixingPrct) );
G2Qmat(1,2) = sum( (g2fd ./ (g2fd + g2sd)) .* (g2sd ./ (1-mixingPrct)) );
G2Qmat(2,2) = sum( (g2sd ./ (g2fd + g2sd)) .* (g2sd ./ (1-mixingPrct)) );

%% Determine three cell correlations

%Number of iterations for bootstrap
numItr =1000;


clear allT11 allT22 allsizemCorr2 mCorr3 AllmCorr3

ctrStore = 0;
N = 6;

sizeBdry=bdry;
AllinputPos=unique(colonyIDs);
%Probability of selection in Bootstrap is proportional to the bdry size
selProb = sizeBdry ./ sum(sizeBdry);
cumProb = cumsum(selProb);

for ctrItr = 1:numItr
    
    PercentofAnalysisComplete = 100*(ctrItr/numItr)
    %Initialize
    ctrxx = 0;
    
    %Ctr of cells thrown out
    ctrThrownOut = 0;
    
    %Initialize
    ngen_max = 10;
    for u = 1
        for v = 1:ngen_max
            mCorr3{u,v} = zeros(N,N,N);
        end
    end
    
    %Initialize the one-point distribution
    ssDist = zeros(1,N);
    
    %Bootstrap
    tmpPos = randi(length(AllinputPos),1,length(AllinputPos));
    for k = 1:length(AllinputPos)
        tmpPos(k) = find( rand < cumProb , 1 , 'first');
    end
    inputPos = sort(AllinputPos((tmpPos)));
    
    clear Alldistu numBdryCells AllposAct AllFISHcounts LinDist cellID posID
    
    for fctr = 1:length(inputPos)
        
        clear genVec idVec genSig totNumCells cellSig cellSize cellArea cellFrame cellParent cellAncestorFrame cellPos bdryCells allParents
        
        %Load analyzed tracked
        load(['Colony_' num2str(inputPos(fctr)) '.mat']);
        
        bdryCells=colony.CellID;
        
        
        %Compute the correlation functions
        
        
        numBdryCells(fctr) = length(bdryCells);
        
        G1thr=10;
        Megthr=100;
        MThr=150;
        P1Hthr=75;
        %States= [Meg G1/2H G2H LES P1H Mphage];
        for k=1:length(bdryCells)-2
            for i = k+1:length(bdryCells)-1
                for j = i+1:length(bdryCells)
                    %State of the kCell
                    
                    g2=colony.gata2(k)+1;
                    
                    kG2prob=[g2fd(g2) g2sd(g2)]./(g2fd(g2)+g2sd(g2));
                    
                    if colony.gata1(k) >= Megthr
                        kState = [1 0 0 0 0 0];
                    elseif colony.pu1(k) >= MThr
                        kState = [0 0 0 0 0 1];
                    elseif (colony.pu1(k) > P1Hthr)
                        kState = [0 0 0 0 1 0]; %P1H
                    elseif( colony.gata1(k) >= G1thr)
                        kState = [0 1 0 0 0 0]; %G1H
                    else
                        kState2=[0 0 kG2prob(2) kG2prob(1) 0 0];
                        kState=kState2./(sum(kState2));
                    end
                    
                    %Get the state of iCell
                    
                    g2=colony.gata2(i)+1;
                    
                    iG2prob=[g2fd(g2) g2sd(g2)]./(g2fd(g2)+g2sd(g2));
                    
                    if colony.gata1(i) >= Megthr
                        iState = [1 0 0 0 0 0];
                    elseif colony.pu1(i) >= MThr
                        iState = [0 0 0 0 0 1];
                    elseif (colony.pu1(i) > P1Hthr)
                        iState = [0 0 0 0 1 0]; %P1H
                    elseif( colony.gata1(i) >= G1thr)
                        iState = [0 1 0 0 0 0]; %G1H
                    else
                        iState2=[0 0 iG2prob(2) iG2prob(1) 0 0];
                        iState=iState2./(sum(iState2));
                    end
                    
                    %get the State of the jCell
                    %get probabilistic assignment for G2 and P1
                    g2=colony.gata2(j)+1;
                    
                    jG2prob=[g2fd(g2) g2sd(g2)]./(g2fd(g2)+g2sd(g2));
                    
                    if colony.gata1(j) >= Megthr
                        jState = [1 0 0 0 0 0];
                    elseif colony.pu1(j) >= MThr
                        jState = [0 0 0 0 0 1];
                    elseif (colony.pu1(j) > P1Hthr)
                        jState = [0 0 0 0 1 0]; %P1H
                    elseif( colony.gata1(j) >= G1thr)
                        jState = [0 1 0 0 0 0]; %G1H
                    else
                        jState2=[0 0 jG2prob(2) jG2prob(1) 0 0];
                        jState=jState2./(sum(jState2));
                    end
                    
                    
                    %Distance between ij ik and jk
                    distu=[];
                    
                    distu(1)=colony.deltag(i,j);
                    distu(2)=colony.deltag(i,k);
                    distu(3)=colony.deltag(k,j);
                    
                    
                    %u number of generations to the CA of two more
                    %closely related nodes
                    %v number of generations to CA ofall three nodes
                    thisu = min(distu);
                    thisv = max(distu);
                    
                    
                    %Reorder the states such that k state is always
                    %the one furthest away
                    if (distu(1) > distu(2) )
                        if ( distu(2) > distu(3) ) %I is furthest away
                            tmpState = iState;
                            iState = kState;
                            kState = tmpState;
                        elseif ( distu(3) > distu(2) ) %J is furthest away
                            tmpState = jState;
                            jState = kState;
                            kState = tmpState;
                        end
                    end
                    
                    
                    %Get the triple multiple for the three
                    if thisu == 1 && thisv>0
                        %point correlatoin functions
                        thismCorr3 = [];
                        for x = 1:N
                            for y = 1:N
                                for z=1:N
                                    thismCorr3(x,y,z) = iState(x)*jState(y)*kState(z);
                                end
                            end
                        end
                        
                        mCorr3{thisu,thisv} = mCorr3{thisu,thisv} + thismCorr3;
                    end
                    
                end
            end
        end
    end
    
    
    %Normalize the mCorr3s
    for u = 1
        for v= 1:ngen_max
            for k = 1:N
                tmpMat3 = mCorr3{u,v}(:,:,k);
                mCorr3{u,v}(:,:,k) = tmpMat3 + tmpMat3';
            end
            
            
            mCorr3Norm{u,v} = mCorr3{u,v} ./ sum(sum(sum(mCorr3{u,v})));
            
        end
    end
    %Store all the correlation matrices
    AllmCorr3{ctrItr} = mCorr3Norm;
    
end
%% Save all the AllmCorr3s

save('C:\Users\Justin Wheat\Desktop\KCA_resubmission\AllmCorrs3.mat','AllmCorr3')

%% Calculate means and std of experimental 3-states at u=1 and v=2:5
clear mCorr3Mean mCorr3Std
%Get the mean and std of each element
for u = 1       %restrict u to sister cells
    for v = 2:5     %restrict to v = 5
        for i = 2:5 %only want states 2:5
            for j = 2:5 %only want states 2:5
                for k = 2:5 %only want states 2:5
                    
                    alltmpVals = [];
                    for ctrItr = 1:numItr
                        
                        alltmpVals = [alltmpVals AllmCorr3{ctrItr}{u,v}(i,j,k)];
                        
                    end
                    %this will generate a 4x4x4 matrix for 4 values of v
                    mCorr3Mean{u,v-1}(i-1,j-1,k-1) = mean(alltmpVals);
                    mCorr3Std{u,v-1}(i-1,j-1,k-1) = std(alltmpVals);    %the v-1 comes from the fact that we are looking at cousins (v=2) but need to store into the first index
                    
                    
                end
            end
        end
    end
end
% %% Visualize
%
% figure(2)
% vd = 2;
%     for k = 1:4
%     subplot(3,4,k)
%     %figure
%
%     imagesc(mCorr3Mean{1,vd}(:,:,k))
%     if (k==1), title('Actual'); end
%
%     %Fractional error
%     %imagesc(mCorr3Std{1,2}(:,:,k)./AllmCorr3{1}{1,2}(:,:,k))
%
%     colorbar
%
%     set(gca,'FontSize',20)
%     set(gca,'LineWidth',3)
%
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
%
%     colormap(parula)
%     %colormap(bone)
%
%     %caxis(my_Caxis(k,:))
%
%
%     set(gcf,'color','w')
%
%
%
%     end

%     vd = 3;
%     for k = 1:4
%     subplot(3,4,k+4)
%     %figure
%
%     imagesc(mCorr3Mean{1,vd}(:,:,k))
%     if (k==1), title('Actual'); end
%
%     %Fractional error
%     %imagesc(mCorr3Std{1,2}(:,:,k)./AllmCorr3{1}{1,2}(:,:,k))
%
%     colorbar
%
%     set(gca,'FontSize',20)
%     set(gca,'LineWidth',3)
%
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
%
%     colormap(parula)
%     %colormap(bone)
%
%     %caxis(my_Caxis(k,:))
%
%
%     set(gcf,'color','w')
%
%
%
%     end
%
%     vd = 4;
%     for k = 1:4
%     subplot(3,4,k+8)
%     %figure
%
%     imagesc(mCorr3Mean{1,vd}(:,:,k))
%     if (k==1), title('Actual'); end
%
%     %Fractional error
%     %imagesc(mCorr3Std{1,2}(:,:,k)./AllmCorr3{1}{1,2}(:,:,k))
%
%     colorbar
%
%     set(gca,'FontSize',20)
%     set(gca,'LineWidth',3)
%
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
%
%     colormap(parula)
%     %colormap(bone)
%
%     %caxis(my_Caxis(k,:))
%
%
%     set(gcf,'color','w')
%
%
%
%     end


%% Generate 3 cell correlations from each dynamical model

clear th1Corr3 th2Corr3 th3Corr3 th4Corr3
tMat = combinedTransMat;    %derived from KCA, fully connected model
tMat(5,:)=[];
tMat(:,5)=[];

p0 = mssDist(1:4);   %starting population frequencies
N = 4; %4 states in order G12H, G2H, LES, P1H
ngen = 4; %go to v = 5
eyeMat = eye(N); %Identity matrix

%Calculate the theoretical three-point correlation functions for all four
%possible chain based models

%Th4 fully Connected Transition Matrix
for u = 1 : ngen-1
    for v = 1 : ngen-u  %v =1 is actually a cousin cell, storage purposes only
        thisCorr3 = zeros(N,N,N); %(Na,Nb,Nc) NaNb are the states of the closer two points Nc is the state of the furthest relative
        
%         Compute the initial distribution
        pI = tMat^(ngen-u-v)*p0';
        for i = 1:N
            
            %Get the probability of Nc, i.e. the last common ancestor
            pNc = tMat^(u+v) * eyeMat(:,i);
            
            %Gen the probability of Nd (parent of Na and Nb)
            pNd = tMat^u * eyeMat(:,i);
            
            
            %Initialize
            sumNaNb = zeros(N);
            
            for j = 1:N
                sumNaNb = sumNaNb + pNd(j) * ( (tMat^u * eyeMat(:,j)) * (tMat^u * eyeMat(:,j))' );
            end
            
            for j = 1:N
                thisCorr3(:,:,j) = thisCorr3(:,:,j) + sumNaNb * pNc(j) * pI(i);
            end
            
        end
        
        th4Corr3{u,v} = thisCorr3;
        
    end
    
    
end

%now for the Markov Chain

tMat(3:end,1) = 0;
tMat(1:2,4) = 0;

for k =1:N
    tMat(:,k) = tMat(:,k) / sum(tMat(:,k));
end

for u = 1 : ngen-1
    for v = 1 : ngen-u  %v =1 is actually a cousin cell, storage purposes only
        thisCorr3 = zeros(N,N,N); %(Na,Nb,Nc) NaNb are the states of the closer two points Nc is the state of the furthest relative
        
%         Compute the initial distribution
        pI = tMat^(ngen-u-v)*p0';
        for i = 1:N
            
            %Get the probability of Nc, i.e. the last common ancestor
            pNc = tMat^(u+v) * eyeMat(:,i);
            
            %Gen the probability of Nd (parent of Na and Nb)
            pNd = tMat^u * eyeMat(:,i);
            
            
            %Initialize
            sumNaNb = zeros(N);
            
            for j = 1:N
                sumNaNb = sumNaNb + pNd(j) * ( (tMat^u * eyeMat(:,j)) * (tMat^u * eyeMat(:,j))' );
            end
            
            for j = 1:N
                thisCorr3(:,:,j) = thisCorr3(:,:,j) + sumNaNb * pNc(j) * pI(i);
            end
            
        end
        
        th2Corr3{u,v} = thisCorr3;
        
    end
    
    
end

%now for the Partially Connected Model

tMat(:,1) = [1;0;0;;0];
tMat(:,4) = [0;0;0;1];


for k =1:N
    tMat(:,k) = tMat(:,k) / sum(tMat(:,k));
end

for u = 1 : ngen-1
    for v = 1 : ngen-u  %v =1 is actually a cousin cell, storage purposes only
        thisCorr3 = zeros(N,N,N); %(Na,Nb,Nc) NaNb are the states of the closer two points Nc is the state of the furthest relative
        
%         Compute the initial distribution
        pI = tMat^(ngen-u-v)*p0';
        for i = 1:N
            
            %Get the probability of Nc, i.e. the last common ancestor
            pNc = tMat^(u+v) * eyeMat(:,i);
            
            %Gen the probability of Nd (parent of Na and Nb)
            pNd = tMat^u * eyeMat(:,i);
            
            
            %Initialize
            sumNaNb = zeros(N);
            
            for j = 1:N
                sumNaNb = sumNaNb + pNd(j) * ( (tMat^u * eyeMat(:,j)) * (tMat^u * eyeMat(:,j))' );
            end
            
            for j = 1:N
                thisCorr3(:,:,j) = thisCorr3(:,:,j) + sumNaNb * pNc(j) * pI(i);
            end
            
        end
        
        th3Corr3{u,v} = thisCorr3;
        
    end
    
    
end


%Finally the fully irreversible model

tMat(3:4,2) = [0;0];


for k =1:N
    tMat(:,k) = tMat(:,k) / sum(tMat(:,k));
end

for u = 1 : ngen-1
    for v = 1 : ngen-u  %v =1 is actually a cousin cell, storage purposes only
        thisCorr3 = zeros(N,N,N); %(Na,Nb,Nc) NaNb are the states of the closer two points Nc is the state of the furthest relative
        
%         Compute the initial distribution
        pI = tMat^(ngen-u-v)*p0';
        for i = 1:N
            
            %Get the probability of Nc
            pNc = tMat^(u+v) * eyeMat(:,i);
            
            %Gen the probability of Nd (parent of Na and Nb)
            pNd = tMat^u * eyeMat(:,i);
            
            
            %Initialize
            sumNaNb = zeros(N);
            
            for j = 1:N
                sumNaNb = sumNaNb + pNd(j) * ( (tMat^u * eyeMat(:,j)) * (tMat^u * eyeMat(:,j))' );
            end
            
            for j = 1:N
                thisCorr3(:,:,j) = thisCorr3(:,:,j) + sumNaNb * pNc(j) * pI(i);
            end
            
        end
        
        th1Corr3{u,v} = thisCorr3;
        
    end
    
    
end


%% Graph the results comparing transitions

% Plot the three point freqs for v 2-4
clc
clear AllCorr3Th AllCorr3Obs AllCorr3ObsErr
colors=['r','y','k','b'];
shapes=['o','hexagram','diamond','*'];
figure(3)
%note that N below is max 4! we have converted to a four HSPC state system now
%to exclude megakaryoctes. This was done above when generating mCorr3Mean and mCorr3Std by going from
%2:6 rather 1:N and indexing in 1:5.
for t=1:3       %V value. All ij are u=1
    ctrx = 0;
    for k = 1:4 %states of kth cell
        for i =1:N%sister cell i
            for j = i:N %sister cell j
                ctrx = ctrx + 1;
                All1Corr3Th(ctrx) = th1Corr3{1,t}(i,j,k);
                All2Corr3Th(ctrx) = th2Corr3{1,t}(i,j,k);
                All3Corr3Th(ctrx) = th3Corr3{1,t}(i,j,k);
                All4Corr3Th(ctrx) = th4Corr3{1,t}(i,j,k);
                
                
                AllCorr3Obs(ctrx) = mCorr3Mean{1,t}(i,j,k);
                AllCorr3ObsErr(ctrx) = mCorr3Std{1,t}(i,j,k);
            end
        end
        
        tmpY = log10(AllCorr3Obs);
        tmpX1 = log10(All1Corr3Th);
        tmpX2 = log10(All2Corr3Th);
        tmpX3 = log10(All3Corr3Th);
        tmpX4 = log10(All4Corr3Th);
       
        
        tmpYErr = AllCorr3ObsErr./AllCorr3Obs;
        
        errorbar(tmpX3,tmpY,tmpYErr,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0.5 0.5 0.5],'Color','k','MarkerSize',10)
        set(gca,'FontSize',36)
        set(gca,'LineWidth',3)
        hold on
        
        errorbar(tmpX1,tmpY,tmpYErr,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','Color','k','MarkerSize',10)
        set(gca,'FontSize',36)
        set(gca,'LineWidth',3)
        hold on
        errorbar(tmpX2,tmpY,tmpYErr,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0.25 0.25 0.25],'Color','k','MarkerSize',10)
        set(gca,'FontSize',36)
        set(gca,'LineWidth',3)
        hold on
        
        errorbar(tmpX4,tmpY,tmpYErr,'o','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','Color','k','MarkerSize',10)
        set(gca,'FontSize',36)
        hold on
        
        
        clear All1Corr3Th All2Corr3Th All3Corr3Th All4Corr3Th AllCorr3Obs AllCorr3ObsErr
        
    end
    hold on
    
end
hold on


plot([-5 0],[-5 0],'k:','LineWidth',2,'MarkerSize',30)
xlim([-5.0 0])
ylim([-5.0 0])
%
xlabel('Log(predicted value)', 'Fontsize', 30)
ylabel('Log(actual value)', 'Fontsize', 30)
%
hold off


%% Calculate total error for each model over all lineage distances up to n=5
clear MAE

for t=1:3   %this calls into the
    a=sum(sum(sum(abs(mCorr3Mean{1,t}-th1Corr3{1,t})))); %MAE for each k state at v
    MAE(1,t)=a;
    
    a=sum(sum(sum(abs(mCorr3Mean{1,t}-th2Corr3{1,t})))); %MAE for each k state at v
    MAE(2,t)=a;
    
    a=sum(sum(sum(abs(mCorr3Mean{1,t}-th3Corr3{1,t})))); %MAE for each k state at v
    MAE(3,t)=a;
    
    a=sum(sum(sum(abs(mCorr3Mean{1,t}-th4Corr3{1,t})))); %MAE for each k state at v
    MAE(4,t)=a;
    

end
figure(4)
imagesc(MAE)
