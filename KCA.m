%% KCA
clear all
load(['colonyIDs.mat']); 
load('KCA_datamatrix_withStateAssignments.mat');
%%

% %Gata2 fitting and Q matrix calculation
clear fd sd pdf_nbinommixture x mixingPrct

x=KCA.gata2(KCA.gata1<250);

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
%
%% Now assemble the combined Qmat for all transitions

N=4;
Qmat=eye(N);
Qmat(2,2)=G2Qmat(2,2);
Qmat(2,3)=G2Qmat(2,1);
Qmat(3,2)=G2Qmat(1,2);

Qmat(3,3)=1-(sum(Qmat(2,3)+Qmat(4,3)));

for jj = 1:N
    Qmat(:,jj) = Qmat(:,jj) ./ sum(Qmat(:,jj));
end

AllinputPos=colonyIDs;
for c = 1:length(AllinputPos)
            load(['Colony_' num2str(AllinputPos(c)) '.mat']);

    sizeBdry(c) = length(colony.CellID);
end

%Number of iterations for bootstrap
numItr =5000;
clear allT11 allT22 allsizemCorr2 

ctrStore = 0;


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
    for u = 1:ngen_max
        mCorr2{u} = zeros(N);
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
        P1Hthr=75;
        %States= [G1/2H G2H LES P1H];
        for i = 1:length(bdryCells)-1
            if( colony.gata1(i) <= Megthr )
                for j = i+1:length(bdryCells)
                    if( colony.gata1(j) <= Megthr )
                        %Get the state of iCell
                        
                        
                        g2=colony.gata2(i)+1;
                        iG2prob=[g2fd(g2) g2sd(g2)]./(g2fd(g2)+g2sd(g2));
                        
                        if (colony.pu1(i) > P1Hthr(1))
                            iState = [0 0 0 1]; %P1H
                        elseif( colony.gata1(i) >= G1thr)
                            iState = [1 0 0 0]; %G1H
                        else
                            iState2=[ 0 iG2prob(2) iG2prob(1) 0];
                            iState=iState2./(sum(iState2));
                        end
                        
                        %get the State of the jCell
                        
                        g2=colony.gata2(j)+1;
                        
                        iG2prob=[g2fd(g2) g2sd(g2)]./(g2fd(g2)+g2sd(g2));
                        
                        if (colony.pu1(j) > P1Hthr(1))
                            jState = [0 0 0 1]; %P1H
                        elseif( colony.gata1(j) >= G1thr)
                            jState = [1 0 0 0]; %G1H
                        else
                            iState2=[ 0 iG2prob(2) iG2prob(1) 0];
                            jState=iState2./(sum(iState2));
                        end
                        
                        
                        
                        distu=colony.deltag(i,j);
                        
                        if ( distu > 0 )
                            mCorr2{distu} = mCorr2{distu} + iState' * jState;
                        end
                    end
                end
            end
        end
        %Get the FISH counts at the boundary
        for i = 1:length(bdryCells)
            if( colony.gata1(i) <= Megthr )
                
                iCell = bdryCells(i);
                cellThrownOut = 0;
                p1=colony.pu1(i)+1;
                g2=colony.gata2(i)+1;
                
                iP1prob=[fd(p1) sd(p1)]./(fd(p1)+sd(p1));
                iG2prob=[g2fd(g2) g2sd(g2)]./(g2fd(g2)+g2sd(g2));
                
                if (colony.pu1(i) > P1Hthr(1))
                    iState = [0 0 0 1]; %P1H
                elseif( colony.gata1(i) >= G1thr)
                    iState = [1 0 0 0]; %G1H
                else
                    iState2=[ 0 iG2prob(2) iG2prob(1) 0];
                    iState=iState2./(sum(iState2));
                end
                ssDist = ssDist + iState;
                
                ctrxx = ctrxx + 1;
                AllFISHcounts(ctrxx) = colony.pu1(i);
                AllFISHcounts(ctrxx) = colony.gata1(i);
                AllFISHcounts(ctrxx) = colony.gata2(i);
                AllFISHstate(ctrxx,:) = iState;
            end
        end
        
    end
    %Normalize
    ssDist = ssDist ./ sum(ssDist);
    emp_ssDist = ssDist;
    
    %Normalize the mCorr2s
    for u = 1:ngen_max
        mCorr2{u} = mCorr2{u} + mCorr2{u}';
        mCorr2Norm{u} = mCorr2{u} ./ sum(sum(mCorr2{u}));
    end
    
    
    clear AllEigs AllStayRates sizemCorr2 AllmTrans AllmTransSignle avgmTrans
    %Infer the rate from u = 1
    intvPlot = 1;
    checkStore = 1;
    for i =1:intvPlot:6
        
        %Infer the rates
        ssDist = (inv(Qmat) * sum(mCorr2Norm{i})')';
        mTmp =(inv(Qmat)*mCorr2Norm{i}*inv(Qmat')) ./ (sqrt(ssDist)' * sqrt(ssDist));
        
        if ( sum(sum(isnan(mTmp))) == 0 )
            
            [V D] = eig(mTmp);
            mTrans = (V *D.^(1/((2*i))) * V') .* (sqrt(ssDist)'* (1./sqrt(ssDist)));
            
            %Normalize
            for jj = 1:N
                mTrans(:,jj) = mTrans(:,jj) ./ sum(mTrans(:,jj));
            end
            
            %Store the eigenvalues of T
            [V D] = eig(mTrans);
            AllEigs(:,i) = sort(diag(D));
            
            
            
            
            %Store the results
            AllmTrans{i} = mTrans;
        else
            checkStore = 0; %do not store
        end
        
        
        
    end
    
    
    %store the results
    if ( checkStore )
        ctrStore = ctrStore + 1;
        allpredmTrans{ctrStore} = AllmTrans;
        allssDist(ctrStore,:) = emp_ssDist;
    end
    
    %Store all the correlation matrices
    AllmCorr{ctrItr} = mCorr2Norm;
    
    
    
end
%% Get out of source transitions probabilities
for ii=1:numItr
    for jj=1:6
        AllG2trans(ii,1,jj) = allpredmTrans{ii}{jj}(1,2);
        AllG2trans(ii,2,jj) = allpredmTrans{ii}{jj}(3,2);
        AllG2trans(ii,3,jj) = allpredmTrans{ii}{jj}(4,2);
    end
end


AllG2trans(imag(AllG2trans)>0.001 | AllG2trans > 1 ) =NaN;
AllG2trans(AllG2trans<0)=0;

AllG2trans=real(AllG2trans);

for ii=1:numItr
    for jj=1:6
        AllLEStrans(ii,1,jj) = allpredmTrans{ii}{jj}(1,3);
        AllLEStrans(ii,2,jj) = allpredmTrans{ii}{jj}(2,3);
        AllLEStrans(ii,3,jj) = allpredmTrans{ii}{jj}(4,3);
    end
end


AllLEStrans(imag(AllLEStrans)>0.001 | AllLEStrans > 1 ) =NaN;
AllLEStrans(AllLEStrans<0)=0;

AllLEStrans=real(AllLEStrans);

for ii=1:numItr
    for jj=1:6
        AllG1trans(ii,1,jj) = allpredmTrans{ii}{jj}(2,1);
        AllG1trans(ii,2,jj) = allpredmTrans{ii}{jj}(3,1);
        AllG1trans(ii,3,jj) = allpredmTrans{ii}{jj}(4,1);
    end
end


AllG1trans(imag(AllG1trans)>0.001 | AllG1trans > 1 ) =NaN;
AllG1trans(AllG1trans<0)=0;

AllG1trans=real(AllG1trans);

for ii=1:numItr
    for jj=1:6
        AllP1trans(ii,1,jj) = allpredmTrans{ii}{jj}(1,4);
        AllP1trans(ii,2,jj) = allpredmTrans{ii}{jj}(2,4);
        AllP1trans(ii,3,jj) = allpredmTrans{ii}{jj}(3,4);
    end
end


AllP1trans(imag(AllP1trans)>0.001 | AllP1trans > 1 ) =NaN;
AllP1trans(AllP1trans<0)=0;

AllP1trans=real(AllP1trans);
figure
subplot(1,4,2)
hold on
colors=['r','k','b']
inc = 1;
incLim = 6;
for k=1:3
    errorbar((1:inc:incLim)',reshape(nanmean(AllG2trans(:,k,1:inc:incLim)),1,6),reshape(nanstd(AllG2trans(:,k,1:inc:incLim)),1,6),'MarkerFaceColor',colors(k),'MarkerEdgeColor','k', 'Marker','o','LineWidth',2,'MarkerSize',10,'Color',colors(k))
    set(gca,'FontSize',24)
    set(gca,'LineWidth',3)
    hold on
    xlim([0 7])
end
hold off

subplot(1,4,3)
hold on
colors=['r','y','b']
inc = 1;
incLim = 6;
for k=1:3
    errorbar((1:inc:incLim)',reshape(nanmean(AllLEStrans(:,k,1:inc:incLim)),1,6),reshape(nanstd(AllLEStrans(:,k,1:inc:incLim)),1,6),'MarkerFaceColor',colors(k),'MarkerEdgeColor','k','Marker','o','LineWidth',2,'MarkerSize',10,'Color',colors(k))
    set(gca,'FontSize',24)
    set(gca,'LineWidth',3)
    hold on
    xlim([0 7])
end
hold off

subplot(1,4,1)
hold on
colors=['y','k','b']
inc = 1;
incLim = 6;
for k=1:3
    errorbar((1:inc:incLim)',reshape(nanmean(AllG1trans(:,k,1:inc:incLim)),1,6),reshape(nanstd(AllG1trans(:,k,1:inc:incLim)),1,6),'MarkerFaceColor',colors(k),'MarkerEdgeColor','k', 'Marker','o','LineWidth',2,'MarkerSize',10,'Color',colors(k))
    set(gca,'FontSize',24)
    set(gca,'LineWidth',3)
    hold on
    xlim([0 7])
end
hold off

subplot(1,4,4)
hold on
colors=['r','y','k']
inc = 1;
incLim = 6;
for k=1:3
    errorbar((1:inc:incLim)',reshape(nanmean(AllP1trans(:,k,1:inc:incLim)),1,6),reshape(nanstd(AllP1trans(:,k,1:inc:incLim)),1,6),'MarkerFaceColor',colors(k),'MarkerEdgeColor','k','Marker','o','LineWidth',2,'MarkerSize',10,'Color',colors(k))
    set(gca,'FontSize',24)
    set(gca,'LineWidth',3)
    hold on
    xlim([0 7])
end
hold off

%%
%Plot the inferred rates for retaining the state of the parent cell
%(diagonal entries) of the transition matrix.

%Compute trans Mat
clear predTransMat stdTransMat
for u = 1:6
    
    for iC = 1:N
        for jC = 1:N
            %Get the element from store matrices
            for kk = 1:length(allpredmTrans)
                thisElem(kk) = allpredmTrans{kk}{u}(iC,jC);
                if ( imag(thisElem(kk)) > 0.001 || thisElem(kk) > 1)
                    thisElem(kk) = nan;
                end
            end
            predTransMat{u}(iC,jC) = nanmean(thisElem);
            stdTransMat{u}(iC,jC) = nanstd(thisElem);
        end
    end
    
    
end

figure
hold on

inc = 1;
incLim = 6;

clear diagTMat diagTstd

for u = 1:incLim %distance
    for k = 1:N %state
        diagTMat{k}(u) = real(predTransMat{u}(k,k));
        diagTstd{k}(u) = real(stdTransMat{u}(k,k));
    end
end

colors = {'r' 'g' [0.5 0.5 0.5] 'b','y'};

ctrSubplot = 0;
for kk = 1:N
    
    ctrSubplot = ctrSubplot + 1;
    
    errorbar([1:inc:incLim],diagTMat{kk},diagTstd{kk},'o', 'MarkerFaceColor' , colors{ctrSubplot} , 'LineWidth',2,'MarkerSize',12,'MarkerEdgeColor','k','Color','k')
    hold on
    
    set(gca,'FontSize',36)
    set(gca,'LineWidth',3)
    set(gca,'XTickLabel',[]);
    
    
    xlim([0 incLim+1])
    set(gca,'XTick',[1:incLim])
    
    line([0,6],[diagTMat{kk}(1),diagTMat{kk}(1)])
    hold on
    
    
end

%Combine the first two distances u = 1 nad u=2  to get the rates with
%lowest std
for iC = 1:N
    for jC= 1:N
        [tmpVal tmpPos] = min([stdTransMat{1}(iC,jC) stdTransMat{2}(iC,jC)]);
        combinedTransMat(iC,jC) = predTransMat{tmpPos}(iC,jC);
        combinedStdMat(iC,jC) = stdTransMat{tmpPos}(iC,jC);
    end
end


% Show pairs at different distances
%Plot the two-cell correlation matrix of sisters
N=4;
for u=1:6
    plotU =u; %Which u to plot
    %Change this to plot the two-cell correlation matrix of first cousins
    %plotU=2 and that of second cousins, plotU=3  (as in Fig S5)
    
    sumMCorr = zeros(N);
    
    for ctrItr = 1:length(AllmCorr)
        
        sumMCorr = sumMCorr + AllmCorr{ctrItr}{plotU};
        
        
    end
    
    avgMCorr = sumMCorr ./ length(AllmCorr);
    s=sum(avgMCorr,1);
    
    inMat=bsxfun(@rdivide,avgMCorr,s)
    %Plot the correlation functions with the correct ordering
    
    
    myMat = (inMat)
    figure (1)
    subplot(2,3,u)
    bar(myMat','stacked')
    hold on
end

hold off