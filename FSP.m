%This is the  file for the FSP analysis for each gene in each
%state. You need to make sure the correct gene is selected, that the
%correct table subset is selected (fn2) based on the state of the cell in
%question, and that you update the names for the number of TS and the burst
%sizes. Also make sure the half life is set correctly (td).

clear all
load('KLdatamatrix.mat'); 

%% Subset on the group you want to analyze
fn2 = fn;
idsim=find(fn2.CellType =='CMP');% & fn.p1TS <=2);
dat=table2array(fn2(idsim,[27,29,30])); %[mrna, nascent total, TS] Gata2 = 27,29,30, Gata1 = 19,21,22, PU1 = 11,13,14

activTS = find(dat(:,3)>=1);
b= nanmean(dat(activTS,2)./dat(activTS,3));
nu = mean(dat(:,3)./2);
mx = max(dat(:,1))*2;    %This will set the FSP truncation to 2* the max mRNA observed
td = 20; %Must update this for the half life of the mRNA you want to analyze

P0 = zeros(mx*4,1);
P0(1)=1;

tau =  0.1:0.1:60;%Set range for parameter sweep for tau
kini = 0.1:0.1:10;  %initiation rates


Tend =360;
paramOptimazation = zeros(numel(tau),3);

for ii =1:numel(tau)
    ii
    tic
    tau0 = tau(ii);
    paramOptimazation(ii,1) = tau0;
    [Pm, G]=timePropagate(tau0, mx,nu,b,td,P0,Tend);

    Obs = zeros(size(Pm,1),1);
    h = histcounts(dat(:,1),'Normalization','pdf','BinMethod','integers');
    Obs(1:numel(h)) = h;

    TSf = histcounts(dat(:,3),'Normalization','pdf','BinMethod','integers');
    ObsTS = TSf;
    paramOptimazation(ii,2) = getKLD(Obs, Pm);    %KLD for mRNA distribution
    paramOptimazation(ii,3) = getKLD(ObsTS,G);    %KLD for TS distribution
    toc
end
% 
%% 
figure(1)
subplot(1,2,1)
scatter(paramOptimazation(:,1),paramOptimazation(:,2))
subplot(1,2,2)
scatter(paramOptimazation(:,1),paramOptimazation(:,3))


[r,c] = find(paramOptimazation(:,2) == max(paramOptimazation(:,2)));
% [r,c] = find(paramOptimazation(:,2) == max(paramOptimazation(:,2)));

tID = paramOptimazation(r,1)

kon = nu/tID
koff = (1-nu)/tID
kini = b*koff
kd = log(2)/td
kr = nu*kini

figure(2)
[Pm, G]=timePropagate(tID,mx,nu,b,td,P0,Tend);
subplot(1,2,1)
histogram(dat(:,1),'Normalization','pdf','BinMethod','sturges','FaceAlpha',0.8,'FaceColor','k')
hold on
plot(Pm,'Color','r','LineWidth',2)
subplot(1,2,2)
histogram(dat(:,3),'Normalization','pdf','BinMethod','integers','FaceAlpha',0.8,'FaceColor','k')
hold on
plot([0 1 2], G, 'Color','r','LineWidth',2)
hold off

