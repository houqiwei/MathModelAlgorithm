close all;
clear all;
popSize=60;%种群的数量
n=10;%城市的数量
nSalemen=2;%旅行商的数量
minTour=3;%最少的访问的城市
numIter=1;%当前的迭代次数
maxIter=2000;%最大的迭代次数
c1=0.5;%认知系数
c2=0.7;%社会学习系数
w=0.96-numIter/maxIter;%惯性权重
city=[ 263.1000  728.7000;
    218.3333  612.0000;
    168.4000  542.6000;
    182.6000  462.6000;
    117.0000  403.2857;
    92.8333  309.5000;
    12.8333  451.3333;
    161.2000  658.6000;
    110.2000  564.8000;
    105.0000  473.8000];%随机生成的城市的坐标
for i=1:n
    for j=1:n
        cityDist(i,j)=sqrt((city(i,1)-city(j,1)).^2+(city(i,2)-city(i,2)).^2);
    end
end%计算城市之间的距离
%初始化路径和断点
popRoute(1,:)=(1:n);
popBreak(1,:)=rand_breaks(nSalemen,minTour,n);
for i=2:popSize
    popRoute(i,:)=randperm(n);
    popBreak(i,:)=rand_breaks(nSalemen,minTour,n);
end
%初始化速度
for i=1:popSize
    velocity(i,:)=round(rand(1,n)*n);
end

figure('Name','MTSP_PSO | 过程','Numbertitle','off');
subplot(2,2,1);
pclr=~get(0,'DefaultAxesColor');
for i=1:n
    x(1,i)=city(i,1);
    y(1,i)=city(i,2);
end
plot(x,y,'*','Color',pclr);
axis([-100,500,0,850]);
set(gca,'XTick',-100:100:500);
set(gca,'YTick',0:100:850);

%xlabel('m');
%ylabel('m');
title('坐标点');
for j=1:popSize
    pRoute=popRoute(j,:);
    pBreak=popBreak(j,:);
    d=0;
    rngs=[[1 pBreak+1];[pBreak n]]';
    for k=1:nSalemen
        for p=rngs(k,1):rngs(k,2)-1
            d=d+cityDist(pRoute(p),pRoute(p+1));
        end
        d=d+cityDist(pRoute(rngs(k,1)),pRoute(rngs(k,2)));
    end
    eachPathDis(j)=d;
end
indivdualBestFitness=eachPathDis;
[minDist,index]=min(eachPathDis);
globalBestFitness=minDist;
for i=1:popSize
    globalBestRoute(i,:)=popRoute(index,:);
    globalBestBreak(i,:)=popBreak(index,:);
end
indivdualBestRoute=popRoute;
indivdualBestBreak=popBreak;
%迭代开始
for i=1:maxIter
    globalBest(numIter)=globalBestFitness;
    numIter=numIter+1;
    %开始PSO的过程
    pij_xij=GenerateChangeNums(popRoute,indivdualBestRoute);
    pij_xij=HoldByOdds(pij_xij,c1);
    pgj_xij=GenerateChangeNums(popRoute,globalBestRoute);
    pgj_xij=HoldByOdds(pgj_xij,c2);
    
    velocity=HoldByOdds(velocity,w);
    
    popRoute=PathExchange(popRoute,velocity);
    popRoute=PathExchange(popRoute,pij_xij);
    popRoute=PathExchange(popRoute,pgj_xij);
    %计算粒子群优化后各路径的距离
    for j=1:popSize
        pRoute=popRoute(j,:);
        pBreak=popBreak(j,:);
        d=0;
        rngs=[[1 pBreak+1];[pBreak n]]';
        for k=1:nSalemen
            for p=rngs(k,1):rngs(k,2)-1
                d=d+cityDist(pRoute(p),pRoute(p+1));
            end
            d=d+cityDist(pRoute(rngs(k,1)),pRoute(rngs(k,2)));
        end
        eachPathDis(j)=d;
    end
    
    %加入flip和swap过程，此时要更新断点
    newPopRoute=zeros(popSize,n);
    newPopBreak=zeros(popSize,nSalemen-1);
    randomOrder=randperm(popSize);
    for j=4:4:popSize
        %chooseIdx=randperm(popSize);
        rts=popRoute(randomOrder(j-3:j),:);
        
        brs=popBreak(randomOrder(j-3:j),:);
        
        dists=eachPathDis(randomOrder(j-3:j));
        [ignore,index]=min(dists);
        bestOf4Route=rts(index,:);
        bestOf4Break=brs(index,:);
        routeInsertionPoints = sort(ceil(n*rand(1,2)));
        I = routeInsertionPoints(1);
        J = routeInsertionPoints(2);
        for k=1:4
            tmpPopRoute(k,:)=bestOf4Route;
            tmpPopBreak(k,:)=bestOf4Break;
        end
        for k=1:4
            switch k
                case 2
                    tmpPopRoute(k,I:J)=tmpPopRoute(k,J:-1:I);
                case 3
                    tmpPopBreak(k,:)=rand_breaks(nSalemen,minTour,n);
                case 4
                    tmpPopRoute(k,I:J)=tmpPopRoute(k,J:-1:I);
                    tmpPopBreak(k,:)=rand_breaks(nSalemen,minTour,n);%
                otherwise
            end
        end
        newPopRoute(j-3:j,:)=tmpPopRoute;
        newPopBreak(j-3:j,:)=tmpPopBreak;
    end
    popRoute=newPopRoute;
    popBreak=newPopBreak;
    %改进部分完成，但是好像还有点问题
    %计算改进后的路径距离
    for j=1:popSize
        pRoute=popRoute(j,:);
        pBreak=popBreak(j,:);
        d=0;
        rngs=[[1 pBreak+1];[pBreak n]]';
        for k=1:nSalemen
            for p=rngs(k,1):rngs(k,2)-1
                d=d+cityDist(pRoute(p),pRoute(p+1));
            end
            d=d+cityDist(pRoute(rngs(k,1)),pRoute(rngs(k,2)));
        end
        newEachPathDis(j)=d;
    end
    %然后更新全局和个体最优
    [minNewPath,idx]=min(newEachPathDis);
    if minNewPath<globalBestFitness
        globalBestFitness=minNewPath;
        for k=1:popSize
            globalBestRoute(k,:)=popRoute(idx,:);
            globalBestBreak(k,:)=popBreak(idx,:);
        end
        
    end
    IsChange=newEachPathDis<indivdualBestFitness;
    indivdualBestRoute(find(IsChange),:)=popRoute(find(IsChange),:);
    indivdualBestBreak(find(IsChange),:)=popBreak(find(IsChange),:);
end
subplot(2,2,2);

plot((1:maxIter),globalBest(1:maxIter));
title(sprintf('迭代后的最优结果是:%1.4f',globalBest(maxIter)));
pRoute=globalBestRoute(1,:);
pBreak=globalBestBreak(1,:);
subplot(2,2,3);
pathPlot(pRoute,pBreak,city,nSalemen,n);
title(sprintf('最后的路由'));

%迭代结束

function Hold=HoldByOdds(Hold,odds)
[x y]=size(Hold);
for i=1:x
    for j=1:y
        if rand>odds
            Hold(i,j)=0;
        end
    end
end

end

function changeNums=GenerateChangeNums(popRoute,BestRoute)
[x y]=size(popRoute);
changeNums=zeros(x,y);
for i=1:x
    pop=BestRoute(i,:);
    pop1=popRoute(i,:);
    for j=1:y
        NoFromBestVar=pop(j);%y选择的是最优的某个位置
        for k=1:y
            NoFromRoute=pop1(k);%k是该粒子的某个位置
            if (NoFromRoute==NoFromBestVar)&&(j~=k)
                changeNums(i,j)=k;%i，j表示最优粒子的某行某列，里面的值是粒子的某列
                pop1(k)=pop1(j);
                pop1(j)=NoFromRoute;
            end
        end
    end
end
end

function Route=PathExchange(popRoute,Hold)
[x y]=size(popRoute);
for i=1:x
    pop=popRoute(i,:);
    a=Hold(i,:);
    for j=1:y
        if a(1,j)~=0
            pop1=pop(1,j);
            pop(1,j)=pop(1,a(1,j));
            pop(1,a(1,j))=pop1;
        end
    end
    Route(i,:)=pop;
end
end

function pathPlot(pRoute,pBreak,city,nSalemen,n)
clr=[1 0 0;0 0 1;0.67 0 1;0 1 0;1 0.5 0];
rngs=[[1 pBreak+1];[pBreak n]]';
for i=1:nSalemen
    rtes=pRoute([rngs(i,1):rngs(i,2) rngs(i,1)]);
    plot(city(rtes,1),city(rtes,2),'o-','Color',clr(i,:));
    hold on;
end
end

function breaks=rand_breaks(nSalemen,minTour,n)
nBreak=nSalemen-1;%4
dof=n-minTour*nSalemen;
addto=ones(1,dof+1);
for k=2:nBreak
    addto=cumsum(addto);
end
cumProb=cumsum(addto)/sum(addto);
if minTour==1
    tmpBreaks=randperm(n-1);
    breaks=sort(tmpBreaks(1:nBreak));
else
    nAdjust=find(rand<cumProb,1)-1;
    spaces=ceil(nBreak*rand(1,nAdjust));
    adjust=zeros(1,nBreak);
    for kk=1:nBreak
        adjust(kk)=sum(spaces==kk);
    end
    breaks=minTour*(1:nBreak)+cumsum(adjust);
end
end