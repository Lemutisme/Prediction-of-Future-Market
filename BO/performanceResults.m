clear all, clc;

%变量初始化
dataFile='BO-5min.csv';
slpg=70;
PV=42000;

Sample=[datenum('01/1/2005'),datenum('04/14/2020')];

barsBack=17001;


Length=8410;
StopPct=0.021;


resultLabel={'Profit','WorstDrawDown','StDev','#trades'};

tradeTable=zeros(0,6); % trade#, Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position

resultSample=zeros(length(Length),length(StopPct),length(resultLabel));
resultOutSample=zeros(length(Length),length(StopPct),length(resultLabel));
i=1;j=1;
E0=100000;

trades=[];
position=0;

%上载资料
d=ezread(dataFile);
d.numTime=datenum(d.Time);
d.numTime=datenum(d.Date)+d.numTime-floor(d.numTime);
d.N=length(d.numTime);
d.M=5;

figure(1); clf(1); plot(d.numTime,d.Close,'b'); datetick('x');
%返回;

%样本中的索引:
indSample1=max(sum(d.numTime<Sample(1))+1,barsBack);
indSample2=max(sum(d.numTime<(Sample(2)+1)),barsBack);




%计算统计数据

    L=Length;
    S=StopPct;

    disp(['calculating performance for Length = ' num2str(L) ', StopPct = ' num2str(S)]);
    
    %计算 HH and LL
    HH=zeros(length(d.numTime),1);
    LL=zeros(length(d.numTime),1);
    
    for k=(barsBack+1):length(d.numTime)
        HH(k)=max(d.High((k-L):(k-1)));
        LL(k)=min(d.Low((k-L):(k-1)));
    end
    
        %设定初始条件:
        
        position=0;
        E=zeros(length(d.numTime),1)+E0;
        DD=zeros(length(d.numTime),1);
        trades=zeros(length(d.numTime),1);
        Emax=E0;
        
        %运行时间和交易:
        for k=(barsBack+1):length(d.numTime)
            traded=false;
            delta=PV*(d.Close(k)-d.Close(k-1))*position;
            
            if (position== 0)
                 
                buy=d.High(k)>=HH(k);
                sell=d.Low(k)<=LL(k);
                
                if (buy && sell)
                   delta = -slpg+PV*(LL(k)-HH(k));
                   tmp=[-1,LL(k),k,PV,slpg/2,-1;+1,HH(k),k,PV,slpg/2,0]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position

                   tradeTable=[tradeTable;tmp];
                   trades(k)=1;
                else
                    if(buy)
                        delta = -slpg/2 + PV*(d.Close(k)-HH(k));
                        position= 1;
                        traded=true;
                        benchmarkLong=d.High(k);
                        tmp=[+1,HH(k),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                   
                        trades(k)=0.5;
                    end
                    if(sell)
                        delta = -slpg/2 - PV*(d.Close(k)-LL(k));
                        position=-1;
                        traded=true;
                        benchmarkShort=d.Low(k);
                        tmp=[-1,LL(k),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                        trades(k)=0.5;
                    end
                end
                
            end
            
            if (position== 1 && ~traded) 
                sellShort=d.Low(k)<=LL(k);
                sell=d.Low(k)<=(benchmarkLong*(1-S));
                
                if(sellShort && sell)
                    %做空的拷贝
                    if(sellShort)
                        delta=delta-slpg-2*PV*(d.Close(k)-LL(k));
                        position=-1;
                        benchmarkShort=d.Low(k);
                        tmp=[-1,LL(k),k,PV,slpg/2,0;-1,LL(k),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                        trades(k)=1;
                    end
                
                else
                    if(sell)
                        delta=delta-slpg/2-PV*(d.Close(k)-(benchmarkLong*(1-S)));%min(Open,stopPrice)
                        position=0;
                        tmp=[-1,(benchmarkLong*(1-S)),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                        trades(k)=0.5;
                    end
                    
                    if(sellShort)
                        delta=delta-slpg-2*PV*(d.Close(k)-LL(k));%min(Open,LL(k))
                        position=-1;
                        benchmarkShort=d.Low(k);
                        tmp=[-1,LL(k),k,PV,slpg/2,0;-1,LL(k),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                        trades(k)=1;
                    end
                end
                
                benchmarkLong=max(d.High(k),benchmarkLong);
                
            end
            
            if (position==-1 && ~traded)
                buyLong=d.High(k)>=HH(k);
                buy=d.High(k)>=(benchmarkShort*(1+S));
                
                if(buyLong && buy)
                    %做多的拷贝
                    if(buyLong)
                        delta=delta-slpg+2*PV*(d.Close(k)-HH(k));
                        position=1;
                        benchmarkLong=d.High(k);
                        tmp=[+1,HH(k),k,PV,slpg/2,0;+1,HH(k),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                        trades(k)=1;
                    end
                
                else
                    if(buy)
                        delta=delta-slpg/2+PV*(d.Close(k)-(benchmarkShort*(1+S)));
                        position=0;                        
                        tmp=[+1,(benchmarkShort*(1+S)),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                        trades(k)=0.5;
                    end
                    
                    if(buyLong)
                        delta=delta-slpg+2*PV*(d.Close(k)-HH(k));
                        position=1;
                        benchmarkLong=d.High(k);
                        tmp=[+1,HH(k),k,PV,slpg/2,0;+1,HH(k),k,PV,slpg/2,position]; % Buy+trade#, 1/Sell-1 , Price, Time, PV, slpg, position
                        tradeTable=[tradeTable;tmp];
                        trades(k)=1;
                    end
                end
                
                benchmarkShort=min(d.Low(k),benchmarkShort);
            end
            
            if (position== 0 && traded),  end
            
            if (position== 1 && traded),  end
            
            if (position==-1 && traded),  end
            
            %更新权益
            E(k)=E(k-1)+delta; 
            %计算亏损
            Emax=max(Emax, E(k));
            DD(k)=E(k)-Emax;
        end
        
        %计算股票曲线的统计数据
        %收益损失计算
        PnL=[zeros(barsBack,1); E((barsBack+1):end)-E(barsBack:(end-1))];
        %采样时间段:
        resultSample(i,j,:)=[E(indSample2)-E(indSample1),min(DD(indSample1:indSample2)),std(PnL(indSample1:indSample2)),sum(trades(indSample1:indSample2))]; %{'Profit','WorstDrawDown','StDev','#trades'}

        

%Buy+1/Sell-1 , Price, Time, PV, slpg, position
disp('Trade, Price, TimeStamp, PointValue, slippage, position');
for i = 1: length(tradeTable(:,1))
  if(d.numTime(tradeTable(i,3))>Sample(1) && d.numTime(tradeTable(i,3))<=Sample(2))
    if(tradeTable(i,1)>0), tmp='Buy , '; else tmp='Sell, '; end;
    disp([tmp num2str(tradeTable(i,2)) ', ' datestr(d.numTime(tradeTable(i,3)),'yyyy-mm-dd HH:MM') ', ' num2str(tradeTable(i,4)) ', ' num2str(tradeTable(i,5)) ', ' num2str(tradeTable(i,6))]);
  end;
end;

%daily equity
ind=1:length(d.numTime);
ind=[ind(floor(d.numTime(1:end-1))~=floor(d.numTime(2:end))),ind(end)];
dates=floor(d.numTime(ind));
dE=E(ind);
dPNL=[0;dE(2:end)-dE(1:end-1)];


disp(' ');
disp('Performance for Sample Period')
disp(['Profit             : ' num2str(resultSample(1,1,1))]);
disp(['WorstDrawDown      : ' num2str(resultSample(1,1,2))]);
disp(['NPTWDD             : ' num2str(-resultSample(1,1,1)/resultSample(1,1,2)*100) ' %']);
disp(['#trades            : ' num2str(resultSample(1,1,4))]);
disp(['Avg Rate Of Return : ' num2str(mean(dPNL(d.numTime(ind)>Sample(1) & d.numTime(ind)<(Sample(2)+1)))/E0*253*100) ' %']);
disp(['Standard Deviation : ' num2str(std(dPNL(d.numTime(ind)>Sample(1) & d.numTime(ind)<(Sample(2)+1)))/E0*(253^0.5)*100) ' %']);
disp(['Sharpe Ratio       : ' num2str((mean(dPNL(d.numTime(ind)>Sample(1) & d.numTime(ind)<(Sample(2)+1)))/E0*253*100)/(std(dPNL(d.numTime(ind)>Sample(1) & d.numTime(ind)<(Sample(2)+1)))/E0*(253^0.5)*100))]);
disp(' ');
disp('you will need to process trade table to calculate: % of Winners, Average Winner, Average Loser, # of Winners,# of Losers, and some other characteristics');
 
figure(2);plot(d.numTime,E);
hold on;plot(d.numTime(ind),dE);hold off;
datetick('x');
figure(3);plot(d.numTime(ind(2000:2500)),dE(2000:2500)); datetick('x','mm/yyyy');

diary N2.txt