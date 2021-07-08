%% Convert Sm/Nd ratios into SiO2 contents %%

%%
clear all

minSiO2 = 45;
maxSiO2 = 75;

%% DATA FOR IGNEOUS ROCKS

% Import data from excel
dataSiO2 = xlsread('Data_IgneousRocks');
A = length(dataSiO2);

% Outlier screening in the igneous rock database
for i=1:A
    if dataSiO2(i,1) < minSiO2 || dataSiO2(i,1) >= maxSiO2
       dataSiO2(i,1) = NaN;
    end
    
    if dataSiO2(i,2) < 0.05 || dataSiO2(i,2) >= 0.4
       dataSiO2(i,2) = NaN;
    end
end

dataSiO2(any(isnan(dataSiO2),2),:) = [];


%% Calculate mean for every 1wt.% SiO2 %%

dataSiO2 = sortrows(dataSiO2);
A = length(dataSiO2);

AvdataSiO2 = NaN(maxSiO2-minSiO2,4);
ii = minSiO2;

x = 1;

for i=2:A
 % Calculate averages for each bin age
        if dataSiO2(i,1)>= ii+1
        y = i-1;
        AvdataSiO2(ii-minSiO2+1,1) = ii;
        AvdataSiO2(ii-minSiO2+1,2) = nanmean(dataSiO2(x:y,2));
        AvdataSiO2(ii-minSiO2+1,3) = 2*nanstd(dataSiO2(x:y,2))/sqrt(y-x);

 % Remove outlier above or below the mean +/-2SD in each bin age
        for k=x:y
    
            if dataSiO2(k,2) > AvdataSiO2(ii-minSiO2+1,2) +  AvdataSiO2(ii-minSiO2+1,3)*sqrt(y-x) || dataSiO2(k,2) < AvdataSiO2(ii-minSiO2+1,2) -  AvdataSiO2(ii-minSiO2+1,3)*sqrt(y-x)
            dataSiO2(k,2) = NaN;
            end
        end
                     
        AvdataSiO2(ii-minSiO2+1,4) = sum(~isnan(dataSiO2(x:y,2)));  
        ii=ii+1;
        x = i;
        end  
end 

y = i;
AvdataSiO2(ii-minSiO2+1,1) = ii;
AvdataSiO2(ii-minSiO2+1,2) = nanmean(dataSiO2(x:y,2));
AvdataSiO2(ii-minSiO2+1,3) = 2*nanstd(dataSiO2(x:y,2))/sqrt(y-x);

for k=x:y
    if dataSiO2(k,2) > AvdataSiO2(ii-minSiO2+1,2) + AvdataSiO2(ii-minSiO2+1,3)*sqrt(y-x) || dataSiO2(k,2) < AvdataSiO2(ii-minSiO2+1,2) -  AvdataSiO2(ii-minSiO2+1,3)*sqrt(y-x)
    dataSiO2(k,2) = NaN;
    end
end

AvdataSiO2(ii-minSiO2+1,2) = nanmean(dataSiO2(x:y,2));
AvdataSiO2(ii-minSiO2+1,3) = 2*nanstd(dataSiO2(x:y,2))/sqrt(sum(~isnan(dataSiO2(x:y,2))));
AvdataSiO2(ii-minSiO2+1,4) = sum(~isnan(dataSiO2(x:y,2)));  

Nbdata = sum(~isnan(dataSiO2(:,2)));

disp('_________________________')
disp(['Number of data taken for igneous rock regression: ',sprintf('%1.0f',Nbdata)])

figure(1)

hold on
errorbar(AvdataSiO2(:,2),AvdataSiO2(:,1),AvdataSiO2(:,3),'horizontal','om');
disp('Regression for igneous data - SiO2 contents vs. 147Sm/144Nd');
fitresult_SiO2 = fit(AvdataSiO2(:,2),AvdataSiO2(:,1),'poly3');
plot(fitresult_SiO2,'-m');
xlabel('147Sm/144Nd');
ylabel('SiO2 (wt%)');
plot(dataSiO2(:,2),dataSiO2(:,1),'.b');
hold off
disp('_________________________')



%% Estimate the standard errors for the different coefficients in regression by bootstrapping the residuals

B = length(AvdataSiO2);
SiO2_Modelled = NaN(B,3);
SiO2_Modelled(:,1) = AvdataSiO2(:,2); % Sm/Nd
SiO2_Modelled(:,2) = AvdataSiO2(:,1); % SiO2
p1 = fitresult_SiO2.p1; p2 = fitresult_SiO2.p2; p3 = fitresult_SiO2.p3; p4 = fitresult_SiO2.p4; 
SiO2_bootstrap = NaN(B,1);

for i=1:B
    % Calculate SiO2 predicted
    SiO2_Modelled(i,3) = p1 * SiO2_Modelled(i,1)^3 + p2 * SiO2_Modelled(i,1)^2 + p3 * SiO2_Modelled(i,1) + p4;
     % Calculate residuals
    SiO2_Modelled(i,4) = SiO2_Modelled(i,2) - SiO2_Modelled(i,3);
end

N = 1000;
coef_fit_SiO2 = NaN(N,4);

for i=1:N
    for ii=1:B
        if SiO2_Modelled(ii,3) + SiO2_Modelled(ii,4) > SiO2_Modelled(ii,3)
        SiO2_bootstrap(ii,1) = random('uniform',SiO2_Modelled(ii,3),SiO2_Modelled(ii,3) + SiO2_Modelled(ii,4));
        end
        
        if SiO2_Modelled(ii,3) + SiO2_Modelled(ii,4) < SiO2_Modelled(ii,3)
        SiO2_bootstrap(ii,1) = random('uniform',SiO2_Modelled(ii,3) + SiO2_Modelled(ii,4), SiO2_Modelled(ii,3));
        end
    end
    fitresult_SiO2 = fit(SiO2_Modelled(:,1),SiO2_bootstrap(:,1),'poly3');
    coef_fit_SiO2(i,1) = fitresult_SiO2.p1; coef_fit_SiO2(i,2) = fitresult_SiO2.p2; coef_fit_SiO2(i,3) = fitresult_SiO2.p3; coef_fit_SiO2(i,4) = fitresult_SiO2.p4;
end

Average_p1 = mean(coef_fit_SiO2(:,1));
Average_p2 = mean(coef_fit_SiO2(:,2));
Average_p3 = mean(coef_fit_SiO2(:,3));
Average_p4 = mean(coef_fit_SiO2(:,4));

Average_fit_SiO2 = NaN(B,1);
 
for ii=1:B
    Average_fit_SiO2(ii,1) = Average_p1 * SiO2_Modelled(ii,1)^3 + Average_p2  * SiO2_Modelled(ii,1)^2 + Average_p3 * SiO2_Modelled(ii,1) + Average_p4;
end

fit_SiO2 = NaN(B,N);

for i=1:N
    for ii=1:B
        fit_SiO2(ii,i) = coef_fit_SiO2(i,1) * SiO2_Modelled(ii,1)^3 + coef_fit_SiO2(i,2) * SiO2_Modelled(ii,1)^2 + coef_fit_SiO2(i,3) * SiO2_Modelled(ii,1) + coef_fit_SiO2(i,4);
    end
end


figure(1)
hold on
for i=1:N
    plot(SiO2_Modelled(:,1),fit_SiO2(:,i),'-c');
end
plot(SiO2_Modelled(:,1),Average_fit_SiO2(:,1),'-k');
hold off


%% MONTE-CARLO SIMULATION TO ESTIMATE THE AVERAGE SiO2 CONTENT OF THE CRUST FROM 147Sm/144Nd RATIOS OF SEDIMENTS

% Import sedimentary rock database
data = xlsread('Data_Sed.xls');
A = length(data);
data = sortrows(data);

%% Calculate mean Sm/Nd ratios every 200 Ma

% Outlier screening in the sedimentary rock database
AgeStrati = data (:,1);
Sm147Nd144 = data (:,2);
Nd143Nd144_0 = data (:,3);
Nd143Nd144_T = data (:,4);
EpsNd_T = data (:,5);
TDM = data (:,6);

for i=2:A
    if TDM(i,1)> 4.56 || TDM(i,1) < 0 
       TDM(i,1) = NaN;
    end
    
    if Sm147Nd144(i,1) < 0.05 || Sm147Nd144(i,1) >= 0.4
       Sm147Nd144(i,1) = NaN;
    end
    
    if TDM(i,1)< AgeStrati(i,1)
        if TDM(i,1) < AgeStrati(i,1)-0.1
        TDM(i,1) = NaN;
        else
        TDM(i,1) = AgeStrati(i,1);
        end
    end
    
    if isnan(TDM(i,1))
      AgeStrati(i,1) = NaN;
      Sm147Nd144(i,1) = NaN;
      Nd143Nd144_0(i,1) = NaN;
      Nd143Nd144_T(i,1) = NaN;
      EpsNd_T(i,1) = NaN;
    end
end

AgeStrati(any(isnan(AgeStrati),2),:) = [];
Sm147Nd144(any(isnan(Sm147Nd144),2),:) = [];
Nd143Nd144_0(any(isnan(Nd143Nd144_0),2),:) = [];
Nd143Nd144_T(any(isnan(Nd143Nd144_T),2),:) = [];
EpsNd_T(any(isnan(EpsNd_T),2),:) = [];
TDM(any(isnan(TDM),2),:) = [];

A = length(TDM);
Summary_average = NaN(19,12);
ii = 1;
x = NaN(10,1);
y = NaN(10,1);
x(1,1) = 1;
z = 1;

for i=2:A
 % Calculate averages for each bin age
    limit = 0.2+ 0.2*(ii-1);
    if (AgeStrati(i,1) > limit-0.0000000000001)       
        y(z,1) = i-1;
        Summary_average(z,1) = limit - 0.1;
        Summary_average(z,2) = sum(~isnan(Sm147Nd144(x(z,1):y(z,1)))); 
        Summary_average(z,3) = nanmean(Sm147Nd144(x(z,1):y(z,1),1));
        Summary_average(z,4) = 2*std(Sm147Nd144(x(z,1):y(z,1),1))/sqrt(y(z,1)-x(z,1)+1);
        
        for k=x(z,1):y(z,1)
            if Sm147Nd144(k,1)  > Summary_average(z,3) +  Summary_average(z,3) * sqrt(y(z,1)-x(z,1)+1) || Sm147Nd144(k,1)  < Summary_average(z,3) - Summary_average(z,3) * sqrt(y(z,1)-x(z,1)+1)
            Sm147Nd144(k,1) = NaN;
            end
        end
   
        Summary_average(z,3) = nanmean(Sm147Nd144(x(z,1):y(z,1),1));
        Summary_average(z,4) = 2*nanstd(Sm147Nd144(x(z,1):y(z,1),1))/sqrt(sum(~isnan(Sm147Nd144(x(z,1):y(z,1)))));       
        x(z+1,1) = i;
        z = z+1;
        ii=ii+1;
    end
end 

% Calculate averages for the last bin
limit = 0.2+ 0.2*(ii-1);
y(z,1) = A;  
Summary_average(z,1) = limit - 0.1;
Summary_average(z,2) = sum(~isnan(Sm147Nd144(x(z,1):y(z,1)))); 
Summary_average(z,3) = nanmedian(Sm147Nd144(x(z,1):y(z,1),1));
Summary_average(z,4) = 2*mad(Sm147Nd144(x(z,1):y(z,1),1))/sqrt(y(z,1)-x(z,1)+1);

for k=x(z,1):y(z,1)
   if Sm147Nd144(i,1)  > Summary_average(z,3) +  Summary_average(z,3) * sqrt(y(z,1)-x(z,1)+1) || Sm147Nd144(i,1)  < Summary_average(z,3) - Summary_average(z,3) * sqrt(y(z,1)-x(z,1)+1)
      Sm147Nd144(i,1) = NaN;
   end
end
        
Summary_average(z,3) = nanmean(Sm147Nd144(x(z,1):y(z,1),1));
Summary_average(z,4) = 2*nanstd(Sm147Nd144(x(z,1):y(z,1),1))/sqrt(sum(~isnan(Sm147Nd144(x(z,1):y(z,1)))));
       
Summary_average(20,:) = [];
z = length(Summary_average);

MaxNbPoints = max(Summary_average(:,2));
Sm147Nd144_sorted = NaN(MaxNbPoints,z);

for i=1:z
    for ii=1:(y(i,1)-x(i,1)+1)
    
        Sm147Nd144_sorted(ii,i) = Sm147Nd144(x(i,1)+ii-1,1);
            
    end
end

% % Moving average with a kernel of 200 Ma
MovingAv_SampleTime = NaN(A,4);
yy = NaN(A,1);
zz = 1;
yy(1,1) = 1;
limit = 0.2;

for i=2:A
    while (AgeStrati(i,1) > limit-0.00000000000001)
       limit = AgeStrati(zz,1) + 0.2;
       yy(zz,1) = i-1;
       MovingAv_SampleTime(zz,1) = nanmean(AgeStrati(zz:yy(zz,1)));
       MovingAv_SampleTime(zz,2) = sum(~isnan(Sm147Nd144(zz:yy(zz,1))));
       MovingAv_SampleTime(zz,3) = nanmean(Sm147Nd144(zz:yy(zz,1)));
       MovingAv_SampleTime(zz,4) = 2*nanstd(Sm147Nd144(zz:yy(zz,1),1))/sqrt(yy(zz,1)-zz+1);               
       zz = zz+1;   
    end
end
%-------------------------------

figure(3)
hold on
plot(AgeStrati,Sm147Nd144,'.k');
plot(Summary_average(:,1),Summary_average(:,3),'r');
errorbar(Summary_average(:,1),Summary_average(:,3),Summary_average(:,4),'or');
xlabel('Age strati (Ga) ');
ylabel('147Sm/144Nd');

figure(4)
hold on
plot(MovingAv_SampleTime(:,1),MovingAv_SampleTime(:,3),'b')
errorbar(MovingAv_SampleTime(:,1),MovingAv_SampleTime(:,3),MovingAv_SampleTime(:,4),'b')
plot(Summary_average(:,1),Summary_average(:,3),'-r');
errorbar(Summary_average(:,1),Summary_average(:,3),Summary_average(:,4),'or');
legend('data','Moving average 200Ma','Moving average 200Ma','Average by bin of 200Ma','Average by bin of 200Ma');
xlabel('Age strati (Ga)');
ylabel('147Sm/144Nd');
ylim([0.08 0.16]);


%% Calculation of SiO2 contents from Sm/Nd ratios of sediments

Modeled_SiO2_average = NaN(z,N);
Modeled_SiO2_average_random = NaN(z,N);
Average_modeled_SiO2_average = NaN(z,2);

figure(5)
hold on
M = 10000;

for i=1:M
    for ii = 1:z
        SmNd_average = random('normal',Summary_average(ii,3),Summary_average(ii,4)./2);
        W = randi(N);
        Modeled_SiO2_average_random(ii,i) = coef_fit_SiO2(W,1) * SmNd_average^3 + coef_fit_SiO2(W,2) * SmNd_average^2 + coef_fit_SiO2(W,3) * SmNd_average + coef_fit_SiO2(W,4);
    end
end

for ii = 1:z
    Average_modeled_SiO2_average(ii,1) = nanmean(Modeled_SiO2_average_random(ii,:));
    Average_modeled_SiO2_average(ii,2) = 2*nanstd(Modeled_SiO2_average_random(ii,:));  
end

plot(Summary_average(:,1),Average_modeled_SiO2_average(:,1),'-r');
errorbar(Summary_average(:,1),Average_modeled_SiO2_average(:,1),Average_modeled_SiO2_average(:,2),'or');
ylim([50 80]);
xlabel('Age (Ga)');
ylabel('SiO2 content of the crust (.wt%)');
hold off







