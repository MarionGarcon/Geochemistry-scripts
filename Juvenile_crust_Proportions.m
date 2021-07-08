%% Calulcate juvenile crust proportions

clear all

%% Import Sedimentary database
data = xlsread('Data_Sed.xls');
A = length(data);
data = sortrows(data);

AgeStrati = data (:,1);
Sm147Nd144 = data (:,2);
Nd143Nd144_0 = data (:,3);
Nd143Nd144_T = data (:,4);
EpsNd_T = data (:,5);
TDM = data (:,6);
CrustAge = data (:,7); % Crust age is the difference between TDM and Age Stati

%% Outlier screening in the Sedimentary database
for i=2:A
    if TDM(i,1)> 4.56 || TDM(i,1) < 0 
       TDM(i,1) = NaN;
    end
    if TDM(i,1)< AgeStrati(i,1)
        if TDM(i,1) < (AgeStrati(i,1)-0.1)
        TDM(i,1) = NaN;
        else
        TDM(i,1) = AgeStrati(i,1);
        CrustAge(i,1) = 0;
        end
    end
    
    if isnan(TDM(i,1))
      AgeStrati(i,1) = NaN;
      Sm147Nd144(i,1) = NaN;
      Nd143Nd144_0(i,1) = NaN;
      Nd143Nd144_T(i,1) = NaN;
      EpsNd_T(i,1) = NaN;
      CrustAge(i,1) = NaN;
    end 
end

AgeStrati(any(isnan(AgeStrati),2),:) = [];
Sm147Nd144(any(isnan(Sm147Nd144),2),:) = [];
Nd143Nd144_0(any(isnan(Nd143Nd144_0),2),:) = [];
Nd143Nd144_T(any(isnan(Nd143Nd144_T),2),:) = [];
EpsNd_T(any(isnan(EpsNd_T),2),:) = [];
CrustAge(any(isnan(CrustAge),2),:) = [];
TDM(any(isnan(TDM),2),:) = [];

A = length(TDM);
Summary_average = NaN(19,12);
ii = 1;
x = NaN(10,1);
y = NaN(10,1);
x(1,1) = 1;
z = 1;

%% Calculate averages for each age bin
for i=2:A    
    limit = 0.2+ 0.2*(ii-1);
   
       if AgeStrati(i,1) > (limit-0.0000000000001)     
        y(z,1) = i-1;
        Summary_average(z,1) = limit - 0.1;
        Summary_average(z,2) = sum(~isnan(Sm147Nd144(x(z,1):y(z,1)))); 
        % 147Sm/144Nd ratios
        Summary_average(z,3) = nanmean(Sm147Nd144(x(z,1):y(z,1),1));
        Summary_average(z,4) = 2*std(Sm147Nd144(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1));
        % EpsNd(T)         
        Summary_average(z,5) = nanmean(EpsNd_T(x(z,1):y(z,1),1));
        Summary_average(z,6) = 2*std(EpsNd_T(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1));
        % TDM 
        Summary_average(z,7) = nanmean(TDM(x(z,1):y(z,1),1));
        Summary_average(z,8) = 2*std(TDM(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1));
        % Age of the eroded crust 
        Summary_average(z,9) = nanmean(CrustAge(x(z,1):y(z,1),1));
        Summary_average(z,10) = 2*std(CrustAge(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1));
        % 143Nd/144Nd ratios
        Summary_average(z,11) = nanmean(Nd143Nd144_T(x(z,1):y(z,1),1));
        Summary_average(z,12) = 2*std(Nd143Nd144_T(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1));
        
        x(z+1,1) = i;
        z = z+1;
        ii=ii+1; 
       end
end 
   
limit = 0.2+ 0.2*(ii-1);
y(z,1) = A;  

% Calculate averages for the last bin
Summary_average(z,1) = limit - 0.1;
Summary_average(z,2) = y(z,1)-x(z,1)+1; 
Summary_average(z,3) = nanmean(Sm147Nd144(x(z,1):y(z,1),1));
Summary_average(z,4) = 2*std(Sm147Nd144(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1),1);       
Summary_average(z,5) = nanmean(EpsNd_T(x(z,1):y(z,1),1));
Summary_average(z,6) = 2*std(EpsNd_T(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1),1);
Summary_average(z,7) = nanmean(TDM(x(z,1):y(z,1),1));
Summary_average(z,8) = 2*std(TDM(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1),1);
Summary_average(z,9) = nanmean(CrustAge(x(z,1):y(z,1),1));
Summary_average(z,10) = 2*std(CrustAge(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1),1);
Summary_average(z,11) = nanmean(Nd143Nd144_T(x(z,1):y(z,1),1));
Summary_average(z,12) = 2*std(Nd143Nd144_T(x(z,1):y(z,1),1)/sqrt(y(z,1)-x(z,1)+1));
        
MaxNbPoints = max(Summary_average(:,2));

Sm147Nd144_sorted = NaN(MaxNbPoints,z);
EpsNd_T_sorted = NaN(MaxNbPoints,z);
TDM_sorted = NaN(MaxNbPoints,z);
CrustAge_sorted = NaN(MaxNbPoints,z);

for i=1:z
    for ii=1:(y(i,1)-x(i,1)+1)
        Sm147Nd144_sorted(ii,i) = Sm147Nd144(x(i,1)+ii-1,1);
        EpsNd_T_sorted(ii,i) = EpsNd_T(x(i,1)+ii-1,1);
        TDM_sorted(ii,i) = TDM(x(i,1)+ii-1,1);
        CrustAge_sorted(ii,i) = CrustAge(x(i,1)+ii-1,1);
    end
end

%% Show the variations of various ratios through time
figure(1)
hold on
subplot(2,2,1)
boxplot(Sm147Nd144_sorted,Summary_average(:,1));
xlabel('Age strati (Ga) ');
ylabel('147Sm/144Nd');

subplot(2,2,2)
boxplot(EpsNd_T_sorted,Summary_average(:,1));
xlabel('Age strati (Ga) ');
ylabel('Epsilon Nd (T)');

subplot(2,2,3)
boxplot(TDM_sorted,Summary_average(:,1));
xlabel('Age strati (Ga) ');
ylabel('Nd model age (Ga)');

subplot(2,2,4)
boxplot(CrustAge_sorted,Summary_average(:,1));
xlabel('Age strati (Ga) ');
ylabel('Age of the eroded crust (Ga)');
hold off

for i=1:z
    if Summary_average(i,8) == 0
       Summary_average(i,1) = NaN;
    end
end
Summary_average(any(isnan(Summary_average),2),:) = [];

z = length(Summary_average);

figure(2)
hold on
plot(AgeStrati(:,1),EpsNd_T(:,1),'.k');
errorbar(Summary_average(:,1),Summary_average(:,5),Summary_average(:,6),'-or');
ylim([-20 10]);
xlabel('Age strati (Ga) ');
ylabel('Eps Nd (T)');
hold off

figure(3)
hold on
plot(AgeStrati,TDM,'.k');
errorbar(Summary_average(:,1),Summary_average(:,7),Summary_average(:,8),'-or');
xlabel('Age strati (Ga) ');
ylabel('Nd model age i.e. TDM (Ga)');
ylim([0 4]);
hold off

figure(5)
hold on
plot(AgeStrati(:,1),CrustAge(:,1),'.k');
errorbar(Summary_average(:,1),Summary_average(:,9),Summary_average(:,10),'-or');
xlabel('Age strati (Ga) ');
ylabel('Average age of the eroded crust (Ga)');
hold off

%% Calculation of the relative proportion of juvenile crust through time

%Initial conditions and definition of variables
Age_init = 4.4; %in Ga. Can be changed to 3.8 Ga

x_NewCrust = NaN(z,1);
t_old = NaN(z,1);
t_new = NaN(z,1);
t_new(:,1)= Summary_average(:,1);
t_old(z,1) = Age_init - Summary_average(z,1);
t_newPerc(:,1) = t_new(:,1);
t_newPerc(z+1,1) = Age_init;

for i=0:(z-1)
    x_NewCrust(z-i,1) = 1 - (Summary_average(z-i,9) / t_old(z-i,1));
        if x_NewCrust(z-i,1)<0 
        x_NewCrust(z-i,1)=0;
        end
    ii = i+1;
        if i < z-1
        t_old(z-ii,1) = (1-x_NewCrust(z-i,1)) * t_old(z-i,1) + t_new(z-i,1) - t_new(z-ii,1);
        end
end

PropNewCrust(:,1) = x_NewCrust(:,1) .*100;

figure(6)
barProp = [PropNewCrust(:,1)];
bar(t_new(:,1),barProp(:,:));
xlabel('Age (Ga)');
ylabel('Relative percentage of juvenile crust created through time (%)');

% Monte-Carlo simulation to estimate errors on XJuv (error propagation)
N=10000;
x_NewCrust_error = NaN(z,N);
t_old_error = NaN(z,N);
t_old_error(z,:) = 0.1;
PropNewCrust_error = NaN(z,N);

for j=1:N
    for i=0:(z-1)
        RandomAge = random('Normal',Summary_average(z-i,9),Summary_average(z-i,10)/2);
        x_NewCrust_error(z-i,j) = 1 - (RandomAge / t_old_error(z-i,j));
        if x_NewCrust_error(z-i,j)<0 
           x_NewCrust_error(z-i,j)=0;
        end
        
        ii = i+1;
  
        if i < z-1
           t_old_error(z-ii,j) = (1-x_NewCrust_error(z-i,j)) * t_old_error(z-i,j) + t_new(z-i,1) - t_new(z-ii,1);
        end
    end
        PropNewCrust_error(:,j) = x_NewCrust_error(:,j) .*100;    
end 

AveragePropNewCrust_error = NaN(z,1);
SDPropNewCrust_error = NaN(z,1);

for i=1:z
    AveragePropNewCrust_error(i,1)= mean(PropNewCrust_error(i,:));
    SDPropNewCrust_error(i,1)= 2*std(PropNewCrust_error(i,:));
end

figure(7)
hold on
plot(t_new(:,1),PropNewCrust(:,1),'-b');
errorbar(t_new(:,1),AveragePropNewCrust_error(:,1),SDPropNewCrust_error(:,1),'ob');
xlabel('Age (Ga)');
ylabel('Relative percentage of juvenile crust created through time (%)');
ylim([0 100]);
hold off
