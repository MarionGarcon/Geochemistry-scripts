%% Estimate volume of crust for different beta based on the relative proportions of juvenile crust deduced from Sediments

clear all

% Import the relative proportions of juvenile crust deduced from Nd isotopic composition of sediments
load x_NewCrust.mat
A = length(x_NewCrust);

% Different beta values
B1 = 0; % No recycling
B2 = 1; % Recycling = New crust
B3 = 1.5; % Recycling = 1.5*New crust
B4 = 0.5; % New crust = 2*Recycling

Volume_CC = NaN(A,5);
Volume_CC(A,:) = 1; % The total volume of crust present at 3.7Ga is initially defined to be equal to 1. Will be renormalized to present later in the script. 
Volume_CC(:,2) = 1; 

Volume_NewCrust = NaN(A,5); %New crust means juvenile here
Volume_NewCrust(A,:) = x_NewCrust(A,2); 

Volume_RecyCrust = NaN(A,4);
Volume_RecyCrust(:,1) = 0; 
Volume_RecyCrust(A,2) = -x_NewCrust(A,2); 
Volume_RecyCrust(A,3) = -1.5*x_NewCrust(A,2); 
Volume_RecyCrust(A,4) = -0.5*x_NewCrust(A,2); 
Volume_RecyCrust(A,5) = x_NewCrust(A,2)-1;

x_RecyCrust = NaN(A,5);
x_RecyCrust(A,5) = 1-x_NewCrust(A,2);

for i=1:(A-1)
    
    %Case for which B1=0; No recycling
    Volume_CC(A-i,1) = Volume_CC(A-i+1,1)/(1-x_NewCrust(A-i,2));
    Volume_NewCrust(A-i,1) = x_NewCrust(A-i,2)*Volume_CC(A-i,1);
    
    %Case for which B2=1; Recycling = New crust
    Volume_NewCrust(A-i,2) = x_NewCrust(A-i,2)* Volume_CC(A-i,2);
    Volume_RecyCrust(A-i,2) = -x_NewCrust(A-i,2)* Volume_CC(A-i,2);
    
    %Case for which B3=1.5; Recycling = 1.5*New crust
    Volume_CC(A-i,3) = Volume_CC(A-i+1,3)/(1+0.5*x_NewCrust(A-i,2));
    Volume_NewCrust(A-i,3) = x_NewCrust(A-i,2)* Volume_CC(A-i,3);
    Volume_RecyCrust(A-i,3) = -1.5*x_NewCrust(A-i,2)* Volume_CC(A-i,3);
    
    %Case for which B4=0.5; New crust = 2*Recycling
    Volume_CC(A-i,4) = Volume_CC(A-i+1,4)/(1-x_NewCrust(A-i,2)*0.5);
    Volume_NewCrust(A-i,4) = x_NewCrust(A-i,2)* Volume_CC(A-i,4);
    Volume_RecyCrust(A-i,4) = -0.5*x_NewCrust(A-i,2)* Volume_CC(A-i,4);
    
end

% Normalization to present day volume of continental crust
Volume_CC_Today = NaN(A,5);
Volume_NewCrust_Today = NaN(A,5);
Volume_NewCrust_Today_bis = NaN(A,5);
Volume_Recycled_Today = NaN(A,5);
Volume_Recycled_Today_bis = NaN(A,5);

for i=1:A
    %Case for which B1=0; No recycling
    Volume_CC_Today(i,1) = Volume_CC(i,1)/Volume_CC(1,1);
    Volume_NewCrust_Today(i,1) =  x_NewCrust(i,2)* Volume_CC_Today(i,1);
    Volume_Recycled_Today(i,1) =  0;
    
    %Case for which B2=1; Recycling = New crust
    Volume_CC_Today(i,2) = Volume_CC(i,2)/Volume_CC(1,2);
    Volume_NewCrust_Today(i,2) =  x_NewCrust(i,2)* Volume_CC_Today(i,2);
    Volume_Recycled_Today(i,2) =  -x_NewCrust(i,2)* Volume_CC_Today(i,2);
    
    %Case for which B3=2; Recycling = 2*New crust
    Volume_CC_Today(i,3) = Volume_CC(i,3)/Volume_CC(1,3);
    Volume_NewCrust_Today(i,3) =  x_NewCrust(i,2)* Volume_CC_Today(i,3);
    Volume_Recycled_Today(i,3) =  -1.5*x_NewCrust(i,2)* Volume_CC_Today(i,3);
    
    %Case for which B4=0.5; New crust = 2*Recycling
    Volume_CC_Today(i,4) = Volume_CC(i,4)/Volume_CC(1,4);
    Volume_NewCrust_Today(i,4) =  x_NewCrust(i,2)* Volume_CC_Today(i,4);
    Volume_Recycled_Today(i,4) =  -0.5*x_NewCrust(i,2)* Volume_CC_Today(i,4);
end

figure(1)
hold on
stairs(x_NewCrust(:,1),Volume_CC_Today(:,1),'-r');
stairs(x_NewCrust(:,1),Volume_NewCrust_Today(:,1),'-b');
stairs(x_NewCrust(:,1),Volume_Recycled_Today(:,1),'-c');
axis([0 4 -0.2 1]);
xlabel('Age (Ga) ');
ylabel('Volume of crust normalized to present-day');
legend('Total crust','Juvenile crust','Recycled crust');
hold off

figure(2)
hold on
stairs(x_NewCrust(:,1),Volume_CC_Today(:,2),'-r');
stairs(x_NewCrust(:,1),Volume_NewCrust_Today(:,2),'-b');
stairs(x_NewCrust(:,1),Volume_Recycled_Today(:,2),'-c');
axis([0 4 -1 1]);
xlabel('Age (Ga) ');
ylabel('Volume of crust normalized to present-day');
legend('Total crust','Juvenile crust','Recycled crust');
hold off

figure(3)
hold on
stairs(x_NewCrust(:,1),Volume_CC_Today(:,3),'-r');
stairs(x_NewCrust(:,1),Volume_NewCrust_Today(:,3),'-b');
stairs(x_NewCrust(:,1),Volume_Recycled_Today(:,3),'-c');
axis([0 4 -5 5]);
xlabel('Age (Ga) ');
ylabel('Volume of crust normalized to present-day');
legend('Total crust','Juvenile crust','Recycled crust');
hold off

figure(4)
hold on
stairs(x_NewCrust(:,1),Volume_CC_Today(:,4),'-r');
stairs(x_NewCrust(:,1),Volume_NewCrust_Today(:,4),'-b');
stairs(x_NewCrust(:,1),Volume_Recycled_Today(:,4),'-c');
axis([0 4 -1 1]);
xlabel('Age (Ga) ');
ylabel('Volume of crust normalized to present-day');
legend('Total crust','Juvenile crust','Recycled crust');
hold off

