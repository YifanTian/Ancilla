

E2_1 = load('ancillat_FE.txt'); 
%E2_2 = load('ancillat_FE_2.txt'); 
E1_1 = load('FiniteTE.txt'); 
%E1_2 = load('FiniteTE_2.txt');

%E2_1 = [E2_1(:,1),E2_1(:,2)];
%E2_2 = [E2_2(:,1),E2_2(:,2)];
E3 = ones(40,1)*(-4.2580352);
%p1 = polyfit();

plot(E1_1(:,1),[E1_1(:,2)-E2_1(:,2)],'-','LineWidth',1.5)
legend('E(ancilla) - E(exact)')
ylabel('Energy')
xlabel('Temperature(beta)')

plot(E1_1(:,1),[E1_1(:,2),E3(:,1)],'-','LineWidth',1.8)
hold on
plot(E2_1(:,1),E2_1(:,2),'o','LineWidth',1.8)
legend('exact finite','exact ground','ancilla')
title('Finite temperature');
xlabel('Temperature(beta)')
ylabel('Energy')

%E20 = load('myfile_1_20.txt'); 
%figure(1); plot(E2_1(:,1),[E1(:,2),E2_1(:,2),E3(:,1)],'-o','LineWidth',1.5);
figure(1); plotfit(E1(:,1),[E1(:,2),E3(:,1)],'-o','LineWidth',1.5);
hold on
figure(1); plot(E2_1(:,1),[E2_1(:,2),E3(:,1)],'-o','LineWidth',1.5);
%legend('termal density matrix','ancilla','exact')
%ylim([-45,-35])
%ylim([-0.5,2])
grid on
%legend('energy')
title('Finite temperature');
xlabel('temperature')
ylabel('energy')


% Sn = load('myfile_Sn.txt');
% %C = load('myfile_C.txt');
% SS = load('myfile_SS.txt');
% figure(2); plot(Sn(1:99,1),Sn(1:99,2:3),'-','LineWidth',1.5);
% legend('Sz','Sx')
% hold on
% figure(2); plot(Sn(1:99,1),SS(:,2),'o','LineWidth',1.5);

%legend('SS')