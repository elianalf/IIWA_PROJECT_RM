clear all
close all

%% Definzione Traiettoria
 f=800; 
 N=8;
%p = [  x    y    z  ];
p(1,:) = [0.6 0 0.178];
p(2,:) = [0.7 0.1 0.18];  pos_c=2;   %numero del punto terminale della traiettoria circolare
p(3,:) = [0.65 0.2 0.19];
p(4,:) = [0.6 0.3 0.20]; 
p(5,:) = [0.5 0.25 0.21];
p(6,:) = [0.45 0.4 0.22];
p(7,:) = [0.35 0.3 0.23];
p(8,:) = [0.3 0.5 0.24];  
for i=1:N
    test(i)= sqrt(p(i,1)^2+p(i,2)^2);
    if(test(i)>0.75)
        error =  test(i)
        i;
    end
end
N=size(p,1);
%fisica realizzabilità


centro = [ 0.6,  0.1, 0.178]; %centroentro arco circonferenza da percorrere
ro = norm(p(1,:) - centro);
sci = ro*3*pi/2;        %valore iniziale della s
scf=ro*( 2*pi  );           %valore finale della s

T = [0 0.5 1.5 2 2.5 3 3.5 4]';

Ts=T*f; %istanti di tempo in numero campioni 
T_size=size(T,1);
delta_t = [0 0 200 200 200 200 150 150];

anticipo=sum(delta_t);
tempo = (T(1):(1/f):(T(end)- anticipo/f)); %tempo totale percorso in ms
tempo=tempo(1:end-1)'; 

% DELTA_T
DELTA_T(1) = 0;
for j=2:1:N
    DELTA_T(j) = DELTA_T(j-1) + delta_t(j);
end

%Calcolo s'(t)=s1(t), s1d(t), s1dd(t)
for j=2:1:N
    differenza =  p(j,:) - p(j-1,:);
    campioni = (T(j)-T(j-1))*f;
%     campioni=length(T(j-1):camp_al_sec:T(j));
    
    if (j~=pos_c) %diverso da 
        [sp ,spd ,spdd ] = lspb(0, norm(differenza), campioni);
        s1(1:size(sp),j-1)=sp;
        s1d(1:size(sp),j-1)=spd;
        s1dd(1:size(sp),j-1)=spdd;
    else 
        [sp ,spd ,spdd ] = lspb(sci, scf, campioni);
        s1c = sp;
        s1cd = spd;
        s1cdd = spdd;
        s1(1:size(sp),j-1)=sp;
        s1d(1:size(sp),j-1)=spd;
        s1dd(1:size(sp),j-1)=spdd;
    end
    
end
%Calcolo s(t), sd(t), sdd(t)
for j=2:1:N
    differenza =  p(j,:) - p(j-1,:);
    
    a = 0;
    b = Ts(j-1)-DELTA_T(j);
    c = Ts(j)-DELTA_T(j);
    d = Ts(T_size)-DELTA_T(N);
    
    bracket_camp = c - b;
    
    % S
    s(a+1:b+1,j-1) = 0;
    s(b+1:c,j-1) = s1(1:bracket_camp,j-1);
    
    if (j~=pos_c)
        s(c:d,j-1) = norm(differenza);
    else
        s(c:d,j-1) = scf;
    end 
    
    % SD
    sd(a+1:b+1,j-1) = 0;
    sd(b+1:c,j-1) = s1d(1:bracket_camp,j-1);
    sd(c:d,j-1) = 0;
    
    % SDD
    sdd(a+1:b+1,j-1) = 0;
    sdd(b+1:c,j-1) = s1dd(1:bracket_camp,j-1);
    sdd(c:d,j-1) = 0;
end


% Calcola pe,ped,pedd
dimp=size(tempo,1);
pe =zeros(dimp,3);
ped = zeros(dimp,3);
pedd = zeros(dimp,3);
pe(1,:) = p(1,:);

for j=2:1:N
    differenza =  p(j,:) - p(j-1,:);
    if (j~=pos_c)
        pe = pe + ((s(:,j-1)*differenza)/norm(differenza)); 
        ped = ped + ((sd(:,j-1)*differenza)/norm(differenza));
        pedd = pedd + ((sdd(:,j-1)*differenza)/norm(differenza));
        
    else
                                 
        [pe, ped, pedd] = circolare (s1c, s1cd, s1cdd, centro, ro, pe, ped, pedd, Ts(pos_c), s(:,j-1));  
    end
end

% 
% figure(1)
% hold on
% subplot(3,1,1);   
% hold on; grid on;
% plot (tempo, pe(:,1));
% xlabel('tempo[s]');ylabel('x');
% title('Posizione end-effector');
% 
% subplot(3,1,2);   
% hold on; grid on;
% plot (tempo, pe(:,2));
% xlabel('tempo[s]');ylabel('y');
% 
% subplot(3,1,3);   
% hold on;  grid on;
% plot (tempo, pe(:,3));
% xlabel('tempo[s]');ylabel('z');
% 
% figure(3)
% hold on
% subplot(3,1,1);   
% hold on; grid on;
% plot (tempo, ped(:,1));
% xlabel('tempo[s]');ylabel('Vx');
% title('Velocità end-effector');
% 
% subplot(3,1,2);   
% hold on; grid on;
% plot (tempo, ped(:,2));
% xlabel('tempo[s]');ylabel('Vy');
% 
% subplot(3,1,3);   
% hold on;  grid on;
% plot (tempo, ped(:,3));
% xlabel('tempo[s]');ylabel('Vz');
% 
% figure(4)
% hold on
% subplot(3,1,1);   
% hold on; grid on;
% plot (tempo, pedd(:,1));
% xlabel('tempo[s]');ylabel('Ax');
% title('Accelerazione end-effector');
% 
% subplot(3,1,2);   
% hold on; grid on;
% plot (tempo, pedd(:,2));
% xlabel('tempo[s]');ylabel('Ay');
% 
% subplot(3,1,3);   
% hold on;  grid on;
% plot (tempo, pedd(:,3));
% xlabel('tempo[s]');ylabel('Az');
% 

figure(4)
hold on
plot3(pe(:,1),pe(:,2),pe(:,3),'b');
grid on;
xlabel('x');ylabel('y');zlabel('z');