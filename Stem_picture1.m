load('0_result.mat')
x1=Event_t1;
x2=Event_k1;
d0=6034;
s0=0;
for i=1:1:d0
    s0=s0+Event_k1(i);
end
a0=s0/d0
figure();
subplot(2,1,1)
stem(x1,x2);

% load('1_result.mat')
% y1=Event_t2;
% y2=Event_k2;
% d1=6672;
% s1=0;
% for i=1:1:d1
%     s1=s1+Event_k2(i);
% end
% a1=s1/d1
% subplot(5,1,2)
% stem(y1,y2);

% load('3_result.mat')
% z1=Event_t3;
% z2=Event_k3;
% d3=975;
% s3=0;
% for i=1:1:d3
%     s3=s3+Event_k3(i);
% end
% a3=s3/d3
% subplot(5,1,3)
% stem(z1,z2);

load('4_result.mat')
u1=Event_t4;
u2=Event_k4;
d4=845;
s4=0;
for i=1:1:d4
    s4=s4+Event_k4(i);
end
a4=s4/d4
subplot(2,1,2)
stem(u1,u2);

% load('5_result.mat')
% p1=Event_t5;
% p2=Event_k5;
% d5=1042;
% s5=0;
% for i=1:1:d5
%     s5=s5+Event_k5(i);
% end
% a5=s5/d5
% subplot(5,1,5)
% stem(p1,p2);
% 
% A=(a0+a1+a3+a4+a5)/5