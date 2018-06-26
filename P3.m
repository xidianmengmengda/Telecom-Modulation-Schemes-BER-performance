%% 3.1
clear;
close all;
clc;

load('A2BPart3.mat')
fc = 20E6;
Ns = 1024;  
Rs = 1.5e6;
fs = Rs*Ns;

Nb = length(x);

t = linspace(0,1/Rs * ceil(Nb/2),Ns*ceil(Nb/2)+1);
t = t(1:end-1);
f = linspace(-fs/2,fs/2,length(t));


%% 3.2
xi = x(1:2:end);
xq = x(2:2:end);


%% 3.3
ci = cos(2*pi*fc*t);
cq = sin(2*pi*fc*t);
yi = kron(xi,ones(1,Ns)) .* ci;
yq = [kron(xq,ones(1,Ns)) zeros(1,Ns)] .* cq;

y = yi - yq;

%% 3.4
ri = sign(sum(reshape(y .* ci,Ns,ceil(Nb/2))));
rq = sign(sum(reshape(y .* -cq,Ns,ceil(Nb/2))));
rq_plot = rq;
rq = rq(1:end-1);
r = zeros(1,Nb);
r(1:2:end) = ri;
r(2:2:end) = rq;


error = sum(r ~= x);
%% 3.7
PE_A = qfunc((2*10^(-5/10))^0.5);
PE_B = qfunc((2*10^(5/10))^0.5);
PE_C = qfunc((2*10^(10/10))^0.5);


%% 3.8
Tsym = 1/Rs;
Tsam = 1/fs;

SNRA = 10*log10(2*10^(-5/10)*Tsam/Tsym);
SNRB = 10*log10(2*10^(5/10)*Tsam/Tsym);
SNRC = 10*log10(2*10^(10/10)*Tsam/Tsym);
ya = awgn(y,SNRA,'measured');
yb = awgn(y,SNRB,'measured');
yc = awgn(y,SNRC,'measured');

rai = sum(reshape(ya .* ci,Ns,ceil(Nb/2)));
raq = sum(reshape(ya .* -cq,Ns,ceil(Nb/2)));
rbi = sum(reshape(yb .* ci,Ns,ceil(Nb/2)));
rbq = sum(reshape(yb .* -cq,Ns,ceil(Nb/2)));
rci = sum(reshape(yc .* ci,Ns,ceil(Nb/2)));
rcq = sum(reshape(yc .* -cq,Ns,ceil(Nb/2)));
raq_plot = raq;
rbq_plot = rbq;
rcq_plot = rcq;
raq = raq(1:end-1);
rbq = rbq(1:end-1);
rcq = rcq(1:end-1);


ra = zeros(1,Nb);
ra(1:2:end) = rai;
ra(2:2:end) = raq;
ra = sign(ra);
errorA = sum(ra~= x);

rb = zeros(1,Nb);
rb(1:2:end) = rbi;
rb(2:2:end) = rbq;
rb = sign(rb);
errorB = sum(rb~= x);

rc = zeros(1,Nb);
rc(1:2:end) = rci;
rc(2:2:end) = rcq;
rc = sign(rc);
errorC = sum(rc~=x);
%% 3.7
figure;
subplot(2,2,1)
plot(ri,rq_plot);
subplot(2,2,2)
plot(rai,raq_plot);
subplot(2,2,3)
plot(rbi,rbq_plot);
subplot(2,2,4)
plot(rci,rcq_plot);


%% 3.11

test_x = string_to_bits(test_msg_str);
test_Nb = numel(test_x);

test_t = linspace(0,1/Rs * ceil(test_Nb/2),Ns*ceil(test_Nb/2)+1);
test_t = test_t(1:end-1);



test_xi = test_x(1:2:end);
test_xq = test_x(2:2:end);



test_ci = cos(2*pi*fc*test_t);
test_cq = sin(2*pi*fc*test_t);
test_yi = kron(test_xi,ones(1,Ns)) .* test_ci;
test_yq = kron(test_xq,ones(1,Ns)) .* test_cq;
test_y = test_yi - test_yq;



%% 3.12

qpsk_Nb = (length(qpsk_msg)/Ns);


qpsk_t = linspace(0,length(qpsk_msg)/fs,length(qpsk_msg)+1);
qpsk_t=qpsk_t(1:end-1);

qpsk_ci = cos(2*pi*fc*qpsk_t);
qpsk_cq = sin(2*pi*fc*qpsk_t);

qpsk_ri = sign(sum(reshape(qpsk_msg.*qpsk_ci,[Ns,qpsk_Nb])));
qpsk_rq = sign(sum(reshape(qpsk_msg.*-qpsk_cq,[Ns,qpsk_Nb])));

qpsk_r = zeros(1,2*qpsk_Nb);
qpsk_r(1:2:end) = qpsk_ri;
qpsk_r(2:2:end) = qpsk_rq;

qpsk_received_msg = bits_to_string(qpsk_r)




















% qpsk_Nb = 2*numel(qpsk_msg)/Ns;
% 
% qpsk_t = linspace(0,1/Rs * ceil(qpsk_Nb/2),length(qpsk_msg)+1);
% qpsk_t = qpsk_t(1:end-1);
% 
% qpsk_ci = cos(2*pi*fc*qpsk_t);
% qpsk_cq = sin(2*pi*fc*qpsk_t);
% 
% qpsk_ri = sum(reshape(qpsk_msg .* qpsk_ci,Ns,ceil(qpsk_Nb/2)));
% qpsk_rq = sum(reshape(qpsk_msg .* -qpsk_cq,Ns,ceil(qpsk_Nb/2)));
% 
% 
% qpsk_r = zeros(1,qpsk_Nb);
% qpsk_r(1:2:end) = qpsk_ri;
% qpsk_r(2:2:end) = qpsk_rq;
% qpsk_r = sign(qpsk_r);
% 
% qpsk_string_received = bits_to_string(qpsk_r);
% 
