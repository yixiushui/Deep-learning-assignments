function [sys,x0,str,ts] = spacemodel(t,x,u,flag)
%主函数
%主函数包含四个输出：
%                 sys数组包含某个子函数返回的值
%                 x0为所有状态的初始化向量
%                 str是保留参数，总是一个空矩阵
%                 Ts返回系统采样时间
%函数的四个输入分别为采样时间t、状态x、输入u和仿真流程控制标志变量flag
%输入参数后面还可以接续一系列的附带参数simStateCompliance
switch flag,
  case 0,
      [sys,x0,str,ts]=mdlInitializeSizes;
  case 1,
    sys=mdlDerivatives(t,x,u);

  case 3,
    sys=mdlOutputs(t,x,u);

  case {2,4,9}
    sys=[];
  otherwise
    error(['Unhandled flag = ', num2str(flag)]);
 
end
%主函数结束
%下面是各个子函数，即各个回调过程
function [sys,x0,str,ts]=mdlInitializeSizes
global node c b Fai
node = 7;
c = 0.1*[-1.5,-1,-0.5,0,0.5,1.0,1.5;
        -1.5,-1,-0.5,0,0.5,1.0,1.5;
        -1.5,-1,-0.5,0,0.5,1.0,1.5;
        -1.5,-1,-0.5,0,0.5,1.0,1.5;
        -1.5,-1,-0.5,0,0.5,1.0,1.5];
b = 10;
Fai = 5*eye(2);
%初始化回调子函数
%提供状态、输入、输出、采样时间数目和初始状态的值
%初始化阶段，标志变量flag首先被置为0，S-function首次被调用时
%该子函数首先被调用，且为S-function模块提供下面信息
%该子函数必须存在
sizes = simsizes;
%生成sizes数据结构，信息被包含在其中
sizes.NumContStates  = 2*node;
%连续状态数，缺省为0
sizes.NumDiscStates  = 0;
%离散状态数，缺省为0
sizes.NumOutputs     = 3;
%输出个数，缺省为0
sizes.NumInputs      = 11;
%输入个数，缺省为0
sizes.DirFeedthrough = 1;
%是否存在直馈通道，1存在，0不存在
sizes.NumSampleTimes = 0;
%采样时间个数，至少是一个
sys = simsizes(sizes);
%返回size数据结构所包含的信息
x0  = 0.1*ones(1,2*node);

%设置初始状态
str = [];
%保留变量置空
ts  = [];
%设置采样时间
function sys=mdlDerivatives(t,x,u)
global node c b Fai
qd1 = u(1);
d_qd1 = u(2);
dd_qd1 = u(3);
qd2 = u(4);
d_qd2 = u(5);
dd_qd2 = u(6);

q1 = u(7);
d_q1 = u(8);
q2 = u(9);
d_q2 = (10);

e1 = qd1-q1;
e2 = qd2-q2;
de1 = d_qd1-d_q1;
de2 = d_qd2-d_q2;
e = [e1;e2];
de = [de1;de2];
r = de+Fai*e;

qd = [qd1;qd2];
dqd = [d_qd1;d_qd2];
ddqd = [dd_qd1;dd_qd2];

z1 = [e(1);de(1);qd(1);dqd(1);ddqd(1)];
z2 = [e(2);de(2);qd(2);dqd(2);ddqd(2)];

for j = 1:1:node
    h1(j) = exp(-norm(z1-c(:,j))^2/(b*b));
    h2(j) = exp(-norm(z2-c(:,j))^2/(b*b));
end
F = 1.5*eye(node);
for i = 1:1:node
    sys(i) = 1.5*h1(i)*r(1);
    sys(i+node) = 1.5*h2(i)*r(2);
end
%计算导数回调子函数
%给定t,x,u计算连续状态的导数，可以在此给出系统的连续状态方程
%该子函数可以不存在
%sys表示状态导数，即dx

function sys=mdlOutputs(t,x,u)
global node c b Fai
qd1 = u(1);
d_qd1 = u(2);
dd_qd1 = u(3);
qd2 = u(4);
d_qd2 = u(5);
dd_qd2 = u(6);

q1 = u(7);
d_q1 = u(8);
q2 = u(9);
d_q2 = (10);

e1 = qd1-q1;
e2 = qd2-q2;
de1 = d_qd1-d_q1;
de2 = d_qd2-d_q2;
e = [e1;e2];
de = [de1;de2];
r = de+Fai*e;

qd = [qd1;qd2];
dqd = [d_qd1;d_qd2];
ddqd = [dd_qd1;dd_qd2];

W_f1 = [x(1:node)]';
W_f2 = [x(node+1:2*node)]';
z1 = [e(1);de(1);qd(1);dqd(1);ddqd(1)];
z2 = [e(2);de(2);qd(2);dqd(2);ddqd(2)];

for j = 1:1:node
    h1(j) = exp(-norm(z1-c(:,j))^2/(b*b));
    h2(j) = exp(-norm(z2-c(:,j))^2/(b*b));
end

fn = [W_f1*h1';W_f2*h2'];
Kv = 20*eye(2);
bd = 0.1;
epn = 0.20;
v = -(bd+epn)*sign(r);
fn_norm = norm(fn);
tol = fn+Kv*r-v;

sys(1) = tol(1);
sys(2) = tol(2);
sys(3) = fn_norm;


%计算输出回调函数
%给定t,x,u计算输出，可以在此描述系统的输出方程
%该子函数必须存在
%sys表示输出，即y
