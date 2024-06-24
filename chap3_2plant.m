function [sys,x0,str,ts] = sfuntmpl(t,x,u,flag)
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
global p g
%初始化回调子函数
%提供状态、输入、输出、采样时间数目和初始状态的值
%初始化阶段，标志变量flag首先被置为0，S-function首次被调用时
%该子函数首先被调用，且为S-function模块提供下面信息
%该子函数必须存在
sizes = simsizes;
%生成sizes数据结构，信息被包含在其中
sizes.NumContStates  = 4;
%连续状态数，缺省为0
sizes.NumDiscStates  = 0;
%离散状态数，缺省为0
sizes.NumOutputs     = 5;
%输出个数，缺省为0
sizes.NumInputs      = 3;
%输入个数，缺省为0
sizes.DirFeedthrough = 0;
%是否存在直馈通道，1存在，0不存在
sizes.NumSampleTimes = 0;
%采样时间个数，至少是一个
sys = simsizes(sizes);
%返回size数据结构所包含的信息
x0  = [0.09,0,-0.09,0];
%设置初始状态
str = [];
%保留变量置空
ts  = [];
%设置采样时间
p = [2.9,0.76,0.87,3.04,0.87];
g = 9.8;
function sys=mdlDerivatives(t,x,u)
global p g
M = [p(1)+p(2)+2*p(3)*cos(x(3)),p(2)+p(3)*cos(x(3));
    p(2)+p(3)*cos(x(3)),p(3)];
V = [-p(3)*x(4)*sin(x(3)),-p(3)*(x(2)+x(4))*sin(x(3));
    p(3)*x(2)*sin(x(3)),0];
G = [p(4)*g*cos(x(1))+p(5)*g*cos(x(1)+x(3));
    p(5)*g*cos(x(1)+x(3))];
dq = [x(2);x(4)];
F = 0.2*sign(dq);
told = [0.1*sin(t);0.1*sin(t)];
tol = u(1:2);
S = inv(M)*(tol-V*dq-G-F-told);
sys(1) = x(2);
sys(2) = S(1);
sys(3) = x(4);
sys(4) = S(2);
%计算导数回调子函数
%给定t,x,u计算连续状态的导数，可以在此给出系统的连续状态方程
%该子函数可以不存在

%sys表示状态导数，即dx

function sys=mdlOutputs(t,x,u)
%计算输出回调函数
%给定t,x,u计算输出，可以在此描述系统的输出方程
%该子函数必须存在
global p g
M = [p(1)+p(2)+2*p(3)*cos(x(3)),p(2)+p(3)*cos(x(3));
    p(2)+p(3)*cos(x(3)),p(3)];
V = [-p(3)*x(4)*sin(x(3)),-p(3)*(x(2)+x(4))*sin(x(3));
    p(3)*x(2)*sin(x(3)),0];
G = [p(4)*g*cos(x(1))+p(5)*g*cos(x(1)+x(3));
    p(5)*g*cos(x(1)+x(3))];
dq = [x(2);x(4)];
F = 0.2*sign(dq);
told = [0.1*sin(t);0.1*sin(t)];

qd1 = 0.1*sin(t);
d_qd1 = 0.1*cos(t);
dd_qd1 = -0.1*sin(t);
qd2 = 0.1*sin(t);
d_qd2 = 0.1*cos(t);
dd_qd2 = -0.1*sin(t);

q1 = x(1);
d_q1 = dq(1);
q2 = x(3);
d_q2 = dq(2);
q = [q1;q2];

e1 = qd1-q1;
e2 = qd2-q2;
de1 = d_qd1-d_q1;
de2 = d_qd2-d_q2;
e = [e1;e2];
de = [de1;de2];
Fai = 5*eye(2);
dqd = [d_qd1;d_qd2];
ddqd = [dd_qd1;dd_qd2];
dqr = dqd+Fai*e;
ddqr = ddqd+Fai*de;
f = M*ddqr+V*dqr+G+F;
f_form = norm(f);

sys(1) = x(1);
sys(2) = x(2);
sys(3) = x(3);
sys(4) = x(4);
sys(5) = f_form;
%sys表示输出，即y


