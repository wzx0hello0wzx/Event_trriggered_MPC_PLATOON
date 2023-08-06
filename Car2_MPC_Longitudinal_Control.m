function [sys,x0,str,ts] = Car2_MPC_Longitudinal_Control(t,x,u,flag)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%      SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                           restoring the complete simulation state of the
%                           model. The allowed values are: 'DefaultSimState',
%                           'HasNoSimState' or 'DisallowSimState'. If this value
%                           is not speficified, then the block's compliance with
%                           simState feature is set to 'UknownSimState'.


%   Copyright 1990-2010 The MathWorks, Inc.

%
% The following outlines the general structure of an S-function.
%
switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes;

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3
    sys=mdlOutputs(t,x,u);

  case{1,4,9}%Unused Flags
    sys=[];
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
%     DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
    error(['unhandled flag = ',num2str(flag)]); % Error handling
end

% end main

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;
%都是相对于本模块而言的，不是整个控制系统
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 6;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [30;3;0;40;3;0];
global U2;   %U 控制量aref
U2=[0];    %初始化U
global k2;%初始化U
k2=1;
%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0.05 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state

% end mdlInitializeSizes


%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = x;

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)

 global a b u_piao2;  %a、b分别表示线性化后的状态矩阵和输入矩阵 u_piao是计算出来的值
 global U2 k2;   %控制输入 速度和转角
 global kesi2;
 tic
%     profile on;
    Nx=3;   %状态量的个数
    Nu =1;  %控制量的个数
    Np =60; %预测步长
    Nc=30;  %控制步长
    Row=10; %松弛因子
    fprintf('Update start2, t=%6.3f\n',t)
    %t_d =u(3)*3.1415926/180;%CarSim输出的为角度，角度转换为弧度，u(3)为phi
    
      r(3)=u(6);  %参考加速度
      r(2)=u(5);%参考速度
      r(1)=u(4)-10;%车辆2的参考位置
      vd1=u(6);   %期望控制量(期望加速度)
    
    %矩阵的处理，先分配空间，再向空间幅值
    kesi2=zeros(Nx+Nu,1);    %组合后的状态量
    kesi2(1)=u(1)-r(1);  %u(1)==X(1)=x 当前时刻的状态
    kesi2(2)=u(2)-r(2);  %u(2)==X(2)=y
    kesi2(3)=u(3)-r(3);   %u(3)==X(3)=phi
    kesi2(4)=U2(1);   %上一时刻的控制输入v
    fprintf('kesi2=%4.2f\n',kesi2);
    %kesi(5)=U(2);   %上一时刻的控制输入delta_f
    fprintf('Update start2, u(1)=%4.2f\n',U2(1))
    %fprintf('Update start, u(2)=%4.2f\n',U(2))
    
    T=0.05; %采样时间
    %T_all=40;%临时设定，总的仿真时间，主要功能是防止计算期望轨迹越界
    % Mobile Robot Parameters
    %L = 2.6;
    % Mobile Robot variable
    
    % ==========================================================
    %矩阵初始化
    % ==========================================================
    u_piao2=zeros(Nx,Nu);  %u_piao是计算出来的值
     persistent u_real2 V2 Event_t2 Event_k2 eve2;
 if(isempty(u_real2))
     u_real2=zeros(Nc,Nu);
     V2=zeros(Nc,Nu);
     Event_k2=[];
     Event_t2=[];
     eve2=1;
 end

w=0.002;
% w=0;%DMPC

  if(k2==1||k2>30||abs(0-V2(k2))>w)

      Event_t2(eve2)=t;
      Event_k2(eve2)=k2;
      eve2=eve2+1;

      k2=1;
%     Q=100*eye(Nx*Np,Nx*Np); 
    q=[100 0 0;0 70 0;0 0 700];
    S_q=size(q,1);
    for m=0:1:Np-1
        Q(m*S_q+1:(m+1)*S_q,m*S_q+1:(m+1)*S_q)=q;
    end
    R=30*eye(Nu*Nc);
    a=[1, T,    0;
       0, 1,    T;
       0, 0, 1-T/0.5];
    b=[0;0;T/0.5];
    A_cell=cell(2,2);   %组合后的Ak
    B_cell=cell(2,1);
    A_cell{1,1}=a;
    A_cell{1,2}=b;
    A_cell{2,1}=zeros(Nu,Nx);
    A_cell{2,2}=eye(Nu);
    B_cell{1,1}=b;      %组合后的Bk
    B_cell{2,1}=eye(Nu);
    A=cell2mat(A_cell); %转化成矩阵形式
    B=cell2mat(B_cell);
    C=[1 0 0 0;0 1 0 0;0 0 1 0]; %组合后的C，值输出状态量X
    PHI_cell=cell(Np,1);    %预测公式的PHI
    THETA_cell=cell(Np,Nc); %预测公式的THETA
    for j=1:1:Np	
        PHI_cell{j,1}=C*A^j;
        for z2=1:1:Nc	%在给THETA赋值的时候，书上的THATA矩阵表述与这里不同。
            if z2<=j
                THETA_cell{j,z2}=C*A^(j-z2)*B;
            else 
                THETA_cell{j,z2}=zeros(Nx,Nu);
            end
        end
    end
    PHI=cell2mat(PHI_cell);%size(PHI)=[Nx*Np Nx+Nu]
    THETA=cell2mat(THETA_cell);%size(THETA)=[Nx*Np Nu*(Nc+1)]
    H_cell=cell(2,2);
    H_cell{1,1}=THETA'*Q*THETA+R;
    H_cell{1,2}=zeros(Nu*Nc,1);
    H_cell{2,1}=zeros(1,Nu*Nc);
    H_cell{2,2}=Row;            %What's this meanning?
    H=cell2mat(H_cell);

    error=PHI*kesi2;
    f_cell=cell(1,2);
    f_cell{1,1}=2*error'*Q*THETA;
    f_cell{1,2}=0;
%     f=(cell2mat(f_cell))';
    f=cell2mat(f_cell);
    
 %% 以下为约束生成区域
 %不等式约束
    A_t=zeros(Nc,Nc);%见falcone论文 P181，下三角全为1
    for p=1:1:Nc
        for q=1:1:Nc
            if q<=p  
                A_t(p,q)=1;
            else 
                A_t(p,q)=0;
            end
        end 
    end 
    A_I=kron(A_t,eye(Nu));%对应于falcone论文约束处理的矩阵A,求克罗内克积
    Ut=kron(ones(Nc,1),U2);%此处感觉论文里的克罗内科积有问题,暂时交换顺序
    umin=[-2];%维数与控制变量的个数相同
    umax=[2];
    delta_umin=[-0.2];
    delta_umax=[0.2];
    Umin=kron(ones(Nc,1),umin);
    Umax=kron(ones(Nc,1),umax);
    A_cons_cell={A_I zeros(Nu*Nc,1);-A_I zeros(Nu*Nc,1)};
    b_cons_cell={Umax-Ut;-Umin+Ut};
    A_cons=cell2mat(A_cons_cell);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围
    b_cons=cell2mat(b_cons_cell);%（求解方程）状态量不等式约束的取值
   % 状态量约束
    M=10; 
    delta_Umin=kron(ones(Nc,1),delta_umin);
    delta_Umax=kron(ones(Nc,1),delta_umax);
    lb=[delta_Umin;0];%（求解方程）状态量下界，包含控制时域内控制增量和松弛因子
    ub=[delta_Umax;M];%（求解方程）状态量上界，包含控制时域内控制增量和松弛因子
    
    %% 开始求解过程
%     options = optimset('Algorithm','active-set');
    options = optimset('Algorithm','interior-point-convex'); 
    [X,fval,exitflag]=quadprog(H,f,A_cons,b_cons,[],[],lb,ub,[],options);
      fprintf('exitflag2=%d\n',exitflag);
      fprintf('H2=%4.2f\n',H(1,1));
      fprintf('f2=%4.2f\n',f(1,1));
    %% 计算输出，反馈
    for l=1:1:Nc
    u_piao2(l)=X(l); %计算出的控制增量v
    %u_piao(2)=X(2); %计算出的控制增量delta
    U2(1)=kesi2(4)+u_piao2(1);%用于存储上一个时刻的控制量v
    u_real2(l)=kesi2(4)+u_piao2(l)+vd1;%u(k)
    end
    k2=k2+1;
    fprintf('k2=%4.2f\n',k2);
    sys = u_real2(1);
    Y=zeros(Nx*Np,1);
    XX=zeros(Nc,Nu);
    for i=1:1:Nc
    XX(i)=X(i);
    end
    Y=error+THETA*XX;%预测的输出序列
    for d=1:1:Nc
        V2(d)=Y(2+(d-1)*3);%获取速度偏差值序列
    end

  elseif (k2>1&&k2<=30&&abs(0-V2(k2))<=w)
        sys = u_real2(k2);
        k2=k2+1;
        fprintf('k2=%4.2f\n',k2);
    end
toc


fileDir = 'C:\Users\lenovo\Desktop\Paper1_simulinkModel\'; % 保存文件的路径
savePath = strcat(fileDir, num2str(1), '_result.mat'); % 拼接路径和文件名
% num2str(n) 在循环中可以生成1_result.mat, 2_result.mat....
% 拼接后的文件为 ： 文件路径/1_result.mat
save(savePath, 'Event_t2', 'Event_k2'); % 保存变量val1,val2到1_result.mat中
% end mdlOutputs