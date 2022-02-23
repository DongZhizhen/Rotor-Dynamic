% The source code comes from https://wenku.baidu.com/view/22a34940be1e650e52ea99a4.html

%求解转子系统前三个临界转速和主振型的传递矩阵法
clc
clear
%等截面轴参数
l1=0.12;
d=0.04;
A=pi*d*d/4;
%轮盘参数
D=0.5;
h=0.025;
%盘轴材料参数（忽略轴的质量）
a=1;
u=0.3;
rou=7800;
E=2.0e11;
G=E/(2*(1+u));
I=pi*(d^4)/64;
K1=2.0e7;
v1=6*E*I/(a*G*A*l1*l1);
mi=rou*pi*D^2/4;     %轮盘的集中质量
Jp=mi*D^2/8;
Jd=Jp/2;
Ji=Jp-Jd;
%参数的数组形式
L=[l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 0 0];
M=[0 mi mi mi mi mi mi 0 0 0 0 0 mi  mi 0];
K=[K1 0 0 0 0 0 0 K1 0 0 0 K1 0 0 0];
v=[v1 v1 v1 v1 v1 v1 v1 v1 v1 v1 v1 v1 v1 0 0];
J=[0 Ji Ji Ji Ji Ji Ji 0 0 0 0 0 Ji Ji 0];
k=0;
Tit=['第一阶频率的振型和弯矩图';'第二阶频率的振型和弯矩图';'第三阶频率的振型和弯矩图'];
for w=0:0.01:4000
    for i=1:15
        T(:,:,i)=[1+(L(i)^3)*(1-v(i))*(M(i)*w^2-K(i))/(6*E*I) L(i)+L(i)^2*J(i)*w^2/(2*E*I) L(i)^2/(2*E*I) L(i)^3*(1-v(i))/(6*E*I);
            (L(i)^2)*(M(i)*w^2-K(i))/(2*E*I) 1+L(i)*J(i)*w^2/(E*I) L(i)/(E*I) L(i)^2/(2*E*I);
            L(i)*(M(i)*w^2-K(i)) J(i)*w^2 1 L(i);
            M(i)*w^2-K(i) 0 0 1];
    end
    H=T(:,:,1);
    for i2=2:15
        H=T(:,:,i2)*H;
    end
    F=H(3,1)*H(4,2)-H(3,2)*H(4,1);
    if F*(-1)^k<0 %求解临界转速
        k=k+1;
        wi(k)=w;
        w=wi(k);
        ni(k)=wi(k)*30/pi;
    end
end
for i1=1:3
    w=wi(i1);
    for j=1:14
        T(:,:,j)=[1+(L(j)^3)*(1-v(j))*(M(j)*w^2-K(j))/(6*E*I) L(j)+L(j)^2*J(j)*w^2/(2*E*I) L(j)^2/(2*E*I) L(j)^3*(1-v(j))/(6*E*I);
            (L(j)^2)*(M(j)*w^2-K(j))/(2*E*I) 1+L(j)*J(j)*w^2/(E*I) L(j)/(E*I) L(j)^2/(2*E*I);
            L(j)*(M(j)*w^2-K(j)) J(j)*w^2 1 L(j);
            M(j)*w^2-K(j) 0 0 1];
    end
    H=T(:,:,1);
    for j=2:15
        H=T(:,:,j)*H;
    end
    b=-H(4,1)/H(4,2);
    X(:,1)=([1 b 0 0]');
    for n=2:16
        X(:,n)=T(:,:,n-1)*X(:,n-1);  %相邻两质点右边的传递关系
    end
    for j1=1:15
        y(j1)=X(1,j1);
        z(j1)=X(3,j1);
        x(j1)=(j1-1)*l1;
    end
    y(16)=X(1,16);
    x(16)=1.56;
    z(16)=X(3,16);
    y=y/max(abs(y)); %归一化
    z=z/max(abs(z));
    subplot(3,1,i1)
    plot(x,y,'b-',x,z,'r:')
    title(Tit(i1,:))
    xlabel('轴长'), ylabel('不平衡值')
    axis([0,1.56,-1.2,1.2])
    grid on
    z;
end
legend('振型','弯矩')
ni
wi


%转子系统的不平衡响应
clc
clear
ww=[153.68 216.23 550.22 1449 2488.1 2495.6 3250.1 3714.7];  %前8阶固有频率
n=ww*30/pi;                                                  %前8阶转频
wi=[0.9*ww(1) (ww(1)+ww(2))/2];                %0.9*ww(1)和(ww(1)+ww(2))/2
%等截面轴参数
l1=0.12;
d=0.04;
A=pi*d*d/4;
%轮盘参数
D=0.5;
h=0.025;
%盘轴材料参数（忽略轴的质量）
rou=7800;
E=2.0e11;
I=pi*(d^4)/64;
K1=2.0e7;
m=rou*pi*D^2/8;           %轮盘的集中质量
Jp=m*D^2/8;  Jd=Jp/2;
J1=Jp-Jd;
u1=0.8e-4;
%参数的数组形式
L=[l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1 l1];
M=[0 m m m m m m 0 0 0 0 0 m m];
K=[K1 0 0 0 0 0 0 K1 0 0 0 K1 0 0];
J=[0 J1 J1 J1 J1 J1 J1 0 0 0 0 0 J1 J1];
Tit=['wi(1)时的振动响应图';'wi(2)时的振动响应图'];
U=[0 0 0 0 0 0 0 0 0 0 0 0 0 u1];
for i=1:2
    w=wi(i);
    n(i)=w*30/pi;
    for j=1:14
        T(:,:,j)=[1+(L(j)^3)*(M(j)*w^2-K(j))/(6*E*I) L(j)+L(j)^2*J(j)*w^2/(2*E*I) L(j)^2/(2*E*I) L(j)^3/(6*E*I) L(j)^3/(6*E*I)*U(j)*w^2;
            (L(j)^2)*(M(j)*w^2-K(j))/(2*E*I) 1+L(j)*J(j)*w^2/(E*I) L(j)/(E*I) L(j)^2/(2*E*I) L(j)^2/(2*E*I)*U(j)*w^2;
            L(j)*(M(j)*w^2-K(j)) J(j)*w^2 1 L(j) L(j)*U(j)*w^2;
            M(j)*w^2-K(j) 0 0 1 U(j)*w^2;
            0 0 0 0 1];
    end
    G=T(:,:,1);
    for j1=2:14
        H=T(:,:,j1)*G;
        G=H;
    end
    D1=H([3 4],[1 2]);
    B=H([3 4],[5 2]);B(:,1)=-B(:,1);
    C=H([3 4],[1 5]);C(:,2)=-C(:,2);
    b=det(B)/det(D1);c=det(C)/det(D1);
    X(:,1)=([b c 0 0 1]');
    for n=2:14
        X(:,n)=T(:,:,n-1)*X(:,n-1);    %相邻两质点右边的传递关系
    end

    y(1)=X(1,1);
    x(1)=0;
    for j2=2:13
        y(j2)=X(1,j2);
        x(j2)=x(j2-1)+L(j2-1);
    end
    y(14)=X(1,14);
    x(14)=1.56;
    xi=0:0.05:1.56;
    yi=interp1(x,y,xi,'spline');
    subplot(2,1,i)
    plot(xi,yi,'b-o','LineWidth',1.5)
    title(Tit(i,:))
    xlabel('轴长'), ylabel('不平衡值')
    grid on  
end
