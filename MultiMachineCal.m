clear
close all

global N2 S tset

P0 = [0.92 0.96 0.94];
P1 = [0.74 0.79 0.81 0.71];
P2 = [0.71 0.87 0.79 0.81 0.82];
N0 = [8 5];
N1 = [5 4 2 3];
N2 = [8 2 4 6 3];
beta = [1 2];
gamma = [-2 0];
tset = 1;
t = 300;

PR = zeros(2,t+1);
CR = zeros(1,t+1);
WIP0 = zeros(size(N0,2),t+1);
WIP1 = zeros(size(N1,2),t+1);
WIP2 = zeros(size(N2,2),t+1);
ST0 = zeros(size(P0,2),t+1);
ST1 = zeros(size(P1,2),t+1);
ST2 = zeros(size(P2,2),t+1);
BL0 = zeros(size(P0,2)-1,t+1);
BL1 = zeros(size(P1,2)-1,t+1);
BL2 = zeros(size(P2,2)-1,t+1);

pf0 = zeros(1,size(P0,2));
pb0 = zeros(1,size(P0,2));
pf1 = zeros(1,size(P1,2));
pb1 = zeros(1,size(P1,2));
pf2 = zeros(1,size(P2,2));
pb2 = zeros(1,size(P2,2));

X0 = zeros(max(N0)+1,size(N0,2));
X0(1,:) = 1;
X1 = zeros(max(N1)+1,size(N1,2));
X1(1,:) = 1;
X2 = zeros(max(N2)+1,size(N2,2));
X2(1,:) = 1;

S = (N1(1)+1)*(N2(1)+1)*(tset(1)+1)*2;
X = zeros(S,1);
X(1) = 1;

for i = 1:t+1
    pf0(1) = P0(1);
    for j = 2:size(P0,2)
        pf0(j) = P0(j)*(1-X0(1,j-1));
    end

    pb1(size(P1,2)) = P1(size(P1,2));
    pb2(size(P2,2)) = P2(size(P2,2));
    for j = size(P1,2)-1:-1:1
        pb1(j) = P1(j)*(1-(1-pb1(j+1))*X1(N1(j+1)+1,j+1));
    end
    for j = size(P2,2)-1:-1:1
        pb2(j) = P2(j)*(1-(1-pb2(j+1))*X2(N2(j+1)+1,j+1));
    end
    
    A = zeros(S,S);
    for n1 = 0:N1(1)
        for n2 = 0:N2(1)
            for sta = 1:2
                for setup = 0:tset(1)
                    g1 = calnum(sta,[n1,n2],setup);
                    if (sta == 1)&&(setup == 0)
                        if n2 < beta(1,1)*n1+gamma(1,1)
                            if (n1 ~= 0)&&(n2 ~= 0)
                                g2 = calnum(2,[n1,n2],1);
                                A(g2,g1) = (1-pb1(1))*(1-pb2(1));
                                g2 = calnum(2,[n1,n2-1],1);
                                A(g2,g1) = (1-pb1(1))*pb2(1);
                                g2 = calnum(2,[n1-1,n2],1);
                                A(g2,g1) = pb1(1)*(1-pb2(1));
                                g2 = calnum(2,[n1-1,n2-1],1);
                                A(g2,g1) = pb1(1)*pb2(1);
                            elseif (n1 == 0)&&(n2 ~= 0)
                                g2 = calnum(2,[n1,n2],1);
                                A(g2,g1) = 1-pb2(1);
                                g2 = calnum(2,[n1,n2-1],1);
                                A(g2,g1) = pb2(1);
                            elseif (n1 ~= 0)&&(n2 == 0)
                                g2 = calnum(2,[n1,n2],1);
                                A(g2,g1) = 1-pb1(1);
                                g2 = calnum(2,[n1-1,n2],1);
                                A(g2,g1) = pb1(1);
                            elseif (n1 == 0)&&(n2 == 0)
                                g2 = calnum(2,[n1,n2],1);
                                A(g2,g1) = 1;
                            end
                        else
                            g3 = calnum(1,[n1,n2],0);
                            if (n1 == N1(1))&&(n2 ~= 0)
                                g4 = calnum(1,[n1,n2],0);
                                A(g4,g3) = (1-pb2(1))*(pf0(size(P0,2))*pb1(1)+1-pb1(1));
                                g4 = calnum(1,[n1-1,n2],0);
                                A(g4,g3) = (1-pb2(1))*((1-pf0(size(P0,2)))*pb1(1));
                                g4 = calnum(1,[n1,n2-1],0);
                                A(g4,g3) = pb2(1)*(pf0(size(P0,2))*pb1(1)+1-pb1(1));
                                g4 = calnum(1,[n1-1,n2-1],0);
                                A(g4,g3) = pb2(1)*((1-pf0(size(P0,2)))*pb1(1));
                            elseif (n1 == N1(1))&&(n2 == 0)
                                g4 = calnum(1,[n1,n2],0);
                                A(g4,g3) = 1-(1-pf0(size(P0,2)))*pb1(1);
                                g4 = calnum(1,[n1-1,n2],0);
                                A(g4,g3) = (1-pf0(size(P0,2)))*pb1(1);
                            elseif (n2 == 0)&&(n1 == 0)
                                g4 = calnum(1,[n1,n2],0);
                                A(g4,g3) = (1-pf0(size(P0,2)))*(1-pb1(1))+(1-pf0(size(P0,2)))*pb1(1);
                                g4 = calnum(1,[n1+1,n2],0);
                                A(g4,g3) = pf0(size(P0,2))*pb1(1)+pf0(size(P0,2))*(1-pb1(1));
                            elseif (n1 ~= 0)&&(n1 ~=N1(1))&&(n2 == 0)
                                g4 = calnum(1,[n1,n2],0);
                                A(g4,g3) = pf0(size(P0,2))*pb1(1)+(1-pf0(size(P0,2)))*(1-pb1(1));
                                g4 = calnum(1,[n1+1,n2],0);
                                A(g4,g3) = pf0(size(P0,2))*(1-pb1(1));
                                g4 = calnum(1,[n1-1,n2],0);
                                A(g4,g3) = (1-pf0(size(P0,2)))*pb1(1);
                            elseif (n1 == 0)&&(n2~=0)
                                g4 = calnum(1,[n1,n2],0);
                                A(g4,g3) = (1-pb2(1))*(1-pf0(size(P0,2)));
                                g4 = calnum(1,[n1+1,n2],0);
                                A(g4,g3) = (1-pb2(1))*pf0(size(P0,2));
                                g4 = calnum(1,[n1,n2-1],0);
                                A(g4,g3) = pb2(1)*(1-pf0(size(P0,2)));
                                g4 = calnum(1,[n1+1,n2-1],0);
                                A(g4,g3) = pb2(1)*pf0(size(P0,2));
                            else
                                g4 = calnum(1,[n1,n2],0);
                                A(g4,g3) = (1-pb2(1))*(pf0(size(P0,2))*pb1(1)+(1-pf0(size(P0,2)))*(1-pb1(1)));
                                g4 = calnum(1,[n1+1,n2],0);
                                A(g4,g3) = (1-pb2(1))*pf0(size(P0,2))*(1-pb1(1));
                                g4 = calnum(1,[n1-1,n2],0);
                                A(g4,g3) = (1-pb2(1))*(1-pf0(size(P0,2)))*pb1(1);
                                g4 = calnum(1,[n1,n2-1],0);
                                A(g4,g3) = pb2(1)*(pf0(size(P0,2))*pb1(1)+(1-pf0(size(P0,2)))*(1-pb1(1)));
                                g4 = calnum(1,[n1+1,n2-1],0);
                                A(g4,g3) = pb2(1)*pf0(size(P0,2))*(1-pb1(1));
                                g4 = calnum(1,[n1-1,n2-1],0);
                                A(g4,g3) = pb2(1)*(1-pf0(size(P0,2)))*pb1(1);
                            end
                        end
                    elseif (sta == 2)&&(setup == 0)
                        if n2 > beta(1,2)*n1+gamma(1,2)
                            if (n1 ~= 0)&&(n2 ~= 0)
                                g2 = calnum(1,[n1,n2],1);
                                A(g2,g1) = (1-pb1(1))*(1-pb2(1));
                                g2 = calnum(1,[n1,n2-1],1);
                                A(g2,g1) = (1-pb1(1))*pb2(1);
                                g2 = calnum(1,[n1-1,n2],1);
                                A(g2,g1) = pb1(1)*(1-pb2(1));
                                g2 = calnum(1,[n1-1,n2-1],1);
                                A(g2,g1) = pb1(1)*pb2(1);
                            elseif (n1 == 0)&&(n2 ~= 0)
                                g2 = calnum(1,[n1,n2],1);
                                A(g2,g1) = 1-pb2(1);
                                g2 = calnum(1,[n1,n2-1],1);
                                A(g2,g1) = pb2(1);
                            elseif (n1 ~= 0)&&(n2 == 0)
                                g2 = calnum(1,[n1,n2],1);
                                A(g2,g1) = 1-pb1(1);
                                g2 = calnum(1,[n1-1,n2],1);
                                A(g2,g1) = pb1(1);
                            elseif (n1 == 0)&&(n2 == 0)
                                g2 = calnum(1,[n1,n2],1);
                                A(g2,g1) = 1;
                            end
                        else
                            g3 = calnum(2,[n1,n2],0);
                            if (n2 == N2(1))&&(n1 ~= 0)
                                g4 = calnum(2,[n1,n2],0);
                                A(g4,g3) = (1-pb1(1))*(pf0(size(P0,2))*pb2(1)+1-pb2(1));
                                g4 = calnum(2,[n1,n2-1],0);
                                A(g4,g3) = (1-pb1(1))*((1-pf0(size(P0,2)))*pb2(1));
                                g4 = calnum(2,[n1-1,n2],0);
                                A(g4,g3) = pb1(1)*(pf0(size(P0,2))*pb2(1)+1-pb2(1));
                                g4 = calnum(2,[n1-1,n2-1],0);
                                A(g4,g3) = pb1(1)*((1-pf0(size(P0,2)))*pb2(1));
                            elseif (n2 == N2(1))&&(n1 == 0)
                                g4 = calnum(2,[n1,n2],0);
                                A(g4,g3) = 1-(1-pf0(size(P0,2)))*pb2(1);
                                g4 = calnum(2,[n1,n2-1],0);
                                A(g4,g3) = (1-pf0(size(P0,2)))*pb2(1);
                            elseif (n2 == 0)&&(n1 == 0)
                                g4 = calnum(2,[n1,n2],0);
                                A(g4,g3) = (1-pf0(size(P0,2)))*(1-pb2(1))+(1-pf0(size(P0,2)))*pb2(1);
                                g4 = calnum(2,[n1,n2+1],0);
                                A(g4,g3) = pf0(size(P0,2))*pb2(1)+pf0(size(P0,2))*(1-pb2(1));
                            elseif (n2 ~= 0)&&(n2 ~=N2(1))&&(n1 == 0)
                                g4 = calnum(2,[n1,n2],0);
                                A(g4,g3) = pf0(size(P0,2))*pb2(1)+(1-pf0(size(P0,2)))*(1-pb2(1));
                                g4 = calnum(2,[n1,n2+1],0);
                                A(g4,g3) = pf0(size(P0,2))*(1-pb2(1));
                                g4 = calnum(2,[n1,n2-1],0);
                                A(g4,g3) = (1-pf0(size(P0,2)))*pb2(1);
                            elseif (n2 == 0)&&(n1~=0)
                                g4 = calnum(2,[n1,n2],0);
                                A(g4,g3) = (1-pb1(1))*(1-pf0(size(P0,2)));
                                g4 = calnum(2,[n1,n2+1],0);
                                A(g4,g3) = (1-pb1(1))*pf0(size(P0,2));
                                g4 = calnum(2,[n1-1,n2],0);
                                A(g4,g3) = pb1(1)*(1-pf0(size(P0,2)));
                                g4 = calnum(2,[n1-1,n2+1],0);
                                A(g4,g3) = pb1(1)*pf0(size(P0,2));
                            else
                                g4 = calnum(2,[n1,n2],0);
                                A(g4,g3) = (1-pb1(1))*(pf0(size(P0,2))*pb2(1)+(1-pf0(size(P0,2)))*(1-pb2(1)));
                                g4 = calnum(2,[n1,n2+1],0);
                                A(g4,g3) = (1-pb1(1))*pf0(size(P0,2))*(1-pb2(1));
                                g4 = calnum(2,[n1,n2-1],0);
                                A(g4,g3) = (1-pb1(1))*(1-pf0(size(P0,2)))*pb2(1);
                                g4 = calnum(2,[n1-1,n2],0);
                                A(g4,g3) = pb1(1)*(pf0(size(P0,2))*pb2(1)+(1-pf0(size(P0,2)))*(1-pb2(1)));
                                g4 = calnum(2,[n1-1,n2+1],0);
                                A(g4,g3) = pb1(1)*pf0(size(P0,2))*(1-pb2(1));
                                g4 = calnum(2,[n1-1,n2-1],0);
                                A(g4,g3) = pb1(1)*(1-pf0(size(P0,2)))*pb2(1);
                            end
                        end
                    elseif (setup > 0)&&(setup < tset(1))
                        if (n1 ~= 0)&&(n2 ~= 0)
                            g2 = calnum(sta,[n1,n2],setup+1);
                            A(g2,g1) = (1-pb1(1))*(1-pb2(1));
                            g2 = calnum(sta,[n1,n2-1],setup+1);
                            A(g2,g1) = (1-pb1(1))*pb2(1);
                            g2 = calnum(sta,[n1-1,n2],setup+1);
                            A(g2,g1) = pb1(1)*(1-pb2(1));
                            g2 = calnum(sta,[n1-1,n2-1],setup+1);
                            A(g2,g1) = pb1(1)*pb2(1);
                        elseif (n1 == 0)&&(n2 ~= 0)
                            g2 = calnum(sta,[n1,n2],setup+1);
                            A(g2,g1) = 1-pb2(1);
                            g2 = calnum(sta,[n1,n2-1],setup+1);
                            A(g2,g1) = pb2(1);
                        elseif (n1 ~= 0)&&(n2 == 0)
                            g2 = calnum(sta,[n1,n2],setup+1);
                            A(g2,g1) = 1-pb1(1);
                            g2 = calnum(sta,[n1-1,n2],setup+1);
                            A(g2,g1) = pb1(1);
                        elseif (n1 == 0)&&(n2 == 0)
                            g2 = calnum(sta,[n1,n2],setup+1);
                            A(g2,g1) = 1;
                        end
                    elseif setup == tset(1)
                        if sta ==1
                            if (n1 == N1(1))&&(n2 ~= 0)
                                g2 = calnum(1,[n1,n2],0);
                                A(g2,g1) = (1-pb2(1))*(pf0(size(P0,2))*pb1(1)+1-pb1(1));
                                g2 = calnum(1,[n1-1,n2],0);
                                A(g2,g1) = (1-pb2(1))*((1-pf0(size(P0,2)))*pb1(1));
                                g2 = calnum(1,[n1,n2-1],0);
                                A(g2,g1) = pb2(1)*(pf0(size(P0,2))*pb1(1)+1-pb1(1));
                                g2 = calnum(1,[n1-1,n2-1],0);
                                A(g2,g1) = pb2(1)*((1-pf0(size(P0,2)))*pb1(1));
                            elseif (n1 == N1(1))&&(n2 == 0)
                                g2 = calnum(1,[n1,n2],0);
                                A(g2,g1) = 1-(1-pf0(size(P0,2)))*pb1(1);
                                g2 = calnum(1,[n1-1,n2],0);
                                A(g2,g1) = (1-pf0(size(P0,2)))*pb1(1);
                            elseif (n2 == 0)&&(n1 == 0)
                                g2 = calnum(1,[n1,n2],0);
                                A(g2,g1) = (1-pf0(size(P0,2)))*(1-pb1(1))+(1-pf0(size(P0,2)))*pb1(1);
                                g2 = calnum(1,[n1+1,n2],0);
                                A(g2,g1) = pf0(size(P0,2))*pb1(1)+pf0(size(P0,2))*(1-pb1(1));
                            elseif (n1 ~= 0)&&(n1 ~=N1(1))&&(n2 == 0)
                                g2 = calnum(1,[n1,n2],0);
                                A(g2,g1) = pf0(size(P0,2))*pb1(1)+(1-pf0(size(P0,2)))*(1-pb1(1));
                                g2 = calnum(1,[n1+1,n2],0);
                                A(g2,g1) = pf0(size(P0,2))*(1-pb1(1));
                                g2 = calnum(1,[n1-1,n2],0);
                                A(g2,g1) = (1-pf0(size(P0,2)))*pb1(1);
                            elseif (n1 == 0)&&(n2~=0)
                                g2 = calnum(1,[n1,n2],0);
                                A(g2,g1) = (1-pb2(1))*(1-pf0(size(P0,2)));
                                g2 = calnum(1,[n1+1,n2],0);
                                A(g2,g1) = (1-pb2(1))*pf0(size(P0,2));
                                g2 = calnum(1,[n1,n2-1],0);
                                A(g2,g1) = pb2(1)*(1-pf0(size(P0,2)));
                                g2 = calnum(1,[n1+1,n2-1],0);
                                A(g2,g1) = pb2(1)*pf0(size(P0,2));
                            else
                                g2 = calnum(1,[n1,n2],0);
                                A(g2,g1) = (1-pb2(1))*(pf0(size(P0,2))*pb1(1)+(1-pf0(size(P0,2)))*(1-pb1(1)));
                                g2 = calnum(1,[n1+1,n2],0);
                                A(g2,g1) = (1-pb2(1))*pf0(size(P0,2))*(1-pb1(1));
                                g2 = calnum(1,[n1-1,n2],0);
                                A(g2,g1) = (1-pb2(1))*(1-pf0(size(P0,2)))*pb1(1);
                                g2 = calnum(1,[n1,n2-1],0);
                                A(g2,g1) = pb2(1)*(pf0(size(P0,2))*pb1(1)+(1-pf0(size(P0,2)))*(1-pb1(1)));
                                g2 = calnum(1,[n1+1,n2-1],0);
                                A(g2,g1) = pb2(1)*pf0(size(P0,2))*(1-pb1(1));
                                g2 = calnum(1,[n1-1,n2-1],0);
                                A(g2,g1) = pb2(1)*(1-pf0(size(P0,2)))*pb1(1);
                            end
                        else
                            if (n2 == N2(1))&&(n1 ~= 0)
                                g2 = calnum(2,[n1,n2],0);
                                A(g2,g1) = (1-pb1(1))*(pf0(size(P0,2))*pb2(1)+1-pb2(1));
                                g2 = calnum(2,[n1,n2-1],0);
                                A(g2,g1) = (1-pb1(1))*((1-pf0(size(P0,2)))*pb2(1));
                                g2 = calnum(2,[n1-1,n2],0);
                                A(g2,g1) = pb1(1)*(pf0(size(P0,2))*pb2(1)+1-pb2(1));
                                g2 = calnum(2,[n1-1,n2-1],0);
                                A(g2,g1) = pb1(1)*((1-pf0(size(P0,2)))*pb2(1));
                            elseif (n2 == N2(1))&&(n1 == 0)
                                g2 = calnum(2,[n1,n2],0);
                                A(g2,g1) = 1-(1-pf0(size(P0,2)))*pb2(1);
                                g2 = calnum(2,[n1,n2-1],0);
                                A(g2,g1) = (1-pf0(size(P0,2)))*pb2(1);
                            elseif (n2 == 0)&&(n1 == 0)
                                g2 = calnum(2,[n1,n2],0);
                                A(g2,g1) = (1-pf0(size(P0,2)))*(1-pb2(1))+(1-pf0(size(P0,2)))*pb2(1);
                                g2 = calnum(2,[n1,n2+1],0);
                                A(g2,g1) = pf0(size(P0,2))*pb2(1)+pf0(size(P0,2))*(1-pb2(1));
                            elseif (n2 ~= 0)&&(n2 ~=N2(1))&&(n1 == 0)
                                g2 = calnum(2,[n1,n2],0);
                                A(g2,g1) = pf0(size(P0,2))*pb2(1)+(1-pf0(size(P0,2)))*(1-pb2(1));
                                g2 = calnum(2,[n1,n2+1],0);
                                A(g2,g1) = pf0(size(P0,2))*(1-pb2(1));
                                g2 = calnum(2,[n1,n2-1],0);
                                A(g2,g1) = (1-pf0(size(P0,2)))*pb2(1);
                            elseif (n2 == 0)&&(n1~=0)
                                g2 = calnum(2,[n1,n2],0);
                                A(g2,g1) = (1-pb1(1))*(1-pf0(size(P0,2)));
                                g2 = calnum(2,[n1,n2+1],0);
                                A(g2,g1) = (1-pb1(1))*pf0(size(P0,2));
                                g2 = calnum(2,[n1-1,n2],0);
                                A(g2,g1) = pb1(1)*(1-pf0(size(P0,2)));
                                g2 = calnum(2,[n1-1,n2+1],0);
                                A(g2,g1) = pb1(1)*pf0(size(P0,2));
                            else
                                g2 = calnum(2,[n1,n2],0);
                                A(g2,g1) = (1-pb1(1))*(pf0(size(P0,2))*pb2(1)+(1-pf0(size(P0,2)))*(1-pb2(1)));
                                g2 = calnum(2,[n1,n2+1],0);
                                A(g2,g1) = (1-pb1(1))*pf0(size(P0,2))*(1-pb2(1));
                                g2 = calnum(2,[n1,n2-1],0);
                                A(g2,g1) = (1-pb1(1))*(1-pf0(size(P0,2)))*pb2(1);
                                g2 = calnum(2,[n1-1,n2],0);
                                A(g2,g1) = pb1(1)*(pf0(size(P0,2))*pb2(1)+(1-pf0(size(P0,2)))*(1-pb2(1)));
                                g2 = calnum(2,[n1-1,n2+1],0);
                                A(g2,g1) = pb1(1)*pf0(size(P0,2))*(1-pb2(1));
                                g2 = calnum(2,[n1-1,n2-1],0);
                                A(g2,g1) = pb1(1)*(1-pf0(size(P0,2)))*pb2(1);
                            end
                        end
                    end
                end
            end
        end
    end
    
    CPR1M = zeros(1,S);
    CPR2M = zeros(1,S);
    CCRM = zeros(1,S);
    CWIP1M = zeros(1,S);
    CWIP2M = zeros(1,S);
    CST1M = zeros(1,S);
    CST2M = zeros(1,S);
    for n1 = 0:N1(1)
        for n2 = 0:N2(1)
            for sta = 1:2
                for setup = 0:tset(1)
                    g1 = calnum(sta,[n1,n2],setup);
                    if n1 ~= 0
                        CPR1M(g1) = pb1(1);
                    end
                    if n2 ~= 0
                        CPR2M(g1) = pb2(1);
                    end
                    if (sta == 1)&&(n1 ~= N1(1))&&(setup == 0)&&(n2 >= beta(1,1)*n1+gamma(1,1))
                        CCRM(g1) = pf0(size(P0,2));
                    end
                    if (sta == 2)&&(n2 ~= N2(1))&&(setup == 0)&&(n2 <= beta(1,2)*n1+gamma(1,2))
                        CCRM(g1) = pf0(size(P0,2));
                    end
                    if (sta == 1)&&(n1 ~= N1(1))&&(setup == tset(1))
                        CCRM(g1) = pf0(size(P0,2));
                    end
                    if (sta == 2)&&(n1 ~= N2(1))&&(setup == tset(1))
                        CCRM(g1) = pf0(size(P0,2));
                    end

                    if (sta == 1)&&(n1 == N1(1))&&(setup == 0)&&(n2 >= beta(1,1)*n1+gamma(1,1))
                        CCRM(g1) = pf0(size(P0,2))*pb1(1);
                    end
                    if (sta == 2)&&(n2 == N2(1))&&(setup == 0)&&(n2 <= beta(1,2)*n1+gamma(1,2))
                        CCRM(g1) = pf0(size(P0,2))*pb2(1);
                    end
                    if (sta == 1)&&(n1 == N1(1))&&(setup == tset(1))
                        CCRM(g1) = pf0(size(P0,2))*pb1(1);
                    end
                    if (sta == 2)&&(n1 == N2(1))&&(setup == tset(1))
                        CCRM(g1) = pf0(size(P0,2))*pb2(1);
                    end


                    CWIP1M(g1) = n1;
                    CWIP2M(g1) = n2;
                    if n1 == 0
                        CST1M(g1) = P1(1);
                    end
                    if n2 == 0
                        CST2M(g1) = P2(1);
                    end
                end
            end
        end
    end
    pf1(1) = CPR1M*X;
    pf2(1) = CPR2M*X;
    pb0(size(P0,2)) = CCRM*X;
    
    for j = 2:size(P1,2)
        pf1(j) = P1(j)*(1-X1(1,j));
    end
    for j = 2:size(P2,2)
        pf2(j) = P2(j)*(1-X2(1,j));
    end
    for j = size(P0,2)-1:-1:1
        pb0(j) = P0(j)*(1-(1-pb0(j+1))*X0(N0(j)+1,j));
    end

    CPR1 = zeros(1,N1(size(N1,2))+1);
    for j = 2:N1(size(N1,2))+1
        CPR1(j) = P1(size(P1,2));
    end

    CPR2 = zeros(1,N2(size(N2,2))+1);
    for j = 2:N2(size(N2,2))+1
        CPR2(j) = P2(size(P2,2));
    end
    
    CCR = zeros(1,N0(1)+1);
    for j = 1:N0(1)+1
        if j ~= (N0(1)+1)
            CCR(j) = P0(1);
        else
            CCR(j) = P0(1)*pb0(2);
        end
    end
    CWIP0 = zeros(size(N0,2),max(N0)+1);
    CWIP1 = zeros(size(N1,2),max(N1)+1);
    CWIP2 = zeros(size(N2,2),max(N2)+1);
    for j = 1:size(N0,2)
        for k = 1:max(N0)+1
            CWIP0(j,k)=k-1;
        end
    end
    for j = 2:size(N1,2)
        for k = 1:max(N1)+1
            CWIP1(j,k)=k-1;
        end
    end
    for j = 2:size(N2,2)
        for k = 1:max(N2)+1
            CWIP2(j,k)=k-1;
        end
    end
    CBL0 = zeros(size(P0,2)-1,max(N0)+1);
    CBL1 = zeros(size(P1,2)-1,max(N1)+1);
    CBL2 = zeros(size(P2,2)-1,max(N2)+1);
    for j = 1:size(P0,2)-1
        CBL0(j,N0(j)+1) = P0(j)*(1-pb0(j+1));
    end
    for j = 1:size(P1,2)-1
        CBL1(j,N1(j+1)+1) = P1(j)*(1-pb1(j+1));
    end
    for j = 1:size(P2,2)-1
        CBL2(j,N2(j+1)+1) = P2(j)*(1-pb2(j+1));
    end
    CST0 = zeros(size(P0,2),max(N0)+1);
    CST1 = zeros(size(P1,2),max(N1)+1);
    CST2 = zeros(size(P2,2),max(N2)+1);
    for j=2:size(P0,2)
        CST0(j,1) = P0(j);
    end
    for j=1:size(P1,2)
        CST1(j,1) = P1(j);
    end
    for j=1:size(P2,2)
        CST2(j,1) = P2(j);
    end

    PR(1,i) = CPR1*X1(1:N1(size(N1,2))+1,size(N1,2));
    PR(2,i) = CPR2*X2(1:N2(size(N2,2))+1,size(N2,2));
    CR(i) = CCR*X0(1:N0(1)+1,1);
    for j = 1:size(N0,2)
        WIP0(j,i) = CWIP0(j,1:N0(j)+1)*X0(1:N0(j)+1,j);
    end
    WIP1(1,i) = CWIP1M*X;
    for j = 2:size(N1,2)
        WIP1(j,i) = CWIP1(j,1:N1(j)+1)*X1(1:N1(j)+1,j);
    end
    WIP2(1,i) = CWIP2M*X;
    for j = 2:size(N2,2)
        WIP2(j,i) = CWIP2(j,1:N2(j)+1)*X2(1:N2(j)+1,j);
    end
    for j = 1:size(P0,2)-1
        BL0(j,i) = CBL0(j,1:N0(j)+1)*X0(1:N0(j)+1,j);
    end
    for j = 1:size(P1,2)-1
        BL1(j,i) = CBL1(j,1:N1(j+1)+1)*X1(1:N1(j+1)+1,j+1);
    end
    for j = 1:size(P2,2)-1
        BL2(j,i) = CBL2(j,1:N2(j+1)+1)*X2(1:N2(j+1)+1,j+1);
    end
    for j = 2:size(P0,2)
        ST0(j,i) = CST0(j,1:N0(j-1)+1)*X0(1:N0(j-1)+1,j-1);
    end
    ST1(1,i) = CST1M*X;
    for j = 2:size(P1,2)
        ST1(j,i) = CST1(j,1:N1(j)+1)*X1(1:N1(j)+1,j);
    end
    ST2(1,i) = CST2M*X;
    for j = 2:size(P2,2)
        ST2(j,i) = CST2(j,1:N2(j)+1)*X2(1:N2(j)+1,j);
    end

    for k = 1:size(N0,2)
        A2 = zeros(N0(k)+1,N0(k)+1);
        for g=1:N0(k)+1
            for j=1:N0(k)+1
                if g==j
                    if g==1
                        A2(g,j)=1-pf0(k);
                    elseif g~=N0(k)+1
                        A2(g,j)=1-pf0(k)-pb0(k+1)+2*pf0(k)*pb0(k+1);
                    else
                        A2(g,j)=pf0(k)*pb0(k+1)+1-pb0(k+1);
                    end
                end
                if g==j+1
                    if g==2
                        A2(g,j)=pf0(k);
                    else
                        A2(g,j)=pf0(k)*(1-pb0(k+1));
                    end
                end
                if j==g+1
                    A2(g,j)=pb0(k+1)*(1-pf0(k));
                end
            end
        end
        X0(1:N0(k)+1,k) = A2*X0(1:N0(k)+1,k);
    end
    
    for k = 2:size(N1,2)
        A2 = zeros(N1(k)+1,N1(k)+1);
        for g=1:N1(k)+1
            for j=1:N1(k)+1
                if g==j
                    if g==1
                        A2(g,j)=1-pf1(k-1);
                    elseif g~=N1(k)+1
                        A2(g,j)=1-pf1(k-1)-pb1(k)+2*pf1(k-1)*pb1(k);
                    else
                        A2(g,j)=pf1(k-1)*pb1(k)+1-pb1(k);
                    end
                end
                if g==j+1
                    if g==2
                        A2(g,j)=pf1(k-1);
                    else
                        A2(g,j)=pf1(k-1)*(1-pb1(k));
                    end
                end
                if j==g+1
                    A2(g,j)=pb1(k)*(1-pf1(k-1));
                end
            end
        end
        X1(1:N1(k)+1,k) = A2*X1(1:N1(k)+1,k);
    end
    
    for k = 2:size(N2,2)
        A2 = zeros(N2(k)+1,N2(k)+1);
        for g=1:N2(k)+1
            for j=1:N2(k)+1
                if g==j
                    if g==1
                        A2(g,j)=1-pf2(k-1);
                    elseif g~=N2(k)+1
                        A2(g,j)=1-pf2(k-1)-pb2(k)+2*pf2(k-1)*pb2(k);
                    else
                        A2(g,j)=pf2(k-1)*pb2(k)+1-pb2(k);
                    end
                end
                if g==j+1
                    if g==2
                        A2(g,j)=pf2(k-1);
                    else
                        A2(g,j)=pf2(k-1)*(1-pb2(k));
                    end
                end
                if j==g+1
                    A2(g,j)=pb2(k)*(1-pf2(k-1));
                end
            end
        end
        X2(1:N2(k)+1,k) = A2*X2(1:N2(k)+1,k);
    end
    X = A*X;
end

i = 0:t;

plot(i,PR(1,:),'--','color','#0072BD','LineWidth',2);hold on;
plot(i,PR(2,:),'-.','color','#0072BD','LineWidth',2);hold on;
grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'PR$$_i(n)$$'},'Interpreter','latex');
set(gca,'fontsize',25);
axis([0 300 0 0.65]);

% plot(i,CR,'LineWidth',2);hold on;
% grid on;
% xlabel({'$$n$$'},'Interpreter','latex');
% ylabel({'CR$$(n)$$'},'Interpreter','latex');
% set(gca,'fontsize',25);
% axis([0 300 0.75 0.95]);

% plot(i,WIP0(1,:),'-','color','#0072BD','LineWidth',2);hold on;
% plot(i,WIP0(2,:),'-','color','#D95319','LineWidth',2);hold on;
% plot(i,WIP1(1,:),'--','color','#0072BD','LineWidth',2);hold on;
% plot(i,WIP1(2,:),'--','color','#D95319','LineWidth',2);hold on;
% plot(i,WIP1(3,:),'--','color','#EDB120','LineWidth',2);hold on;
% plot(i,WIP1(4,:),'--','color','#7E2F8E','LineWidth',2);hold on;
% plot(i,WIP2(1,:),'-.','color','#0072BD','LineWidth',2);hold on;
% plot(i,WIP2(2,:),'-.','color','#D95319','LineWidth',2);hold on;
% plot(i,WIP2(3,:),'-.','color','#EDB120','LineWidth',2);hold on;
% plot(i,WIP2(4,:),'-.','color','#7E2F8E','LineWidth',2);hold on;
% plot(i,WIP2(5,:),'-.','color','#77AC30','LineWidth',2);hold on;
% grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'WIP$$_{i,j}(n)$$'},'Interpreter','latex');
% set(gca,'fontsize',25);
% axis([0 300 0 8]);

% plot(i,BL0(1,:),'-','color','#0072BD','LineWidth',2);hold on;
% plot(i,BL0(2,:),'-','color','#D95319','LineWidth',2);hold on;
% plot(i,BL1(1,:),'--','color','#D95319','LineWidth',2);hold on;
% plot(i,BL1(2,:),'--','color','#D95319','LineWidth',2);hold on;
% plot(i,BL1(3,:),'--','color','#EDB120','LineWidth',2);hold on;
% plot(i,BL2(1,:),'-.','color','#0072BD','LineWidth',2);hold on;
% plot(i,BL2(2,:),'-.','color','#D95319','LineWidth',2);hold on;
% plot(i,BL2(3,:),'-.','color','#EDB120','LineWidth',2);hold on;
% plot(i,BL2(4,:),'-.','color','#7E2F8E','LineWidth',2);hold on;
% grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'BL$$_{i,j}(n)$$'},'Interpreter','latex');
% set(gca,'fontsize',25);
% axis([0 300 0 0.17]);

% plot(i,ST0(2,:),'-','color','#D95319','LineWidth',2);hold on;
% plot(i,ST0(3,:),'-','color','#EDB120','LineWidth',2);hold on;
% plot(i,ST1(1,:),'--','color','#D95319','LineWidth',2);hold on;
% plot(i,ST1(2,:),'--','color','#D95319','LineWidth',2);hold on;
% plot(i,ST1(3,:),'--','color','#EDB120','LineWidth',2);hold on;
% plot(i,ST1(4,:),'--','color','#7E2F8E','LineWidth',2);hold on;
% plot(i,ST2(1,:),'-.','color','#0072BD','LineWidth',2);
% plot(i,ST2(2,:),'-.','color','#D95319','LineWidth',2);
% plot(i,ST2(3,:),'-.','color','#EDB120','LineWidth',2);
% plot(i,ST2(4,:),'-.','color','#7E2F8E','LineWidth',2);
% plot(i,ST2(5,:),'-.','color','#4DBEEE','LineWidth',2);
% grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'ST$$_{i,j}(n)$$'},'Interpreter','latex');
% set(gca,'fontsize',25);
% axis([0 300 0 1]);

function num = calnum(sta,b,u)
    global N2 S tset
    num = ((b(1)*(N2(1)+1)+b(2))*(tset(1)+1)+u+1)+(sta-1)*(S/2);
end