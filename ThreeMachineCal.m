clear
close all

global N2 S tset
P0 = 0.98;
P1 = 0.75;
P2 = 0.84;
N0 = [];
N1 = 7;
N2 = 8;
ha = [1/2 2];
hb = [-2 3];
tset = 2;
t = 200;

PR = zeros(2,t+1);
CR = zeros(1,t+1);
WIP1 = zeros(1,t+1);
WIP2 = zeros(1,t+1);
ST1 = zeros(1,t+1);
ST2 = zeros(1,t+1);
BL0 = zeros(1,t+1);

S = (N1(1)+1)*(N2(1)+1)*(tset(1)+1)*2;%total state number

A = zeros(S,S);
for n1 = 0:N1(1)
    for n2 = 0:N2(1)
        for sta = 1:2
            for setup = 0:tset(1)
                g1 = calnum(sta,[n1,n2],setup);
                if (sta == 1)&&(setup == 0)
                    if n2 < ha(1,1)*n1+hb(1,1)
                        if (n1 ~= 0)&&(n2 ~= 0)
                            g2 = calnum(2,[n1,n2],1);
                            A(g2,g1) = (1-P1(1))*(1-P2(1));
                            g2 = calnum(2,[n1,n2-1],1);
                            A(g2,g1) = (1-P1(1))*P2(1);
                            g2 = calnum(2,[n1-1,n2],1);
                            A(g2,g1) = P1(1)*(1-P2(1));
                            g2 = calnum(2,[n1-1,n2-1],1);
                            A(g2,g1) = P1(1)*P2(1);
                        elseif (n1 == 0)&&(n2 ~= 0)
                            g2 = calnum(2,[n1,n2],1);
                            A(g2,g1) = 1-P2(1);
                            g2 = calnum(2,[n1,n2-1],1);
                            A(g2,g1) = P2(1);
                        elseif (n1 ~= 0)&&(n2 == 0)
                            g2 = calnum(2,[n1,n2],1);
                            A(g2,g1) = 1-P1(1);
                            g2 = calnum(2,[n1-1,n2],1);
                            A(g2,g1) = P1(1);
                        elseif (n1 == 0)&&(n2 == 0)
                            g2 = calnum(2,[n1,n2],1);
                            A(g2,g1) = 1;
                        end
                    else
                        g3 = calnum(1,[n1,n2],0);
                        if (n1 == N1(1))&&(n2 ~= 0)
                            g4 = calnum(1,[n1,n2],0);
                            A(g4,g3) = (1-P2(1))*(P0(size(P0,2))*P1(1)+1-P1(1));
                            g4 = calnum(1,[n1-1,n2],0);
                            A(g4,g3) = (1-P2(1))*((1-P0(size(P0,2)))*P1(1));
                            g4 = calnum(1,[n1,n2-1],0);
                            A(g4,g3) = P2(1)*(P0(size(P0,2))*P1(1)+1-P1(1));
                            g4 = calnum(1,[n1-1,n2-1],0);
                            A(g4,g3) = P2(1)*((1-P0(size(P0,2)))*P1(1));
                        elseif (n1 == N1(1))&&(n2 == 0)
                            g4 = calnum(1,[n1,n2],0);
                            A(g4,g3) = 1-(1-P0(size(P0,2)))*P1(1);
                            g4 = calnum(1,[n1-1,n2],0);
                            A(g4,g3) = (1-P0(size(P0,2)))*P1(1);
                        elseif (n2 == 0)&&(n1 == 0)
                            g4 = calnum(1,[n1,n2],0);
                            A(g4,g3) = (1-P0(size(P0,2)))*(1-P1(1))+(1-P0(size(P0,2)))*P1(1);
                            g4 = calnum(1,[n1+1,n2],0);
                            A(g4,g3) = P0(size(P0,2))*P1(1)+P0(size(P0,2))*(1-P1(1));
                        elseif (n1 ~= 0)&&(n1 ~=N1(1))&&(n2 == 0)
                            g4 = calnum(1,[n1,n2],0);
                            A(g4,g3) = P0(size(P0,2))*P1(1)+(1-P0(size(P0,2)))*(1-P1(1));
                            g4 = calnum(1,[n1+1,n2],0);
                            A(g4,g3) = P0(size(P0,2))*(1-P1(1));
                            g4 = calnum(1,[n1-1,n2],0);
                            A(g4,g3) = (1-P0(size(P0,2)))*P1(1);
                        elseif (n1 == 0)&&(n2~=0)
                            g4 = calnum(1,[n1,n2],0);
                            A(g4,g3) = (1-P2(1))*(1-P0(size(P0,2)));
                            g4 = calnum(1,[n1+1,n2],0);
                            A(g4,g3) = (1-P2(1))*P0(size(P0,2));
                            g4 = calnum(1,[n1,n2-1],0);
                            A(g4,g3) = P2(1)*(1-P0(size(P0,2)));
                            g4 = calnum(1,[n1+1,n2-1],0);
                            A(g4,g3) = P2(1)*P0(size(P0,2));
                        else
                            g4 = calnum(1,[n1,n2],0);
                            A(g4,g3) = (1-P2(1))*(P0(size(P0,2))*P1(1)+(1-P0(size(P0,2)))*(1-P1(1)));
                            g4 = calnum(1,[n1+1,n2],0);
                            A(g4,g3) = (1-P2(1))*P0(size(P0,2))*(1-P1(1));
                            g4 = calnum(1,[n1-1,n2],0);
                            A(g4,g3) = (1-P2(1))*(1-P0(size(P0,2)))*P1(1);
                            g4 = calnum(1,[n1,n2-1],0);
                            A(g4,g3) = P2(1)*(P0(size(P0,2))*P1(1)+(1-P0(size(P0,2)))*(1-P1(1)));
                            g4 = calnum(1,[n1+1,n2-1],0);
                            A(g4,g3) = P2(1)*P0(size(P0,2))*(1-P1(1));
                            g4 = calnum(1,[n1-1,n2-1],0);
                            A(g4,g3) = P2(1)*(1-P0(size(P0,2)))*P1(1);
                        end
                    end
                elseif (sta == 2)&&(setup == 0)
                    if n2 > ha(1,2)*n1+hb(1,2)
                        if (n1 ~= 0)&&(n2 ~= 0)
                            g2 = calnum(1,[n1,n2],1);
                            A(g2,g1) = (1-P1(1))*(1-P2(1));
                            g2 = calnum(1,[n1,n2-1],1);
                            A(g2,g1) = (1-P1(1))*P2(1);
                            g2 = calnum(1,[n1-1,n2],1);
                            A(g2,g1) = P1(1)*(1-P2(1));
                            g2 = calnum(1,[n1-1,n2-1],1);
                            A(g2,g1) = P1(1)*P2(1);
                        elseif (n1 == 0)&&(n2 ~= 0)
                            g2 = calnum(1,[n1,n2],1);
                            A(g2,g1) = 1-P2(1);
                            g2 = calnum(1,[n1,n2-1],1);
                            A(g2,g1) = P2(1);
                        elseif (n1 ~= 0)&&(n2 == 0)
                            g2 = calnum(1,[n1,n2],1);
                            A(g2,g1) = 1-P1(1);
                            g2 = calnum(1,[n1-1,n2],1);
                            A(g2,g1) = P1(1);
                        elseif (n1 == 0)&&(n2 == 0)
                            g2 = calnum(1,[n1,n2],1);
                            A(g2,g1) = 1;
                        end
                    else
                        g3 = calnum(2,[n1,n2],0);
                        if (n2 == N2(1))&&(n1 ~= 0)
                            g4 = calnum(2,[n1,n2],0);
                            A(g4,g3) = (1-P1(1))*(P0(size(P0,2))*P2(1)+1-P2(1));
                            g4 = calnum(2,[n1,n2-1],0);
                            A(g4,g3) = (1-P1(1))*((1-P0(size(P0,2)))*P2(1));
                            g4 = calnum(2,[n1-1,n2],0);
                            A(g4,g3) = P1(1)*(P0(size(P0,2))*P2(1)+1-P2(1));
                            g4 = calnum(2,[n1-1,n2-1],0);
                            A(g4,g3) = P1(1)*((1-P0(size(P0,2)))*P2(1));
                        elseif (n2 == N2(1))&&(n1 == 0)
                            g4 = calnum(2,[n1,n2],0);
                            A(g4,g3) = 1-(1-P0(size(P0,2)))*P2(1);
                            g4 = calnum(2,[n1,n2-1],0);
                            A(g4,g3) = (1-P0(size(P0,2)))*P2(1);
                        elseif (n2 == 0)&&(n1 == 0)
                            g4 = calnum(2,[n1,n2],0);
                            A(g4,g3) = (1-P0(size(P0,2)))*(1-P2(1))+(1-P0(size(P0,2)))*P2(1);
                            g4 = calnum(2,[n1,n2+1],0);
                            A(g4,g3) = P0(size(P0,2))*P2(1)+P0(size(P0,2))*(1-P2(1));
                        elseif (n2 ~= 0)&&(n2 ~=N2(1))&&(n1 == 0)
                            g4 = calnum(2,[n1,n2],0);
                            A(g4,g3) = P0(size(P0,2))*P2(1)+(1-P0(size(P0,2)))*(1-P2(1));
                            g4 = calnum(2,[n1,n2+1],0);
                            A(g4,g3) = P0(size(P0,2))*(1-P2(1));
                            g4 = calnum(2,[n1,n2-1],0);
                            A(g4,g3) = (1-P0(size(P0,2)))*P2(1);
                        elseif (n2 == 0)&&(n1~=0)
                            g4 = calnum(2,[n1,n2],0);
                            A(g4,g3) = (1-P1(1))*(1-P0(size(P0,2)));
                            g4 = calnum(2,[n1,n2+1],0);
                            A(g4,g3) = (1-P1(1))*P0(size(P0,2));
                            g4 = calnum(2,[n1-1,n2],0);
                            A(g4,g3) = P1(1)*(1-P0(size(P0,2)));
                            g4 = calnum(2,[n1-1,n2+1],0);
                            A(g4,g3) = P1(1)*P0(size(P0,2));
                        else
                            g4 = calnum(2,[n1,n2],0);
                            A(g4,g3) = (1-P1(1))*(P0(size(P0,2))*P2(1)+(1-P0(size(P0,2)))*(1-P2(1)));
                            g4 = calnum(2,[n1,n2+1],0);
                            A(g4,g3) = (1-P1(1))*P0(size(P0,2))*(1-P2(1));
                            g4 = calnum(2,[n1,n2-1],0);
                            A(g4,g3) = (1-P1(1))*(1-P0(size(P0,2)))*P2(1);
                            g4 = calnum(2,[n1-1,n2],0);
                            A(g4,g3) = P1(1)*(P0(size(P0,2))*P2(1)+(1-P0(size(P0,2)))*(1-P2(1)));
                            g4 = calnum(2,[n1-1,n2+1],0);
                            A(g4,g3) = P1(1)*P0(size(P0,2))*(1-P2(1));
                            g4 = calnum(2,[n1-1,n2-1],0);
                            A(g4,g3) = P1(1)*(1-P0(size(P0,2)))*P2(1);
                        end
                    end
                elseif (setup > 0)&&(setup < tset(1))
                    if (n1 ~= 0)&&(n2 ~= 0)
                        g2 = calnum(sta,[n1,n2],setup+1);
                        A(g2,g1) = (1-P1(1))*(1-P2(1));
                        g2 = calnum(sta,[n1,n2-1],setup+1);
                        A(g2,g1) = (1-P1(1))*P2(1);
                        g2 = calnum(sta,[n1-1,n2],setup+1);
                        A(g2,g1) = P1(1)*(1-P2(1));
                        g2 = calnum(sta,[n1-1,n2-1],setup+1);
                        A(g2,g1) = P1(1)*P2(1);
                    elseif (n1 == 0)&&(n2 ~= 0)
                        g2 = calnum(sta,[n1,n2],setup+1);
                        A(g2,g1) = 1-P2(1);
                        g2 = calnum(sta,[n1,n2-1],setup+1);
                        A(g2,g1) = P2(1);
                    elseif (n1 ~= 0)&&(n2 == 0)
                        g2 = calnum(sta,[n1,n2],setup+1);
                        A(g2,g1) = 1-P1(1);
                        g2 = calnum(sta,[n1-1,n2],setup+1);
                        A(g2,g1) = P1(1);
                    elseif (n1 == 0)&&(n2 == 0)
                        g2 = calnum(sta,[n1,n2],setup+1);
                        A(g2,g1) = 1;
                    end
                elseif setup == tset(1)
                    if sta ==1
                        if (n1 == N1(1))&&(n2 ~= 0)
                            g2 = calnum(1,[n1,n2],0);
                            A(g2,g1) = (1-P2(1))*(P0(size(P0,2))*P1(1)+1-P1(1));
                            g2 = calnum(1,[n1-1,n2],0);
                            A(g2,g1) = (1-P2(1))*((1-P0(size(P0,2)))*P1(1));
                            g2 = calnum(1,[n1,n2-1],0);
                            A(g2,g1) = P2(1)*(P0(size(P0,2))*P1(1)+1-P1(1));
                            g2 = calnum(1,[n1-1,n2-1],0);
                            A(g2,g1) = P2(1)*((1-P0(size(P0,2)))*P1(1));
                        elseif (n1 == N1(1))&&(n2 == 0)
                            g2 = calnum(1,[n1,n2],0);
                            A(g2,g1) = 1-(1-P0(size(P0,2)))*P1(1);
                            g2 = calnum(1,[n1-1,n2],0);
                            A(g2,g1) = (1-P0(size(P0,2)))*P1(1);
                        elseif (n2 == 0)&&(n1 == 0)
                            g2 = calnum(1,[n1,n2],0);
                            A(g2,g1) = (1-P0(size(P0,2)))*(1-P1(1))+(1-P0(size(P0,2)))*P1(1);
                            g2 = calnum(1,[n1+1,n2],0);
                            A(g2,g1) = P0(size(P0,2))*P1(1)+P0(size(P0,2))*(1-P1(1));
                        elseif (n1 ~= 0)&&(n1 ~=N1(1))&&(n2 == 0)
                            g2 = calnum(1,[n1,n2],0);
                            A(g2,g1) = P0(size(P0,2))*P1(1)+(1-P0(size(P0,2)))*(1-P1(1));
                            g2 = calnum(1,[n1+1,n2],0);
                            A(g2,g1) = P0(size(P0,2))*(1-P1(1));
                            g2 = calnum(1,[n1-1,n2],0);
                            A(g2,g1) = (1-P0(size(P0,2)))*P1(1);
                        elseif (n1 == 0)&&(n2~=0)
                            g2 = calnum(1,[n1,n2],0);
                            A(g2,g1) = (1-P2(1))*(1-P0(size(P0,2)));
                            g2 = calnum(1,[n1+1,n2],0);
                            A(g2,g1) = (1-P2(1))*P0(size(P0,2));
                            g2 = calnum(1,[n1,n2-1],0);
                            A(g2,g1) = P2(1)*(1-P0(size(P0,2)));
                            g2 = calnum(1,[n1+1,n2-1],0);
                            A(g2,g1) = P2(1)*P0(size(P0,2));
                        else
                            g2 = calnum(1,[n1,n2],0);
                            A(g2,g1) = (1-P2(1))*(P0(size(P0,2))*P1(1)+(1-P0(size(P0,2)))*(1-P1(1)));
                            g2 = calnum(1,[n1+1,n2],0);
                            A(g2,g1) = (1-P2(1))*P0(size(P0,2))*(1-P1(1));
                            g2 = calnum(1,[n1-1,n2],0);
                            A(g2,g1) = (1-P2(1))*(1-P0(size(P0,2)))*P1(1);
                            g2 = calnum(1,[n1,n2-1],0);
                            A(g2,g1) = P2(1)*(P0(size(P0,2))*P1(1)+(1-P0(size(P0,2)))*(1-P1(1)));
                            g2 = calnum(1,[n1+1,n2-1],0);
                            A(g2,g1) = P2(1)*P0(size(P0,2))*(1-P1(1));
                            g2 = calnum(1,[n1-1,n2-1],0);
                            A(g2,g1) = P2(1)*(1-P0(size(P0,2)))*P1(1);
                        end
                    else
                        if (n2 == N2(1))&&(n1 ~= 0)
                            g2 = calnum(2,[n1,n2],0);
                            A(g2,g1) = (1-P1(1))*(P0(size(P0,2))*P2(1)+1-P2(1));
                            g2 = calnum(2,[n1,n2-1],0);
                            A(g2,g1) = (1-P1(1))*((1-P0(size(P0,2)))*P2(1));
                            g2 = calnum(2,[n1-1,n2],0);
                            A(g2,g1) = P1(1)*(P0(size(P0,2))*P2(1)+1-P2(1));
                            g2 = calnum(2,[n1-1,n2-1],0);
                            A(g2,g1) = P1(1)*((1-P0(size(P0,2)))*P2(1));
                        elseif (n2 == N2(1))&&(n1 == 0)
                            g2 = calnum(2,[n1,n2],0);
                            A(g2,g1) = 1-(1-P0(size(P0,2)))*P2(1);
                            g2 = calnum(2,[n1,n2-1],0);
                            A(g2,g1) = (1-P0(size(P0,2)))*P2(1);
                        elseif (n2 == 0)&&(n1 == 0)
                            g2 = calnum(2,[n1,n2],0);
                            A(g2,g1) = (1-P0(size(P0,2)))*(1-P2(1))+(1-P0(size(P0,2)))*P2(1);
                            g2 = calnum(2,[n1,n2+1],0);
                            A(g2,g1) = P0(size(P0,2))*P2(1)+P0(size(P0,2))*(1-P2(1));
                        elseif (n2 ~= 0)&&(n2 ~=N2(1))&&(n1 == 0)
                            g2 = calnum(2,[n1,n2],0);
                            A(g2,g1) = P0(size(P0,2))*P2(1)+(1-P0(size(P0,2)))*(1-P2(1));
                            g2 = calnum(2,[n1,n2+1],0);
                            A(g2,g1) = P0(size(P0,2))*(1-P2(1));
                            g2 = calnum(2,[n1,n2-1],0);
                            A(g2,g1) = (1-P0(size(P0,2)))*P2(1);
                        elseif (n2 == 0)&&(n1~=0)
                            g2 = calnum(2,[n1,n2],0);
                            A(g2,g1) = (1-P1(1))*(1-P0(size(P0,2)));
                            g2 = calnum(2,[n1,n2+1],0);
                            A(g2,g1) = (1-P1(1))*P0(size(P0,2));
                            g2 = calnum(2,[n1-1,n2],0);
                            A(g2,g1) = P1(1)*(1-P0(size(P0,2)));
                            g2 = calnum(2,[n1-1,n2+1],0);
                            A(g2,g1) = P1(1)*P0(size(P0,2));
                        else
                            g2 = calnum(2,[n1,n2],0);
                            A(g2,g1) = (1-P1(1))*(P0(size(P0,2))*P2(1)+(1-P0(size(P0,2)))*(1-P2(1)));
                            g2 = calnum(2,[n1,n2+1],0);
                            A(g2,g1) = (1-P1(1))*P0(size(P0,2))*(1-P2(1));
                            g2 = calnum(2,[n1,n2-1],0);
                            A(g2,g1) = (1-P1(1))*(1-P0(size(P0,2)))*P2(1);
                            g2 = calnum(2,[n1-1,n2],0);
                            A(g2,g1) = P1(1)*(P0(size(P0,2))*P2(1)+(1-P0(size(P0,2)))*(1-P2(1)));
                            g2 = calnum(2,[n1-1,n2+1],0);
                            A(g2,g1) = P1(1)*P0(size(P0,2))*(1-P2(1));
                            g2 = calnum(2,[n1-1,n2-1],0);
                            A(g2,g1) = P1(1)*(1-P0(size(P0,2)))*P2(1);
                        end
                    end
                end
            end
        end
    end
end

X = zeros(S,1);
X(1,1) = 1;
CPR1 = zeros(1,S);
CPR2 = zeros(1,S);
CCR = zeros(1,S);
CWIP = zeros(2,S);
CST = zeros(2,S);
CBL = zeros(1,S);
for n1 = 0:N1(1)
    for n2 = 0:N2(1)
        for sta = 1:2
            for setup = 0:tset(1)
                g1 = calnum(sta,[n1,n2],setup);
                if n1 ~= 0
                    CPR1(g1) = P1(1);
                end
                if n2 ~= 0
                    CPR2(g1) = P2(1);
                end
                if (sta == 1)&&(n1 ~= N1(1))&&(setup == 0)&&(n2 >= ha(1,1)*n1+hb(1,1))
                    CCR(g1) = P0(size(P0,2));
                end
                if (sta == 2)&&(n2 ~= N2(1))&&(setup == 0)&&(n2 <= ha(1,2)*n1+hb(1,2))
                    CCR(g1) = P0(size(P0,2));
                end
                if (sta == 1)&&(n1 ~= N1(1))&&(setup == tset(1))
                    CCR(g1) = P0(size(P0,2));
                end
                if (sta == 2)&&(n1 ~= N2(1))&&(setup == tset(1))
                    CCR(g1) = P0(size(P0,2));
                end

                if (sta == 1)&&(n1 == N1(1))&&(setup == 0)&&(n2 >= ha(1,1)*n1+hb(1,1))
                    CCR(g1) = P0(size(P0,2))*P1(1);
                end
                if (sta == 2)&&(n2 == N2(1))&&(setup == 0)&&(n2 <= ha(1,2)*n1+hb(1,2))
                    CCR(g1) = P0(size(P0,2))*P2(1);
                end
                if (sta == 1)&&(n1 == N1(1))&&(setup == tset(1))
                    CCR(g1) = P0(size(P0,2))*P1(1);
                end
                if (sta == 2)&&(n1 == N2(1))&&(setup == tset(1))
                    CCR(g1) = P0(size(P0,2))*P2(1);
                end

                if n1 == 0
                    CST(1,g1) = P1(1);
                end
                if n2 == 0
                    CST(2,g1) = P2(1);
                end
                if (sta == 1)&&(n1 == N1(1))&&(setup == 0)&&(n2 >= ha(1,1)*n1+hb(1,1))
                    CBL(g1) = P0(size(P0,2))*(1-P1(1));
                end
                if (sta == 1)&&(n1 == N1(1))&&(setup == tset(1))
                    CBL(g1) = P0(size(P0,2))*(1-P1(1));
                end
                if (sta == 2)&&(n2 == N2(1))&&(setup == 0)&&(n2 <= ha(1,2)*n1+hb(1,2))
                    CBL(g1) = P0(size(P0,2))*(1-P2(1));
                end
                if (sta == 2)&&(n2 == N2(1))&&(setup == tset(1))
                    CBL(g1) = P0(size(P0,2))*(1-P2(1));
                end
                CWIP(1,g1) = n1;
                CWIP(2,g1) = n2;
            end
        end
    end
end

for i = 1:t+1
    PR(1,i) = CPR1*X;
    PR(2,i) = CPR2*X;
    CR(i) = CCR*X;
    WIP1(i) = CWIP(1,:)*X;
    WIP2(i) = CWIP(2,:)*X;
    ST1(i) = CST(1,:)*X;
    ST2(i) = CST(2,:)*X;
    BL0(i) = CBL*X;
    X=A*X;
end

i = 0:t;
plot(i,PR(1,:),'--','LineWidth',3);hold on;
plot(i,PR(2,:),'--','LineWidth',3);hold on;
grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'PR$$_i(n)$$'},'Interpreter','latex');
legend({'PR$$_1(n)$$','PR$$_2(n)$$'},'Interpreter','latex');
set(gca,'fontsize',25);
axis([0 200 0 0.75]);

% plot(i,CR,'--','LineWidth',3);hold on;
% grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'CR$$(n)$$'},'Interpreter','latex');
% legend({'CR$$(n)$$'},'Interpreter','latex');
% set(gca,'fontsize',25);
% axis([0 200 0.8 1]);

% plot(i,WIP1,'--','LineWidth',3);hold on;
% plot(i,WIP2,'--','LineWidth',3);hold on;
% grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'WIP$$_{i,j}(n)$$'},'Interpreter','latex');
% legend({'WIP$$_{1,0}(n)$$','WIP$$_{2,0}(n)$$'},'Interpreter','latex');
% set(gca,'fontsize',25);
% axis([0 200 0 3.5]);

% plot(i,ST1,'--','LineWidth',3);hold on;
% plot(i,ST2,'--','LineWidth',3);hold on;
% grid on;xlabel({'$$n$$'},'Interpreter','latex');ylabel({'ST$$_{i,j}(n)$$'},'Interpreter','latex');
% legend({'ST$$_{1,1}(n)$$','ST$$_{2,1}(n)$$'},'Interpreter','latex');
% set(gca,'fontsize',25);
% axis([0 200 0 0.9]);

function num = calnum(sta,b,u)
    global N2 S tset
    num = ((b(1)*(N2(1)+1)+b(2))*(tset(1)+1)+u+1)+(sta-1)*(S/2);
end