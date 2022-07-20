  

%This code implements the generalized sensing function 
% presented in the accepted paper at the ACM\IEEE IPSN2022 

% The title of the paper is:"Eliminating Design Effort: A Reconfigurable 
% Sensing Framework For Chipless, Backscatter Tag"

close all
clear
clc
%A set of equivalent circuit parameters is derived from 
% the IoT devices and sensed objects present in the scene. 
% In the code below, we sense soil with 5% moisture(dielectric constant of 5) 
% using a Wi-Fi signal(4.9-5.8GHz) .

ep=5;  % Dielectric constant of soil with 5% moisture
a=0.2; %Parameters of submodel A
b=7;   %Parameters of submodel A
c=4.5; %Parameters of submodel A

i=1;j=1;

for Lc = 0.56:0.01:13.7                   %Search range for circuit parameters
    for c1 = 0.1:0.1:1.73                 %Search range for circuit parameters
        for cc = 0.001:0.001:0.05         %Search range for circuit parameters
            %submodel A : Relationship between user requirements (sensing target type and frequency band) and label equivalent circuit
            f(1,i) = 1/(2*pi*sqrt(Lc*1e-9*(0.9+c/ep)* ...
                     (c1*1e-12*ep*a+cc*1e-12*ep+8.854187817e-12*ep/1000*b)));   
            if f(1,i)<5.8*10^9&&f(1,i)>4.9*10^9  %wifi frequency band
                Lc_valid(1,j) = Lc;              %Collection of legal inductors
                c1_valid(1,j) = c1;              %Collection of legal capacitors
                cc_valid(1,j) = cc;              %Collection of legal capacitors
                j = j+1;
            end
            i = i+1;
        end
    end
end
%% submodel B : Relationship between equivalent circuit parameters and DGS size
% get [a,w] through Lc_fitting relation

%  the relationship bwteen [a,w] and equivlent inductor
%  Coefficients (with 95% confidence bounds): Lc
       p00_Lc =       112.2  ;%(15.16, 209.3)
       p10_Lc =       -71.31 ;%  (-131.3, -11.25)
       p01_Lc  =      -3.273 ;%   (-72.86, 66.32)
       p20_Lc  =       19.42 ;% (4.771, 34.07)
       p11_Lc  =       -22.2 ;% (-51.14, 6.747)
       p02_Lc  =       89.87 ;%  (22.87, 156.9)
       p30_Lc  =      -2.636 ;%  (-4.399, -0.8731)
       p21_Lc  =       5.895 ;%  (1.175, 10.62)
       p12_Lc  =      -12.63 ;%  (-27.34, 2.081)
       p03_Lc  =      -59.99 ;% (-119.2, -0.8136)
       p40_Lc  =      0.1755 ;%  (0.07091, 0.2802)
       p31_Lc  =     -0.5732 ;%   (-0.9251, -0.2214)
       p22_Lc  =       1.322 ;%  (-0.02941, 2.674)
       p13_Lc  =     -0.3111 ;%  (-6.058, 5.436)
       p04_Lc  =        34.9 ;%  (2.31, 67.49)
       p50_Lc  =   -0.004509 ;% (-0.006962, -0.002057)
       p41_Lc  =      0.0167 ;%  (0.006618, 0.02679)
       p32_Lc  =   -0.002057 ;%  (-0.04991, 0.0458)
       p23_Lc  =     -0.4036 ;%   (-0.6429, -0.1643)
       p14_Lc  =       1.946 ;%  (0.6852, 3.206)
       p05_Lc  =      -10.89 ;%  (-18.55, -3.221)

%  the relationship bwteen [a,w] and equivlent capacitor
       p00_c1 =      -3.479 ;%  (-26.3, 19.34)
       p10_c1 =       4.199 ;%  (-8.918, 17.31)
       p01_c1 =      -13.85 ;%  (-25.23, -2.458)
       p20_c1 =      -1.563 ;%  (-4.55, 1.425)
       p11_c1 =        8.05 ;% (3.349, 12.75)
       p02_c1 =       -5.61 ;% (-13.2, 1.979)
       p30_c1 =      0.2535 ;% (-0.08354, 0.5905)
       p21_c1 =      -1.495 ;% (-2.24, -0.7503)
       p12_c1 =      0.7597 ;% (-1.075, 2.594)
       p03_c1 =       3.628 ;% (-1.68, 8.936)
       p40_c1 =    -0.01859 ;% (-0.03744, 0.0002486)
       p31_c1 =      0.1219 ;%(0.06855, 0.1753)
       p22_c1 =     -0.1387 ;% (-0.31, 0.0326)
       p13_c1 =      0.4517 ;% (-0.1394, 1.043)
       p04_c1 =      -3.186 ;% (-5.85, -0.5225)
       p50_c1 =   0.0005053 ;% (8.759e-05, 0.000923)
       p41_c1 =   -0.003559 ;% (-0.005016, -0.002102)
       p32_c1 =    0.005382 ;% (-0.0005351, 0.0113)
       p23_c1 =    -0.01164 ;% (-0.0371, 0.01382)
       p14_c1 =    -0.05703 ;% (-0.1729, 0.05889)
       p05_c1 =      0.7871 ;% (0.1768, 1.397)

m=1;n=1;

for x=5.04:0.05:15                     %Search range for DGS size parameters
    for y=0.1:0.02:1.5                 %Search range for DGS size parameters
        % Fitting equation of equivalent inductance of DGS structure
        Lc_predict1(m,1) = p00_Lc + p10_Lc*x + p01_Lc*y + p20_Lc*x.^2 + ...
                           p11_Lc*x.*y + p02_Lc*y.^2 + p30_Lc*x.^3 + ...
                           p21_Lc*x.^2.*y + p12_Lc*x.*y.^2 + p03_Lc*y.^3 + ...
                           p40_Lc*x.^4 + p31_Lc*x.^3.*y + p22_Lc*x.^2.*y.^2 + ...
                           p13_Lc*x.*y.^3 + p04_Lc*y.^4 + p50_Lc*x.^5 + ...
                           p41_Lc*x.^4.*y + p32_Lc*x.^3.*y.^2 + ...
                           p23_Lc*x.^2.*y.^3 + p14_Lc*x.*y.^4 + p05_Lc*y.^5;

        % Fitting equation of equivalent inductance of DGS structure
        c1_predict1(m,1) = p00_c1 + p10_c1*x + p01_c1*y + p20_c1*x.^2 + ... 
                           p11_c1*x.*y + p02_c1*y.^2 + p30_c1*x.^3 + p21_c1*x.^2.*y + ...
                           p12_c1*x.*y.^2 + p03_c1*y.^3 + p40_c1*x.^4 + ...
                           p31_c1*x.^3.*y + p22_c1*x.^2.*y.^2 + ...
                           p13_c1*x.*y.^3 + p04_c1*y.^4 + p50_c1*x.^5 + ...
                           p41_c1*x.^4.*y + p32_c1*x.^3.*y.^2 + ...
                           p23_c1*x.^2.*y.^3 + p14_c1*x.*y.^4 + p05_c1*y.^5;
        
        if Lc_predict1(m,1)<= max(Lc_valid)&&Lc_predict1(m,1)>= min(Lc_valid)&& c1_predict1(m,1)<= max(c1_valid)&&c1_predict1(m,1) >= min(c1_valid)
            
            a_valid(1,n)=x;              %Collection of legal length a
            w_valid(1,n)=y;              %Collection of legal width w
            n=n+1; 
        end
        m=m+1;
    end
end

%By combining the above two sections, we can calculate the DGS size 
% that meets theusers' requirments, so as to realize the rapid design of passive tags


%% Draw a scatter diagram to show the legal LC, CC, C1 found in last section

%Generate mesh
[x,y]=meshgrid(5:0.005:12,0.01:0.002:1.5);

% Fitting equation of equivalent inductance of DGS structure
Lc_predict = p00_Lc + p10_Lc*x + p01_Lc*y + p20_Lc*x.^2 + p11_Lc*x.*y + ...
             p02_Lc*y.^2 + p30_Lc*x.^3 + p21_Lc*x.^2.*y + p12_Lc*x.*y.^2 + ...
             p03_Lc*y.^3 + p40_Lc*x.^4 + p31_Lc*x.^3.*y + p22_Lc*x.^2.*y.^2 + ...
             p13_Lc*x.*y.^3 + p04_Lc*y.^4 + p50_Lc*x.^5 + p41_Lc*x.^4.*y + ...
             p32_Lc*x.^3.*y.^2+p23_Lc*x.^2.*y.^3 + p14_Lc*x.*y.^4 + p05_Lc*y.^5;

% Fitting equation of equivalent inductance of DGS structure
c1_predict = p00_c1 + p10_c1*x + p01_c1*y + p20_c1*x.^2 + p11_c1*x.*y + ...
             p02_c1*y.^2 + p30_c1*x.^3 + p21_c1*x.^2.*y + p12_c1*x.*y.^2 + ...
             p03_c1*y.^3 + p40_c1*x.^4 + p31_c1*x.^3.*y + p22_c1*x.^2.*y.^2 + ...
             p13_c1*x.*y.^3 + p04_c1*y.^4 + p50_c1*x.^5 + p41_c1*x.^4.*y + ...
             p32_c1*x.^3.*y.^2 + p23_c1*x.^2.*y.^3 + p14_c1*x.*y.^4 + p05_c1*y.^5;

%Upper and lower boundaries of legal equivalent circuit parameters 
%calculated by submodel A
Lc_valid_matrix1 = max(Lc_valid)*ones(size(x,1),size(x,2));
Lc_valid_matrix2 = min(Lc_valid)*ones(size(x,1),size(x,2));

c1_valid_matrix3 = max(c1_valid)*ones(size(x,1),size(x,2));
c1_valid_matrix4 = min(c1_valid)*ones(size(x,1),size(x,2));

%The intersection of the fitting surface of the equivalent circuit 
% and the legal upper and lower boundaries
[B,I] = find(abs(Lc_predict-max(Lc_valid))<=0.002 );
for k=1:length(B)
    xx(k,:) = x(B(k),I(k));
    yy(k,:) = y(B(k),I(k));
    Lc_inter(k,:) = Lc_predict(B(k),I(k));
end

[B1,I1] = find(abs(Lc_predict-min(Lc_valid))<= 0.002);
for k=1:length(B1)
    xx1(k,:) = x(B1(k),I1(k));
    yy1(k,:) = y(B1(k),I1(k));
    Lc_inter1(k,:) = Lc_predict(B1(k),I1(k));
end

[B2,I2] = find(abs(c1_predict-max(c1_valid))<= 0.002);
for k=1:length(B2)
    xx2(k,:) = x(B2(k),I2(k));
    yy2(k,:) = y(B2(k),I2(k));
    c1_inter(k,:) = c1_predict(B2(k),I2(k));
end

[B3,I3] = find(abs(c1_predict-min(c1_valid))<= 0.002);
for k=1:length(B3)
    xx3(k,:) = x(B3(k),I3(k));
    yy3(k,:) = y(B3(k),I3(k));
    c1_inter1(k,:) = c1_predict(B3(k),I3(k));
end
%% Draw Lc
figure('color',[1 1 1 ])
% Draw curved surface of Lc
s=surf(x,y,Lc_predict);
map=colormap(flipud(bone));
s.EdgeColor='none';
s.LineStyle='--';
s.FaceAlpha=1;
s.FaceLighting='gouraud';
hold on

% Draw the upper boundary of Lc
s1=surf(x,y,Lc_valid_matrix1);
s1.EdgeColor='none';
s1.LineStyle='--';
s1.FaceColor=[213 41 65]/255;
s1.FaceAlpha=0.3;
hold on
% Draw the lower boundary of Lc
s2=surf(x,y,Lc_valid_matrix2);
s2.EdgeColor='none';
s2.LineStyle='--';
s2.FaceColor=[128 128 128]/255;
s2.FaceAlpha=0.6;
hold on

% Draw an intersection line 0f Lc
scatter3(xx,yy,Lc_inter,10,'filled','MarkerFaceColor',[213 41 65]/255);
hold on
scatter3(xx1,yy1,Lc_inter1,10,'filled','MarkerFaceColor',[213 41 65]/255);
xlim([5,12]);
legend('Lc = f(a,w)','z = Lc_{max}','z = Lc_{min}','FontSize',8,'FontWeight','normal','FontName','Arial')
legend('boxoff')
xlabel('a (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
ylabel('w (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
zlabel('Lc (nH)','FontSize',10,'FontWeight','bold','FontName','Arial');
set(gca,'Position',[0.2 0.2 0.4 0.4],'FontSize',10,'FontWeight','bold','FontName','Arial');
%% Draw cc 

cc = cal_cc(5,0.01,12,1,0.1,15)*1e12;

% Generate mesh
[x_cc,y_cc] = meshgrid(5:0.01:12,0.1:0.01:1.5);

% Calculate the matrix of upper and lower boundaries
cc_valid_matrix1 = max(cc_valid)*ones(size(x_cc,1),size(x_cc,2));
cc_valid_matrix2 = min(cc_valid)*ones(size(x_cc,1),size(x_cc,2));

figure('color',[1 1 1 ])

% Draw curved surface of cc
s_cc=surf(x_cc,y_cc,cc');
map=colormap(flipud(bone));
s_cc.EdgeColor='none';
s_cc.LineStyle='--';
s_cc.FaceAlpha=1;
s_cc.FaceLighting='gouraud';
hold on

% Draw the upper boundary of cc
s1_cc=surf(x_cc,y_cc,cc_valid_matrix1);
s1_cc.EdgeColor='none';
s1_cc.LineStyle='--';
s1_cc.FaceColor=[213 41 65]/255;
s1_cc.FaceAlpha=0.3;
hold on

% Draw the lower boundary of cc
s2_cc=surf(x_cc,y_cc,cc_valid_matrix2);
s2_cc.EdgeColor='none';
s2_cc.LineStyle='--';
s2_cc.FaceColor=[128 128 128]/255;
s2_cc.FaceAlpha=0.6;
xlim([5,12])
ylim([0.1,1.5])
xlabel('a (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
ylabel('w (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
zlabel('cc (pF)','FontSize',10,'FontWeight','bold','FontName','Arial');
legend('cc = f(a,w)','cc_m_a_x','cc_m_i_n','FontSize',8,'FontWeight','bold','FontName','Arial');
set(gca,'Position',[0.2 0.2 0.4 0.4],'FontSize',10,'FontWeight','bold','FontName','Arial');

%% Draw c1
figure('color',[1 1 1 ])

% Draw curved surface of cc
s=surf(x,y,c1_predict);
xlim([5,12]);
map=colormap(flipud(bone));
s.EdgeColor='none';
s.LineStyle='--';
s.FaceAlpha=1;
s.FaceLighting='gouraud';
hold on

% Draw the upper boundary of c1
s4 = surf(x,y,c1_valid_matrix3);
hold on

% Draw the lower boundary of c1
s5 = surf(x,y,c1_valid_matrix4);
hold on

% Draw the intersection line between the surface and the upper boundary
scatter3(xx2,yy2,c1_inter,10,'filled','MarkerFaceColor',[49 97 173]/255);
hold on
% Draw the intersection line between the surface and the lower boundary
scatter3(xx3,yy3,c1_inter1,10,'filled','MarkerFaceColor',[49 97 173]/255);

s2.EdgeColor='none';
s2.LineStyle='--';
s2.FaceColor=[49 97 153]/255;
s2.FaceAlpha=1;

s3.EdgeColor='none';
s3.LineStyle='--';
s3.FaceColor=[213 41 65]/255;
s3.FaceAlpha=0.4;


s4.EdgeColor='none';
s4.LineStyle='--';
s4.FaceColor=[213 41 65]/255;
s4.FaceAlpha=0.3;

s5.EdgeColor='none';
s5.LineStyle='--';
s5.FaceColor=[128 128 128]/255;
s5.FaceAlpha=0.6;
zlim([0,2.5])
legend('c1 = f(a,w)','z = c1_{max}','z = c1_{min}','FontSize',8,'FontWeight','normal','FontName','Arial')
legend('boxoff')
xlabel('a (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
ylabel('w (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
zlabel('c1 (pF)','FontSize',10,'FontWeight','bold','FontName','Arial');
set(gca,'Position',[0.2 0.2 0.4 0.4],'FontSize',10,'FontWeight','bold','FontName','Arial');

%% Draw the intersection part of the equivalent circuit parameter set (meet the users' requirment)

save('yy1.mat','yy1');
save('xx1.mat','xx1');

save('yy2.mat','yy2');
save('xx2.mat','xx2');

save('yy3.mat','yy3');
save('xx3.mat','xx3');

z1 = zeros(1,length(a_valid));

load xx1.mat
load yy1.mat

load xx2.mat
load yy2.mat

load xx3.mat
load yy3.mat

figure('color',[1 1 1])

plot(xx,yy,'linewidth',2,'color',[213 41 65]/255);
hold on
plot(xx1,yy1,'linewidth',2,'color',[213 41 65]/255);
hold on
plot(xx2,yy2,'linewidth',2,'color',[49 97 153]/255);
hold on

plot(xx3,yy3,'linewidth',2,'color',[49 97 153]/255);
hold on

% Draw the part of the set consisting of all legal A and w that meet users'
% requirement
scatter(a_valid,w_valid,4,'filled','MarkerFaceColor',[128 128 128]/255);

xlim([5,12])
xlabel('a (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
ylabel('w (mm)','FontSize',10,'FontWeight','bold','FontName','Arial');
zlabel('Lc (nH)','FontSize',10,'FontWeight','bold','FontName','Arial');
legend('Lc_m_a_x','Lc_m_i_n','c1_m_a_x','c1_m_i_n','FontSize',8,'FontWeight','normal','FontName','Arial');
set(gca,'Position',[0.2 0.2 0.4 0.4],'FontSize',10,'FontWeight','bold','FontName','Arial');