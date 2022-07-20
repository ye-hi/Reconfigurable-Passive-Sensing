%% Calculate the capacitance of cc according to the structure of parallel plate capacitor
function cc = fun1(a_s,a_d,a_e,w_s,w_d,w_e)
h = 18*10^-6;
k = 8.99*10^9;
c = 2e-3;
m = 1;

for j = a_s:a_d:a_e
    n=1;
    for i = w_s:w_d:w_e
        a = j*1e-3;
        w = i*1e-4;
        s = h*3*(a-2*w)/(4*k*pi*w);
        s1 = 3*h*w/(4*k*pi*a);
        s2 = 2*h*w/(4*k*pi*(a-c)/2);
        cc(m,n) = s+s1+s2;
        n=n+1;
    end
    m = m+1;
end
end