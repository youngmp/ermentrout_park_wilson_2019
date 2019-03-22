function dydt = izkF(t,y,a,b)
    dydt = [0.04*y(1)^2 + 5*y(1) + 140-y(2) + 10;
            a*(b*y(1)-y(2))];
end