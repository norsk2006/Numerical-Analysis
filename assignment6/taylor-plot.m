n = 0.2;

for i = 1:4
    
    n = n*10;

    [t,w,y,ge] = taylorbrine(n);

    figure

    plot(t,w,'red')
    
    figure

    plot(t,ge,'blue')
    
end

% set appropriate value to decfactor here - replace the 0 value with correct value
decfactor = -3;

% complete the function, taylorbrine(), below

function [t,w,y,ge] = taylorbrine(n)

    h = 10/n; %step size
    t = linspace(0, 10, 10/h + 1); %vector of values of t
    y = 500*(50-t).*(log(50./(50-t))+1);
    w(1) = y(1);

    for i = 1:n
        w(i+1) = w(i) + h*(500+(w(i)/(t(i)-50))) + (h^2/2)*(500/(t(i)-50)) - (h^3/6)*(500/(t(i)-50)^2);
    end

   ge = abs(w-y);


end