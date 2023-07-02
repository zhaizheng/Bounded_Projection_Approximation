x = linspace(-10,10,1000);
y = 100000*((min(x-0,0)).^2+(min(5-x,0)).^2);
plot(x,y,'linewidth',4);
z = max(y)*ones(size(x));
z(logical((x<5).*(x>0))) = 0;
hold on 
plot(x,z,'linewidth',4);