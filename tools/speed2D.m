function v = speed2D(x,y,t)
v = zeros(length(t),1);

for i = 2:length(t)-1
    v(i) = sqrt((x(i+1)-x(i-1))^2+(y(i+1)-y(i-1))^2)/(t(i+1)-t(i-1));
end
v(1) = v(2);
v(end) = v(end-1);
end