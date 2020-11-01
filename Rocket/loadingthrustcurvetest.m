data = xlsread('/home/rahul/Documents/school/Current/EPQ/Rocket/F10_Thrust.xls','Sheet1');

Time = data(:,1);
Thrust = data(:,2);

interp1(Time, Thrust, 5)

xq = 0:0.001:10;
plot(Time, Thrust, 'o', xq, interp1(Time, Thrust, xq), ':.');
xlim([0 10]);