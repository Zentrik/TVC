q = [0.8985957085822401 -5.112254282927898E-7 -0.4387775661053514 2.1784904076780705E-9];
coord = [0 0 10.223155596222691];

w = q(1);
x = q(2);
y = q(3);
z = q(4);

coordx = coord(1);
coordy = coord(2);
coordz = coord(3);

a = -x * coordx - y * coordy - z * coordz;
b = w * coordx + y * coordz - z * coordy;
c = w * coordy - x * coordz + z * coordx;
d = w * coordz + x * coordy - y * coordx;

[a b c d];
quatmultiply(q, [0 coord]);

[-a * x + b * w - c * z + d * y, -a * y + b * z + c * w - d * x, -a * z - b * y + c * x + d * w]

quatmultiply(quatmultiply(q, [0 coord]), quatinv(q));
quatrotate(quatinv(q), coord);