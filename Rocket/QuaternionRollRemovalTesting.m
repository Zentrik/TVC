local_roll = quatexp([0 0 0 pi/4]);
quat = quatmultiply([0.71 0.71 0 0], local_roll); %current orientation, first quat then rotation about local z axis
% quatmultiply(changeofbasis([0.71 0.71 0 0], local_roll), [0.71 0.71 0 0]) %testing change of basis by working out quat using local roll

k = [0 0 1];
desired = k;
vec = quatrot(quat, k); %vector representing current rocket orientation, should be invariant to local_roll
error = quatnormalize([norm(vec) * norm(desired) + dot(vec, desired), cross(vec, desired)]); % rotates vec to desired
% alpha = 2 * atan(quat(4)/quat(1));
% error = quatnormalize(quatmultiply([cos(alpha/2) 0 0 sin(alpha/2)], quatinv(quat))); % need to add rotation to desired if desired not k

% quatmultiply(o, quat); should give local_roll
local_error = changeofbasis(quatinv(quat), error) % if o is in new basis, then old basis must be local

% quat = sym('quat', [1 4], 'real');
% assumeAlso(sum(quat.^2) == 1)
% e = sym('e', [1 4], 'real');
% assumeAlso(sum(e.^2) == 1)
% w = sym('w',  [1 3], 'real');
% extra_derivative = simplify(.5 * quatmultiply(quat, quatmultiply([0 w], quatmultiply(e, quatinv(quat)))) + .5 * quatmultiply(quat, quatmultiply(e, .5 * quatmultiply(quatinv(quat), [0 w]))));
% test = subs(extra_derivative, quat, [1 0 0 0])
subs(test, , [1 0 0 0])

function y = quatrot(quat, vec)
    y = quatrotate(quatinv(quat), vec);
end

function quot = changeofbasis(rotation_from_old_basis_to_new, rotation_in_new_basis)
    quot = quatmultiply(rotation_from_old_basis_to_new, quatmultiply(rotation_in_new_basis, quatinv(rotation_from_old_basis_to_new))); % rotate from new basis to old, rotate, then rotate back to new basis
end