using LinearAlgebra

# Quaternion kinematics for the error-state Kalman filter
# Joan Sola

export quatL, quatR, skew, wexp, wexp_w, conjugate, rotate, slerp, slerp_quat, qualLog, to_matrix

function quatL(quat)
    S = zeros(4, 4)
    S += quat[1] * I(4)
    
    S[2:4, 1] = quat[2:4]
    S[1, 2:4] = -S[2:4, 1]
    S[2:4, 2:4] += skew(quat[2:4])

    return S
end

function quatR(quat)
    S = zeros(4, 4)
    S += quat[1] * I(4)
    
    S[2:4, 1] = quat[2:4]
    S[1, 2:4] = -S[2:4, 1]
    S[2:4, 2:4] -= skew(quat[2:4])

    return S
end

function skew(w)
    [0    -w[3]  w[2];
     w[3]  0    -w[1];
    -w[2]  w[1]  0]
end

function wexp(w, approx=false)
    if approx
        t = Taylor1(Float64, 20)
        return evaluate([cos(t / 2); [1; 1; 1] * eps() * sin(t / 2) / t], norm(w))
    end

    # Written in terms of θ² = w ⋅ w rather than θ = norm(w) so that it stays
    # smooth (and ForwardDiff differentiable) at w = 0. `norm` is not
    # differentiable there — it returns a NaN partial — and the old `theta <
    # eps()` early return produced a *constant*, so `ForwardDiff.derivative` of
    # anything wrapping `wexp` silently evaluated to zero whenever it was
    # evaluated at w = 0. That is exactly the case hit by the coast-time
    # Jacobian in the guidance problem, whose reference value of t_coast is 0.
    thetasq = dot(w, w)

    if thetasq < 1e-12 # series expansions of cos(θ/2) and sin(θ/2)/θ in θ²
        return [1 - thetasq / 8 + thetasq^2 / 384; w * (1/2 - thetasq / 48 + thetasq^2 / 3840)]
    end

    theta = sqrt(thetasq)

    return [cos(theta / 2); w * (sin(theta / 2) / theta)]
end

function wexp_w(w)
    theta = norm(w)
    
    if theta < eps()
        return I(3)
    end

    J = I(3) - (1 - cos(theta))/theta^2 * skew(w) + (theta - sin(theta)) / theta^3 * skew(w)^2

    return J #joan sola
end

function conjugate(quat)
    return [quat[1]; - quat[2:4]]
end

function rotate(quat, vector)
    # tmp = quatL(quat) * quatL([0; vector]) * conjugate(quat)
    # return tmp[2:4]
    return to_matrix(quat) * vector
end

#interpolates from v to w by frac ∈ [0, 1]
function slerp(v, w, frac) 
    # axis = cross(v, w) / norm(cross(v, w))
    # angle = atan(norm(cross(v, w)), v' * w) / 2
    # # tmp = [v' * w / 2; cross(v, w)]
    # # quat = quat / norm(tmp)

    # interpolator = [cos(frac * angle); sin(frac * angle) * axis]

    # return rotate(interpolator, v)

    angle = acos(v' * w)
    return (sin((1 - frac) * angle) * v + sin(frac * angle) * w) / sin(angle)
end

#interpolates from v to w by frac ∈ [0, 1]
function slerp_quat(q0, q1, frac)
    Δq = quatL(conjugate(q0)) * q1 # quaternion rotating q0 to q1
    # Δq_t = wexp(frac * quatLog(Δq)[2:4])

    axis, angle = quatLogAxisAngle(Δq) # angle is the half angle, i.e. Δq = [cos(angle); sin(angle) * axis]

    # `frac` used to be missing here, so this returned normalize(q1) for every
    # value of `frac`, i.e. it was not an interpolation at all.
    Δq_t = [cos(frac * angle); sin(frac * angle) * axis]
    qt = quatL(q0) * Δq_t

    return qt
end

function quatLogAxisAngle(quat)
    if norm(quat[2:4]) ≈ 0
        axis = [0; 0; 1]
    else
        axis = quat[2:4] / norm(quat[2:4])
    end
    angle = atan(norm(quat[2:4]), quat[1])

    return axis, angle
end

function quatLog(quat)
    if norm(quat[2:4]) ≈ 0
        axis = [0; 0; 1]
    else
        axis = quat[2:4] / norm(quat[2:4])
    end
    angle = atan(norm(quat[2:4]), quat[1])

    return [log(norm(quat)); axis * angle]
end

function to_matrix(quat)
    #return quatR(conjugate(quat)) * quatL(quat) [2:4, 2:4]
    return (quat[1]^2 - quat[2:4] ⋅ quat[2:4]) * I(3) + 2 * quat[2:4] * quat[2:4]' + 2 * quat[1] * skew(quat[2:4]) # note, quat[2:4] * quat[2:4]' gives a matrix
end