using LinearAlgebra, TaylorSeries

# Quaternion kinematics for the error-state Kalman filter
# Joan Sola

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
    theta = norm(w)

    if approx
        t = Taylor1(Float64, 20)
        return evaluate([cos(t / 2); [1; 1; 1] * eps() * sin(t / 2) / t], theta)
    end

    if theta < eps()
        return [1; 0; 0; 0] # code will break here
    end

    return [cos(theta / 2); w / theta * sin(theta / 2)]
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
    
    # from quatLog.
    axis = Δq[2:4] / norm(Δq[2:4])
    angle = frac * atan(norm(Δq[2:4]), Δq[1]) / 2

    Δq_t = [cos(angle); sin(angle) * axis]
    qt = quatL(q0) * Δq_t

    return qt
end

function quatLog(quat)
    axis = quat[2:4] / norm(quat[2:4])
    angle = atan(norm(quat[2:4]), quat[1])

    return [log(norm(quat)); axis * angle]
end

function to_matrix(quat)
    #return quatR(conjugate(quat)) * quatL(quat) [2:4, 2:4]
    return (quat[1]^2 - quat[2:4] ⋅ quat[2:4]) * I(3) + 2 * quat[2:4] * quat[2:4]' + 2 * quat[1] * skew(quat[2:4]) # note, quat[2:4] * quat[2:4]' gives a matrix
end