function vor2dl(gj, gj1, x, z, xj, zj, xj1, zj1)

    # x, z is the target

    #l = sqrt((x2 - x1)^2 + (z2 - z1)^2)
    rj = sqrt((x - xj)^2 + (z - zj)^2)
    rj1 = sqrt((x - xj1)^2 + (z - zj1)^2)
    r = sqrt((xj - xj1)^2 + (zj - zj1)^2)

    thetaj = acos(minimum([(r^2 + rj^2 - rj1^2)/(2*rj*r) 1.0]))
    thetaj1 = -acos(minimum([(r^2 + rj1^2 - rj^2)/(2*rj1*r) 1.0]))

    alpha = atan((zj1 - zj)/(xj1 - xj))

    trans = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]

    up = (z/(2*pi)) * ((gj1 - gj)/(xj1 - xj)) * log(rj1/rj) + (gj*(xj1 - xj) + (gj1 - gj)*(x - xj))*(thetaj1 - thetaj)/(2*pi*(xj1 - xj))
    wp = -(gj*(xj1 - xj) + (gj1 - gj)*(x - xj))*log(rj/rj1)/(2*pi*(xj1 - xj)) + (z/(2*pi)) * ((gj1 - gj)/(xj1 - xj)) * ((xj1 - xj)/z + thetaj1 - thetaj)

    upa = (z/(2*pi)) * ((-gj)/(xj1 - xj)) * log(rj1/rj) + (gj*(xj1 - xj) - gj*(x - xj))*(thetaj1 - thetaj)/(2*pi*(xj1 - xj))
    upb = (z/(2*pi)) * ((gj1)/(xj1 - xj)) * log(rj1/rj) + (gj1*(x - xj))*(thetaj1 - thetaj)/(2*pi*(xj1 - xj))
    wpa = -(gj*(xj1 - xj) + (-gj)*(x - xj))*log(rj/rj1)/(2*pi*(xj1 - xj)) + (z/(2*pi)) * ((- gj)/(xj1 - xj)) * ((xj1 - xj)/z + thetaj1 - thetaj)
    wpb = -((gj1)*(x - xj))*log(rj/rj1)/(2*pi*(xj1 - xj)) + (z/(2*pi)) * ((gj1)/(xj1 - xj)) * ((xj1 - xj)/z + thetaj1 - thetaj)

    (u, w) = trans*[up; wp]
    (ua, wa) = trans*[upa; wpa]
    (ub, wb) = trans*[upb; wpb]

    return u, w, ua, wa, ub, wb
end
