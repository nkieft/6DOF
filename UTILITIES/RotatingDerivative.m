function ddt_r = RotatingDerivative(ddt_i, omega, r)

    ddt_r = ddt_i-cross(omega, r);

end
