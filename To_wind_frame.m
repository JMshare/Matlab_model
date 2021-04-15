function [L, D, M] = To_wind_frame(F, M, Alpha)
    [R_V] = Wind_frame_matrix(Alpha);
    FM = [F;M];
    FMw = R_V*FM;
    L = FMw(1);
    D = FMw(3);
    M = M(2);
end