%computes a new step size such that the error estimate is equal to some
%tolerance value tol, using the PI controller
function hnew = newstep(tol, err, errold, hold, k)
hnew = (tol/err).^(2/(3*k))*(tol/errold).^(-1/(3*k))*hold;
end