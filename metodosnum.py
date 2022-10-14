import numpy

def error_bound(a, b, err):
    n = numpy.log((b - a) / err) / numpy.log(2)
    return int(numpy.ceil(n))
    
def validate_interval(f, x0, x1):
    var1 = f(x0)
    var2 = f(x1)
    return  var1*var2 < 0
# solve for root using bisection method
def bisection(f, interval, tol):
    """
    param f: find root for function
    param interval: within range
    param tol: root accuracy tolerance
    """

    # extract interval start and end points
    x0, x1 = interval[0], interval[1]

    if x0>x1:
        xaux = x0
        x0=x1
        x1=xaux

    # check interval can be used to solve for root
    if not validate_interval(f, x0, x1):
        return

    # iterations required to find the root within a specified error bound
    n = error_bound(x0, x1, tol)

    counter = 1

    # iterate over error bound
    while True:

        # calculate root approximation
        root_approx = x0 + ((x1 - x0) / 2)

        # evaluate y at current estimate
        y = f(root_approx)

        # check tolerance condition
        if -tol < y < tol:
            # check that error bound actually worked
            print(counter, n)

            # return root approximation
            return root_approx

        # check if next segment is left of bisection
        if validate_interval(f, x0, root_approx):
            x1 = root_approx
        else:
            x0 = root_approx

        # increment counter
        counter += 1