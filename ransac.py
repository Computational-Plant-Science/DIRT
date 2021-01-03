import numpy
import scipy  # use numpy if scipy unavailable
import scipy.linalg  # use numpy if scipy unavailable


# Copyright (c) 2004-2007, Andrew D. Straw. All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.

#     * Neither the name of the Andrew D. Straw nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

def ransac(data, model, n, k, t, d, debug=False, return_all=False):
    """
    Fit model parameters to data using the RANSAC algorithm.
    
    This implementation written from pseudocode found at
    http://en.wikipedia.org/w/index.php?title=RANSAC&oldid=116358182
    and was adapted to the DIRT pipeline by Alexander Bucksch

    {{{
    Given:
        data - a set of observed data points
        model - a model that can be fitted to data points
        n - the minimum number of data values required to fit the model
        k - the maximum number of iterations allowed in the algorithm
        t - a threshold value for determining when a data point fits a model
        d - the number of close data values required to assert that a model fits well to data
    Return:
        bestfit - model parameters which best fit the data (or nil if no good model is found)
    iterations = 0
    bestfit = nil
    besterr = something really large
    while iterations < k {
        maybeinliers = n randomly selected values from data
        maybemodel = model parameters fitted to maybeinliers
        alsoinliers = empty set
        for every point in data not in maybeinliers {
            if point fits maybemodel with an error smaller than t
                 add point to alsoinliers
        }
        if the number of elements in alsoinliers is > d {
            % this implies that we may have found a good model
            % now test how good it is
            bettermodel = model parameters fitted to all points in maybeinliers and alsoinliers
            thiserr = a measure of how well model fits these points
            if thiserr < besterr {
                bestfit = bettermodel
                besterr = thiserr
            }
        }
        increment iterations
    }
    return bestfit
    }}}
    """

    iterations = 0
    best_fit = None
    best_err = numpy.inf
    best_inlier_idxs = None

    while iterations < k:
        maybe_idxs, test_idxs = random_partition(n, data.shape[0])
        maybe_inliers = data[maybe_idxs, :]
        test_points = data[test_idxs]
        maybe_model = model.fit(maybe_inliers)
        test_err = model.get_error(test_points, maybe_model)
        also_idxs = test_idxs[test_err < t]  # select indices of rows with accepted points
        also_inliers = data[also_idxs, :]

        if debug:
            print('test_err.min()', test_err.min())
            print('test_err.max()', test_err.max())
            print('numpy.mean(test_err)', numpy.mean(test_err))
            print('iteration %d:len(alsoinliers) = %d' % (
                iterations, len(also_inliers)))

        if len(also_inliers) > d:
            betterdata = numpy.concatenate((maybe_inliers, also_inliers))
            bettermodel = model.fit(betterdata)
            better_errs = model.get_error(betterdata, bettermodel)
            # print besterr
            thiserr = numpy.mean(better_errs)
            if thiserr < best_err:
                best_fit = bettermodel
                best_err = thiserr
                best_inlier_idxs = numpy.concatenate((maybe_idxs, also_idxs))

        iterations += 1
    if best_fit is None:
        return 'nan'
    if return_all:
        return best_fit, {'inliers': best_inlier_idxs}
    else:
        return best_fit


def random_partition(n, n_data):
    """return n random rows of data (and also the other len(data)-n rows)"""
    all_idxs = numpy.arange(n_data)
    numpy.random.shuffle(all_idxs)
    idxs1 = all_idxs[:n]
    idxs2 = all_idxs[n:]
    return idxs1, idxs2


class LinearLeastSquaresModel:
    """linear system solved using linear least squares

    This class serves as an example that fulfills the model interface
    needed by the ransac() function.
    """

    def __init__(self, input_columns, output_columns, debug=False):
        self.input_columns = input_columns
        self.output_columns = output_columns
        self.debug = debug

    def fit(self, data):
        A = numpy.vstack([data[:, i] for i in self.input_columns]).T
        B = numpy.vstack([data[:, i] for i in self.output_columns]).T
        x, resids, rank, s = scipy.linalg.lstsq(A, B)
        return x

    def get_error(self, data, model):
        A = numpy.vstack([data[:, i] for i in self.input_columns]).T
        B = numpy.vstack([data[:, i] for i in self.output_columns]).T
        B_fit = scipy.dot(A, model)
        err_per_point = numpy.sum((B - B_fit) ** 2, axis=1)  # sum squared error per row
        return err_per_point


def ransac_fit(X, Y):
    # setup model
    n_inputs = 1
    n_outputs = 1
    Xnew = []
    Ynew = []
    for idx, i in enumerate(X):
        Xnew.append([])
        Xnew[idx].append(i)
    for idx, i in enumerate(Y):
        Ynew.append([])
        Ynew[idx].append(i)
    all_data = numpy.hstack((Xnew, Ynew))
    # all_data=zip(X, Y)
    input_columns = range(n_inputs)  # the first columns of the array
    output_columns = [n_inputs + i for i in range(n_outputs)]  # the last columns of the array
    debug = False
    model = LinearLeastSquaresModel(input_columns, output_columns, debug=debug)

    linear_fit, resids, rank, s = scipy.linalg.lstsq(all_data[:, input_columns],
                                                     all_data[:, output_columns])

    # run RANSAC algorithm
    ransac_fit, ransac_data = ransac(all_data, model,
                                     20, 1000, 7e3, 30,  # misc. parameters
                                     debug=debug, return_all=True)

    return X, numpy.dot(Xnew, ransac_fit)[:, 0]
