
from typing import Dict, Union, List, Collection

import geopandas as gpd
import numpy as np
import pandas as pd
from pyinterpolate import Deconvolution
from tqdm import trange
import matplotlib.pyplot as plt

from pyinterpolate.processing.checks import check_limits
from pyinterpolate.processing.preprocessing.blocks import Blocks, PointSupport
from pyinterpolate.processing.transform.transform import get_areal_centroids_from_agg
from pyinterpolate.variogram import build_experimental_variogram, TheoreticalVariogram, ExperimentalVariogram
from pyinterpolate.variogram.regularization.aggregated import regularize
def plot_variograms(deconvolutionObj: Deconvolution, figPath: str,titleStr: str):
    """
    Function shows experimental semivariogram, theoretical semivariogram and regularized semivariogram after
    semivariogram regularization with ``transform()`` method.
    """
    lags = deconvolutionObj.initial_experimental_variogram.lags

    plt.figure(figsize=(12, 6))
    plt.plot(lags,
             deconvolutionObj.initial_experimental_variogram.experimental_semivariances, 'bo')
    plt.plot(lags,
             deconvolutionObj.initial_theoretical_agg_model.predict(lags), color='r',
             linestyle='--')

    if deconvolutionObj.final_optimal_variogram is not None:
        plt.plot(lags, deconvolutionObj.final_optimal_variogram[:, 1], 'go')

        plt.plot(lags,
                 deconvolutionObj.final_theoretical_model.predict(lags), color='black', linestyle='dotted')
        plt.legend(['Experimental semivariogram of areal data', 'Initial Semivariogram of areal data',
                    'Regularized data points, iteration {}'.format(deconvolutionObj.iter),
                    'Optimized theoretical point support model'])
        plt.title(titleStr+'Semivariograms comparison. Deviation value: {}'.format(deconvolutionObj.optimal_deviation))
    else:
        plt.legend(['Experimental semivariogram of areal data', 'Initial Semivariogram of areal data'])
        plt.title(titleStr+ 'initial_fit: Semivariograms comparison')
    plt.xlabel('Distance')
    plt.ylabel('Semivariance')

    plt.savefig(fname=figPath)
    # plt.show()

def plot_deviations(deconvolutionObj: Deconvolution, figPath: str,titleStr: str):
    plt.figure(figsize=(12, 6))
    plt.plot(
        np.arange(len(deconvolutionObj.deviations)),
        [x / deconvolutionObj.initial_deviation for x in deconvolutionObj.deviations]
    )
    plt.xlabel('Iteration')
    plt.ylabel('Deviation')
    plt.title(titleStr)
    plt.savefig(fname=figPath)
    # plt.show()

def plot_weights(deconvolutionObj: Deconvolution, figPath: str,titleStr: str):
    plt.figure(figsize=(12, 6))
    plt.plot(
        np.arange(len(deconvolutionObj.weights)),
        [np.mean(weight) for weight in deconvolutionObj.weights]
    )
    plt.xlabel('Iteration')
    plt.ylabel('Average weight')
    plt.title(titleStr)
    plt.savefig(fname=figPath)
    # plt.show()
