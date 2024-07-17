import numpy as np
from sklearn.mixture import GaussianMixture
import pandas as pd

def calc_bimodality_index(data):
    # Convert the DataFrame column to a numpy array if input is a pandas DataFrame
    if isinstance(data, pd.DataFrame):
        data = data.values.flatten()
    
    # Reshape the data to 2D array for GaussianMixture model
    data = data.reshape(-1, 1)
    
    # Fit a mixture of two normal distributions
    gmm = GaussianMixture(n_components=2, covariance_type='tied')
    gmm.fit(data)
    
    # Extract parameters
    means = gmm.means_.flatten()
    std_dev = np.sqrt(gmm.covariances_).flatten()[0]
    weights = gmm.weights_
    
    # Sort the means and weights to ensure mean1 < mean2
    if means[0] > means[1]:
        means = means[::-1]
        weights = weights[::-1]
    
    mean1, mean2 = means
    weight1, weight2 = weights
    
    # Calculate the standardized distance
    delta = abs(mean1 - mean2) / std_dev
    
    # Calculate the bimodality index
    pi = min(weight1, weight2)
    bi = delta * np.sqrt(pi * (1 - pi))
    
    return bi
