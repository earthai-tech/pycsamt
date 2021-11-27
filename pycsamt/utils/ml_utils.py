# -*- coding: utf-8 -*-
#       This module is a part of pycsamt utils packages related to the ML
#       online module to predict the layer name.

import numpy as np 
import matplotlib.pyplot as plt 
from random import seed
from random import randrange
from random import random
from csv import reader
from sklearn.metrics import accuracy_score 
from sklearn.datasets import make_blobs 
from sklearn.base import BaseEstimator, TransformerMixin    


__all__= ['neuron', 'decision_boundary', 'ReduceImage']

###################TEST DATA ##################################################
X, y = make_blobs (n_samples = 100 , n_features =2 ,
                   centers =2 , random_state =0 )
y= y.reshape((y.shape[0], 1))
###############################################################################

def quick_view (X, c=y, cmap='summer'): 
    """ quick view datasets X"""
    fig =plt.figure(figsize =(7, 3))
    ax= fig.add_subplot(111)
    ax.scatter(X[:, 0], X[:, 1], c=y, cmap =cmap)
 
    return ax, fig

def initialize_model (X):
    """ Initialize the model parameters `W` and `b` with random values"""
    # W shape = number of features in X including the bias 
    W  = np.random.rand(X.shape[1], 1)
    # b is a real value 
    b = np.random.randn(1)
    
    return W , b 

def model (X, W, b, optimizer ='sigmoid'): 
    """ Create model and set the activation `A` """
    # Z= X.W + b
    Z = X.dot(W) +b 
    if optimizer =='sigmoid': 
        A = 1 / ( 1+ np.exp(-Z))
    return  A 

def log_loss (A, y): 
    """ Compute the model loss error and return a value."""
    # Loss = -1/m SUM (y * log(A) +(1-y)* log(1-A) )
    return  1/ len(y) * np.sum ( -y * np.log(A) - (1- y)* np.log(1- A))

def weights (A, X,  y): 
    """Compute the gradients <Jacobian> using the activation `A` , `y` 
    applied on `X`."""
    dw = 1/ len(y) * np.dot( X.T, A-y) # dw = 1/m XT .(A-y)
    db = 1/ len(y) * np.sum(A-y)  # db = 1/m (A-y)
    
    return dw, db 

def update_Ws(dw, db, W, b, alpha): 
    """ Update the weights at every `n_epochs iteration. `
    iterations."""
    
    W= W - alpha * dw 
    b= b - alpha * db 
    return W, b 

def predict (X, W, b): 
    """ Predict new value from model X using the log loss  and return ``True``
    if probability >=0.5 and ``False `` otherwise."""
    #  The predicted class is binary class 1 and 0.
    A= model(X, W, b)
    # print('Score =', A )
    return A >=0.5

def neuron (X, y, alpha =0.1, n_epochs=100, optimizer ='sigmoid',
            scoring ='accuracy', plot ='yes', **kws ): 
    """ Build the first neuron and return the updated weights and the 
    history duiring the learnings composed of weights `W`, `b` and  """
    # initialize weight values W and b
    W, b = initialize_model(X) 
    
    history , Loss =list(), list()
    # loop for the n_epochs 
    for i in range(n_epochs): 
        A = model(X, W, b, optimizer) # return activation 
        Loss.append(log_loss(A, y)) # keep the loss value 
        dw, db = weights(A, X, y) # compute jacobian 
        W, b = update_Ws(dw, db, W, b, alpha) # update weights 
        history.append([W, b, Loss, i]) # save all history during the learnings 
    # we can make the predictions and compute the probabilities using 
    #`accuracy_score `of scikit-learn 
    y_pred =predict(X, W, b) # W and b are the new updated weights 
    proba = accuracy_score (y, y_pred) 
    print(f"Performance < {scoring}>=", proba)
    
    # can visualize the plot after the n_epochs 
    if plot =='yes':
        plt.plot(Loss, **kws)
        plt.xlabel('number of iterations')
        plt.ylabel('Error activation')
        plt.show()
  
    return W, b, history

def decision_boundary (W,b,  X, y): 
    """ compute the decision boundary using the default abscisses"""
    # specifity the plot abscisses from -1 to 4  (z(x1, x2)=0)
    x0 = np.linspace(-1, 4 , 100) # with 100 points
    # in decision boundary Z= 0  then Z= x1w1 + x2w2 +b  =0 --> w2 = -(x1w1)/w2
    x1 = (- x0 * W[0]  -b )/ W[1]
    
    ax, fig = quick_view(X, c=y)
    ax.plot(x0, x1, c= 'orange', lw =3)
    plt.show()
    
       
class ReduceImage(BaseEstimator, TransformerMixin):
    def __init__(self, normalize = True): 
        self.normalize = normalize
    def fit (self, X): 
        """Fit data and return object """
        return self 
    def transform (self, X): 
        """Transform data, normalize data and 
        converted 64*64 image into a flatterned array. 
        :param X: ndarray dim =2 or 3
        :return: flattened array (n-examples, n_images 64*64) 
        64 x 64 flatterned correspond to 4096 features
        """
         # normalize data [0-1] using linalg norm 
        def normalize_using_linalg_norm(X): 
            return X/ np.linalg.norm(X)
        def flatten_variables(X):  # flatten ndarray 
            if X.ndim ==3: 
                X_= np.zeros((len(X), X.shape[1]*X.shape[2] ))
                for k in range(len(X)): 
                    X_[k]= X[k, :, :].flatten()
            else :X_=  X.flatten()
            return  X_
        if self.normalize: 
            X = normalize_using_linalg_norm(X)
        X = flatten_variables(X)
        return X
  

# https://machinelearningmastery.com/implement-backpropagation-algorithm-scratch-python/  
# Load a CSV file
def load_csv(filename):
    dataset = list()
    with open(filename, 'r') as file:
    	csv_reader = reader(file)
    	for row in csv_reader:
    		if not row:
    			continue
    		dataset.append(row)
    return dataset

# Convert string column to float
def str_column_to_float(dataset, column):
	for row in dataset:
		row[column] = float(row[column].strip())

# Convert string column to integer
def str_column_to_int(dataset, column):
	class_values = [row[column] for row in dataset]
	unique = set(class_values)
	lookup = dict()
	for i, value in enumerate(unique):
		lookup[value] = i
	for row in dataset:
		row[column] = lookup[row[column]]
	return lookup

# Find the min and max values for each column
def dataset_minmax(dataset):
	minmax = list()
	stats = [[min(column), max(column)] for column in zip(*dataset)]
	return stats

# Rescale dataset columns to the range 0-1
def normalize_dataset(dataset, minmax):
	for row in dataset:
		for i in range(len(row)-1):
			row[i] = (row[i] - minmax[i][0]) / (minmax[i][1] - minmax[i][0])

# Split a dataset into k folds
def cross_validation_split(dataset, n_folds):
	dataset_split = list()
	dataset_copy = list(dataset)
	fold_size = int(len(dataset) / n_folds)
	for i in range(n_folds):
		fold = list()
		while len(fold) < fold_size:
			index = randrange(len(dataset_copy))
			fold.append(dataset_copy.pop(index))
		dataset_split.append(fold)
	return dataset_split

# Calculate accuracy percentage
def accuracy_metric(actual, predicted):
	correct = 0
	for i in range(len(actual)):
		if actual[i] == predicted[i]:
			correct += 1
	return correct / float(len(actual)) * 100.0

# Evaluate an algorithm using a cross validation split
def evaluate_algorithm(dataset, algorithm, n_folds, *args):
	folds = cross_validation_split(dataset, n_folds)
	scores = list()
	for fold in folds:
		train_set = list(folds)
		train_set.remove(fold)
		train_set = sum(train_set, [])
		test_set = list()
		for row in fold:
			row_copy = list(row)
			test_set.append(row_copy)
			row_copy[-1] = None
		predicted = algorithm(train_set, test_set, *args)
		actual = [row[-1] for row in fold]
		accuracy = accuracy_metric(actual, predicted)
		scores.append(accuracy)
	return scores

# Calculate neuron activation for an input
def activate(weights, inputs):
	activation = weights[-1]
	for i in range(len(weights)-1):
		activation += weights[i] * inputs[i]
	return activation

# Transfer neuron activation
def transfer(activation):
	return 1.0 / (1.0 + np.exp(-activation))

# Forward propagate input to a network output
def forward_propagate(network, row):
	inputs = row
	for layer in network:
		new_inputs = []
		for neuron in layer:
			activation = activate(neuron['weights'], inputs)
			neuron['output'] = transfer(activation)
			new_inputs.append(neuron['output'])
		inputs = new_inputs
	return inputs

# Calculate the derivative of an neuron output
def transfer_derivative(output):
	return output * (1.0 - output)

# Backpropagate error and store in neurons
def backward_propagate_error(network, expected):
	for i in reversed(range(len(network))):
		layer = network[i]
		errors = list()
		if i != len(network)-1:
			for j in range(len(layer)):
				error = 0.0
				for neuron in network[i + 1]:
					error += (neuron['weights'][j] * neuron['delta'])
				errors.append(error)
		else:
			for j in range(len(layer)):
				neuron = layer[j]
				errors.append(expected[j] - neuron['output'])
		for j in range(len(layer)):
			neuron = layer[j]
			neuron['delta'] = errors[j] * transfer_derivative(neuron['output'])

# Update network weights with error
def update_weights(network, row, l_rate):
	for i in range(len(network)):
		inputs = row[:-1]
		if i != 0:
			inputs = [neuron['output'] for neuron in network[i - 1]]
		for neuron in network[i]:
			for j in range(len(inputs)):
				neuron['weights'][j] += l_rate * neuron['delta'] * inputs[j]
			neuron['weights'][-1] += l_rate * neuron['delta']

# Train a network for a fixed number of epochs
def train_network(network, train, l_rate, n_epoch, n_outputs):
	for epoch in range(n_epoch):
		for row in train:
			outputs = forward_propagate(network, row)
			expected = [0 for i in range(n_outputs)]
			expected[row[-1]] = 1
			backward_propagate_error(network, expected)
			update_weights(network, row, l_rate)

# Initialize a network
def initialize_network(n_inputs, n_hidden, n_outputs):
	network = list()
	hidden_layer = [{'weights':[random() 
                             for i in range(n_inputs + 1)]} 
                 for i in range(n_hidden)]
	network.append(hidden_layer)
	output_layer = [{'weights':[random() 
                             for i in range(n_hidden + 1)]} 
                 for i in range(n_outputs)]
	network.append(output_layer)
	return network

# Make a prediction with a network
def predict(network, row):
	outputs = forward_propagate(network, row)
	return outputs.index(max(outputs))

# Backpropagation Algorithm With Stochastic Gradient Descent
def back_propagation(train, test, l_rate, n_epoch, n_hidden):
	n_inputs = len(train[0]) - 1
	n_outputs = len(set([row[-1] for row in train]))
	network = initialize_network(n_inputs, n_hidden, n_outputs)
	train_network(network, train, l_rate, n_epoch, n_outputs)
	predictions = list()
	for row in test:
		prediction = predict(network, row)
		predictions.append(prediction)
	return(predictions)

if __name__=='main':
    # Test Backprop on Seeds dataset
    seed(1)
    # load and prepare data
    filename = 'seeds_dataset.csv'
    dataset = load_csv(filename)
    for i in range(len(dataset[0])-1):
    	str_column_to_float(dataset, i)
    # convert class column to integers
    str_column_to_int(dataset, len(dataset[0])-1)
    # normalize input variables
    minmax = dataset_minmax(dataset)
    normalize_dataset(dataset, minmax)
    # evaluate algorithm
    n_folds = 5
    l_rate = 0.3
    n_epoch = 500
    n_hidden = 5
    scores = evaluate_algorithm(dataset, back_propagation,
                                n_folds, l_rate, n_epoch, n_hidden)
    print('Scores: %s' % scores)
    print('Mean Accuracy: %.3f%%' % (sum(scores)/float(len(scores)))) 
    
if __name__ =='__main__':
    plot_kws ={'lw':3, 'ls':'--'}
    W, b, *_= neuron(X, y,**plot_kws )
    
    decision_boundary (W, b, X, y)

