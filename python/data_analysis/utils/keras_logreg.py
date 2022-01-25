import numpy as np

from sklearn.linear_model._base import BaseEstimator
from sklearn.model_selection import GridSearchCV

from keras.models import Sequential
from keras.layers import Dense
from keras.regularizers import L1L2
from keras.wrappers.scikit_learn import KerasClassifier
from keras.utils import np_utils

# from tensorflow.config import run_functions_eagerly
# run_functions_eagerly(True)        

def get_model(alpha=0.5, activation='sigmoid', input_dim=10, optimizer='sgd', loss='binary_crossentropy', metrics='binary_accuracy'):
    model = Sequential() 
    
    model.add( Dense(1,  # output dim is 2, one score per each class
                     activation=activation,
                     kernel_regularizer=L1L2(l1=alpha, l2=1.0-alpha),
                     input_dim=input_dim)  # input dimension = number of features your data has 
    )
    
    model.compile(optimizer=optimizer,
                  loss=loss,
                  metrics=metrics
    )
        
    return model

class keras_logreg(BaseEstimator): 
    
    def __init__(self, alpha=0.5, activation='sigmoid', input_dim=10, optimizer='sgd', loss='binary_crossentropy', metrics='binary_accuracy', epochs=10, verbose=True): 
        
        self.alpha = alpha 
        self.activation = activation 
        self.optimizer = optimizer 
        self.loss = loss 
        self.metrics = metrics 
        self.epochs = epochs
        self.input_dim = input_dim
        self.verbose = verbose 
        self.model = get_model(alpha, activation, input_dim, optimizer, loss, metrics) 
    
    def fit(self, X, y):
        self.model.fit(X, y, epochs=self.epochs, shuffle=True, verbose=self.verbose) 
        return self 
    
    def predict(self, X):         
        y_pred = self.model.predict(X) 
        return y_pred
    
    def predict_proba(self, X):
        y_prob = self.model.predict_proba(X) 
        return y_prob 
    
    def score(self, X, y): 
        score_ = self.model.evaluate(X, y, verbose=self.verbose)[-1] 
        return score_ 
    
    def get_coefs(self): 
        self.coef_ = self.model.get_weights() 
        return self 
    
class keras_logreg_cv(BaseEstimator):

    def __init__(self, cv=5, alphas=[0, 1, 100], alpha=0.5, activation='sigmoid', input_dim=10, optimizer='sgd', loss='binary_crossentropy', metrics='binary_accuracy', epochs=10, verbose=False, n_jobs=None): 
                
        self.cv = cv 
        self.alphas = alphas 
        self.n_jobs = n_jobs 
        
    def fit(self, X, y):
        activations = ['softmax', 'softplus', 'softsign', 'relu', 'tanh', 'sigmoid', 'hard_sigmoid', 'linear'] 
        optimizers = ['rmsprop', 'adam'] 
        alphas = np.linspace(0,1,10) 
        
        param_grid = dict(alpha=alphas) 

        clf = KerasClassifier(build_fn=get_model, epochs=10, batch_size=1, verbose=0) 
        grid = GridSearchCV(clf, param_grid=param_grid, cv=self.cv, n_jobs=self.n_jobs) 
        
        Y = np_utils.to_categorical(y, X.shape[-1]) 
        print(X.shape, Y.shape) 
        
        grid.fit(X, Y) 
        
        # self.model = grid.best_estimator_  
        # self.model.fit(X, y) # maybe useless to fit again 
        
        return self 
    
