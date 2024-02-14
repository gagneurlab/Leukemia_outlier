import os
import sys
import time
import pickle
import numpy as np
import pandas as pd

from torch.utils.tensorboard import SummaryWriter
import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader

import math
import random
import numpy as np



class driver_dataset(Dataset):
    def __init__(self, features, labels):
        self.features = torch.Tensor(features.values)
        self.labels = torch.Tensor(labels)

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        feature = self.features[idx]
        label = self.labels[idx]        
        return feature, label


class EarlyStopping:
    def __init__(self, patience=5, min_delta=0):

        self.patience = patience
        self.min_delta = min_delta
        self.counter = 0
        self.early_stop = False
        self.best_loss = 10**6
        
        self.best_acc = 0
        self.acc_counter = 0

    def __call__(self, validation_loss, validation_acc):
        if (self.best_loss - validation_loss ) >= self.min_delta:
            self.counter +=1
            self.best_loss = validation_loss
            if self.counter >= self.patience:  
                self.early_stop = True
        else:
            self.counter = 0
                
        if (self.best_acc <= validation_acc):
            self.acc_counter += 1 
            self.best_acc = validation_acc
            if self.acc_counter >= self.patience:  
                self.early_stop = True
        else:
            self.early_stop = 0
    
    
class nn_model(nn.Module):
    def __init__(self, input_dim, hidden_dim, lr):
        
        super().__init__()   
        self.lr = lr
        
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.dropout1 = nn.Dropout(0.1)
        
        self.input_to_hidden = nn.Linear(input_dim, hidden_dim[0])
        
        layers = []
        for i in range(len(hidden_dim) - 1):
            layers.append(nn.Linear(hidden_dim[i], hidden_dim[i + 1]))
            layers.append(nn.ReLU())
            #layers.append(nn.Dropout(0.1))
        self.hidden_layers = nn.Sequential(*layers)
        
        self.last_layer = nn.Linear(hidden_dim[-1], 1)
        
    def forward(self, x):
        x = self.relu(self.input_to_hidden(x))
        #x = self.dropout1(x)
        x = self.hidden_layers(x)
        out = self.sigmoid(self.last_layer(x))
        return out
    
    def set_writer(self, roundom_round, hidden_dim, lr, train_index_num):
        # create tensorboard logger
        log_dir = "./tensorboard_logs" + "/random_" + str(roundom_round) + "_" + str(hidden_dim) + "_" + str(train_index_num) + "_" + str(lr)
        self.writer = SummaryWriter(log_dir=log_dir)
    
    def optimizer(self):
        return torch.optim.Adam(self.parameters(), lr = self.lr, weight_decay=1e-5)

    def reset_weights(self, model):
        '''
        Try resetting model weights to avoid
        weight leakage after each cross validaton training
        '''
        for layer in model.children():
            if hasattr(layer, 'reset_parameters'):
                print(f'Reset trainable parameters of layer = {layer}')
                layer.reset_parameters()
            self.reset_weights(layer)
    
    def fit(self, train_loader, test_loader, epochs):    
        
        optimizer = self.optimizer()
        #early_stopping = EarlyStopping(patience = 5)
        
        for epoch in range(epochs):
            self.train()
            for X_train, y_train in train_loader:
                y_train = y_train.unsqueeze(1)
           
                # Clear gradients
                optimizer.zero_grad()

                y_hat = self.forward(X_train)
                
                # get loss for the predicted output
                loss_fn = nn.BCELoss()
                loss = loss_fn(y_hat, y_train.float())
                
                loss.backward()
                optimizer.step()

                train_acc = (y_hat.round() == y_train).float().mean()                

                print('epoch {}, train loss {}, train acc {}'.format(epoch, loss.item(), train_acc))
                print("Number of positive genes predicted: {}, and actual positives in training data: {}".format(sum(y_hat.round()), sum(y_train)))
                #self.writer.add_scalar('Loss/Train', loss.item(),  epoch)
                
            # Validation step
            print('Validating')
            self.eval()
            total = 0
            correct = 0
            validation_loss = 0 
            with torch.no_grad():
                for X_test, y_test in test_loader:
                    y_test = y_test.unsqueeze(1)
                    y_test_hat = self.forward(X_test)
                    total += y_test.size(0)
                    correct += (y_test_hat.round() == y_test).sum().item()
                    batch_loss = loss_fn(y_test_hat, y_test.float()).item()
                    validation_loss += batch_loss
                validation_loss /= len(test_loader)
                validation_acc = correct / total
                
                print("maximum validation set prediction: " , max(y_test_hat))
                #self.writer.add_scalar('Loss/Validation', validation_loss,  epoch)
                print('epoch {}, validation loss {}, validation acc {}\n\n'.format(epoch, validation_loss, validation_acc))
            
        self.writer.flush()
        self.writer.close()
            
