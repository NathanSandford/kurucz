import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch.autograd import Variable
import matplotlib.pyplot as plt

description = "Trains NN on normalized spectra"

'''
Defaults
'''
kbase = Path('/global/scratch/nathan_sandford/kurucz')
krun_base = 'kurucz_run'
kout = 'kurucz_out'
spec_dir = 'synthetic_spectra'
NN_dir = 'NN_results'
element_labels = ['Fe', 'Ca', 'Ni', 'Si', 'Ti', 'Co', 'Mg', 'Cr', 'Na', 'K', 'C', 'N', 'Sc', 'Ba',
                  'Al', 'Y']
label_names = ['Teff', 'logg', 'v_micro'] + element_labels
frac_train = 0.75

'''
Parse Args
'''
parser = argparse.ArgumentParser(description=description)
parser.add_argument("spec_name", metavar="spectra", help="File containing spectra")
parser.add_argument("--kbase", "-kb",
                    help=f"Kurucz Base directory (default: {kbase})")
parser.add_argument("--kout", "-ko",
                    help=f"Kurucz Output Directory (default: {kout})")
parser.add_argument("--spec_dir",
                    help=f"Spectra Directory (default: {spec_dir})")
parser.add_argument("--NN_dir", "-nn",
                    help="Neural Network Directory (default: {NN_dir})")
parser.add_argument("--train_fraction", "-fraction",
                    help="Fraction of spectra to train on (default: {frac_train})")
args = parser.parse_args()

if args.kbase:
    kbase = Path(args.kbase)
assert kbase.is_dir(), f'Base directory {kbase} does not exist'

if args.kout:
    kout = kbase.joinpath(args.kout)
else:
    kout = kbase.joinpath(kout)
assert kout.is_dir(), f'Output directory {kout} does not exist'

if args.spec_dir:
    spec_dir = kout.joinpath(args.spec_dir)
else:
    spec_dir = kout.joinpath(spec_dir)
assert spec_dir.is_dir(), f'Spectra directory {spec_dir} does not exist'

if args.NN_dir:
    NN_dir = kout.joinpath(args.NN_dir)
else:
    NN_dir = kout.joinpath(NN_dir)
assert NN_dir.is_dir(), f'Spectra directory {NN_dir} does not exist'

if args.spec_name[-3:] == '.h5':
    input_file = spec_dir.joinpath(args.spec_name)
    output_file = NN_dir.joinpath(args.spec_name[:-3] + '_NN.npz')
else:
    input_file = spec_dir.joinpath(args.spec_name + '.h5')
    output_file = NN_dir.joinpath(args.spec_name + '_NN.npz')

if args.training_fraction:
    frac_train = args.training_fraction

'''
Load Spectra
'''
print('Restoring Spectra')
spectra = pd.read_hdf(input_file, 'spectra').values.T
labels = pd.read_hdf(input_file, 'labels').loc[label_names].values.T

n_train = int(frac_train * spectra.shape[0])
training_spectra = spectra[:n_train,:]
training_labels = labels[:n_train,:]
validation_spectra = spectra[n_train:,:]
validation_labels = labels[n_train:,:]


'''
Define NN
'''
def neural_net(training_labels, training_spectra, validation_labels, validation_spectra,\
             num_neurons = 300, num_steps=1e5, learning_rate=0.001):

    '''
    Training neural networks to emulate spectral models
    
    training_labels has the dimension of [# training spectra, # stellar labels]
    training_spectra has the dimension of [# training spectra, # wavelength pixels]
    The validation set is used to independently evaluate how well the neural networks
    are emulating the spectra. If the networks overfit the spectral variation, while 
    the loss function will continue to improve for the training set, but the validation 
    set should show a worsen loss function.
    The training is designed in a way that it always returns the best neural networks
    before the networks start to overfit (gauged by the validation set).
    num_neurons = number of neurons per hidden layer in the neural networks. 
    We assume a 2 hidden-layer neural networks.
    Increasing this number increases the complexity of the network, which can 
    capture a more subtle variation of flux as a function of stellar labels, but
    increasing the complexity could also lead to overfitting. And it is also slower 
    to train with a larger network.
    num_steps = how many steps to train until convergence. 
    1e5 is good for the specific NN architecture and learning I used by default, 
    but bigger networks take more steps, and decreasing the learning rate will 
    also change this. You can get a sense of how many steps are needed for a new 
    NN architecture by plotting the loss function evaluated on both the training set 
    and a validation set as a function of step number. It should plateau once the NN 
    is converged.  
    learning_rate = step size to take for gradient descent
    This is also tunable, but 0.001 seems to work well for most use cases. Again, 
    diagnose with a validation set if you change this. 
    
    returns:
        training loss and validation loss (per 1000 steps)
        the codes also outputs a numpy saved array ""NN_normalized_spectra.npz" 
        which can be imported and substitute the default neural networks (see tutorial)
    '''
    
    # run on cuda
    dtype = torch.cuda.FloatTensor
    torch.set_default_tensor_type('torch.cuda.FloatTensor')

    # scale the labels, optimizing neural networks is easier if the labels are more normalized
    x_max = np.max(training_labels, axis = 0)
    x_min = np.min(training_labels, axis = 0)
    x = (training_labels - x_min)/(x_max - x_min) - 0.5
    x_valid = (validation_labels-x_min)/(x_max-x_min) - 0.5

    # dimension of the input
    dim_in = x.shape[1]

    # dimension of the output
    num_pixel = training_spectra.shape[1]

    # define neural networks
    model = torch.nn.Sequential(
        torch.nn.Linear(dim_in, num_neurons),
        torch.nn.Sigmoid(),
        torch.nn.Linear(num_neurons, num_neurons),
        torch.nn.Sigmoid(),
        torch.nn.Linear(num_neurons, num_pixel)
    )
    model.cuda()

    # assume L2 loss
    loss_fn = torch.nn.MSELoss(reduction='mean')
    
    # make pytorch variables
    x = Variable(torch.from_numpy(x)).type(dtype)
    y = Variable(torch.from_numpy(training_spectra), requires_grad=False).type(dtype)
    x_valid = Variable(torch.from_numpy(x_valid)).type(dtype)
    y_valid = Variable(torch.from_numpy(validation_spectra), requires_grad=False).type(dtype)

    # weight_decay is for regularization. Not required, but one can play with it. 
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay = 0)

    # train the neural networks
    t = 0
    current_loss = np.inf
    training_loss =[]
    validation_loss = []
    while t < num_steps:
        y_pred = model(x)
        loss = loss_fn(y_pred, y)*1e4
        y_pred_valid = model(x_valid)
        loss_valid = loss_fn(y_pred_valid, y_valid)*1e4
    
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        t += 1

        # print loss function to monitor
        if t % 1000 == 0:
            loss_data = loss.data.item()
            loss_valid_data = loss_valid.data.item()
            training_loss.append(loss_data)
            validation_loss.append(loss_valid_data)
            print('Step ' + str(t) \
                  + ': Training set loss = ' + str(int(loss_data*1000.)/1000.) \
                  + ' / Validation set loss = ' + str(int(loss_valid_data*1000.)/1000.))
 
            # record the weights and biases if the validation loss improves
            if loss_valid_data < current_loss:
                current_loss = loss_valid
                model_numpy = []
                for param in model.parameters():
                    model_numpy.append(param.data.cpu().numpy())
 
    # extract the weights and biases
    w_array_0 = model_numpy[0]
    b_array_0 = model_numpy[1]
    w_array_1 = model_numpy[2]
    b_array_1 = model_numpy[3]
    w_array_2 = model_numpy[4]
    b_array_2 = model_numpy[5]

    # save parameters and remember how we scaled the labels
    np.savez(output_file,\
             w_array_0 = w_array_0,\
             w_array_1 = w_array_1,\
             w_array_2 = w_array_2,\
             b_array_0 = b_array_0,\
             b_array_1 = b_array_1,\
             b_array_2 = b_array_2,\
             x_max=x_max,\
             x_min=x_min,\
             training_loss = training_loss,\
             validation_loss = validation_loss)

    return training_loss, validation_loss

print('Training Spectra')
training_loss, validation_loss = neural_net(training_labels, training_spectra,\
                                            validation_labels, validation_spectra,\
                                            num_neurons = 300, num_steps=1e5, learning_rate=0.001)
print('Training Complete!')
