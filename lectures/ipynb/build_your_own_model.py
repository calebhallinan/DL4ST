# Test your model building skills by building your own Neural Network!

# You will need to download the packages below. The best way to do this is to create a conda env or jupyter env

# 1. Create a new conda environment named pytorch4st with Python 3.10 (recommended for PyTorch)
!conda create -n pytorch4st python=3.10

# 2. Activate the environment
!conda activate pytorch4st

# 3. Install PyTorch
!pip install pytorch torchvision torchaudio

# 4. Install matplotlib
!pip install matplotlib


####################################################################################################################



# read in packages
import torch
from torch import nn
from torch.utils.data import DataLoader # wraps an iterable around the Dataset
from torchvision import datasets, transforms #  stores the samples and their corresponding labels
from torch.utils.data import Dataset
import matplotlib.pyplot as plt


####################################################################################################################


# We will use the MNIST dataset

# 1. Define transform (convert images to tensors and normalize)
transform = transforms.Compose([
    transforms.ToTensor(),  # convert image to tensor
    transforms.Normalize((0.5,), (0.5,))  # normalize to mean=0.5, std=0.5
])

# 2. Download and load the training and test datasets
train_dataset = datasets.MNIST(root="data", train=True, download=True, transform=transform)
test_dataset = datasets.MNIST(root="data", train=False, download=True, transform=transform)

# 3. Create data loaders
train_loader = DataLoader(train_dataset, batch_size=64, shuffle=False)
test_loader = DataLoader(test_dataset, batch_size=64, shuffle=False)

# 4. Inspect a batch
images, labels = next(iter(train_loader))
print(f"Batch shape: {images.shape}, Labels shape: {labels.shape}")



# plot an example to confirm
plt.imshow(images[0].squeeze(), cmap="gray")
plt.title(f"Label: {labels[0]}")
plt.show()



####################################################################################################################


### Now it's your turn! create a class called "MyNeuralNetwork" and train it to predict the MNIST dataset labels. 
# Feel free to use a basic NN, CNN, bring in a pretrained mdoel, or whatever you want to try! 
# Your goal is to try and get the best accuracy you can before class ends. I would also recommend 
# trying to use techniques we have discussed in class such as dropout, batch normalization, number 
# of layers and nodes, etc. 

# Remember to use a training/validation loss curve to see if your model is overfitting.

# You can use previous lecture code to help build the model. However, for this exercise 
# please do not use ChatGPT! The goal is to get you famaliar with the syntax of pytorch, 
# so make sure to write it out yourself.


