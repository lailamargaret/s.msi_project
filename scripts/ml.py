#This script does the machine learning process and saves its results to multiple files
#Code modified from Google's Machine Learning Crash Course
from constants import _MSI_LOCI
import math

import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics
import tensorflow as tf
from tensorflow.python.data import Dataset

tf.logging.set_verbosity(tf.logging.ERROR)
pd.options.display.max_rows = 10
pd.options.display.float_format = '{:.1f}'.format

def feature_list(infile):
	"""
	Creates a list of features given an input file
	"""
        feature_list = []
        with open(infile, 'r') as f:
                for line in f:
                        feature_list.append(line.replace('\n', ''))
        return feature_list

def preprocess_features(input_df, feature_list):
  """Prepares input features from a data set.
    Args:
    input_df: A Pandas DataFrame expected to contain features
    Returns:
    A DataFrame that contains the features to be used for the model, including
    synthetic features (if created).
  """
  
  selected_features = input_df[feature_list]

  processed_features = selected_features.copy()
 
  # Create a synthetic feature if desired.
  
  return processed_features

def preprocess_targets(input_df):
  """Prepares target features (i.e., labels) from data set.

  Args:
    input_df: A Pandas DataFrame expected to contain data
  Returns:
    A DataFrame that contains the target feature.
  """
  output_targets = pd.DataFrame()
   
  # median_house_value is above a set threshold.
  output_targets["msi_status"] = input_df["msi_status"]
  return output_targets


def construct_feature_columns(input_features):
  """Construct the TensorFlow Feature Columns.

  Args:
    input_features: The names of the numerical input features to use.
  Returns:
    A set of feature columns
  """
  return set([tf.feature_column.numeric_column(my_feature)
              for my_feature in input_features])

def my_input_fn(features, targets, batch_size=10, shuffle=True, num_epochs=None):
    """Trains a linear classification model.
  
    Args:
      features: pandas DataFrame of features
      targets: pandas DataFrame of targets
      batch_size: Size of batches to be passed to the model
      shuffle: True or False. Whether to shuffle the data.
      num_epochs: Number of epochs for which data should be repeated. None = repeat indefinitely
    Returns:
      Tuple of (features, labels) for next data batch
    """
    
    # Convert pandas data into a dict of np arrays.
    features = {key:np.array(value) for key,value in dict(features).items()}                                            
 
    # Construct a dataset, and configure batching/repeating.
    ds = Dataset.from_tensor_slices((features,targets)) # warning: 2GB limit
    ds = ds.batch(batch_size).repeat(num_epochs)
    
    # Shuffle the data, if specified.
    if shuffle:
      ds = ds.shuffle(10000)
    
    # Return the next batch of data.
    features, labels = ds.make_one_shot_iterator().get_next()
    return features, labels

def train_linear_classifier_model(
    learning_rate,
    steps,
    batch_size,
    regularization_strength,
    training_examples,
    training_targets,
    validation_examples,
    validation_targets):
  """Trains a linear classification model.
  
  In addition to training, this function also prints training progress information,
  as well as a plot of the training and validation loss over time.
  
  Args:
    learning_rate: A `float`, the learning rate.
    steps: A non-zero `int`, the total number of training steps. A training step
      consists of a forward and backward pass using a single batch.
    batch_size: A non-zero `int`, the batch size.
    regularization_strength: A 'float', the strength of l1 regularization (currently does not work?)
    training_examples: A `DataFrame` containing one or more columns from
      dataframe to use as input features for training.
    training_targets: A `DataFrame` containing exactly one column from
      dataframe to use as target for training.
    validation_examples: A `DataFrame` containing one or more columns from
      dataframe to use as input features for validation.
    validation_targets: A `DataFrame` containing exactly one column from
      dataframe use as target for validation.
      
  Returns:
    A `LinearClassifier` object trained on the training data.
  """

  periods = 20
  steps_per_period = steps / periods
  
  # Create a linear classifier object.
  my_optimizer = tf.train.FtrlOptimizer(learning_rate=learning_rate, l1_regularization_strength=regularization_strength)
  my_optimizer = tf.contrib.estimator.clip_gradients_by_norm(my_optimizer, 5.0)
  linear_classifier = tf.estimator.LinearClassifier(
      feature_columns=construct_feature_columns(training_examples),
      optimizer=my_optimizer
  )
  
  # Create input functions.
  training_input_fn = lambda: my_input_fn(training_examples, 
                                          training_targets["msi_status"], 
                                          batch_size=batch_size)
  predict_training_input_fn = lambda: my_input_fn(training_examples, 
                                                  training_targets["msi_status"], 
                                                  num_epochs=1, 
                                                  shuffle=False)
  predict_validation_input_fn = lambda: my_input_fn(validation_examples, 
                                                    validation_targets["msi_status"], 
                                                    num_epochs=1, 
                                                    shuffle=False)
  
  # Train the model, but do so inside a loop so that we can periodically assess
  # loss metrics.
  print("Training model...")
  print("LogLoss (on training data):")
  training_log_losses = []
  validation_log_losses = []
  for period in range (0, periods):
    # Train the model, starting from the prior state.
    linear_classifier.train(
        input_fn=training_input_fn,
        steps=steps_per_period
    )
    # Take a break and compute predictions.    
    training_probabilities = linear_classifier.predict(input_fn=predict_training_input_fn)
    training_probabilities = np.array([item['probabilities'] for item in training_probabilities])
    
    validation_probabilities = linear_classifier.predict(input_fn=predict_validation_input_fn)
    validation_probabilities = np.array([item['probabilities'] for item in validation_probabilities])
    
    training_log_loss = metrics.log_loss(training_targets, training_probabilities)
    validation_log_loss = metrics.log_loss(validation_targets, validation_probabilities)
    # Occasionally print the current loss.
    print("  period %02d : %0.9f" % (period, training_log_loss))
    # Add the loss metrics from this period to our list.
    training_log_losses.append(training_log_loss)
    validation_log_losses.append(validation_log_loss)
  print("Model training finished.")
  
  # Output a graph of loss metrics over periods.
  plt.ylabel("LogLoss")
  plt.xlabel("Periods")
  plt.title("LogLoss vs. Periods")
  plt.tight_layout()
  plt.plot(training_log_losses, label="training")
  plt.plot(validation_log_losses, label="validation")
  plt.legend()
  plt.savefig('/home/upload/msi_project/ML/histogram_features/top_9/%d_%f_%d_%f_loss.png' % (steps, learning_rate, batch_size, regularization_strength))
  plt.clf()
  return linear_classifier

def run_test(params):
  """
  Calls train_linear_classifier given the list that contains the input
  Creates ROC plot
  """
  linear_classifier = train_linear_classifier_model(
            learning_rate=params[0],
            steps=params[1],
            batch_size=params[2],
            regularization_strength=params[3],
            training_examples=training_examples,
            training_targets=training_targets,
            validation_examples=validation_examples,
            validation_targets=validation_targets)


  predict_validation_input_fn = lambda: my_input_fn(validation_examples, 
                                                  validation_targets["msi_status"], 
                                                  num_epochs=1, 
                                                  shuffle=False)


  evaluation_metrics = linear_classifier.evaluate(input_fn=predict_validation_input_fn)

  print("AUC on the validation set: %0.2f" % evaluation_metrics['auc'])
  print("Accuracy on the validation set: %0.2f" % evaluation_metrics['accuracy'])

  validation_probabilities = linear_classifier.predict(input_fn=predict_validation_input_fn)

  # Get just the probabilities for the positive class.
  validation_probabilities = np.array([item['probabilities'][1] for item in validation_probabilities])

  false_positive_rate, true_positive_rate, thresholds = metrics.roc_curve(
      validation_targets, validation_probabilities)
  plt.plot(false_positive_rate, true_positive_rate, 'b', label = 'AUC = %0.2f' % evaluation_metrics['auc'])
  plt.plot([0, 1], [0, 1], 'r--')
  plt.xlim([0, 1])
  plt.ylim([0, 1])
  plt.ylabel('True Positive Rate (Sensitivity)')
  plt.xlabel('False Positive Rate (1 - Specificity)')
  _ = plt.legend(loc=2)
  plt.savefig('/home/upload/msi_project/ML/histogram_features/top_9/%d_%f_%d_%f_roc.png' % (params[1], params[0], params[2], params[3]))
  plt.clf()
  with open('/home/upload/msi_project/ML/histogram_features/top_9/hyperparameters.txt', 'a') as f:
    f.write(str(in_learning_rate) + '\t' +
            str(in_steps) + '\t' +		 
            str(in_batch_size) + '\t' +
            str(regularization_strength) + '\t' +
            str(evaluation_metrics['auc']) + '\t' + 
            str(evaluation_metrics['accuracy']) + '\n'
            )

  names = linear_classifier.get_variable_names()
  values = [linear_classifier.get_variable_value(x) for x in names]


  return names, values

#--------------------------------------------------------------------- Main -------------------------------------------------------------#

#Produce feature list for use in preprocessing step
feature_list = feature_list('/home/upload/msi_project/ML/histogram_features/top_9/top_9_histogram_features.txt')
#sets for preprocessing steps
training_set = pd.read_csv("/home/upload/msi_project/ML/histogram_features/top_9/training_set_top_9_lengths_full.txt", sep="\t")
validation_set = pd.read_csv("/home/upload/msi_project/ML/histogram_features/top_9/validation_set_top_9_lengths_full.txt", sep="\t")
# Preprocess training examples and targets.
training_examples = preprocess_features(training_set, feature_list)
training_targets = preprocess_targets(training_set)
# Preprocess validation examples and targets
validation_examples = preprocess_features(validation_set, feature_list)
validation_targets = preprocess_targets(validation_set)
 
#set the hyperparameters
in_learning_rate = 0.001
in_steps = 200000
in_batch_size = 10
regularization_strength = 1
#consolidate hyperparameters
params=[in_learning_rate, in_steps, in_batch_size, regularization_strength]

#train the model
names, values = run_test(params)

#save weights information
weights = zip(names, values)
outfile = '/home/upload/msi_project/ML/histogram_features/top_9/%d_%f_%d_%d_model_weights.txt' % (in_steps, in_learning_rate, in_batch_size, regularization_strength)
#report the model-derrived weights to a file
with open(outfile, 'w') as f:
  for weight in weights:
    if weight[0].endswith('weights'):
    	f.write(str(weight[0]) + '\t' + str(weight[1]) + '\n')




