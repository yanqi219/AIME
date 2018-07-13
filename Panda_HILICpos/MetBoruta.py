import pandas as pd
import numpy as np
import os as os
from sklearn.ensemble import RandomForestClassifier
from boruta import BorutaPy

def MetBoruta(dir, featuretable, labeltable, max_depth=7, perc=100, iter=100, n_jobs=-1, class_weight='balanced', n_estimators='auto', random_state=1, verbose=2):
    os.chdir(dir)

    features = pd.read_table(featuretable, sep='\t')
    features.head(5)
    label = pd.read_table(labeltable, sep='\t')
    label.head(5)
    print('The dimension of the feature table is:', features.shape)
    print('The dimension of the label table is:', label.shape)

    # Data pre-processing
    i = ["met_"+str(i) for i in range(1, features.__len__()+1)]
    features.index = i

    link = pd.DataFrame(features.loc[:, ['mz', 'time']])
    features = features.drop(['mz', 'time'], axis=1)
    features = features.T

    label.index = label.iloc[:, 0]
    label = label.drop('SampleID', axis=1)

    features_array = np.array(features)
    label_array = np.array(label)

    # Start RF method

    rf = RandomForestClassifier(n_jobs=n_jobs, class_weight=class_weight, max_depth=max_depth)
    # define Boruta feature selection method
    feat_selector = BorutaPy(rf, n_estimators=n_estimators, verbose=verbose, random_state=random_state, perc=perc, max_iter=iter)
    # find all relevant features
    feat_selector.fit(features_array, label_array)

    # check selected features
    features_flag = pd.DataFrame(feat_selector.support_)
    features_flag = features_flag.T
    features_flag.columns = i
    features_flag.index = ["Flag"]

    # check ranking of features
    features_ranking = feat_selector.ranking_

    # call transform() on X to filter it down to selected features
    featurea_filtered = feat_selector.transform(features_array)

    # Get significant features
    features_sig = features.append(features_flag)
    features_sig.to_csv('Boruta_output.csv', sep='\t')

MetBoruta(dir='C:\Users\QiYan\Dropbox\AIME\Panda_C18neg\C18_Controls_ExpoUnexpo\PANDA_input',
          featuretable='C18_residuals_control_expo_unexpo.txt',
          labeltable='C18_residuals_classlabels_control_expo_unexpo.txt',
          max_depth=7,
          perc=95,
          iter=150)