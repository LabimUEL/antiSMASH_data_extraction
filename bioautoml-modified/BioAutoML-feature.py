import warnings
warnings.filterwarnings(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore')
import pandas as pd
import argparse
import subprocess
import shutil
import sys
import os.path
import time
import lightgbm as lgb
from catboost import CatBoostClassifier
from sklearn.metrics import balanced_accuracy_score
# from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import make_scorer
from sklearn.model_selection import cross_val_score
from sklearn.metrics import f1_score
from hyperopt import hp, fmin, tpe, STATUS_OK, Trials

# Testing
# python bioautoml-modified/BioAutoML-feature.py --antismash_train data/processed --n_cpu -1 --output result

def objective_rf(space):

	"""Automated Feature Engineering - Objective Function - Bayesian Optimization"""

	index = list()
	descriptors = {'CDS_motifs': list(range(0, 84)), 'CDS_smCOG': list(range(84, 324)),
				   'PFAM_desc': list(range(324, 1879)), 'PFAM_GO': list(range(1879, 2342)),
				   'PFAM_ID': list(range(2342, 3993))}

	for descriptor, ind in descriptors.items():
		if int(space[descriptor]) == 1:
			index = index + ind

	# print(space)
	# print(index)

	x = df_x.iloc[:, index]

	if int(space['Classifier']) == 0:
		# if len(fasta_label_train) > 2:
		# 	model = AdaBoostClassifier(random_state=63)
		# else:
		model = CatBoostClassifier(n_estimators=500,
									thread_count=n_cpu, nan_mode='Max',
									logging_level='Silent', random_state=63)
	elif int(space['Classifier']) == 1:
		model = RandomForestClassifier(n_estimators=500, n_jobs=n_cpu, random_state=63)
	else:
		model = lgb.LGBMClassifier(n_estimators=500, n_jobs=n_cpu, random_state=63)

	# if len(fasta_label_train) > 2:
	# 	score = make_scorer(f1_score, average='weighted')
	# else:
	score = make_scorer(balanced_accuracy_score)

	kfold = StratifiedKFold(n_splits=10, shuffle=True)

	if not x.empty:
		metric = cross_val_score(model,
								x,
								labels_y,
								cv=kfold,
								scoring=score,
								n_jobs=n_cpu).mean()
	else:
		metric = -10e3
	
	return {'loss': -metric, 'status': STATUS_OK}

def feature_engineering(estimations, train, train_labels, test, foutput):

	"""Automated Feature Engineering - Bayesian Optimization"""

	global df_x, labels_y

	print('Automated Feature Engineering - Bayesian Optimization')

	df_x = pd.read_csv(train)
	labels_y = pd.read_csv(train_labels)

	if test != '':
		df_test = pd.read_csv(test)

	path_bio = foutput + '/best_descriptors'
	if not os.path.exists(path_bio):
		os.mkdir(path_bio)

	param = {'CDS_motifs': [0, 1], 'CDS_smCOG': [0, 1],
			 'PFAM_desc': [0, 1], 'PFAM_GO': [0, 1], 'PFAM_ID': [0, 1],
			 'Classifier': [0, 1, 2]}

	space = {'CDS_motifs': hp.choice('CDS_motifs', [0, 1]),
			 'CDS_smCOG': hp.choice('CDS_smCOG', [0, 1]),
			 'PFAM_desc': hp.choice('PFAM_desc', [0, 1]),
			 'PFAM_GO': hp.choice('PFAM_GO', [0, 1]),
			 'PFAM_ID': hp.choice('PFAM_ID', [0, 1]),
			 'Classifier': hp.choice('Classifier', [0, 1, 2])}

	trials = Trials()
	best_tuning = fmin(fn=objective_rf,
				space=space,
				algo=tpe.suggest,
				max_evals=estimations,
				trials=trials)

	index = list()
	descriptors = {'CDS_motifs': list(range(0, 84)), 'CDS_smCOG': list(range(84, 324)),
				   'PFAM_desc': list(range(324, 1879)), 'PFAM_GO': list(range(1879, 2342)),
				   'PFAM_ID': list(range(2342, 3993))}

	for descriptor, ind in descriptors.items():
		result = param[descriptor][best_tuning[descriptor]]
		if result == 1:
			index = index + ind

	classifier = param['Classifier'][best_tuning['Classifier']]

	btrain = df_x.iloc[:, index]
	path_btrain = path_bio + '/best_train.csv'
	btrain.to_csv(path_btrain, index=False, header=True)

	if test != '':
		btest = df_test.iloc[:, index]
		path_btest = path_bio + '/best_test.csv'
		btest.to_csv(path_btest, index=False, header=True)
	else:
		btest, path_btest = '', ''

	return classifier, path_btrain, path_btest, btrain, btest


def feature_extraction(antismash_train, antismash_test, foutput):

	"""Load features from antiSMASH."""

	path = foutput + '/feat_extraction'
	path_results = foutput

	try:
		shutil.rmtree(path)
		shutil.rmtree(path_results)
	except OSError as e:
		print("Error: %s - %s." % (e.filename, e.strerror))
		print('Creating Directory...')

	if not os.path.exists(path_results):
		os.mkdir(path_results)

	if not os.path.exists(path):
		os.mkdir(path)

	print('Loading antiSMASH features...')

	if antismash_train:
		datasets = {}

		for filename in os.listdir(antismash_train):
			if filename.endswith('.csv'):
				key = filename[:-4]  # Remove the ".csv" extension
				filepath = os.path.join(antismash_train, filename)
				datasets[key] = pd.read_csv(filepath)

		merged_df = pd.DataFrame()
		for key, df in datasets.items():
			if merged_df.empty:
				merged_df = df
			else:
				merged_df = merged_df.join(df.set_index(['nome', 'acc', 'label']), on=['nome', 'acc', 'label'])

		# Reset the index of the merged DataFrame
		X_train = merged_df.reset_index(drop=True)

		y_train = X_train.pop('label')
		X_train.pop('nome')
		X_train.pop('acc')

		ftrain = path + '/ftrain.csv'
		X_train.to_csv(ftrain, index=False)
		flabeltrain = path + '/flabeltrain.csv'
		y_train.to_csv(flabeltrain, index=False, header=True)
		
		ftest, flabeltest = '', ''
		if antismash_test:
			ftest = path + '/ftest.csv'
			flabeltest = path + '/flabeltest.csv'

	return ftrain, flabeltrain, ftest, flabeltest

##########################################################################
##########################################################################


if __name__ == '__main__':
	print('\n')
	print('###################################################################################')
	print('###################################################################################')
	print('##########         BioAutoML- Automated Feature Engineering             ###########')
	print('##########              Author: Robson Parmezan Bonidia                 ###########')
	print('##########         WebPage: https://bonidia.github.io/website/          ###########')
	print('###################################################################################')
	print('###################################################################################')
	print('\n')
	parser = argparse.ArgumentParser()
	parser.add_argument('-antismash_train', '--antismash_train', help='antiSMASH features train directory, e.g., train/')
	parser.add_argument('-antismash_test', '--antismash_test', help='antiSMASH features test directory, e.g., test/')
	parser.add_argument('-estimations', '--estimations', default=50, help='number of estimations - BioAutoML - default = 50')
	parser.add_argument('-n_cpu', '--n_cpu', default=1, help='number of cpus - default = 1')
	parser.add_argument('-output', '--output', help='results directory, e.g., result')

	args = parser.parse_args()
	antismash_train = args.antismash_train
	antismash_test = args.antismash_test
	estimations = int(args.estimations)
	n_cpu = int(args.n_cpu)
	foutput = str(args.output)

	# for train in antismash_train:
	# 	if os.path.exists(train) is True:
	# 		print('Train - %s: Found File' % train)
	# 	else:
	# 		print('Train - %s: File not exists' % train)
	# 		sys.exit()

	# if antismash_test:
	# 	for test in antismash_test:
	# 		if os.path.exists(test) is True:
	# 			print('Test - %s: Found File' % test)
	# 		else:
	# 			print('Test - %s: File not exists' % test)
	# 			sys.exit()

	start_time = time.time()

	ftrain, ftrain_labels, \
		ftest, ftest_labels = feature_extraction(antismash_train, antismash_test, foutput)

	classifier, path_train, path_test, train_best, test_best = \
		feature_engineering(estimations, ftrain, ftrain_labels, ftest, foutput)

	cost = (time.time() - start_time) / 60
	print('Computation time - Pipeline - Automated Feature Engineering: %s minutes' % cost)

	# if len(fasta_label_train) > 2:
	# 	subprocess.run(['python', 'BioAutoML-multiclass.py', '-train', path_train,
	# 					 '-train_label', ftrain_labels, '-test', path_test,
	# 					 '-test_label', ftest_labels, '-test_nameseq',
	# 					 fnameseqtest, '-nf', 'True', '-classifier', str(classifier),
	# 					 '-n_cpu', str(n_cpu), '-output', foutput])
	# else:
	subprocess.run(['python', 'bioautoml-modified/BioAutoML-binary.py', '-train', path_train,
						'-train_label', ftrain_labels, '-test', path_test, '-test_label',
						ftest_labels, '-nf', 'True', '-classifier', str(classifier), 
						'-tuning', 'True', '-n_cpu', str(n_cpu), '-output', foutput])

##########################################################################
##########################################################################
