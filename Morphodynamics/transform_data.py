#from sklearn import preprocessing
import numpy as np

#import features

def import_train_set(train_file_name = 'AllResults.txt'):
	featurelist=[]

	with open(train_file_name,'r') as infile:
		for line in infile:
			featurelist.append(line.strip())

	# so now, featurelist[1] has names of things in form 'Area, MajorAxisLength, ... Class'
	FeatureNames = [x.strip() for x in featurelist[0].split(',')]
	#FeatureNames has form ['Area','MajorAxisLength',....'Class'] which is what I wanted

	AllData = [[float(x.strip()) for x in featurelist[i].split(',')] for i in range(1,len(featurelist))]

	# Data is in form [[1,2,3,....0.0],[3,3,1,...0.0],...[5,3,1,...0.0]], the last input is the class.

	classes = [int(i[-1]) for i in AllData]

	#classes contains the class number from which the data is from

	#want to delete target from AllData.

	X = [i[0:-1] for i in AllData]

	#X has form similar to Data. So when we reshape, we want the output to be
	# X = array([[0,1,2,...]
	#            [1,2,3,...]])

	Data = np.asarray(X,order = 'F')

	#this has the right form, is uses fortran column-major style memory representation vs row major C-style
	# the notation is scientific, where iris data set looks like a float. CHECKED: Both are type numpy.float64
	#both have same indexing calls, so I think we're in business.

	Target = np.asarray(classes) #looks exactly correct, or at least like iris data set target.
	return (Data, Target)

########################################################################
#for training purposes, the number of samples in data must be divisible by 256
def Trim_Train_Data(Data, Target, max_length=None,balance=False):
	####
	#Inputs: Data is numpy array with N samples (rows) and M measures (cols)
	# Target is 1xN samples with ground truth
	# max_length defines maximum length of training data. Should be divisible by 256, might want to code that...
	# balance is boolean if you wish to have same number of samples in each class.
	print('Class lengths are = ',[sum(Target==i) for i in set(Target)])
	if not balance:
		if np.shape(Data)[0] / 256 != np.round(np.shape(Data)[0] / 256) or max_length<np.shape(Data)[0]:
			print("Trimming data for training purposes...")
			if not max_length:
				max_length=256*(np.floor(np.shape(Data)[0] / 256))
			else:
				if max_length / 256 != np.round(max_length / 256):
					#must make it divisible by 256
					max_length = int(np.floor(max_length/256)*256)
					print("Your given max_length was not divisible by 256. New max length is = %d" % max_length)
			#determine percentages of each class.
			cs = np.unique(Target)
			ps=np.zeros(shape=(1,len(cs)))
			ps=ps[0]
			rows_to_take=np.array([])
			for i in range(len(cs)):
				ps[i]=np.sum(Target==cs[i])/len(Target)
				goodrows=np.where(Target==cs[i])[0]
				rows_to_take=np.append(rows_to_take,goodrows[0:int(np.floor(ps[i]*max_length))])

			ad_row=0
			class_ind=0
			while len(rows_to_take)!=max_length:
				#need to supplament.
				goodrows=np.where(Target==cs[class_ind])[0]
				rows_to_take=np.append(rows_to_take,goodrows[int(np.floor(ps[class_ind]*max_length))+1+ad_row])
				class_ind=class_ind+1
				if class_ind>len(cs):
					class_ind=0
					ad_row=ad_row+1
			rows_to_take=rows_to_take.astype(int)
			X_train_scaled = Data[rows_to_take,:]
			Y_train = Target[rows_to_take]
			print("Complete")
		else:
			X_train_scaled = Data
			Y_train = Target
		print("Final training length = %d" % X_train_scaled.shape[0])
		print('Class lengths after trimming are = ',[sum(Y_train==i) for i in set(Y_train)])
		return (X_train_scaled, Y_train)
	else:
		#determine which has the minimum number of cases.
		cs = np.unique(Target)
		lens=np.zeros((len(cs)))
		for i in range(len(cs)):
			lens[i]=sum(Target==cs[i])

		#randomly sample from each class now that number of samples.
		min_len=int(min(lens))
		rows_to_take=np.array([])
		for i in range(len(cs)):
			possiblerows=np.where(Target==cs[i])[0]
			#now sample without replacement.
			rows_to_take=np.append(rows_to_take,np.random.choice(possiblerows,min_len,replace=False))
		if len(rows_to_take)/256 != np.round(len(rows_to_take)/256) or max_length<len(rows_to_take):
			#trim until correct size.
			if not max_length:
				max_length=256*(np.floor(np.shape(Data)[0] / 256))
			else:
				if max_length / 256 != np.round(max_length / 256):
					#must make it divisible by 256
					max_length = int(np.floor(max_length/256)*256)
					print("Your given max_length was not divisible by 256. New max length is = %d" % max_length)
			#use min_len now to delete entries.
			timearound=0
			pheno=len(cs) #start at the end
			while len(rows_to_take)>max_length:
				#entry to delete is
				#first (min_len-round)*range(1,len(np.unique(Target))+1) -1
				#print("%d entry delete" % (((min_len-timearound)*pheno) - 1))
				rows_to_take=np.delete(rows_to_take,((min_len-timearound)*pheno) - 1)
				pheno=pheno-1
				if pheno<1:
					pheno=len(cs)
					timearound=timearound+1
		rows_to_take=rows_to_take.astype(int)
		X_train_scaled=Data[rows_to_take,:]
		Y_train=Target[rows_to_take]
		print("Final training length = %d" % X_train_scaled.shape[0])
		print('Class lengths after trimming are = ',[sum(Y_train==i) for i in set(Y_train)])
		return (X_train_scaled, Y_train)

#############################REMOVE OUTLIER DATA########################
#How? Do this after scaling the data, then compute a z-score. We'll check the data after that.
def Remove_Outliers(Data, Target):
	#for each class, detect outliers.
	#we'll begin by using z-scoring. This assumes data is described by a Guassian
	#which is why it is vital to do this AFTER scaling the data.
	#I plotted the data, it is absolutely not Gaussian.
	#I tried DBSCAN machine learning algorithm but it is really not helpful.
	#However, the data IS perhaps Gaussian after embedding. We can clean the signal AFTER by sending in
	#the emebedded data in 1, 2, or 3 dimensions and removing points that are beyond a standard deviation.
	#Data is TSNE embedded.
	zscores = np.zeros(np.shape(Data))
	for pheno in np.unique(Target):
		#find rows where phenotype is correct.
		prows = np.where(Target==pheno)[0]
		for dim in range(np.shape(Data)[1]):
			#calculate the mean.
			m = np.mean(Data[prows,dim])
			#calculate std.
			s = np.std(Data[prows,dim])
			for example in range(len(prows)):
				zscores[prows[example],dim] = (Data[prows[example],dim] - m) / s

	#now you calculated the zscores for each element. Apply a threshold
	#good "thumb-rule" thresholds can be: 2.5, 3, 3.5, or more.
	zthresh = 2.5

	zscores = zscores > 2.5

	badrows = [i for i in range(np.shape(zscores)[0]) if zscores[i].any() ]

	Data = np.delete(Data, badrows, axis = 0)
	Target = np.delete(Target, badrows, axis = 0)

	return(Data, Target)






##############################POST AUGMENTATION#########################
def Augment_Size(Data, Target, max_copies=0, s=0.2, balance=False, augment_class=None):
	max_copies = int(max_copies)
	#augment only the copies made by scaling the unit based measures.
	#Measures should go: Area, MjrAxis, MnrAxis, Ecc,ConA,EqD,Sol,Ext,Per,conPer,fiber_length,InscribeR,bleb_len

	#first, determine if we desire class balance.
	if balance:
		#determine which class has maximum number of samples.
		cs = np.unique(Target)
		vals = [sum(Target==cs[i]) for i in cs]
		print('Class %d has max number of samples, increasing other classes via size augmentation' % np.argmax(vals))
		for i in range(len(cs)):
			if i!=np.argmax(vals):
				#determine how many samples need to be made.
				to_make=int(vals[np.argmax(vals)] - vals[i])
				#randomly sample rows from Data with the correct phenotype cs[i]
				possible_rows = np.where(Target==cs[i])[0]
				#sample to_make numbers from possible_rows.
				sampled_rows = np.random.choice(possible_rows,to_make,replace=True)
				newrows = Data[sampled_rows,:]
				size_vary = s * np.random.rand(1,to_make)[0]
				#vary size.
				for v in range(to_make):
					if np.random.rand()<0.5:
						newrows[v,0]=newrows[v,0] + newrows[v,0]*size_vary[v]*size_vary[v]
						newrows[v,1]=newrows[v,1] + newrows[v,1]*size_vary[v]
						newrows[v,2]=newrows[v,2] + newrows[v,2]*size_vary[v]
						newrows[v,4]=newrows[v,4] + newrows[v,4]*size_vary[v]*size_vary[v]
						newrows[v,5]=newrows[v,5] + newrows[v,5]*size_vary[v]
						newrows[v,7]=newrows[v,7] + newrows[v,7]*size_vary[v]
						newrows[v,8]=newrows[v,8] + newrows[v,8]*size_vary[v]
						newrows[v,9]=newrows[v,9] + newrows[v,9]*size_vary[v]
						newrows[v,10]=newrows[v,10] + newrows[v,10]*size_vary[v]
						newrows[v,11]=newrows[v,11] + newrows[v,11]*size_vary[v]
					else:
						newrows[v,0]=newrows[v,0] - newrows[v,0]*size_vary[v]*size_vary[v]
						newrows[v,1]=newrows[v,1] - newrows[v,1]*size_vary[v]
						newrows[v,2]=newrows[v,2] - newrows[v,2]*size_vary[v]
						newrows[v,4]=newrows[v,4] - newrows[v,4]*size_vary[v]*size_vary[v]
						newrows[v,5]=newrows[v,5] - newrows[v,5]*size_vary[v]
						newrows[v,7]=newrows[v,7] - newrows[v,7]*size_vary[v]
						newrows[v,8]=newrows[v,8] - newrows[v,8]*size_vary[v]
						newrows[v,9]=newrows[v,9] - newrows[v,9]*size_vary[v]
						newrows[v,10]=newrows[v,10] - newrows[v,10]*size_vary[v]
						newrows[v,11]=newrows[v,11] - newrows[v,11]*size_vary[v]
			Data=np.concatenate((Data,newrows),axis=0)
			yadd = np.ones(to_make) * cs[i]
			Target=np.concatenate((Target,yadd.astype(int)),axis=0)

		Data=Data[np.argsort(Target),:]
		Target = Target[np.argsort(Target)]

	if augment_class is None:
		if max_copies>0:
			print('Augmenting each class with additional %d samples via size augmentation' % max_copies)
			cs = np.unique(Target)
			for i in range(len(cs)):
				#generate n = max_copies of Data.
				possible_rows = np.where(Target==cs[i])[0]
				#sample to_make numbers from possible_rows.
				sampled_rows = np.random.choice(possible_rows,max_copies,replace=True)
				newrows = Data[sampled_rows,:]
				size_vary = s * np.random.rand(1,max_copies)[0]
				#vary size.
				for v in range(max_copies):
					if np.random.rand()<0.5:
						newrows[v,0]=newrows[v,0] + newrows[v,0]*size_vary[v]*size_vary[v]
						newrows[v,1]=newrows[v,1] + newrows[v,1]*size_vary[v]
						newrows[v,2]=newrows[v,2] + newrows[v,2]*size_vary[v]
						newrows[v,4]=newrows[v,4] + newrows[v,4]*size_vary[v]*size_vary[v]
						newrows[v,5]=newrows[v,5] + newrows[v,5]*size_vary[v]
						newrows[v,7]=newrows[v,7] + newrows[v,7]*size_vary[v]
						newrows[v,8]=newrows[v,8] + newrows[v,8]*size_vary[v]
						newrows[v,9]=newrows[v,9] + newrows[v,9]*size_vary[v]
						newrows[v,10]=newrows[v,10] + newrows[v,10]*size_vary[v]
						newrows[v,11]=newrows[v,11] + newrows[v,11]*size_vary[v]
					else:
						newrows[v,0]=newrows[v,0] - newrows[v,0]*size_vary[v]*size_vary[v]
						newrows[v,1]=newrows[v,1] - newrows[v,1]*size_vary[v]
						newrows[v,2]=newrows[v,2] - newrows[v,2]*size_vary[v]
						newrows[v,4]=newrows[v,4] - newrows[v,4]*size_vary[v]*size_vary[v]
						newrows[v,5]=newrows[v,5] - newrows[v,5]*size_vary[v]
						newrows[v,7]=newrows[v,7] - newrows[v,7]*size_vary[v]
						newrows[v,8]=newrows[v,8] - newrows[v,8]*size_vary[v]
						newrows[v,9]=newrows[v,9] - newrows[v,9]*size_vary[v]
						newrows[v,10]=newrows[v,10] - newrows[v,10]*size_vary[v]
						newrows[v,11]=newrows[v,11] - newrows[v,11]*size_vary[v]
				Data=np.concatenate((Data,newrows),axis=0)
				yadd = np.ones(max_copies) * cs[i]
				Target=np.concatenate((Target,yadd.astype(int)),axis=0)

			Data = Data[np.argsort(Target),:]
			Target = Target[np.argsort(Target)]

	else:
		augment_class=int(augment_class)
		if max_copies>0:
			print('Augmenting Class = %d with additional %d samples via size augmentation' % (augment_class,max_copies))
			#generate n = max_copies of Data.
			possible_rows = np.where(Target==augment_class)[0]
			#sample to_make numbers from possible_rows.
			sampled_rows = np.random.choice(possible_rows,max_copies,replace=True)
			newrows = Data[sampled_rows,:]
			size_vary = s * np.random.rand(1,max_copies)[0]
			#vary size.
			for v in range(max_copies):
				if np.random.rand()<0.5:
					newrows[v,0]=newrows[v,0] + newrows[v,0]*size_vary[v]*size_vary[v]
					newrows[v,1]=newrows[v,1] + newrows[v,1]*size_vary[v]
					newrows[v,2]=newrows[v,2] + newrows[v,2]*size_vary[v]
					newrows[v,4]=newrows[v,4] + newrows[v,4]*size_vary[v]*size_vary[v]
					newrows[v,5]=newrows[v,5] + newrows[v,5]*size_vary[v]
					newrows[v,7]=newrows[v,7] + newrows[v,7]*size_vary[v]
					newrows[v,8]=newrows[v,8] + newrows[v,8]*size_vary[v]
					newrows[v,9]=newrows[v,9] + newrows[v,9]*size_vary[v]
					newrows[v,10]=newrows[v,10] + newrows[v,10]*size_vary[v]
					newrows[v,11]=newrows[v,11] + newrows[v,11]*size_vary[v]
				else:
					newrows[v,0]=newrows[v,0] - newrows[v,0]*size_vary[v]*size_vary[v]
					newrows[v,1]=newrows[v,1] - newrows[v,1]*size_vary[v]
					newrows[v,2]=newrows[v,2] - newrows[v,2]*size_vary[v]
					newrows[v,4]=newrows[v,4] - newrows[v,4]*size_vary[v]*size_vary[v]
					newrows[v,5]=newrows[v,5] - newrows[v,5]*size_vary[v]
					newrows[v,7]=newrows[v,7] - newrows[v,7]*size_vary[v]
					newrows[v,8]=newrows[v,8] - newrows[v,8]*size_vary[v]
					newrows[v,9]=newrows[v,9] - newrows[v,9]*size_vary[v]
					newrows[v,10]=newrows[v,10] - newrows[v,10]*size_vary[v]
					newrows[v,11]=newrows[v,11] - newrows[v,11]*size_vary[v]
			Data=np.concatenate((Data,newrows),axis=0)
			yadd = np.ones(max_copies) * augment_class
			Target=np.concatenate((Target,yadd.astype(int)),axis=0)

			Data = Data[np.argsort(Target),:]
			Target = Target[np.argsort(Target)]


	return (Data, Target)


########################################################################
########################################################################
####### IMPORT THE DEV SET #####
########################################################################
########################################################################
def import_dev_set(dev_file_name='DevResults.txt'):
	print('Importing the dev set...')

	#import features
	featurelist=[]

	with open(dev_file_name,'r') as infile:
		for line in infile:
			featurelist.append(line.strip())

	# so now, featurelist[1] has names of things in form 'Area, MajorAxisLength, ... Class'
	FeatureNames = [x.strip() for x in featurelist[0].split(',')]
	#FeatureNames has form ['Area','MajorAxisLength',....'Class'] which is what I wanted

	DevData = [[float(x.strip()) for x in featurelist[i].split(',')] for i in range(1,len(featurelist))]

	# Data is in form [[1,2,3,....0.0],[3,3,1,...0.0],...[5,3,1,...0.0]], the last input is the class.

	Devclasses = [int(i[-1]) for i in DevData]

	#classes contains the class number from which the data is from

	#want to delete target from AllData.

	DevX = [i[0:-1] for i in DevData]

	#X has form similar to Data. So when we reshape, we want the output to be
	# X = array([[0,1,2,...]
	#            [1,2,3,...]])

	X_dev = np.asarray(DevX,order = 'F')

	#add aspect ratio as last column of data
	AR=[]
	for i in range(len(X_dev)):
		AR.append(X_dev[i,1]/X_dev[i,2])

	AR = np.asarray(AR)

	AR = AR.reshape((len(AR),1))

	X_dev = np.append(X_dev,AR,1) #concatenates arrays appropriately.

	#add form factor as last column of data
	#P^2/Area
	FF=[]
	for i in range(len(X_dev)):
		FF.append(X_dev[i,8]*X_dev[i,8] / X_dev[i,0])
	FF = np.asarray(FF)
	FF = FF.reshape((len(FF),1))
	X_dev = np.append(X_dev,FF,1)

	#this has the right form, is uses fortran column-major style memory representation vs row major C-style
	# the notation is scientific, where iris data set looks like a float. CHECKED: Both are type numpy.float64
	#both have same indexing calls, so I think we're in business.

	y_dev = np.asarray(Devclasses) #looks exactly correct, or at least like iris data set target.

	return (X_dev, y_dev, FeatureNames)

########################################################################
#########DATA IS IN THE SAME FORM AS IS FOUND IN IRIS DATASET###########
########################################################################
# Target = Target classes (0-4) for training and validation (type, numpy.int64, array)
# Data = Data for training and validation to be split. (type, numpy.float64, array)
# FeatureNames = Feature names for each column of data. (type, 'str', python list)
########################################################################
#print "Data is now in the same form as that found in Iris Dataset"
#print "Splitting the training dataset into train/val"

def apply_normalization(X_train, max_norm = False, l1_norm = False, l2_norm = False):
	########################################################
	if max_norm:
		print("Normalizing data using l1_norm")
		X_train = X_train / np.max(np.abs(X_train),0)[None,:]
	if l1_norm:
		print("Normalizing data using l1_norm")
		X_train = X_train / np.sum(X_train,0)[None,:]
	if l2_norm:
		print("Normalizing data using l1_norm")
		X_train = X_train / np.sqrt(np.sum(X_train*X_train,0))[None,:]

	return X_train

########################################################################

def preprocess_train_data(X_train, d = 2):

	############### SPLITTING THE DATASET ##################
	#First split the dataset so it is as if we only had a training set then a eval set.
	#X_train, X_test, y_train, y_test = train_test_split(Data, Target, test_size = .3)#.25)#, random_state =
	#default has shuffle = True. test_size sets the proportion of the data set to include in the test, here 25%.
	########################################################
	if d>1:
		print("Increasing dimensionality of dataset using cross terms")
	#################INCREASING FEATURES####################
	poly=preprocessing.PolynomialFeatures(degree = d, interaction_only = True)
	##IN SOME MODELS with 2 polynomial features, we are getting 90% exactly. In some polynomial 3 models,
	# we are getting 90.83%, which is exactly even with deep learning models.

	X_train = poly.fit_transform(X_train)
	#target_feature_names = ['x'.join(['{}^{}'.format(pair[0],pair[1]) for pair in tuple if pair[1]!=0]) for tuple in [zip(FeatureNames,p) for p in poly.powers_]]
	#poly=preprocessing.PolynomialFeatures(degree = 2, interaction_only = True)
	#X_test = poly.fit_transform(X_test)
	#poly=preprocessing.PolynomialFeatures(degree = 2, interaction_only = True)
	#X_dev = poly.fit_transform(X_dev)

	########################################################

	print("Scaling the data")
	################# SCALE THE DATA #######################
	#Scale the data. Each attribute in the dataset must be independently scaled, that is
	# 0 mean, and unit variance. Doing this returns the z-scores of the data
	# Z = (x - mu) / sigma

	scaler = preprocessing.RobustScaler().fit(X_train)#, QuantileTransformer(output_distribution='normal')
	#preprocessing.StandardScaler().fit(X_train) #IMPORTANT NOTE: We are scaling based only on training data!!!!

	X_train_scaled = scaler.fit_transform(X_train)

	#X_test_scaled = scaler.transform(X_test) # will be used later to evaluate the performance.

	#X_dev_scaled = scaler.transform(X_dev)

	##########################################################

	return (X_train_scaled, scaler)#, target_feature_names)

def preprocess_test_data(X_dev, scaler, d = 2):
		############### SPLITTING THE DATASET ##################
		#First split the dataset so it is as if we only had a training set then a eval set.
		#X_train, X_test, y_train, y_test = train_test_split(Data, Target, test_size = .3)#.25)#, random_state =
		#default has shuffle = True. test_size sets the proportion of the data set to include in the test, here 25%.
		########################################################

		print("Increasing dimensionality of dataset using cross terms")
		#################INCREASING FEATURES####################
		poly=preprocessing.PolynomialFeatures(degree = d, interaction_only = True)
		##IN SOME MODELS with 2 polynomial features, we are getting 90% exactly. In some polynomial 3 models,
		# we are getting 90.83%, which is exactly even with deep learning models.

		#X_train = poly.fit_transform(X_train)
		#target_feature_names = ['x'.join(['{}^{}'.format(pair[0],pair[1]) for pair in tuple if pair[1]!=0]) for tuple in [zip(FeatureNames,p) for p in poly.powers_]]
		#poly=preprocessing.PolynomialFeatures(degree = 2, interaction_only = True)
		#X_test = poly.fit_transform(X_test)
		#poly=preprocessing.PolynomialFeatures(degree = 2, interaction_only = True)
		X_dev = poly.fit_transform(X_dev)

		########################################################

		print("Scaling the data")
		################# SCALE THE DATA #######################
		#Scale the data. Each attribute in the dataset must be independently scaled, that is
		# 0 mean, and unit variance. Doing this returns the z-scores of the data
		# Z = (x - mu) / sigma

		#scaler = preprocessing.StandardScaler().fit(X_train) #IMPORTANT NOTE: We are scaling based only on training data!!!!

		#X_train_scaled = scaler.transform(X_train)

		#X_test_scaled = scaler.transform(X_test) # will be used later to evaluate the performance.

		X_dev_scaled = scaler.transform(X_dev)

		##########################################################

		return (X_dev_scaled)


def Add_Measures(Data, FeatureNames=None, add_AR=True, add_FF=True,
					add_convexity=True, add_curl_old=True, add_curl=True,
					add_sphericity=True, add_InscribedArea=True, add_BlebRel=True):
	############### EXPANDING THE DATASET ##################
	#Add measures of Aspect Ratio, Form Factor, Convexity, Curl, and Sphericity
	#Input: Data must be an np array with N (row) examples x M (cols) measures.
	#Measures should go: Area, MjrAxis, MnrAxis, Ecc,ConA,EqD,Sol,Ext,Per,conPer,fiber_length,InscribeR,bleb_len
	########################################################
	if add_AR:
		AR=[]
		for i in range(len(Data)):
			AR.append(Data[i,1]/Data[i,2])

		AR = np.asarray(AR)

		AR = AR.reshape((len(AR),1))

		Data = np.append(Data,AR,1) #concatenates arrays appropriately.
		if FeatureNames is not None:
			FeatureNames.extend(['AR'])

	if add_FF:
		#this measure is really compactness, if you multiply each by 4 pi
		#note this is different from roundness, which would use convex perimeter
		FF=[]
		for i in range(len(Data)):
			FF.append(Data[i,0] / (Data[i,8]*Data[i,8]))
			#FF.append(Data[i,8]*Data[i,8] / Data[i,0])

		FF = np.asarray(FF)
		FF = FF.reshape((len(FF),1))
		Data = np.append(Data,FF,1)
		if FeatureNames is not None:
			FeatureNames.extend(['FF'])

	if add_convexity:
		CC=[]
		for i in range(len(Data)):
			CC.append(Data[i,8] / Data[i,9])

		CC = np.asarray(CC)
		CC = CC.reshape((len(CC),1))
		Data = np.append(Data,CC,1)
		if FeatureNames is not None:
			FeatureNames.extend(['Convexity'])

	if add_curl_old:
		#tells how curled the object is. might help for lamellipodia.
		#curl is length / fiber length. (I assume length here can be major axis length)
		#fiber length definition is (perimeter - sqrt(perimeter^2 - 16*Area)) / 4

		#this definition does not work for a circle. Note that the result will be imaginary.
		#I changed the 16 to a 4Pi. This should be fine.
		cc=[]
		for i in range(len(Data)):
			if (4 * np.pi * Data[i,0]) <= (Data[i,8]*Data[i,8]):
				fiber_length = (Data[i,8] - np.sqrt((Data[i,8]*Data[i,8]) - (4 * np.pi * Data[i,0]))) / np.pi#4
				cc.append(Data[i,1] / fiber_length)
			else:
				fiber_length = Data[i,8] / np.pi#4
				cc.append(Data[i,1] / fiber_length)

		cc = np.asarray(cc)
		cc = cc.reshape((len(cc),1))
		Data = np.append(Data,cc,1)
		if FeatureNames is not None:
			FeatureNames.extend(['Curl_old'])

	if add_curl:
		cc=[]
		for i in range(len(Data)):
			cc.append(Data[i,1] / Data[i,10])

		cc = np.asarray(cc)
		cc = cc.reshape((len(cc),1))
		Data = np.append(Data,cc,1)
		#bound between 0 and 1 if major axis length could be replaced by feret diameter.
		if FeatureNames is not None:
			FeatureNames.extend(['Curl'])

	if add_sphericity:
		ss = []
		for i in range(len(Data)):
			ss.append(Data[i,11] * 2 / Data[i,1])

		ss = np.asarray(ss)
		ss = ss.reshape((len(ss),1))
		Data = np.append(Data,ss,1)
		#bound between 0 and 1 where 1 is a circle, perfectly spherical, and 0 is not at all.
		#would be better if we had feret diameter instead of major axis.
		if FeatureNames is not None:
			FeatureNames.extend(['Sphericity'])

	if add_InscribedArea:
		aa = []
		for i in range(len(Data)):
			aa.append(Data[i,1] * Data[i,1] * np.pi / Data[i,11])

		aa = np.asarray(aa)
		aa = aa.reshape((len(aa),1))
		Data = np.append(Data,aa,1)
		if FeatureNames is not None:
			FeatureNames.extend(['InArea'])

	if add_BlebRel:
		bb = []
		for i in range(len(Data)):
			bb.append(Data[i,12] / Data[i,11])

		bb = np.asarray(bb)
		bb = bb.reshape((len(bb),1))
		Data = np.append(Data,bb,1)
		if FeatureNames is not None:
			FeatureNames.extend(['Bleb_Rel'])

	if FeatureNames is not None:
		return (Data,FeatureNames)
	else:
		return Data

def Exclude_Measures(Data, FeatureNames=None, ex_Area = False, ex_MjrAxis=False, ex_MnrAxis=False, ex_Ecc=False,
 						ex_ConA=False, ex_EqD=False, ex_Sol=False, ex_Ext=False,
						ex_Per=False,ex_conPer=False,ex_FL=False,ex_InR=False,
						ex_bleb=False):
	#Area,MjrAxis,MnrAxis,Ecc,ConA,EqD,Sol,Ext,Per,conPer,FL,InR

	del_cols = []
	if ex_Area:
		del_cols.append(0)
	if ex_MjrAxis:
		del_cols.append(1)
	if ex_MnrAxis:
		del_cols.append(2)
	if ex_Ecc:
		del_cols.append(3)
	if ex_ConA:
		del_cols.append(4)
	if ex_EqD:
		del_cols.append(5)
	if ex_Sol:
		del_cols.append(6)
	if ex_Ext:
		del_cols.append(7)
	if ex_Per:
		del_cols.append(8)
	if ex_conPer:
		del_cols.append(9)
	if ex_FL:
		del_cols.append(10)
	if ex_InR:
		del_cols.append(11)
	if ex_bleb:
		del_cols.append(12)

	Data = np.delete(Data, del_cols, 1)
	if FeatureNames is not None:
		FeatureNames = [i for j, i in enumerate(FeatureNames) if j not in del_cols]
		return (Data,FeatureNames)
	else:
		return Data


def open_and_save_test_data(fpath,csvfilename,txtfilename,ratio):
	#fpath = '/volumes/chris stuff/chemsensing/chemsensing/Y27632_120518/Results/'
	#/Rho_Act_120118/Results_after/'
	#filename = 'FinalResults_after'
	#option to delete certain measures if done so in training.
	#order should go like
	#%frame number%correctedNum%area%centroidx%centroidy%major%minor%eccentricity
	#%orientation%convex area%filledarea%equivDiameter%solidity%extent%perimeter
	#%perimeter old%convex perimeter%fiber length%%max in radii%bleb length%centersx%centersy

	data = np.genfromtxt(fpath+csvfilename+'.csv',delimiter = ',', usecols=[2,5,6,7,9,11,12,13,14,16,17,18,19],skip_header=1)
	#was cols 3,6,7,8,10,12,13,14,15
	frames_cell =np.genfromtxt(fpath+csvfilename+'.csv',delimiter = ',', usecols=[0,1], skip_header=1)
	#add aspect ratio as last column of data

	data[:,0]=data[:,0]*ratio*ratio #area
	data[:,1]=data[:,1]*ratio #mjr
	data[:,2]=data[:,2]*ratio #MnrAxis
	#ecc unitless
	data[:,4]=data[:,4]*ratio*ratio#ConvexArea
	data[:,5]=data[:,5]*ratio#EquivDiameter
	#Solidity
	#Extent
	data[:,8]=data[:,8]*ratio#Perimeter
	data[:,9]=data[:,9]*ratio#conPerim
	data[:,10]=data[:,10]*ratio#FibLen
	data[:,11]=data[:,11]*ratio#max inscribed r
	data[:,12]=data[:,12]*ratio#bleblen

	preds = np.genfromtxt(fpath+'/'+txtfilename+'.txt', delimiter = ' ', usecols=[4,5,6,7],skip_header=1)
	y_target = np.where(np.max(preds,1)>0.7,np.argmax(preds,1),4)
	#y_target = np.reshape(y_target,(len(y_target),1))

	return (data, y_target, frames_cell)
