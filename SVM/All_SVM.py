""""
================================
AUTHOR: CHRIS EDDY
GOAL: Turn CSV file that contains features of each cell at each time input and using SVM
calculate the probabilities for each phenotype of the cell.
Position is converted to micrometers, but no other measure needs to be converted since
training was done without proper micrometer units.
================================
"""

from __future__ import division
import numpy as np
from sklearn import svm, preprocessing, decomposition, linear_model
from sklearn.externals import joblib
from transform_data import *
from inspect import currentframe, getframeinfo


#add gridsearch
#add closer grid search
#add train
#add load model

class Cell_SVM(object):
	"""
	Encapsulates the SVM functionality
	"""
	def __init__(self):
		pass
		#add whatever variables we need

	def load_model(self,model_pkl_path, scaler_pkl_path):
		self.clf = joblib.load(model_pkl_path)
		self.scaler = joblib.load(scaler_pkl_path)

	def set_scale(self, scale):
		self.pixel2dist=scale
		#pixel2dist=0.53763672 #microns/pixel using 20x oil objective lens on confocal
		#pixel2dist=0.656265 #using 20x oil lens, microns/pixel on dark room microscope

	def load_test_data(self, data_path, d = 1, write_txt = True):
		"""
		Encapsulates SVM_csv.py

		Inputs:
			data_path: path to CSV data you wish to be evaluated
			write_txt: write text file with probabilities etc.
			d: dimensionality to which increase features. If d=2, quadratic, cross features are included.

		Outputs:
			Nothing, as of yet.
		"""
		data = np.genfromtxt(data_path+'.csv',delimiter = ',', usecols=[0,1,2,5,6,7,9,11,12,13,14,16,17,18,19],skip_header=1)#note columns differ in these exp
		#usecols=[0,1,3,6,7,8,10,12,13,14,15],skip_header=1)
		locs = np.genfromtxt(data_path+'.csv',delimiter = ',',usecols = [20,21], skip_header=1)#[4,5], skip_header=1)
		#first, reduce the size of the data. We only want data that are around for ~20 or more frames?
		if len(data.shape)==1:
			data=np.reshape(data,(-1,len(data)))
			locs = np.reshape(locs,(-1,len(locs)))


		cellnumbers = data[:,1];
		Times = data[:,0];

		if len(Times)>1:
			if Times[1] % 15 != 0:
				if Times[0]==0:
					Times=Times+1

		if Times[0] % 15 == 0:
			Times = Times / 15;

		#################################################################
		#I ALSO WANT TO PLOT IN MATLAB THE % OCCURENCES OF EACH CELL TYPE.
		# DOING THE NEXT FEW STEPS SLIGHTLY DISTORTS THE DATA
		# I WILL NEED TO ONLY DO CELLS IN MATLAB WITH GREATER THAN 20 FRAMES
		#################################################################

		#uniquecellnumbers, cellcounts = np.unique(cellnumbers, return_counts = True)

		#Good_cellnumbers = uniquecellnumbers[cellcounts>=20]

		#to_delete = []
		#for i in range(len(cellnumbers)):
		#	if cellnumbers[i] not in Good_cellnumbers:
		#		to_delete.append(i)

		#data = np.delete(data, to_delete, axis = 0)

		#Times = np.delete(Times, to_delete)

		#cellnumbers = np.delete(cellnumbers, to_delete)

		data = np.delete(data, [0,1], axis = 1) #delete time column and cellnumber column, since we already have them.


		print('Data is being converted to correct micrometers, assuming measured in pixels')
		#scale columns by following
		data[:,0] =data[:,0]*pixel2dist*pixel2dist#area
		data[:,1] = data[:,1]*pixel2dist#majoraxislength
		data[:,2] = data[:,2]*pixel2dist#minoraxislength
		#3 is eccentricity, unitless
		data[:,4] =data[:,4]*pixel2dist*pixel2dist#4 convex area,
		data[:,5] = data[:,5]*pixel2dist#5 Equiv diam
		#6 solidity, unitless
		#7 extent, unitless
		data[:,8] = data[:,8]*pixel2dist#8 Perimeter
		data[:,9] = data[:,9]*pixel2dist#8 Convex Perimeter
		data[:,10] = data[:,10]*pixel2dist#8 Fiber Length
		data[:,11] = data[:,11]*pixel2dist#8 Max Inscribed Radius
		data[:,12] = data[:,12]*pixel2dist# Bleb_length


		#add form factor as last column of data?
		#Form factor is added inside the Processing data when doing SVM. See "Transform_data.py"

		if len(Measure_Delete)>0:
			data = np.delete(data, Measure_Delete, axis = 1) #delete time column and cellnumber column, since we already have them.

		#so now data should look just like other data used in SVM

		#add aspect ratio as last column of data
		X_data = Add_Measures(data, add_AR=True, add_FF=True, add_convexity=True,
								add_curl_old=True, add_curl=True, add_sphericity=True,
								add_InscribedArea=True, add_BlebRel=True)
		#if you wish to exclude certain measures:
		#Area,MjrAxis,MnrAxis,Ecc,ConA,EqD,Sol,Ext,Per,conPer,FL,InR,bleb_M
		X_data = Exclude_Measures(X_data, ex_Area=False,
								ex_MjrAxis=False, ex_MnrAxis=False, ex_Ecc=False,
		 						ex_ConA=False, ex_EqD=False, ex_Sol=False, ex_Ext=False,
								ex_Per=False,ex_conPer=False,ex_FL=False,ex_InR=False,
								ex_bleb=False)

		####IF THE DATA WAS POLYNOMIAL BEFORE SCALED, DO THAT NOW!
		#frameinfo = getframeinfo(currentframe())
		#print("IF YOUR SCALER IS A POLYNOMIAL, YOU NEED TO EDIT THE POLYNOMIAL FEATURES, LINE %d CODE" % (frameinfo.lineno + 2))
		#d = 1
		if d==2:
			print("Expanding feature set to include quadratic, cross terms.")
			poly=preprocessing.PolynomialFeatures(degree = d, interaction_only = True)
			X_data_exp = poly.fit_transform(X_data)

			#FIRST, SCALE THE DATA USING THE SCALER
			X_data_scaled = scaler.transform(X_data_exp)
		else:
			X_data_scaled = scaler.transform(X_data)

		#GATHER PROBABILITIES
		Probs = clf.predict_proba(X_data_scaled)

		#Gather Predictions
		Predictions = clf.predict(X_data_scaled)

		#write to file
		#frames, cell numbers, Probs
		frames = Times.reshape((len(Times),1))
		cellnumbers = cellnumbers.reshape((len(cellnumbers),1))

		if write_txt:
			output = np.append(locs,Probs,1)
			output = np.append(cellnumbers,output,1)
			output = np.append(frames, output,1)

			Descriptors = ['frame', 'cellnumber','x-cent','y-cent','actinedge','filopodia','bleb','lamellipodia']

			thefile = open(data_path+'.txt', 'w')
			for item in Descriptors:
			  thefile.write("%s " % item)
			thefile.write("\n")
			thefile.close()

			File = file(data_path+'.txt','a')
			np.savetxt(File, output)
			File.close()

	def _load_train_data(self, training_path):
		"""
		training_path: path to training txt file, produced by Matlab code
		"""
		#import features
		print "Importing the training set..."
		featurelist=[]

		with open(training_path,'r') as infile:
			for line in infile:
				featurelist.append(line.strip())

		# so now, featurelist[1] has names of things in form 'Area, MajorAxisLength, ... Class'
		FeatureNames = [x.strip() for x in featurelist[0].split(',')]

		del FeatureNames[-1]

		#FeatureNames has form ['Area','MajorAxisLength',....'Class'] which is what I wanted

		AllData = [[float(x.strip()) for x in featurelist[i].split(',')] for i in range(1,len(featurelist))]

		# Data is in form [[1,2,3,....0.0],[3,3,1,...0.0],...[5,3,1,...0.0]], the last input is the class.

		classes = [int(i[-1]) for i in AllData]

		#classes contains the class number from which the data is from

		#want to delete target from AllData.

		X = [i[0:-1] for i in AllData]

		Data_og = np.asarray(X,order = 'F')
		Target_og = np.asarray(classes) #looks exactly correct, or at least like iris data set target.

		Target_og[Target_og==4]=2 #combine hemisphere and smallbleb

		#don't recommend doing augment_size. The only thing it scales is the distance per pixel ratio. The problem here is augmenting the area
		Data, Target = Augment_Size(Data_og, Target_og, max_copies=0, s=0.2, balance=False, augment_class=None)

		Data, FeatureNames = Add_Measures(Data, FeatureNames, add_AR=True, add_FF=True, add_convexity=True,
								add_curl_old=True, add_curl=True, add_sphericity=True,
								add_InscribedArea=True, add_BlebRel=True)

		#if you wish to exclude certain measures:
		#Area,MjrAxis,MnrAxis,Ecc,ConA,EqD,Sol,Ext,Per,conPer,FL,InR,bleb
		Data, FeatureNames = Exclude_Measures(Data, FeatureNames, ex_Area=False,
								ex_MjrAxis=False, ex_MnrAxis=False, ex_Ecc=False,
		 						ex_ConA=False, ex_EqD=False, ex_Sol=False, ex_Ext=False,
								ex_Per=False,ex_conPer=False,ex_FL=False,ex_InR=False,
								ex_bleb=False)

		self.X_train = Data
		self.FeatureNames = FeatureNames
		self.Y_train = Target

	def _load_dev_set_data(self, dev_path):
		########################################################################
		########################################################################
		####### IMPORT THE DEV SET #####
		########################################################################
		########################################################################

		print 'Importing the dev set...'

		#import features
		featurelist=[]
		with open(dev_path,'r') as infile:
			for line in infile:
				featurelist.append(line.strip())

		# so now, featurelist[1] has names of things in form 'Area, MajorAxisLength, ... Class'
		DevFeatureNames = [x.strip() for x in featurelist[0].split(',')]
		#FeatureNames has form ['Area','MajorAxisLength',....'Class'] which is what I wanted
		del DevFeatureNames[-1]

		DevData = [[float(x.strip()) for x in featurelist[i].split(',')] for i in range(1,len(featurelist))]

		# Data is in form [[1,2,3,....0.0],[3,3,1,...0.0],...[5,3,1,...0.0]], the last input is the class.

		Devclasses = [int(i[-1]) for i in DevData]

		#classes contains the class number from which the data is from

		#want to delete target from AllData.

		DevX = [i[0:-1] for i in DevData]

		X_dev = np.asarray(DevX,order = 'F')
		X_dev, DevFeatureNames = Add_Measures(X_dev, DevFeatureNames, add_AR=True, add_FF=True, add_convexity=True,
								add_curl_old=True, add_curl=True, add_sphericity=True,
								add_InscribedArea=True, add_BlebRel=True)

		#if you wish to exclude certain measures:
		#Area,MjrAxis,MnrAxis,Ecc,ConA,EqD,Sol,Ext,Per,conPer,FL,InR
		X_dev, DevFeatureNames = Exclude_Measures(X_dev, DevFeatureNames, ex_Area=False,
								ex_MjrAxis=False, ex_MnrAxis=False, ex_Ecc=False,
		 						ex_ConA=False, ex_EqD=False, ex_Sol=False, ex_Ext=False,
								ex_Per=False,ex_conPer=False,ex_FL=False,ex_InR=False,
								ex_bleb=False)
		#this has the right form, is uses fortran column-major style memory representation vs row major C-style
		# the notation is scientific, where iris data set looks like a float. CHECKED: Both are type numpy.float64
		#both have same indexing calls, so I think we're in business.

		y_dev = np.asarray(Devclasses) #looks exactly correct, or at least like iris data set target.

		y_dev[y_dev==4]=2 #combine classes

		self.X_dev = X_dev
		self.Y_dev = y_dev

	def train(self, training_path, dev_path, d = 1, C=1000, gamma = 0.01):
		"""
		Encapsulates CellShape_SVM_CrossTerms.py


		d: Default dimensionality of feature set. If d=2, this will expand feature set to include cross-terms

		C: SVM regularization parameter. For large values of C, the optimization will choose a smaller-margin hyperplane
		   if that hyperplane does a better job of getting all the training points classified correctly. Conversely, a
		   very small value of C will cause the optimizer to look for a larger-margin separating hyperplane, even if that
		   hyperplane misclassifies more points.
		   Used in both RBF and Linear kernel

		gamma: defines how far the influence of a single training example reaches, with low values meaning 'far' and high
		       values meaning 'close'. ... The gamma parameters can be seen as the inverse of the radius of influence of
			   samples selected by the model as support vectors.


		returns:
			models as attribute to Cell_SVM.
		"""
		#load training data
		self._load_train_data(training_path)
		#load test set data
		self._load_dev_set_data(dev_path)
		print "Data is now in the same form as that found in Iris Dataset"
		print "Splitting the training dataset into train/val"

		############### SPLITTING THE DATASET ##################
		#First split the dataset so it is as if we only had a training set then a eval set.
		X_train, X_test, y_train, y_test = train_test_split(self.X_train, self.Y_train, test_size = .20, stratify=Target)
		#default has shuffle = True. test_size sets the proportion of the data set to include in the test, here 25%.
		########################################################


		#################INCREASING FEATURES####################
		if d==2:
			print "Increasing dimensionality of dataset using cross terms"
			poly=preprocessing.PolynomialFeatures(degree = d, interaction_only = True)
			##IN SOME MODELS with 2 polynomial features, we are getting 90% exactly. In some polynomial 3 models,
			# we are getting 90.83%, which is exactly even with deep learning models.

			X_train = poly.fit_transform(X_train)
			temp = ['x'.join(['{}^{}'.format(pair[0],pair[1]) for pair in tuple if pair[1]!=0]) for tuple in [zip(self.FeatureNames,p) for p in poly.powers_]]
			target_feature_names=['1']
			target_feature_names.extend(temp)
			poly=preprocessing.PolynomialFeatures(degree = 2, interaction_only = True)
			X_test = poly.fit_transform(X_test)
			poly=preprocessing.PolynomialFeatures(degree = 2, interaction_only = True)
			X_dev = poly.fit_transform(X_dev)
			del target_feature_names[1]
		else:
			target_feature_names=[]
			target_feature_names.extend(self.FeatureNames)

		#will result in N(n,d) features where n is the number of features, d is the degree of the polynomial.
		#N(n,d) = (n+d)! / (n!d!)
		#since we are choosing to do the interaction terms only (ignore square terms),
		#We get N - len(Data[0,:]) features.

		#The transformation works like so: If your dataset is [a,b,c] features...
		# Output: [1, a, b, c, ab, bc, ca]
		########################################################


		print "Scaling the data"
		################# SCALE THE DATA #######################
		#Scale the data. Each attribute in the dataset must be independently scaled, that is
		# 0 mean, and unit variance. Doing this returns the z-scores of the data
		# Z = (x - mu) / sigma

		scaler = preprocessing.RobustScaler().fit(X_train)#preprocessing.StandardScaler().fit(X_train) #IMPORTANT NOTE: We are scaling based only on training data!!!!

		self.scaler = scaler

		X_train_scaled = scaler.transform(X_train)

		X_test_scaled = self.scaler.transform(X_test) # will be used later to evaluate the performance.

		X_dev_scaled = scaler.transform(self.X_dev)

		##########################################################

		# there are a couple of multilabel classifications in SVM
		# 1) One-vs-All Classifier consists in fitting one classifier per class, the class is
		#    fitted against all other classes. (most efficient and interpretable)
		# 2) One-vs-One Classifier constructs one classifier per pair of classes.
		# several other choices. One-vs-all is usually the default, so lets proceed with that to see how we do.

		# http://scikit-learn.org/stable/modules/multiclass.html

		#NOTE: there are multilabel classifiers as well.

		#http://scikit-learn.org/stable/auto_examples/svm/plot_iris.html

		print "Training the models"

		# we create an instance of SVM and fit out data.
		#if we wanted to plot the support vectors, We do not scale our
		#data. However, even with scaling, we can unscale the points.

		#C = 1000 # SVM regularization parameter

		#For large values of C, the optimization will choose a smaller-margin hyperplane
		#if that hyperplane does a better job of getting all the training points
		#classified correctly. Conversely, a very small value of C will cause the
		#optimizer to look for a larger-margin separating hyperplane, even if that
		#hyperplane misclassifies more points.

		#models = svm.SVC(kernel='linear',C=C)
		models = (svm.SVC(kernel='linear', C=C, probability = True, class_weight='balanced'),
		          svm.LinearSVC(loss = 'hinge',C=C),
		          linear_model.SGDClassifier(loss='hinge'),
		          svm.SVC(kernel='rbf', gamma=0.01, C=1000, probability = True, class_weight='balanced'),
		          svm.SVC(kernel='poly', degree=2, C=C, class_weight='balanced'))

		#SVC and NuSVC implement the "one against one" approach (Knerr et al., 1990) for
		#multiclass classification. If n_class is the number of classes, then
		#n_class * (n_class - 1) / 2 classifiers are constructed and each one trains data from two classes.
		#On the other hand, LinearSVC implements "one vs the rest" multiclass strategy,
		#thus training n_class models. If there are only two classes, only one model is trained
		#In the binary case, the probabilities are calibrated using Platt scaling:
		#logistic regression on the SVM's scores, fit by an additional cross validation
		#on the training set. In the multiclass case, this is extended as per Wu et al. (2004).

		#Wu, Lin and Weng, "Probability estimates for multiclass classification by pairwise coupling", JMLR 5:975-1005, 2004.
		#Platt "Probabilistic outputs for SVMs and comparisons to regularized likelihood methods".
		#models = models.fit(X_train_scaled,y_train)
		models = [clf.fit(X_train_scaled, y_train) for clf in models]
		#Scores = models.score(X_train_scaled,y_train)
		SelfScores = [clf.score(X_train_scaled,y_train) for clf in models]
		ValScores = [clf.score(X_test_scaled,y_test) for clf in models]
		DevScores = [clf.score(X_dev_scaled,self.y_dev) for clf in models]
		ps = [clf.predict(X_dev_scaled) for clf in models]

		print 'Results are as follows: '
		print 'Self Scores = ', SelfScores
		print 'Validation Scores = ', ValScores
		print 'Dev Set Scores = ', DevScores

		#Write performance into txt file

		titles = ('SVC with linear kernel',
		          'LinearSVC (linear kernel)',
				  'SGDClassifier (linear kernel)',
		          'SVC with RBF kernel',
		          'SVC with polynomial (degree 2) kernel')


		now=datetime.datetime.now()
		dtime=str(now.month)+'-'+str(now.day)+'-'+str(now.year)+' at '+str(now.hour)+'-'+str(now.minute)+'-'+str(now.second)
		file1=open("Model_Notes.txt","a")
		file1.write("This has been automatically generated to record models from CellShape_SVM_CrossTerms.py \n")
		file1.write("Generated on %s \n" % dtime)
		file1.write("Using Training data from %s \n" % trainingfname)
		file1.write("Using dev data from %s \n" % dev_file_name)
		file1.write("Model uses degree = %d polynomial features \n" % d)
		file1.write("with the following features: \n")
		file1.write("[")
		for m in range(len(FeatureNames)):
			file1.write("%s " % FeatureNames[m])
		file1.write("]\n")
		for m in range(len(titles)):
			file1.write('Model: %s, ' % titles[m])
			file1.write('Train Score: %.4f, ' % SelfScores[m])
			file1.write('Val Score: %.4f, ' % ValScores[m])
			file1.write('Dev Score: %.4f \n' % DevScores[m])
		file1.write('\n')
		file1.close()

		self.models = models

	def save_model(self, model, scaler, model_name, scaler_name):
		"""
		Although arguments should already be included in self, I placed so user would not forget.
		model_name: 'SVC_linear_kernel_date'
		scaler_name: 'SVC_linear_scaler_date'

		from sklearn.externals import joblib
		joblib.dump(models[0],'SVC_linear_kernel.pkl')
		DON'T FORGET TO SAVE THE SCALER AS WELL
		joblib.dump(scaler,'SVC_scaler.pkl')
		"""
		from sklearn.externals import joblib
		joblib.dump(model, model_name+'.pkl')
		joblib.dump(scaler,scaler_name+'.pkl')

		print "Saved model and scaler to current working directory as:"
		print("   "+model_name+'.pkl')
		print("   "+scaler_name+'.pkl')


	#Still need to add GridSearch and Closer GridSearch
	


#for load model
model_pkl_path = "/users/czeddy/documents/svm/SVC_rbf_010820_16942.pkl"
scaler_pkl_path = "/users/czeddy/documents/svm/SVC_rbf_scaler_010820_16942.pkl"
#test data
data_path = '/users/czeddy/documents/workingfolder/Drug_Treated/BEST_CELLS/New/Control/cell9/cell9'
