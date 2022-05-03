""""
================================
AUTHOR: CHRIS EDDY
DATE: 05/03/22
GOAL: We want to do PCA transform of a set of test data using a common basis
	calculated from a randomly selected, very large set of training data.
	We also want to pass this data along to Matlab for plotting and data
	analysis purposes. Here is an application using object-oriented programming.
================================
"""
from __future__ import division
import numpy as np
from transform_data import *
from inspect import currentframe, getframeinfo
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import scipy.io
import os

class PCA_Shape_Analysis(object):
	def __init__(self):
		pass

	def import_data(self, fpath, pixel2dist=0.5376372):
		"""
		INPUTS
		------------------------------------
		fpath: string, filepath to .csv file analyzed using Matlab code
				Previous_Cell_v2.m, which measures cell shape.
				ex. '/users/czeddy/Documents/MATLAB/Random_Walk/25R_15D.csv'

		pixel2dist: float, microns per pixel. Used to convert data recorded in
				pixels to micron size scale where appropriate.

		OUTPUTS
		-------------------------------------
		X_data: 2D numpy float array, shape N instances x M shape measures
		"""
		#check if file exists.
		if not os.path.isfile(fpath):
			raise ValueError("'fpath' argument does not point to an existing .csv file.")

		data = np.genfromtxt(fpath, delimiter = ',', usecols=[0,1,2,5,6,7,9,11,12,13,14,16,17,18,19],skip_header=1)#note columns differ in these exp
		#usecols=[0,1,3,6,7,8,10,12,13,14,15],skip_header=1)
		locs = np.genfromtxt(fpath, delimiter = ',',usecols = [20,21], skip_header=1)#[4,5], skip_header=1)
		#locs = np.genfromtxt(fpath+filename+'.csv',delimiter = ',',usecols = [3,4], skip_header=1)#[4,5], skip_header=1)

		#convert locs from pixels to micrometers
		#pixel2dist=0.53763672 #microns/pixel using 20x oil objective lens on confocal
		#pixel2dist=0.656265 #using 20x oil lens, microns/pixel, #0.538 for confocal?

		#Note, no other measure is scaled appropriately, they are all in units of pixels. However, because they are only
		#used for classification, this is not a problem--they are all linear transformations

		#names are as follows
		#Column, Measure
		#0, Time ; 1, cellnumner ; 3, Area ; 6, MajorAxisLength ; 7, MinorAxisLength; 8, Eccentricity ;
		#10, Convex Area ; 12, Equivalent Diameter ; 13, Solidity ; 14, Extent ; 15, Perimeter

		cellnumbers = data[:,1];
		Times = data[:,0];

		if Times[1] % 15 != 0:
			if Times[0]==0:
				Times=Times+1

		if Times[0] % 15 == 0:
			Times = Times / 15;


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

		#See "Transform_data.py"
		#so now data should look just like other data used in SVM
		#add measures
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
		return X_data, cellnumbers, Times, locs

	def reduce_training(self, X_data, cellnums, frames, locs, max_samples=15000):
		"""
		INPUTS
		------------------------------------
		X_data: 2D numpy float array, shape n_instances x n_features

		cellnums: 1D numpy float array, shape n_instances

		frames: 1D numpy float array, shape n_instances

		locs: 2D numpy float array, shape n_instances x 2

		max_samples: integer, specifies the maximum first dimension size of X_data

		OUTPUTS
		-------------------------------------
		X_data: 2D numpy float array, shape min(n_instances, max_samples) x n_features
			the reduced form of X_data after randomly selecting N = max_samples
			number of rows of X_data.

		cellnums: 1D numpy float array, shape min(n_instances, max_samples)
			the reduced from of cellnums, corresponding to the same rows taken
			from X_data

		frames: 1D numpy float array, shape min(n_instances, max_samples)
			the reduced from of frames, corresponding to the same rows taken
			from X_data

		locs: 2D numpy float array, shape min(n_instances, max_samples) x 2
			the reduced from of locs, corresponding to the same rows taken
			from X_data

		"""
		if X_data.shape[0]>max_samples:
			#randomly select 15,000 indeces without replacement.
			inds = np.random.choice(list(range(X_data.shape[0])), max_samples, replace=False)
			return X_data[inds,:], cellnums[inds], frames[inds], locs[inds,:]
		else:
			print("data has {} examples, smaller than the given 'max_samples' of {}".format(X_data.shape[0],max_samples))
			return X_data, cellnums, frames, locs

	def scaler_train(self, X_data):
		"""
		INPUTS
		------------------------------------
		X_data: 2D numpy float array, shape N instances x M shape measures;
				used for training purposes of PCA.

		OUTPUTS
		-------------------------------------
		m: 1D numpy float array, shape n_features
				the means calculated from each column of X_data

		s: 1D numpy float array, shape n_features
				the standard deviations calculated from each column of X_data
		"""
		#calculate means
		m = np.mean(X_data,axis=0) #returns an n dim array
		#calculate standard deviation
		s = np.std(X_data,axis=0)
		return m, s

	def apply_scaler(self, X_data, mean, std):
		"""
		INPUTS
		------------------------------------
		X_data: 2D numpy float array, shape n_instance x n_features;
				used for training purposes of PCA.

		mean: 1D numpy float array, n_features vector, of the means of each feature

		std: 1D numpy float array, n_features vector, of the std of each feature

		OUTPUTS
		-------------------------------------
		X_sub_mean: 2D numpy float array, shape n_instance x n_features;
				normalized following the scaling of data from which the mean and std were derived.
		"""
		assert len(mean)==X_data.shape[1], "number of features in data ({}) does not match length of mu vector ({}).".format(X_data.shape[1], len(mean))
		assert len(std)==X_data.shape[1], "number of features in data ({}) does not match length of std vector ({}).".format(X_data.shape[1], len(std))

		#convert to zero mean, unit variance.
		X_sub_mean = (X_data - mean[None,:]) / std[None,:] #subtracts and divides correctly.

		return X_sub_mean

	def PCA_fit(self, X_sub_mean, n_dim_keep=3):
		"""
		INPUTS
		------------------------------------
		X_sub_mean: 2D numpy float array, shape n_instance x n_features;
				normalized data.

		OUTPUTS
		-------------------------------------
		PCA_vecs: 2D numpy float array, shape n_features x n_dim_keep
				contains the first n_dim_keep PCA vectors used to transform
				data into PCA space.
		"""
		#calculate covariance matrix.
		CoVar = np.zeros(shape=(np.shape(X_sub_mean)[1],np.shape(X_sub_mean)[1]))
		n_samples = np.shape(X_sub_mean)[0]
		for measure_i in range(np.shape(X_sub_mean)[1]):
			for measure_j in range(np.shape(X_sub_mean)[1]):
				CoVar[measure_i,measure_j]= (1/(n_samples-1))*np.sum(X_sub_mean[:,measure_i] * X_sub_mean[:,measure_j])

		#calculate eigenvectors of covariance matrix
		vals, vecs = np.linalg.eig(CoVar)
		#import pdb;pdb.set_trace()
		#sort
		i_order = np.argsort(-vals) #argsort works in ascending order. This makes it so the biggest is positive
		vals = vals[i_order]
		vecs = vecs[:,i_order] #the column i of vecs is the eigenvector corresponding to the eigenvalue[i]
		#import pdb;pdb.set_trace()
		#determine how much variance is covered by the first 3 eigenvectors
		print("The first %d PCA vectors describe %.6f of the total variance." %(n_dim_keep, np.sum(vals[:n_dim_keep]) / np.sum(vals)))
		#Now, take only first three vectors.
		PCA_vecs=vecs[:,:n_dim_keep]
		#this is now N x 3 array
		#data is still N x M array
		return PCA_vecs

	def PCA_transform(self, X_sub_mean, PCA_vecs):
		"""
		INPUTS
		------------------------------------
		X_sub_mean: 2D numpy float array, shape n_instance x n_features;
				normalized data.

		PCA_vecs: 2D numpy float array, shape n_features x n_dim_keep
				contains the first n_dim_keep PCA vectors used to transform
				data into PCA space.
		OUTPUTS
		-------------------------------------
		data_PCA: 2D numpy float array, shape n_instance x n_dim_keep
				PCA transformed data.
		"""
		#now, project data by matrix multiplication.
		data_PCA = np.matmul(X_sub_mean,PCA_vecs) #should be M x 3)
		return data_PCA

	def plot_PCA_data(self, data_PCA):
		if data_PCA.shape[1]==3:
			fig = plt.figure()
			ax = plt.axes(projection='3d')
			ax.scatter3D(data_PCA[:,0], data_PCA[:,1], data_PCA[:,2], c='g',linewidth=0.5)
			ax.set_xlabel('V1')
			ax.set_ylabel('V2')
			ax.set_zlabel('V3')
			plt.show()
		elif data_PCA.shape[1]==2:
			fig,ax = plt.subplots(1)
			ax.scatter(data_PCA[:,0], data_PCA[:,1], c='g',linewidth=0.5)
			ax.set_xlabel('V1')
			ax.set_ylabel('V2')
			plt.show()
		else:
			print("dimensionality of PCA data is not correct.")

	def plot_test_data(self, n_traj, data_PCA, cellnums, data_train_PCA):
		inds = np.argwhere(cellnums==n_traj).flatten()
		if data_PCA.shape[1]==3:
			fig = plt.figure()
			ax = plt.axes(projection='3d')
			ax.plot3D(data_PCA[inds,0], data_PCA[inds,1], data_PCA[inds,2], c='r')
			ax.scatter3D(data_PCA[inds,0], data_PCA[inds,1], data_PCA[inds,2], c='g',linewidth=0.5)
			ax.set_xlabel('V1')
			ax.set_ylabel('V2')
			ax.set_zlabel('V3')
			ax.set_xlim(np.min(data_train_PCA[:,0]), np.max(data_train_PCA[:,0]))
			ax.set_ylim(np.min(data_train_PCA[:,1]), np.max(data_train_PCA[:,1]))
			ax.set_zlim(np.min(data_train_PCA[:,2]), np.max(data_train_PCA[:,2]))
			plt.show()
		elif data_PCA.shape[1]==2:
			fig,ax = plt.subplots(1)
			ax.plot(data_PCA[inds,0], data_PCA[inds,1], c='r')
			ax.scatter(data_PCA[inds,0], data_PCA[inds,1], c='g',linewidth=0.5)
			ax.set_xlabel('V1')
			ax.set_ylabel('V2')
			ax.set_xlim(np.min(data_train_PCA[:,0]), np.max(data_train_PCA[:,0]))
			ax.set_ylim(np.min(data_train_PCA[:,1]), np.max(data_train_PCA[:,1]))
			plt.show()
		else:
			print("dimensionality of PCA data is not correct.")


	def save_data_to_mat(self, data_PCA, cellnums, locs, frames, fpath):
		scipy.io.savemat(fpath, dict(data = data_PCA, frames = frames, locs = locs, cellnums = cellnums))

	def help(self):
		print("The order you should process your data is as follows: \n")
		print("P = PCA_Shape_Analysis()")
		print("train_data, train_frames = P.import_data(fpath = '/users/czeddy/documents/WorkingFolder/All_Combined.csv')")
		print("train_data, train_frames = P.reduce_training(train_data, train_frames, 15000)")
		print("mu, std = P.scaler_train(train_data)")
		print("train_data_scaled = P.apply_scaler(train_data, mu, std)")
		print("PCA_vecs = P.PCA_fit(train_data_scaled)")
		print("train_data_pca = P.PCA_transform(train_data_scaled, PCA_vecs)")
		print("test_data, test_frames = P.import_data(fpath = '/users/czeddy/documents/WorkingFolder/25R_15D_021720/25R_15D_021720_edited.csv')")
		print("test_data_scaled = P.apply_scaler(test_data, mu, std)")
		print("test_data_pca = P.PCA_transform(test_data_scaled, PCA_vecs)")

P = PCA_Shape_Analysis()
train_data, train_cnums, train_frames, train_locs = P.import_data(fpath = '/users/czeddy/documents/WorkingFolder/All_Combined.csv')
train_data, train_cnums, train_frames, train_locs = P.reduce_training(train_data, train_cnums, train_frames, train_locs, 15000)
mu, std = P.scaler_train(train_data)
train_data_scaled = P.apply_scaler(train_data, mu, std)
PCA_vecs = P.PCA_fit(train_data_scaled)
train_data_pca = P.PCA_transform(train_data_scaled, PCA_vecs)
test_data, test_cnums, test_frames, test_locs = P.import_data(fpath = '/users/czeddy/documents/WorkingFolder/25R_15D_021720/25R_15D_021720_edited.csv')
test_data_scaled = P.apply_scaler(test_data, mu, std)
test_data_pca = P.PCA_transform(test_data_scaled, PCA_vecs)
