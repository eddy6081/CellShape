from __future__ import division
import numpy as np
import os
#from sklearn import svm, preprocessing, decomposition, linear_model
#from sklearn.externals import joblib
from transform_data import *
from inspect import currentframe, getframeinfo
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

class Cell_Shape(object):
    """
    Object oriented python scripting of shape analysis.
    Mostly in Matlab code as well, but put here for ease of use.
    """
    def __init__(self,pixel2dist=0.53763672):
        self.pixel2dist=pixel2dist
        #0.53763672 microns/pixel using 20x oil objective lens on confocal
        #0.656265 using 20x oil lens, microns/pixel in dark room
        self.framerate = 15
        #minutes per frame

    def load_data(self,csv_filepath):
        #import the data

        data = np.genfromtxt(fpath+filename+'.csv',delimiter = ',', usecols=[0,1,2,5,6,7,9,11,12,13,14,16,17,18,19],skip_header=1)#note columns differ in these exp
        #usecols=[0,1,3,6,7,8,10,12,13,14,15],skip_header=1)
        locs = np.genfromtxt(fpath+filename+'.csv',delimiter = ',',usecols = [20,21], skip_header=1)#[4,5], skip_header=1)
        #locs = np.genfromtxt(fpath+filename+'.csv',delimiter = ',',usecols = [3,4], skip_header=1)#[4,5], skip_header=1)

        #convert locs from pixels to micrometers
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
        pixel2dist=self.pixel2dist
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

        self.Names = ['Area','MjrAxis','MnrAxis','Ecc','ConA','EqD','Sol',\
                      'Ext','Per','conPer','FL','InR','bleb len','AR','FF',\
                      'Convexity','Curl','Perim Curl', 'Sphericity',\
                      'In Area', 'bleb rel']

        self.X_data = X_data
        self.frames = Times
        self.CellNums = cellnumbers
        self.locs = locs
        print("returned attributes of 'X_data', 'frames', 'CellNums', 'Names', and 'locs'.")

    def load_PCA_train_data(self, txt_filepath, max_use = 15000):
        #import the data

        data = np.loadtxt(txt_filepath+'.txt',delimiter = ',', usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12), skiprows=1)#note columns differ in these exp
        #usecols=[0,1,3,6,7,8,10,12,13,14,15],skip_header=1)
        #randomly select only max_use rows.
        if data.shape[0]>max_use:
            print("Reducing computational load by selecting %d data points of %d for PCA training..." % (max_use, data.shape[0]))
            data = data[np.random.choice(data.shape[0],max_use),:]
        #convert locs from pixels to micrometers
        #names are as follows
        #Column, Measure
        #0, Time ; 1, cellnumner ; 3, Area ; 6, MajorAxisLength ; 7, MinorAxisLength; 8, Eccentricity ;
        #10, Convex Area ; 12, Equivalent Diameter ; 13, Solidity ; 14, Extent ; 15, Perimeter

        print('Data is being converted to correct micrometers, assuming measured in pixels')
        pixel2dist=self.pixel2dist
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

        self.X_PCA = X_data

    def normalize_data(self,data, m=None, s=None):
        # if not hasattr(self,'X_data'):
        #     print("you must first import the data using 'load_data'")
        #     return
        # else:
        if m is None:
            #subtract mean and divide by standard deviation. Removes units from measures.
            m = np.mean(data,axis=0) #returns an n dim array
            s = np.std(data,axis=0)

        X_norm = (data - m[None,:]) / s[None,:]
        return X_norm, m, s


    def PCA_train(self, d_keep=3, use_norm=True):
        if not hasattr(self,'X_PCA'):
            print("you must first import the data using 'load_PCA_train_data'")
            return
        else:
            #calculate means
            if use_norm:
                print("Normalizing data before PCA...")
                X_data, self.PCA_mean_1, self.PCA_std  = self.normalize_data(self.X_PCA)
            else:
                X_data = self.X_PCA

            m = np.mean(X_data,axis=0) #returns an n dim array
            self.PCA_mean_2 = m

            X_sub_mean = X_data - m[None,:] #subtracts correctly.

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
            print("The first %d PCA vectors describe %.6f of the total variance." %(d_keep, np.sum(vals[:d_keep]) / np.sum(vals)))
            #Now, take only first three vectors.
            PCA_vecs=vecs[:,:d_keep]
            #this is now N x 3 array
            #data is still N x M array

            self.PCA_vecs = PCA_vecs
            self.PCA_dim = d_keep
            #now, project data by matrix multiplication.
            data_PCA = np.matmul(X_sub_mean,PCA_vecs) #should be M x 3
            self.PCA_data = data_PCA

    def plot_PCA_train(self):
        if not hasattr(self,'PCA_data'):
            print("you must first run 'PCA_train'.")
            return
        else:
            if self.PCA_dim==3:
                #Plot
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.scatter3D(self.PCA_data[:,0], self.PCA_data[:,1], self.PCA_data[:,2], c='g',linewidth=0.5);
                ax.set_xlabel('V1')
                ax.set_ylabel('V2')
                ax.set_zlabel('V3')
                #
                plt.show()
            elif self.PCA_dim==2:
                #plot
                fig,ax = plt.subplots(1)
                ax.scatter(self.PCA_data[:,0], self.PCA_data[:,1], c='g',linewidth=0.5)
                ax.set_xlabel('V1')
                ax.set_ylabel('V2')
                plt.show()
            else:
                print("The PCA dimensionality is %d and cannot be plotted." %(self.PCA_dim))
                return

    def PCA_transform_data(self, data, normalize = True):
        if not hasattr(self,'PCA_vecs'):
            print("You must first run 'PCA_train'.")
            return
        else:
            if normalize:
                data,_,_ = self.normalize_data(data,self.PCA_mean_1,self.PCA_std)

            data = data - self.PCA_mean_2[None,:]

            data_PCA = np.matmul(data,self.PCA_vecs) #should be M x 3
            return data_PCA

    def plot_PCA_trajectory(self, PCA_data, traj=1):
        if not hasattr(self,'PCA_data'):
            print("you must first run 'PCA_train'.")
            return
        #plot 30% training data with alpha
        #plot trajectory on top.
        #find the inds where traj is
        inds = np.where(self.CellNums.astype(int)==traj,True,False)
        if np.sum(inds)==0:
            print("the desired trajectory has length 0. Cannot plot")
            return
        else:
            #randomly pick like 30% of the training inds to display
            train_inds = np.random.choice(self.PCA_data.shape[0],int(self.PCA_data.shape[0]*0.1))
            if self.PCA_dim==3:
                #Plot
                fig = plt.figure()
                ax = plt.axes(projection='3d')
                ax.scatter3D(self.PCA_data[train_inds,0], self.PCA_data[train_inds,1], self.PCA_data[train_inds,2], c='gray', alpha=0.05, linewidth=0.5)
                ax.set_xlabel('V1')
                ax.set_ylabel('V2')
                ax.set_zlabel('V3')
                #ax.plot3d(PCA_data[inds,0],PCA_data[inds,1],PCA_data[inds,2], 'k', line)
                time = np.arange(int(np.sum(inds)))
                time = time / np.max(time)
                ax.plot3D(PCA_data[inds,0],PCA_data[inds,1],PCA_data[inds,2],'black')
                ax.scatter3D(PCA_data[inds,0],PCA_data[inds,1],PCA_data[inds,2], c=time, cmap='viridis', alpha=1, linewidth=0.7)
                #
                #import pdb;pdb.set_trace()
                plt.show()
            elif self.PCA_dim==2:
                #plot
                fig,ax = plt.subplots(1)
                ax.scatter(self.PCA_data[train_inds,0], self.PCA_data[train_inds,1], c='gray', alpha=0.05, linewidth=0.5)
                ax.set_xlabel('V1')
                ax.set_ylabel('V2')
                time = np.arange(int(np.sum(inds)))
                time = time / np.max(time)
                ax.plot(PCA_data[inds,0],PCA_data[inds,1],'black')
                ax.scatter(PCA_data[inds,0],PCA_data[inds,1], c=time, cmap='viridis', alpha=1, linewidth=0.7)

                plt.show()
            else:
                print("The PCA dimensionality is %d and cannot be plotted." %(self.PCA_dim))
                return

    def MSD_shape(self,measure=0):
        #plot the MSD of whatever measure is passed. use self.X_data
        #need frames and cell numbers.
        print("Calculating MSD of %s" %(self.Names[measure]))
        #determine the number of trajectories
        cellnums = np.unique(self.CellNums)
        #make holder for all cell square displacements.
        all_sq_disp=[[] for i in range(len(cellnums))]

        for i,n_cell in enumerate(cellnums):
            #find which rows n_cell is in self.CellNums.
            inds = np.where(self.CellNums==n_cell,True,False)
            #get time slots for these inds
            frames = self.frames[inds]
            #is frames integer values? No, but it is a whole number. convert.
            frames = frames.astype('int')
            M = self.X_data[inds,measure]
            #Need a holder for square displacements. List of lists?
            cell_sq_disp = [[] for x in range(frames[-1]-frames[0])] #empty list of lists.
            #bear in mind that frames may not be continuous
            for t in range(1,len(frames)):
                taus = frames[t:] - frames[:-t]
                taus = taus - 1 #indexing begins at 0.
                #taus=taus.astype('int')
                sq_diff = (M[t:] - M[:-t])**2
                #import pdb;pdb.set_trace()
                for j in range(len(taus)):
                    cell_sq_disp[taus[j]].append(sq_diff[j])

            #append cell_sq_disp
            all_sq_disp[i] = cell_sq_disp #append list of lists to list.
            #we don't need to do this append, it is just to make sure things go right.
            #the more direct way is to append cell_sq_disp to the final holder matrix.

        #get rid of weird indexing
        #all_sq_disp = [all_sq_disp[i][0] if len(all_sq_disp[i])>0 else [] for i in range(len(cellnums))]
        #determine max length.
        ml = int(np.max([len(x) for x in all_sq_disp]))

        #create array to hold
        all_sq_diff_together =[[] for x in range(ml)]
        for traj in range(len(cellnums)):
            for tau in range(len(all_sq_disp[traj])):
                if len(all_sq_disp[traj][tau])>0:
                    all_sq_diff_together[tau].extend(all_sq_disp[traj][tau])

        #calculate means
        msd = [np.mean(all_sq_diff_together[tau]) if len(all_sq_diff_together[tau])>0 else np.nan for tau in range(ml)]
        msd_SEM = [np.std(all_sq_diff_together[tau])/np.sqrt(len(all_sq_diff_together[tau])) if len(all_sq_diff_together[tau])>0 else np.nan for tau in range(ml)]
        msd.insert(0,0.0)
        msd_SEM.insert(0,0.0)
        msd = np.array(msd)
        msd_SEM = np.array(msd_SEM)
        mask = np.isfinite(msd)
        xs = np.arange(len(msd))
        #generate plot
        plt.errorbar(xs[mask], msd[mask], msd_SEM[mask], linestyle='-', marker='o')
        plt.title('MSD %s'%(self.Names[measure]))
        plt.xlabel('tau (frames)')
        plt.ylabel('\sigma^2')
        plt.show()

    def MSD_real_space(self):
        #plot MSD of real space migration. use self.locs which is N x 2
        #need frames and cell numbers.
        print("Calculating MSD of real space")
        #determine the number of trajectories
        cellnums = np.unique(self.CellNums)
        #make holder for all cell square displacements.
        all_sq_disp=[[] for i in range(len(cellnums))]

        for i,n_cell in enumerate(cellnums):
            #find which rows n_cell is in self.CellNums.
            inds = np.where(self.CellNums==n_cell,True,False)
            #get time slots for these inds
            frames = self.frames[inds]
            #is frames integer values? No, but it is a whole number. convert.
            frames = frames.astype('int')
            M = self.locs[inds,:]
            #Need a holder for square displacements. List of lists?
            cell_sq_disp = [[] for x in range(frames[-1]-frames[0])] #empty list of lists.
            #bear in mind that frames may not be continuous
            for t in range(1,len(frames)):
                taus = frames[t:] - frames[:-t]
                taus = taus - 1 #indexing begins at 0.
                #taus=taus.astype('int')
                sq_diff = np.sum((M[t:,:] - M[:-t,:])**2,axis=1) #add deltaX^2 + deltaY^2
                #import pdb;pdb.set_trace()
                for j in range(len(taus)):
                    cell_sq_disp[taus[j]].append(sq_diff[j])

            #append cell_sq_disp
            all_sq_disp[i] = cell_sq_disp #append list of lists to list.
            #we don't need to do this append, it is just to make sure things go right.
            #the more direct way is to append cell_sq_disp to the final holder matrix.

        #get rid of weird indexing
        #all_sq_disp = [all_sq_disp[i][0] if len(all_sq_disp[i])>0 else [] for i in range(len(cellnums))]
        #determine max length.
        ml = int(np.max([len(x) for x in all_sq_disp]))

        #create array to hold
        all_sq_diff_together =[[] for x in range(ml)]
        for traj in range(len(cellnums)):
            for tau in range(len(all_sq_disp[traj])):
                if len(all_sq_disp[traj][tau])>0:
                    all_sq_diff_together[tau].extend(all_sq_disp[traj][tau])

        #calculate means
        msd = [np.mean(all_sq_diff_together[tau]) if len(all_sq_diff_together[tau])>0 else np.nan for tau in range(ml)]
        msd_SEM = [np.std(all_sq_diff_together[tau])/np.sqrt(len(all_sq_diff_together[tau])) if len(all_sq_diff_together[tau])>0 else np.nan for tau in range(ml)]
        msd.insert(0,0.0)
        msd_SEM.insert(0,0.0)
        msd = np.array(msd)
        msd_SEM = np.array(msd_SEM)
        mask = np.isfinite(msd)
        xs = np.arange(len(msd))
        #generate plot
        plt.errorbar(xs[mask], msd[mask], msd_SEM[mask], linestyle='-', marker='o')
        plt.title('MSD Real Space')
        plt.xlabel('tau (frames)')
        plt.ylabel('\sigma^2')
        plt.show()

    def MSD_shape_space(self,PCA_data):
        #plot MSD of real space migration. use self.locs which is N x 2
        #need frames and cell numbers.
        print("Calculating MSD of PCA shape space")
        #determine the number of trajectories
        cellnums = np.unique(self.CellNums)
        #make holder for all cell square displacements.
        all_sq_disp=[[] for i in range(len(cellnums))]

        for i,n_cell in enumerate(cellnums):
            #find which rows n_cell is in self.CellNums.
            inds = np.where(self.CellNums==n_cell,True,False)
            #get time slots for these inds
            frames = self.frames[inds]
            #is frames integer values? No, but it is a whole number. convert.
            frames = frames.astype('int')
            M = PCA_data[inds,:]
            #Need a holder for square displacements. List of lists?
            cell_sq_disp = [[] for x in range(frames[-1]-frames[0])] #empty list of lists.
            #bear in mind that frames may not be continuous
            for t in range(1,len(frames)):
                taus = frames[t:] - frames[:-t]
                taus = taus - 1 #indexing begins at 0.
                #taus=taus.astype('int')
                sq_diff = np.sum((M[t:,:] - M[:-t,:])**2,axis=1) #add deltaX^2 + deltaY^2
                #import pdb;pdb.set_trace()
                for j in range(len(taus)):
                    cell_sq_disp[taus[j]].append(sq_diff[j])

            #append cell_sq_disp
            all_sq_disp[i] = cell_sq_disp #append list of lists to list.
            #we don't need to do this append, it is just to make sure things go right.
            #the more direct way is to append cell_sq_disp to the final holder matrix.

        #get rid of weird indexing
        #all_sq_disp = [all_sq_disp[i][0] if len(all_sq_disp[i])>0 else [] for i in range(len(cellnums))]
        #determine max length.
        ml = int(np.max([len(x) for x in all_sq_disp]))

        #create array to hold
        all_sq_diff_together =[[] for x in range(ml)]
        for traj in range(len(cellnums)):
            for tau in range(len(all_sq_disp[traj])):
                if len(all_sq_disp[traj][tau])>0:
                    all_sq_diff_together[tau].extend(all_sq_disp[traj][tau])

        #calculate means
        msd = [np.mean(all_sq_diff_together[tau]) if len(all_sq_diff_together[tau])>0 else np.nan for tau in range(ml)]
        msd_SEM = [np.std(all_sq_diff_together[tau])/np.sqrt(len(all_sq_diff_together[tau])) if len(all_sq_diff_together[tau])>0 else np.nan for tau in range(ml)]
        msd.insert(0,0.0)
        msd_SEM.insert(0,0.0)
        msd = np.array(msd)
        msd_SEM = np.array(msd_SEM)
        mask = np.isfinite(msd)
        xs = np.arange(len(msd))
        #generate plot
        plt.errorbar(xs[mask], msd[mask], msd_SEM[mask], linestyle='-', marker='o')
        plt.title('MSD Shape Space')
        plt.xlabel('tau (frames)')
        plt.ylabel('\sigma^2')
        plt.show()

    def Autocorr_shape(self,measure=0):
        #calculate autocorrelation for shape measures.
        #do entire trajectory, or just one trajectory.
        print("Calculating Autocorrelation of %s" %(self.Names[measure]))
        #determine the number of trajectories
        cellnums = np.unique(self.CellNums)
        #make holder for all cell square displacements.
        all_autocorr=[[] for i in range(len(cellnums))]

        for traj,n_cell in enumerate(cellnums):
            #find which rows n_cell is in self.CellNums.
            inds = np.where(self.CellNums==n_cell,True,False)
            #get time slots for these inds
            frames = self.frames[inds]
            #is frames integer values? No, but it is a whole number. convert.
            frames = frames.astype('int')
            #get measure
            M = self.X_data[inds,measure]
            #calculate mean and variance
            m = np.mean(M)
            #calculate unnormalized variance for denominator
            if len(frames)>1:
                var_nn = np.sum((M - m)**2)
            else:
                var_nn = np.sum(M**2)
            #cell autocorrelation holder.
            #need a way to distinguish if there are discontinuous lags.
            cell_autocorr = np.empty(((frames[-1]-frames[0])+1,))
            cell_autocorr[:]=np.nan #fill with nan
            for t in range(0,len(frames)):
                if t>0:
                    taus = frames[t:] - frames[:-t]
                    M_t = M[t:]
                    M_tk = M[:-t]
                else:
                    taus = frames - frames
                    M_t = M
                    M_tk = M

                for j in range(len(taus)):
                    if not np.isfinite(cell_autocorr[taus[j]]):
                        if len(frames)>1:
                            cell_autocorr[taus[j]]=(M_t[j]-m)*(M_tk[j]-m)
                        else:
                            cell_autocorr[taus[j]]=M_t[j]*M_tk[j]
                    else:
                        if len(frames)>1:
                            cell_autocorr[taus[j]]+=(M_t[j]-m)*(M_tk[j]-m)
                        else:
                            cell_autocorr[taus[j]]+=M_t[j]*M_tk[j]

            #import pdb;pdb.set_trace()
            cell_autocorr /= var_nn
            all_autocorr[traj]=cell_autocorr
        #now lets compute the average autocorrelation.

        #compute maximum lag time in all trajectories.
        ml = int(np.max([len(x) for x in all_autocorr]))
        #create array to put all autocorrelation trajectories together.
        all_autocorr_together =[[] for x in range(ml)]
        for traj in range(len(cellnums)):
            for tau in range(len(all_autocorr[traj])):
                if np.isfinite(all_autocorr[traj][tau]):
                    all_autocorr_together[tau].append(all_autocorr[traj][tau])

        #compute means
        auto_mean = np.array([np.mean(all_autocorr_together[tau]) if len(all_autocorr_together[tau])>0 else np.nan for tau in range(ml)])
        auto_SEM = np.array([np.std(all_autocorr_together[tau])/np.sqrt(len(all_autocorr_together[tau])) if len(all_autocorr_together[tau])>0 else np.nan for tau in range(ml)])
        mask = np.isfinite(auto_mean) #if any are nan, mask it. This shouldn't be the case though.
        xs = np.arange(len(auto_mean))
        #import pdb;pdb.set_trace()
        #plot autocorr
        #plt.plot(xs[mask],auto_mean[mask],linestyle='--', marker='o')
        plt.errorbar(xs[mask],auto_mean[mask],auto_SEM[mask],linestyle='--', marker='o')
        plt.title('Autocorrelation of %s'%(self.Names[measure]))
        plt.xlabel('delay (frames)')
        plt.ylabel('rho_k')
        plt.grid()
        plt.show()


    def Autocorr_velocity(self):
        #calculate instantaneous velocities (directions too)
        #calculate autocorrelation.
        print("Calculating Autocorrelation of velocity")
        #determine the number of trajectories
        cellnums = np.unique(self.CellNums)
        #make holder for all cell square displacements.
        all_autocorr=[[] for i in range(len(cellnums))]

        for traj,n_cell in enumerate(cellnums):
            #find which rows n_cell is in self.CellNums.
            inds = np.where(self.CellNums==n_cell,True,False)
            #get time slots for these inds
            frames = self.frames[inds]
            #is frames integer values? No, but it is a whole number. convert.
            frames = frames.astype('int')
            #get measure
            M = self.locs[inds,:]
            if M.shape[0]>1:
                #calculate instantaneous velocities
                #what about for non-continuous frames? NaN.
                fr_diff = frames[1:]-frames[:-1]
                V = M[1:,:] - M[:-1,:] / fr_diff[:,None] #pixel per frame
                #convert units
                V = V * self.pixel2dist / self.framerate
                #actually, conversion shouldn't matter.

                #so if discontinuous, we need to make V NaN.
                #V = np.where(fr_diff[:,None]==1.0, V, [np.nan,np.nan])
                V = V[np.where(fr_diff==1.0),:][0]
                #also, add a nan to the end so V is the same length as frames.
                #V = np.concatenate((V,np.array([[np.nan, np.nan]])))
                #no, just delete last frame.
                frames=frames[:-1]
                #easiest method is to delete frames and velocities
                frames=frames[np.where(fr_diff==1.0)]
                #import pdb;pdb.set_trace()

                #calculate mean and variance
                #m = np.nanmean(V,axis=0)
                #calculate unnormalized variance for denominator
                var_nn = np.sum(np.nansum(V**2,axis=0))
                #cell autocorrelation holder.
                #need a way to distinguish if there are discontinuous lags.
                try:
                    cell_autocorr = np.empty(((frames[-1]-frames[0])+1,))
                except:
                    import pdb;pdb.set_trace()
                cell_autocorr[:]=np.nan #fill with nan
                for t in range(0,len(frames)):
                    if t>0:
                        taus = frames[t:] - frames[:-t]
                        M_t = V[t:,:]
                        M_tk = V[:-t,:]
                    else:
                        taus = frames - frames
                        M_t = V
                        M_tk = V

                    for j in range(len(taus)):
                        for dim in range(np.shape(V)[1]):
                            if not np.isfinite(cell_autocorr[taus[j]]):
                                cell_autocorr[taus[j]]=M_t[j,dim]*M_tk[j,dim]
                            else:
                                cell_autocorr[taus[j]]+=M_t[j,dim]*M_tk[j,dim]


                #import pdb;pdb.set_trace()
                cell_autocorr /= var_nn
                #import pdb;pdb.set_trace()
                all_autocorr[traj]=cell_autocorr
        #now lets compute the average autocorrelation.
        #compute maximum lag time in all trajectories.
        ml = int(np.max([len(x) for x in all_autocorr]))
        #create array to put all autocorrelation trajectories together.
        all_autocorr_together =[[] for x in range(ml)]
        for traj in range(len(cellnums)):
            for tau in range(len(all_autocorr[traj])):
                if np.isfinite(all_autocorr[traj][tau]):
                    all_autocorr_together[tau].append(all_autocorr[traj][tau])

        #compute means
        auto_mean = np.array([np.mean(all_autocorr_together[tau]) if len(all_autocorr_together[tau])>0 else np.nan for tau in range(ml)])
        auto_SEM = np.array([np.std(all_autocorr_together[tau])/np.sqrt(len(all_autocorr_together[tau])) if len(all_autocorr_together[tau])>0 else np.nan for tau in range(ml)])
        mask = np.isfinite(auto_mean) #if any are nan, mask it. This shouldn't be the case though.
        xs = np.arange(len(auto_mean))
        #import pdb;pdb.set_trace()
        #plot autocorr
        #plt.plot(xs[mask],auto_mean[mask],linestyle='--', marker='o')
        plt.errorbar(xs[mask],auto_mean[mask],auto_SEM[mask],linestyle='--', marker='o')
        plt.title('Autocorrelation of velocity')
        plt.xlabel('delay (frames)')
        plt.ylabel('rho_k')
        plt.grid()
        plt.show()

    def Autocorr_persist(self):
        #after caluclating from Compute_persist_vecs, do this.
        pass

    def Compute_persist_vecs(self, N, future=0, past=0):
        #calculate persistence vectors using SVD of N trajectories.
        #compute vectors from t to t+1, t+2, ... t+N, ...t+N+future
        #compute vectors from t-N-past, ... t-N, ... t-2, t-1 to t
        #use these vectors in SVD.
        pass

    #need a function for calculation of either mean speed or Instantaneous speed.


fpath = '/users/czeddy/Documents/MATLAB/Random_Walk/'
filename = '25R_15D'
P1 = fpath + filename #use for load_data

PCA_path = '/users/czeddy/Documents/MATLAB/Random_Walk/All_Combined_training'
#the above path is a file that contains all cell states that have been classified,
#in all types of gels.


#Load object
CS = Cell_Shape()

##Load PCA data and plot it!
CS.load_PCA_train_data(PCA_path)
CS.PCA_train()#normalize is default, as is use 3 PCA vectors #works
CS.plot_PCA_train()

#Load a different dataset, project it into PCA
CS.load_data(fpath+filename)
dta = CS.PCA_transform_data(CS.X_data)
CS.plot_PCA_trajectory(dta, traj=2)

#Perform MSD of the PCA transformed data, to see diffusivity in shape space.
CS.MSD_shape_space(dta) #interesting!
#All data shows a constrained random walk in shape space,
#but superdiffusivity in real space.
#Perhaps this can be modeled by a two state system, on or off motility with a
#single barrier with opening width O and barrier width W

#Check the velocity autocorrelation!
#CS.Autocorr_velocity() #works

#Check the autocorrelation of shape to see memory of shape measure!
#CS.Autocorr_shape(measure=0) #works

#See the real space mean square displacement!
CS.MSD_real_space() #works

#See the MSD curve of a measure!
#CS.MSD_shape(measure=0) #works
