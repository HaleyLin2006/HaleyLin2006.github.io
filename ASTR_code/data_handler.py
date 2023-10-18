import pandas as pd
import numpy as np

class data():
    def __init__(self):
        self.df = pd.read_csv("/Users/haleylin/Desktop/ASTR Sphinx/ASTR_code/data.csv")

        # self.df = np.load('data.npz', allow_pickle=True)
        # self.data_param = pd.read_csv('data_param.csv')
    
    def save(self, param=[], model="", radius=[], density=[], v_r=[], v_theta=[], v_phi=[]):
        """
        Saving data into csv. list -> csv

        args:
            param: the parameter used to calculate data
                cs, t_acc, omega, time 
        """
        # ifSave = self.isSaved(param, model)

        # if ifSave == False:
        #     print('Saving data...')
        #     data = [radius, density, v_r, v_theta, v_phi]

        #     # print(len(radius), len(density), len(v_r), len(v_theta), len(v_phi))
        #     # np.savez('data.npz', data)
        #     self.data_param.loc[len(self.data_param)] = [param, model]
        #     self.data_param.to_csv('data_param.csv')

        ifSave = self.isSaved(param, model)

        if ifSave == False:
            print("Saving data...")
            data = [param, model, radius, density, v_r, v_theta, v_phi]
            self.df.loc[len(self.df)] = data            
            self.df.to_csv('data.csv')
    
    def load(self, param=[], model=""):
        tf_param = list(self.df['param'].isin([str(param)]))
        tf_model = list(self.df['model'].isin([str(model)]))
        for i in range(len(tf_param)):
            if tf_param[i] == True and tf_model[i] == True:
                radius  = self.df.iloc[i]['radius']
                density = self.df.iloc[i]['density']
                v_r     = self.df.iloc[i]['v_r']
                v_theta = self.df.iloc[i]['v_theta']
                v_phi   = self.df.iloc[i]['v_phi']

        return radius, density, v_r, v_theta, v_phi
    
    # def delete(self, param=[], model=""):
    
    def isSaved(self, param=[], model=""):
        """
        Check if the combination of parameter and model have been save enough.

        args:
            param: the parameter used to calculate data
                [cs, t_acc, omega, time]
            model: which version of calculation
                "Shu" or "TSC"

        return:
            ifSaved: determine if the combination have save before
                True: yes, the data have been saved before
                False: no, the data need to be save
        """

        # ifSave = False
        # tf_param = list(self.data_param['param'].isin([str(param)]))
        # tf_model = list(self.data_param['model'].isin([str(model)]))
        # for i in range(len(tf_param)):
        #     if tf_param[i] == True and tf_model[i] == True:
        #         ifSave = i
        #         break
        
        # if ifSave != False:
        #     print("Data Exist!")
        
        # return ifSave


        ifSave = False
        tf_param = list(self.df['param'].isin([str(param)]))
        tf_model = list(self.df['model'].isin([str(model)]))
        for i in range(len(tf_param)):
            if tf_param[i] == True and tf_model[i] == True:
                ifSave = True
                break
        
        if ifSave == True:
            print("")
        
        return ifSave
