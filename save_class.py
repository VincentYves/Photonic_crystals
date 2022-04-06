import os 
import numpy as np 

class Data_set():
    def __init__(self,name,root,data_lib = {},data_list = [], presave = 1):
        self.name = name
        self.root = root
        self.directory = self.root + "/" + self.name
        self.data_lib = data_lib
        self.data_list = data_list
        for key, value in self.data_lib.items() :
            print(key, " as key for ", value)
            self.data_list.append(key)
        
        if presave == 1:
            try:
                os.chdir(self.directory)          
            except Exception:
                os.chdir(self.root)
                os.mkdir(self.directory)
                os.chdir(self.directory)
                print("Creating folder ... ", self.directory)
            np.save(self.directory+"/DATALIB_"+self.name+".npy",self.data_lib,allow_pickle = True)
            #np.save(self.directory+"/DATALIST_"+self.name+".npy",self.data_list,allow_pickle = True)

    
    def create_directory(self):
        try:
            os.chdir(self.directory)
        except Exception:
            os.chdir(self.root)
            os.mkdir(self.directory)
            
    def update_data_lib(self,data,name_data):
        try:
            os.chdir(self.directory)          
        except Exception:
            self.create_directory(self.directory)
            os.chdir(self.directory)
            print("Creating folder ... ", self.directory)
        
        self.data_list.append(name_data)
        self.data_lib.update({name_data:data})
        np.save(self.directory+"/DATALIB_"+self.name+".npy",self.data_lib,allow_pickle = True)
        #np.save(self.directory+"/DATALIST_"+self.name+".npy",self.data_list,allow_pickle = True)
        print("====>  Data ",name_data+"_"+self.name, "saved in data library")
        
    def save_data_npy(self,data,name_data):
        try:
            os.chdir(self.directory)          
        except Exception:
            self.create_directory(self.directory)
            os.chdir(self.directory)
            print("Creating folder ... ", self.directory)
        
        np.save(self.directory+"/"+name_data+"_"+self.name+".npy",data,allow_pickle = True)
        print("====>  Data ",name_data+"_"+self.name, "save .npy format")
