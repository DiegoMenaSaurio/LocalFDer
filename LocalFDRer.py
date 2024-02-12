
#import libraries

import pandas as pd
import numpy as np
import re
import os
import sys


#create a class to store a dataframe and to define methods to calculate local and global FDR

class dataframe:
    def __init__(self, table):
        self.table = table
        
#adds a column with the cXcorr of each PSMs

    def xcorr(self):
        c= [1 if i<3 else 1.22 for i in self.table['Charge'].to_list()]
        a = np.log(self.table['XCorr']/c)
        b = np.log(self.table['MH+ [Da]']*2/110)
        self.table['XCorr_corr'] = a/b
        self.table = self.table.sort_values(by=['XCorr_corr'], ascending = False)
        return self.table
    
#filter the dataframe according to the FDR threshold

    def fdr_filter(self,FDR_value):
        decoys_Trues = [decoy_indicator in i for i in self.table['Protein Accessions'].to_list() ]
        self.table['FDR_local'] = decoys_Trues
        decoys = np.cumsum([self.table['FDR_local'].to_list()])
        #totales = np.cumsum([type(i) == str for i in self.table['Protein Accessions']])
        target = np.cumsum([decoy_indicator not in i for i in self.table['Protein Accessions']])
        self.table['FDR']=decoys/target
        self.table = self.table[self.table['FDR']<=FDR_value]
        self.table = self.table.loc[self.table['FDR_local']==False]
        return self.table
    
#filter the dataframe according to the selected modified amino acids

    def aminoacid_filter(self,aminoacidos,aminoacidos_to_exclude):

        if aminoacidos != ['\n']:
            for i in aminoacidos:
                aminoacido = i.replace('\n','')
                sort_by_amino = [bool(re.search(fr'{aminoacido}[0-9]*{modificacion_final}', str(i))) for i in self.table['Modifications'].to_list() ]	#crea lista en funcion de si la PSM tiene el aminoacido modificado(true) o no(false)
                self.table['amino_sort'+aminoacido]=sort_by_amino
        
            vector=np.array(self.table.iloc[:,-(len(aminoacidos)):])
            lista=[]
            for i in range(len(vector)):
                if any(vector[i]):
                    lista.append(i)
                    
            self.table=self.table.reset_index().loc[lista,:]

        else:
            print('Error, you don\'t fill the aminoacid_list in the config file'+'\n')
            print('You will not filter anything!!!')

        for i in aminoacidos_to_exclude:
            if i != '':
                aminoacido_2 = i.replace('\n','')
                sort_by_amino_2 = [bool(re.search(fr'{aminoacido_2}[0-9]*{modificacion_exl_final}', i)) for i in self.table['Modifications'].to_list() ]#crea lista en funcion de si la PSM tiene el aminoacido modificado(true) o no(false)
                self.table['amino_exclude' + aminoacido_2]=sort_by_amino_2
                self.table = self.table.loc[self.table['amino_exclude' + aminoacido_2]==False]
                
                return self.table
            else:
                return self.table

#define the length of the dataframe

    def __len__(self):
        return len(self.table)

#reads the config file and defines variables

try:
    for i in range(len(sys.argv)):
        if '-c' in sys.argv[i]:
            configquotes=sys.argv[i][2:]
    config=configquotes.replace('"','')
    with open(config) as f:
        L=f.readlines()
    f.close()
    for i in L:
        if "infile" in i:
            Path=i[8:-1]
        elif "outfile" in i:
            pathFolder= i[9:-1]
        elif "finalname" in i:
            fileName = i[11:-1]
        elif "modification_to_count" in i:
            modificacion= i[23:-1]
            modificacion_final = "\(" + modificacion + "\)"  
        elif 'amino_acids_to_exclude' in i:
            aminoacidos_to_exclude = i[24:-1].split(',')
        elif "FDR:" in i:
            FDR_value = float(i[5:-1])
        elif 'aminoacids_list' in i:
            aminoacidos = i[17:].split(',')
        elif 'modification_to_exclude' in i:
            modificacion_exl=i[25:-1]
            modificacion_exl_final = "\(" + modificacion_exl + "\)"
        elif 'decoy_prefix' in i:
            decoy_indicator=i[14:-1]

#makes an unique dataframe with all the files PSMs.txt
    try:
        entries = np.array(os.listdir(Path))
        ficheros = []
        psms = []
        for i in entries:
            ficheros.append(Path+"\\"+ i)
        for i in ficheros:
            if 'PSMs.txt' in i :
                psms.append(i)
                
        for i in range(len(psms)):
            if i == 0:
                data_frame_uniq=pd.read_csv(psms[i],sep="\t")
            else:
                data_frame=pd.read_csv(psms[i],sep="\t")
                data_frame_uniq=pd.concat([data_frame_uniq,data_frame],axis=0)

        df=dataframe(data_frame_uniq)
        df_global=dataframe(data_frame_uniq)


#pipeline to calculate FDR local and filter
        df.xcorr()
        df.aminoacid_filter(aminoacidos,aminoacidos_to_exclude)
        df.fdr_filter(FDR_value).to_csv(os.path.join(pathFolder, fileName), sep='\t', index = False)


#pipeline to calculate FDR global and filter
        df_global.xcorr()
        df_global.fdr_filter(FDR_value)
        df_global.aminoacid_filter(aminoacidos,aminoacidos_to_exclude)

#defines the nº of PSMs in each case and makes the resume.txt file
        a=len(df)
        b=len(df_global)

        with open(pathFolder+'\\SUMMARY.txt','w') as w:
            w.write('Nº PSMs with local FDR: '+str(a) +'\n')
            w.write('Nº PSMs with global FDR: '+ str(b)+'\n')
            w.write('Percentage: '+ str(round(a*100/b,2))+'%')
    
    except ValueError or NameError:
        print('Check you config file, maybe something is missing or not well filled')


except FileNotFoundError:
    print('Maybe you don\'t attach a path')