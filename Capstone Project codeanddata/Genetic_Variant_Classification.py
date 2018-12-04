
# coding: utf-8

# In[79]:


# import the necessary packages 
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 
sns.set_style(style='whitegrid')
sns.set(font_scale = 1.5)
import warnings 
warnings.filterwarnings('ignore')
import re
from IPython.display import display


# In[80]:


import ipynb.fs.full.Functions_Variant_Classification as m


# In[81]:


df = m.dataimport("clinvar_conflicting.csv")


# In[82]:


m.datadisplay(df, 10)


# In[83]:


# Checking the number of NAs 
#df.info()
df.isna().sum()


# In[84]:


#It seems that the Nas are related to the MODIFIER IMPACT and down/up stream gene variant  
display(df[(df['INTRON'].isnull()) & (df['EXON'].isnull())].head())


# In[85]:


#Getting a general idea about the correlation 
df.corr()
plt.figure(figsize= (30, 20))
sns.heatmap(df.corr(), annot = True)
plt.tight_layout()
plt.savefig('corrmat.jpg')
plt.show()


# In[86]:


# Checking the target variable distribution
ax = sns.countplot(x='CLASS', data=df)
ax.set(xlabel='CLASS',ylabel='Number of Variants');
plt.savefig('classdist.jpg')
plt.show()


# In[87]:


# Replacing and specifying the types of the variables 
small_df = m.preprocessing1(df)


# In[88]:


small_df.head()


# In[89]:


# Mapping the ordinal variables and specifying the right variable type

small_df = m.preprocessing2(small_df)


# In[90]:


small_df.head()


# In[91]:


# Plotting the correlation matrice of some of the variables in the processed data frame
cor_data = small_df[['AF_ESP', 'AF_EXAC', 'AF_TGP','Protein_position', 'CDS_position', 'cDNA_position','PolyPhen', 'CADD_PHRED', 'POS', 'INTRON', 'EXON', 'CLASS']]
plt.figure(figsize= (15, 8))
sns.heatmap(round(cor_data.corr(), 3), annot = True)
plt.show()


# In[92]:


cor_data2 = small_df[['AF_ESP', 'AF_EXAC', 'AF_TGP','Protein_position', 'CDS_position', 'cDNA_position','CLASS']]
plt.figure(figsize= (10, 6))
sns.heatmap(round(cor_data2.corr(), 3), annot = True)
plt.show()


# In[93]:


cor_data.head()


# In[94]:


# Showing the historgram of some of the important variables 

# cols = ['EXON', 'INTRON', 'CLASS','IMPACT', 'BAM_EDIT', 'LoFtool', 'CADD_PHRED']
# small_df[cols]
cor_data.hist(figsize=(10, 10))
plt.suptitle('Figure 2 - The Histogram plot of some features ',x=0.5, y=1.01, verticalalignment='top', fontsize= 18)
plt.tight_layout()
plt.savefig('histvar.jpg')
plt.show();


# In[95]:


# Checking the variables which have more than 2000 unique value. Some of these variables are just identifiers.
m.check_unique(small_df)


# In[96]:


# Dropping the unnecessary and redundant variables 
#clean_df = small_df.drop(['CLNDISDB', 'CLNDN', 'CLNHGVS', 'CLNVI', 'SYMBOL', 'Feature', 'Codons', 'CADD_RAW', 'MC','Feature_type','Amino_acids','AF_ESP', 'AF_EXAC', 'AF_TGP','Protein_position', 'CDS_position', 'cDNA_position', 'SIFT', 'PolyPhen'],axis=1)
clean_df = small_df.drop(['CLNDISDB', 'CLNDN', 'CLNHGVS', 'CLNVI', 'SYMBOL', 'Feature','AF_TGP','AF_ESP','Protein_position', 'CDS_position',  'Codons', 'CADD_RAW', 'MC','Feature_type','Amino_acids'],axis=1)


# In[97]:


# # clean_df2, clean_df1 for bayesian analysis( I used CHROM variable to create a hierarchical model)
# clean_df1 = clean_df 
# clean_df1 ['CHROM'].replace(['X', 'MT'], ['23', '24'], inplace=True)
# clean_df1['CHROM'] = clean_df1['CHROM'].astype('int64')
# clean_df1 = pd.get_dummies(clean_df1[['AF_EXAC','AF_TGP','CADD_PHRED','POS','LoFtool','cDNA_position','EXON_LEN','AF_ESP',
#  'CDS_position','Protein_position','IMPACT','BLOSUM62','PolyPhen','EXON', 'INTRON','CHROM',
#  'INTRON_LEN', 'BAM_EDIT', 'CLASS']])


# In[98]:


#Finally, applying one-hot encoding

clean_df = pd.get_dummies(clean_df)
print("New  shape after one-hot encoding:" , np.shape(clean_df))


# In[99]:


# Selecting a subset of the dataset for grid search

clean_gs = clean_df.sample(n=20000, random_state=0)


# In[100]:


# splitting and scaling the data
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import scale
from sklearn.metrics import classification_report, accuracy_score, roc_curve, auc, confusion_matrix



# In[101]:


# Selecting the best hyperparameters for Random Forest

from sklearn.ensemble import RandomForestClassifier as rfc


# In[102]:


# top 15 relevent features 

plt.figure(figsize=(20, 20))
    
results = m.get_feature_importance(rfc, clean_df)

m.plt_tpfeat(results[:15], sns)


# In[103]:


len(results[results['Score']> 0.001])


# In[104]:


# subsetting the data for future prediction 

top_feat = results['Feature'][:65].tolist()
top_feat.extend(['CLASS'])
#clean_df1 = clean_df[top_feat]


# In[105]:


# train and test for predition 

X_train, X_test, y_train, y_test = m.spt(clean_df[top_feat] ,train_test_split)
X_train, X_test = m.std(scale, X_train, X_test)


# In[106]:


from sklearn.model_selection import GridSearchCV
#train and test sets for grid search
X_train_gs, X_test_gs, y_train_gs, y_test_gs = m.spt(clean_gs[top_feat] ,train_test_split)
X_train_gs, X_test_gs = m.std(scale, X_train_gs, X_test_gs)

m.sel_rf(rfc,GridSearchCV, X_train_gs, y_train_gs)


# In[107]:


# predicting the target using random forest 

y_pred1 = m.rf_cla(rfc,X_train, X_test, y_train, y_test)

model= {}
model["Random Forest"] = y_pred1

print(confusion_matrix(y_test, y_pred1))
print(classification_report(y_test, y_pred1))
print(accuracy_score(y_test, y_pred1))


# In[108]:


#Selecting the best hyperparameters for XGBClassifier
from xgboost import XGBClassifier as xgb

m.sel_xgb(xgb,GridSearchCV, X_train_gs, y_train_gs)


# In[109]:


# predicting the target using XGBClassifier
y_pred2 = m.xgb_cla(xgb, X_train, X_test, y_train, y_test)

model["XGBoosting"] = y_pred2

print(confusion_matrix(y_test, y_pred2))
print(classification_report(y_test, y_pred2))
print(accuracy_score(y_test, y_pred2))


# In[110]:


# Selecting the best hyperparameters for LogisticRegression

from sklearn.linear_model import LogisticRegression as lr

m.sel_lr(lr,GridSearchCV ,X_train_gs, y_train_gs)


# In[111]:


# predicting the target using LogisticRegression
y_pred3 = m.lr_cla(lr, X_train, X_test, y_train, y_test)


model["Logistic Regression"] = y_pred3
print(confusion_matrix(y_test, y_pred3))
print(classification_report(y_test, y_pred3))
print(accuracy_score(y_test, y_pred3))


# In[112]:


from sklearn.linear_model import SGDClassifier

m.sel_sgd(SGDClassifier,GridSearchCV ,X_train, y_train)


# In[113]:


# predicting the target using SGDClassifier
y_pred4 = m.sgd_cla(SGDClassifier,X_train, X_test, y_train, y_test)

model["SGDClassifier"] = y_pred4

print(confusion_matrix(y_test, y_pred4))
print(classification_report(y_test, y_pred4))
print(accuracy_score(y_test, y_pred4))


# In[114]:


# predicting the target using GaussianNB

from sklearn.naive_bayes import GaussianNB

y_pred5 = m.nb_cla(GaussianNB,X_train, X_test, y_train, y_test)

model["GaussianNB"] = y_pred5

print(confusion_matrix(y_test, y_pred5))
print(classification_report(y_test, y_pred5))
print(accuracy_score(y_test, y_pred5))


# In[115]:


# saving the results to a file
# Plotting the roc_curv of the different models 
plt.figure(figsize=(6, 6))

m.plot_roc_curve(roc_curve,auc, y_test, model)


# In[116]:


model = pd.DataFrame(model)
model.to_csv('results')


# In[117]:


# Selecting the best hyperparameters for SVC
#from sklearn.svm import SVC

#m.sel_svc(SVC,GridSearchCV ,X_train_gs, y_train_gs)


# In[118]:


# predicting the target using SVC

# y_pred = m.svc_cla(SVC, X_train, X_test, y_train, y_test)

# classification_report(y_test, y_pred)
# accuracy_score(y_test, y_pred)
# classification_report(y_test, y_pred)


# In[119]:


# data preparation for the bayesian part 

# X_train1, X_test1, y_train1, y_test1 = m.spt(clean_df1,train_test_split)


# In[120]:


# # Bayesian model 
# import pystan

# model = m.cr_model()
# data = m.model_par(X_train1, X_test1, y_train1)
# y_pred = m.bas_model(model, data, pystan)

# y_pred = y_pre(y_pred)

# print(classification_report(y_test, y_pred))
# print(accuracy_score(y_test, y_pred))


# In[ ]:




