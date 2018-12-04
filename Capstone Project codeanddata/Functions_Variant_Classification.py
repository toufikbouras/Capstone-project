
# coding: utf-8

# In[1]:


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


# In[2]:


# importing the data
def dataimport(data):
    df = pd.read_csv(data)
    return(df)


# In[3]:


# import the data and show a subset of it 
#df = pd.read_csv("clinvar_conflicting.csv")
def datadisplay(df, n):
    
    pd.options.display.max_columns = None
    display(df.tail(n))


# In[4]:


# Some location variables have two values so on this function I chose only one variable of the two.
def choselocvar(column):
    di = []
    column = column.astype('str')
    for i in range(len(column)):
        if (re.search(r'-', column[i])):
            if(column[i].split('-')[0] != '?'):
                di.append(column[i].split('-')[0])
            else:
                di.append(column[i].split('-')[1])
        else:
            di.append(column[i])   
    return(di)


# In[5]:


# Removing the columns with 90% NAs
def selectdata(df):
    df = df.dropna(axis = 1, thresh= len(df) - 0.90* len(df))
    return(df)
    


# In[6]:


# Replacing questionmarks with the mean of the adjacent values 
def replaceqesmark(column):
    column = column.astype('str')
    
    for i in range(len(column)):
        
        if (re.search(r'\?', column[i]) == True):
            column[i] = str((int(column[i+1]) + int(column[i-1]))/2)
                
    return(column)


# In[7]:


# Replacing the NAs with the closest possible value
def replacenadj(column):
    column = column.astype('float')
    k = len(column) - 1
    d = []
    for i in range(k):
    
        if (column[i]== 0):
        
            if((column[i+1]!= 0) & (column[i-1]!= 0)):
                a = column[i+1]
                b = column[i-1]
                d.append((a + b)/2)
                
            elif((column[i + 1]== 0) & (column[i-1]!= 0)):

                b = column[i-1]
                d.append(b + 1)
            elif((column[i + 1]!= 0) & (column[i-1]== 0)):
            
                a = column[i+1]
                d.append(a - 1)
            else:
                d.append(0)
        else:
            d.append(column[i])
        
    d.append(column[65187])
    return(d)


# In[8]:


# This function is prepared to fill NA values for cadd('CADD_PHRED'). I noticed that the CADD
# PHRED is somehow correlated to impact so I used it to fill na.
def caddnareplace(small_df, colname):
    ta = []
    ind = small_df[small_df[colname].isnull()].index
    for i in range(len(small_df[colname])):
    
        if (i not in ind):
            ta.append(small_df[colname][i])
        else:
            
            if (small_df['IMPACT'][i] == 'HIGH'):
                ta.append(small_df[small_df['IMPACT'] == 'HIGH'][colname].mean())
    
            elif(small_df['IMPACT'][i] == 'MODIFIER'):
                ta.append(small_df[small_df['IMPACT'] == 'MODIFIER'][colname].mean())
     
            elif(small_df['IMPACT'][i] == 'MODERATE'):
                ta.append(small_df[small_df['IMPACT'] == 'MODERATE'][colname].mean())
        
            else:
                ta.append(small_df[small_df['IMPACT'] == 'LOW'][colname].mean())
    return(ta)
     


# In[9]:


# The EXON and INTRON columns are written as string ratios so I need to transform them into values 
def transratio(column) :
    di = []
    de = []
    for i in range(len(column)):
        if (column[i] != '0'):
            di.append([int(s) for s in (column[i].split('/'))][0])
            de.append([int(s) for s in (column[i].split('/'))][1])
        else:
            di.append(0)            
            de.append(0)
    return(di, de) 


# In[10]:


# Replace the dash sign in Allele, REF, ALT columns with NO
def transgene(column):
    d = []
    for i in range(len(column)):
        if (re.search(r'-', column[i])):
            d.append('No')
        elif(len(column[i])> 1):
            d.append('Mo') 
        else:
            d.append(column[i])
    return(d)


# In[11]:


def preprocessing1(df):
        
    small_df = selectdata(df)
# Applying the function of chosing location of the varriables
    small_df['CDS_position'] = list(choselocvar(small_df['CDS_position']))
    small_df['Protein_position'] = list(choselocvar(small_df['Protein_position']))
    small_df['cDNA_position'] = list(choselocvar(small_df['cDNA_position']))
    
# Replacing the NAs with 0 values as strings. Later we will replace the 0 with an estimated value
    small_df['cDNA_position'] = small_df['cDNA_position'].replace('nan', '0')
    small_df['Protein_position'] = small_df['CDS_position'].replace('nan', '0')
    small_df['CDS_position'] = small_df['Protein_position'].replace('nan', '0')
    
# Aplying the function of replacing question marks with the mean of the adjuscent values  
    small_df['cDNA_position'] = replaceqesmark(small_df['cDNA_position'])
    small_df['CDS_position'] = replaceqesmark(small_df['CDS_position'])
    small_df['Protein_position'] = replaceqesmark(small_df['Protein_position'])
    
# Replacing some NA values with the mean of the adjuscent values
    small_df['cDNA_position'] = list(replacenadj(small_df['cDNA_position']))
    small_df['CDS_position'] = list(replacenadj(small_df['CDS_position']))
    small_df['Protein_position'] = list(replacenadj(small_df['Protein_position']))
    
# filling the rest of NAs with the forward filling function
    small_df['cDNA_position'] = small_df['cDNA_position'].replace(0, np.nan).fillna(method='ffill') 
    small_df['CDS_position'] = small_df['CDS_position'].replace(0, np.nan).fillna(method='ffill') 
    small_df['Protein_position'] = small_df['Protein_position'].replace(0, np.nan).fillna(method='ffill')

# Replacing NAs with the adequate values
    small_df['CADD_PHRED'] = caddnareplace(small_df, 'CADD_PHRED')
    small_df['EXON'].fillna('0', inplace=True)
    small_df['INTRON'].fillna('0', inplace=True)
   
# transforming the values of EXONs and INTRONs
    a, b= transratio(small_df['EXON'])
    c, d = transratio(small_df['INTRON'])
    small_df['EXON'] = list(a)
    small_df['EXON_LEN'] = list(b)
    small_df['INTRON'] = list(c)
    small_df['INTRON_LEN'] = list(d)
    # Changing three clumns into numbers
    small_df['Allele'] = list(transgene(small_df['Allele']))
    small_df['REF'] = list(transgene(small_df['REF']))
    small_df['ALT'] = list(transgene(small_df['ALT']))
    
# filling NAs with the relevent values and type values 
    small_df['BLOSUM62'].fillna(0, inplace = True)# 
    small_df['BIOTYPE'].fillna('No', inplace= True)# not specified category
    small_df['BAM_EDIT'].fillna('No', inplace= True)#
# Similarly Amino_acids are the proteins themselves  
    small_df['Amino_acids'].fillna('No', inplace = True)
    small_df['Codons'].fillna('No', inplace = True)# this should be dropped
# I will fill the na values with no value
    small_df['PolyPhen'].fillna('No', inplace = True)
    small_df.SIFT = small_df['SIFT'].fillna('No', inplace= True)

## It seems that LOFtool is an important measure in our study, when it is null I think that mean 
#there is no "Lost of function". Therefore, I will be putting "0" value in the NAs
    small_df['LoFtool'].fillna('0', inplace= True)
    small_df['LoFtool'] = small_df['LoFtool'].astype('float')

    small_df['CADD_RAW'] = small_df['CADD_RAW'].astype('float')
    small_df['CADD_PHRED'] = small_df['CADD_PHRED'].astype('float')
    small_df['LoFtool'] = small_df['LoFtool'].astype('float')
# Specifying the type of the different variables and their relevent values
    small_df['STRAND'] = small_df['STRAND'].astype('str')

    small_df['ORIGIN'].fillna(0, inplace= True)
    small_df['ORIGIN'] = small_df['ORIGIN'].astype('str')

    small_df['Allele'] = small_df['Allele'].astype('str')
    small_df['REF'] = small_df['REF'].astype('str')
    small_df['ALT'] = small_df['ALT'].astype('str')

 
    return(small_df)


# In[12]:


def preprocessing2(small_df):    
    
    BIOTYPE_map = {'No':'0','misc_RNA':'1','protein_coding':'2'}
    small_df.BIOTYPE = small_df.BIOTYPE.map(BIOTYPE_map)
    small_df.BIOTYPE = small_df.BIOTYPE.astype('int64')

    BAM_EDIT_map = {'No':'0','OK':'1','FAILED':'2'}
    small_df.BAM_EDIT = small_df.BAM_EDIT.map(BAM_EDIT_map)
    small_df.BAM_EDIT = small_df.BAM_EDIT.astype('int64')

    IMPACT_map = {'MODERATE':'6','LOW':'2','HIGH':'16','MODIFIER':'1'}
    small_df.IMPACT = small_df.IMPACT.map(IMPACT_map)
    small_df.IMPACT = small_df.IMPACT.astype('int64')

    PolyPhen_map = {'benign':'1','probably_damaging':'2','possibly_damaging':'3','unknown':'4', 'No': '0'}
    small_df.PolyPhen = small_df.PolyPhen.map(PolyPhen_map)
    small_df.PolyPhen = small_df.PolyPhen.astype('int64')

    SIFT_map = {'tolerated':'1','deleterious':'3','deleterious_low_confidence':'2', 'No':'0'}
    small_df.SIFT = small_df.SIFT.map(SIFT_map).astype('float')
    small_df.SIFT = small_df.SIFT.fillna(0) 
    small_df.SIFT = small_df.SIFT.astype('int64')
    return small_df


# In[13]:


# Looking into the varibales with a large unique values

def check_unique(small_df):
    cols = small_df.columns
    li = []
    for i in range(len(cols)):
        if(len(small_df[cols[i]].unique()) > 2000):
            li.append(cols[i])
    return(li)


# In[14]:


# splitting the data

def spt(clean_df,train_test_split):
    X = clean_df.drop('CLASS', axis=1)
    y = clean_df['CLASS']
    X_train, X_test, y_train, y_test = train_test_split(X,y,test_size = 0.2, random_state=0)

    return X_train, X_test, y_train, y_test


# In[15]:


# scaling the data

def std(scale,X_train, X_test):
    
    X_train = scale(X_train)
    X_test = scale(X_test)
    return X_train, X_test


# In[16]:


# Selecting the best hyperparameters for Random Forest

def sel_rf(rfc,GridSearchCV, X_train, y_train):
    param_grid = {'n_estimators': np.arange(10, 100, 10),
                  'max_features':['auto', 'log2', None],
                  'criterion': ['gini', 'entropy'],
                  'max_depth': np.arange(2, 20)
                  }

    gr_forest = GridSearchCV(rfc(), param_grid, cv=5, scoring='accuracy', n_jobs= -1)

    gr_forest.fit(X_train, y_train)

    print(gr_forest.best_estimator_)
    print(gr_forest.best_score_)


# In[17]:


# predicting the target using random forest 

def rf_cla(rfc, X_train, X_test, y_train, y_test):
    classifier = rfc(bootstrap=True, class_weight=None, criterion='entropy',
            max_depth=14, max_features=None, max_leaf_nodes=None,
            min_impurity_decrease=0.0, min_impurity_split=None,
            min_samples_leaf=1, min_samples_split=2,
            min_weight_fraction_leaf=0.0, n_estimators=50, n_jobs=-1,
            oob_score=False, random_state=None, verbose=0,
            warm_start=False)

    classifier.fit(X_train, y_train)
    predictions = classifier.predict(X_test)
    
    return predictions

    


# In[18]:


# sort the features according to their importance 

def get_feature_importance(clsf, data):
    
    X = data.drop('CLASS', axis=1)
    y = data['CLASS']
    clsf = clsf(bootstrap=True, class_weight=None, criterion='gini',
            max_depth=16, max_features=None, max_leaf_nodes=None,
            min_impurity_decrease=0.0, min_impurity_split=None,
            min_samples_leaf=1, min_samples_split=2,
            min_weight_fraction_leaf=0.0, n_estimators=70, n_jobs=-1,
            oob_score=False, random_state=0, verbose=0,
            warm_start=False)
    
    clsf.fit(X, y) 
    imp = clsf.feature_importances_.tolist()
    feat = X.columns
    result = pd.DataFrame({'Feature':feat,'Score':imp})
    result = result.sort_values(by=['Score'],ascending=False)
    return result


# In[19]:


# Ploting the top 15 features 
# def plt_tpfeat(results, plt):
    
#     features = results['feat']
#     importances = results['score']

#     plt.figure(figsize=(30, 30))
#     plt.title('Feature Importances')
#     plt.barh(range(len(results)), importances, color='b', align='center')
#     plt.yticks(range(len(results)), features)
#     plt.xlabel('Relative Importance')
#     plt.gca().invert_yaxis()
#     plt.savefig('top_features.jpg')
#     plt.show()

def plt_tpfeat(results, sns):
    
    features = results['Feature']
    importances = results['Score']
    sns.barplot(x="Score", y="Feature", data=results)
    plt.title('Feature Importances')
    plt.xlabel('Relative Importance')
    plt.savefig('topfeatures.jpg')
    plt.show()


# In[20]:


# Selecting the best hyperparameters for xgboost

def sel_xgb(xgb, GridSearchCV, X_train, y_train):
    
    param_grid = {'learning_rate': [0.001, 0.01, 0.1, 1],
                 'max_depth': np.arange(4,10),
                 'n_estimators': np.arange(50, 400, 50)
                  }
    gr_xgb = GridSearchCV(xgb(), param_grid, cv=3, scoring='accuracy', n_jobs= -1)

    gr_xgb.fit(X_train, y_train)
    print(gr_xgb.best_estimator_)
    print(gr_xgb.best_score_)


# In[21]:


# predicting the target using xgboost  

def xgb_cla(xgb, X_train, X_test, y_train, y_test):
    
    classifier = xgb(base_score=0.5, booster='gbtree', colsample_bylevel=1,
       colsample_bytree=1, gamma=0, learning_rate=0.01, max_delta_step=0,
       max_depth=8, min_child_weight=1, missing=None, n_estimators=400,
       n_jobs=1, nthread=None, objective='binary:logistic', random_state=0,
       reg_alpha=0, reg_lambda=1, scale_pos_weight=1, seed=None,
       silent=True, subsample=1)

    classifier.fit(X_train, y_train)

    prediction = classifier.predict(X_test)

    return prediction


# In[22]:


def sel_lr(lr,GridSearchCV ,X_train, y_train):
    
    param_grid = {'C': [0.001, 0.01, 0.1, 1, 10, 100],
                 'tol': [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1],
                 }

    gr_lg = GridSearchCV(lr(), param_grid, cv=5, scoring='accuracy', n_jobs= -1)

    gr_lg.fit(X_train, y_train)

    print(gr_lg.best_estimator_)
    print(gr_lg.best_score_)


# In[23]:


def lr_cla(lr, X_train, X_test, y_train, y_test):
    
    classifier = lr(C=0.001, class_weight=None, dual=False, fit_intercept=True,
          intercept_scaling=1, max_iter=100, multi_class='ovr', n_jobs=1,
          penalty='l2', random_state=None, solver='liblinear', tol=1e-06,
          verbose=0, warm_start=False)

    classifier.fit(X_train, y_train)

    prediction = classifier.predict(X_test)

    return prediction


# In[24]:


# slect the best parameters for svc
def sel_svc(SVC,GridSearchCV ,X_train, y_train):

    param_grid = {'C': np.arange(1,10),
                 'tol': [0.0001, 0.001, 0.01, 0.1, 1],
                 'kernel': ['linear', 'rbf', 'poly', 'sigmoid'],
                 'random_state': [0]}

    gr_SVC = GridSearchCV(SVC(), param_grid, cv=5, scoring='accuracy', n_jobs= -1)

    gr_SVC.fit(X_train, y_train)

    print(gr_SVC.best_estimator_)
    print(gr_SVC.best_score_)
    


# In[25]:


# predict the target variable using svc
###############################################3333
def svc_cla(svc, X_train, X_test, y_train, y_test):
    
    classifier = svc(C=1, cache_size=200, class_weight=None, coef0=0.0,
                      decision_function_shape='ovr', degree=3, gamma='auto', kernel='poly',
                      max_iter=-1, probability=False, random_state=0, shrinking=True,
                      tol=0.0001, verbose=False)

    classifier.fit(X_train, y_train)

    prediction = classifier.predict(X_test)

    return prediction


# In[26]:


def sel_sgd(sgd,GridSearchCV ,X_train, y_train):
    
    param_grid = {
    'loss': ['log'],
    'penalty': ['elasticnet'],
    'alpha': [0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1],
    'l1_ratio': [0, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 1],
    }
    
    gr_sgd = GridSearchCV(sgd(), param_grid, cv=5, scoring='accuracy', n_jobs= -1)

    gr_sgd.fit(X_train, y_train)

    print(gr_sgd.best_estimator_)
    print(gr_sgd.best_score_)


# In[27]:


def sgd_cla(sgd, X_train, X_test, y_train, y_test):
    classifier = sgd(alpha=0.01, average=False, class_weight=None, epsilon=0.1,
               eta0=0.0, fit_intercept=True, l1_ratio=0.2, learning_rate='optimal',
               loss='log', max_iter=None, n_iter=None, n_jobs=-1,
               penalty='elasticnet', power_t=0.5, random_state=0, shuffle=True,
               tol=None, verbose=0, warm_start=False)

    classifier.fit(X_train, y_train)

    prediction = classifier.predict(X_test)

    return prediction


# In[28]:


def nb_cla(gnb, X_train, X_test, y_train, y_test):
    
    classifier = gnb(priors=None)
    
    classifier.fit(X_train, y_train)

    prediction = classifier.predict(X_test)

    return prediction
    


# In[29]:


# plotting ROC-Curves


#from sklearn.metrics import roc_curve, auc

def plot_roc_curve(roc_curve,auc, y_test, model):


    for name, y_pred in model.items():

        fpr, tpr, thresholds = roc_curve(y_test, y_pred)

        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, label='%s: %0.2f' % (name, roc_auc), linewidth=5)

    plt.title('ROC Curve', fontsize=14)
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([-0.1, 1.1])
    plt.ylim([-0.1, 1.1])
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.legend(loc=2, prop={'size': 8})

    plt.savefig('roc_curves.jpg')
    plt.show()


# In[30]:


# Bayesian model specifications 

def cr_model():
    
    model_specification = """
    data {

        int I;// Number of covatriates (exon and intron)
        int H; // Number of covatriates (the other important variables)
        int N; // Number of individuals(training)
        int M; // Number of individuals(testing)
        int<lower=1> JCHR; // numer of groups
        int<lower=1, upper=JCHR> chrom[N];// indv_to_subgroup
        int<lower=1, upper=JCHR> chrom_p[M];// indv_to_subgroup

        int<lower=0, upper=1> y[N];

        row_vector[I] X[N];
        matrix[N, H] X2;

        row_vector[I] X_p[M];
        matrix[M, H] X2_p;

    }
    parameters{

        vector[I] beta[JCHR];
        vector<lower=0>[JCHR] alpha;// one intercept for every chromosome
        vector[H] beta2;
        vector<lower=0>[N] alpha2;

        real mu_alpha;
        real<lower=0> sigma_alpha;
        real mu_beta;
        real<lower=0> sigma_beta;
    }    

    transformed parameters {
      vector[N] y_hat;
      vector[N] m;
      vector[M] mu_pred;
      vector[M] m_p;

      for (i in 1:N) {
        m[i] = X[i]* beta[chrom[i]] + alpha[chrom[i]];
      }
      y_hat = m +  X2*beta2 + alpha2;


      for (p in 1: M){
          m_p[p] = X_p[p]* beta[chrom_p[p]] + alpha[chrom_p[p]];
      }
      mu_pred = m_p + X2_p*beta2 + mean(alpha2);
    }

    model {

        for(l in 1:JCHR){
            beta[l] ~ normal(mu_beta,sigma_beta);
            alpha[l] ~ normal(mu_alpha, sigma_alpha);
        }

        beta2 ~ normal(0 , 10);
        alpha2 ~ cauchy(0, 10);
        mu_alpha ~ normal(0, 10);
        sigma_alpha ~ cauchy(0, 10);
        mu_beta ~ normal(0, 10);
        sigma_beta ~ cauchy(0, 10);

        for (i in 1:N)
            y[i] ~ bernoulli_logit(y_hat[i]);
    }
    generated quantities{
        vector[M] y_pred;
        for (i in 1:M)
            y_pred[i] = bernoulli_logit_rng(mu_pred[i]);
    }


    """
    return model_specification


# In[31]:


# the specification od the data for the model

def model_par(X_train1,X_test1, y_train1):
    
    X_train2 = X_train1.drop(['EXON', 'INTRON'], axis=1)
    X_test2 = X_test1.drop(['EXON', 'INTRON'], axis=1)
    
    details_of_data = {   'I': 2,
                          'H':len(X_train2.columns),
                          'N': len(X_train1),
                          'M':len(X_test1),
                          'JCHR': 24,
                          'chrom': X_train1.CHROM, 
                          'chrom_p':X_test1.CHROM,
                          'X':X_train1[['EXON', 'INTRON']],
                          'X2':X_train2,
                          'X_p':X_test1[['EXON', 'INTRON']],
                          'X2_p': X_test2,
                          'y': y_train1
                    }
    return  details_of_data


# In[32]:


# running the model 
def bas_model(model_specification, details_of_data, pystan):
    
    varying_intercept_fit = pystan.stan(model_code=model_specification, data=details_of_data, iter=1000, chains=2)
    y_pred = varying_intercept_fit['y_pred']
    return y_pred


# In[33]:


def y_pre(y_pred):
    y_pred = [sum(x) for group in zip(y_pred) for x in zip(group)]
    y_pred = [ 0 if (i/1500)<0.5 else 1 for i in sum(y_pred)]
    return y_pred


# In[ ]:



 
    

