############################ Primary data cleaning code #############################
#####################################################################################


################### defining function

def recode_value_Shock(value):
    if value < 1.25:
        return 0
    else:
        return 1
    
def recode_value_AIS(value):
    if value < 3:
        return 0
    else:
        return 1

def recode_value_Age(value):
    if value >=65:
        return 1
    elif value <= 35:
        return 0
    else:
        return 2
    
def recode_value_Severity(value):
    if value < 25:
        return 0
    else:
        return 1
    
def recode_value_GCS_bi(value):
    if value < 9:
        return 1
    else:
        return 0
######################################################
############################################################################### working on omic sub-cohort


address='SWAT_proteins_clinical.csv'


#reading the excel file as a dataframe
df0=pd.read_csv(address)

##################################################################### missing value imputation with mice
    
kds = mf.ImputationKernel(df0.loc[(df0['Sample']==sample_set),'phgcs':'Death_by_Other'],save_all_iterations=True,
                              random_state=100)

## Run the MICE algorithm for 30 iterations
kds.mice(30)
df_imputed = kds.complete_data()
X_imputed=np.array(np.round(df_imputed,3))
    
################################################################# ending imupation and cleaning   

df0.loc[(df0['Sample']==sample_set),'phgcs':'Death_by_Other']=X_imputed

df0['Shock']=np.round(df0['phhr']/df0['phsbp'],2)

#############################################

df0[['Biological_sex','PhTXA','Plasma','LTOWB','PRBC','Platelet','Cryo','Head_AIS','Face_AIS',
     'ISS','TBI','phcpr','phgcs','phsbp','phhr','phrr','Hypotention','Penetrating',
     'Lowest_ED_SBP','ED_intubation','ED_CPR','Ph_ED_Hypotension','ED_mortality',
     'Hospital_days','Ventilator_days','30d_Mortality','ICU_days','Pre_Matched','Time point',
     'Fluids','Total_blood',
     'Plasma_units','PRBC_units','LTOWB_units','Cryo_units','Platelet_units','Platelets']]=\
np.round(df0[['Biological_sex','PhTXA','Plasma','LTOWB','PRBC','Platelet','Cryo','Head_AIS','Face_AIS',
     'ISS','TBI','phcpr','phgcs','phsbp','phhr','phrr','Hypotention','Penetrating',
     'Lowest_ED_SBP','ED_intubation','ED_CPR','Ph_ED_Hypotension','ED_mortality',
     'Hospital_days','Ventilator_days','30d_Mortality','ICU_days','Pre_Matched','Time point',
     'Fluids','Total_blood',
     'Plasma_units','PRBC_units','LTOWB_units','Cryo_units','Platelet_units','Platelets']],0)


df0['Head_AIS_bi']=df0['Head_AIS'].apply(recode_value_AIS)
df0['Shock_bi']=df0['Shock'].apply(recode_value_Shock)
df0['Age_group']=df0['Age'].apply(recode_value_Age)
df0['Severity']=df0['ISS'].apply(recode_value_Severity)
df0['Prehospital_GCS']=df0['phgcs'].apply(recode_value_GCS_bi)
df0['Outcome'] = df0.apply(lambda row: 1 if (row['ICU_days'] > 10) or (row['30d_Mortality'] == 1) else 0, axis=1)
new_cats=['Outcome','Head_AIS_bi','Shock_bi','Age_group','Severity','Prehospital_GCS']



df0['WB_vs_Plasma'] = df0.apply(lambda row: 1 if ((row['LTOWB']==1) and (row['Plasma']==0))  \
                                              else 0 if ((row['LTOWB']==0) and (row['Plasma']==1)) else 2, axis=1)


##################################################
#################################################################################### Working on full SWAT cohort file 

dfmodeling=pd.read_csv('Modeling_CSV.csv')


############# identifying missing values for coagulation parameters before imputation

PT0_misses=dfmodeling[dfmodeling['PT_0h'].isnull()].index.tolist()
INR0_misses=dfmodeling[dfmodeling['INR_0h'].isnull()].index.tolist()
Platelets0_misses=dfmodeling[dfmodeling['Platelets_0h'].isnull()].index.tolist()

################################ missing value imputation (after eliminating individuals with high missing, see associated manuscript methods)

impute_sets=1
imputed_sets=[]
kds_model = mf.ImputationKernel(dfmodeling.loc[:,'30d_Mortality':'External_AIS'],save_all_iterations=True,random_state=100)
    
    
    # Run the MICE algorithm for 30 iterations
kds_model.mice(30)

    # Return the completed dataset.
df_imputed_model = kds_model.complete_data()
Xmodel_imputed=np.array(np.round(df_imputed_model,3))    

dfmodeling.loc[:,'30d_Mortality':'External_AIS']=Xmodel_imputed
imputed_sets.append(dfmodeling)

######################################################### Ending imputation and cleaning


dfmodeling[['Biological_sex','ICU_days','Hospital_days','Ventilator_days','30d_Mortality',
            'Lowest_ED_SBP','ED_intubation','ED_CPR','Ph_ED_Hypotension','ED_mortality',
            'Head_AIS','Face_AIS','Abdomen_AIS','Chest_AIS','Extremity_AIS','External_AIS','TBI',
            'phgcs','phsbp','phhr','phrr','Penetrating',
            'Fluids','Total_blood',
            'Fluids_0h','Fluids_4h','Fluids_24h','Total_blood_0h','Total_blood_4h','Total_blood_24h',
            'Plasma_B_2','Plasma_B_3','Whole_B_2','Whole_B_3','PRBC_B_2','PRBC_B_3',
            'Platelet_B_2','Platelet_B_3','Cryo_B_2','Cryo_B_3','Fluids_B_2','Fluids_B_3',
            'Plasma_0h','PRBC_0h','LTOWB_0h','Cryo_0h','Platelet_0h',
            'Plasma_4h','PRBC_4h','LTOWB_4h','Cryo_4h','Platelet_4h',
            'Plasma_24h','PRBC_24h','LTOWB_24h','Cryo_24h','Platelet_24h',
            'Platelets_0h','Platelets_4h','Platelets_24h']]=\
np.round(dfmodeling[['Biological_sex','ICU_days','Hospital_days','Ventilator_days','30d_Mortality',
            'Lowest_ED_SBP','ED_intubation','ED_CPR','Ph_ED_Hypotension','ED_mortality',
            'Head_AIS','Face_AIS','Abdomen_AIS','Chest_AIS','Extremity_AIS','External_AIS','TBI',
            'phgcs','phsbp','phhr','phrr','Penetrating',
            'Fluids','Total_blood',
            'Fluids_0h','Fluids_4h','Fluids_24h','Total_blood_0h','Total_blood_4h','Total_blood_24h',
            'Plasma_B_2','Plasma_B_3','Whole_B_2','Whole_B_3','PRBC_B_2','PRBC_B_3',
            'Platelet_B_2','Platelet_B_3','Cryo_B_2','Cryo_B_3','Fluids_B_2','Fluids_B_3',
            'Plasma_0h','PRBC_0h','LTOWB_0h','Cryo_0h','Platelet_0h',
            'Plasma_4h','PRBC_4h','LTOWB_4h','Cryo_4h','Platelet_4h',
            'Plasma_24h','PRBC_24h','LTOWB_24h','Cryo_24h','Platelet_24h',
            'Platelets_0h','Platelets_4h','Platelets_24h']],0)

dfmodeling['Shock']=np.round(dfmodeling['phhr']/dfmodeling['phsbp'],2)

def sum_of_squares_of_largest(row):
    largest_three = sorted(row, reverse=True)[:3]
    return sum([i**2 for i in largest_three])

# Apply ISS-making function to the subset of the DataFrame
dfmodeling['ISS'] = dfmodeling[['Head_AIS','Face_AIS','Abdomen_AIS',
                                'Chest_AIS','Extremity_AIS','External_AIS']].apply(sum_of_squares_of_largest, axis=1)

######################################## defining new variables

    
dfmodeling['Severity'] = dfmodeling['ISS'].apply(recode_value_Severity)    
dfmodeling['Head_AIS_bi'] = dfmodeling['Head_AIS'].apply(recode_value_AIS)
dfmodeling['Shock_bi']=dfmodeling['Shock'].apply(recode_value_Shock)
dfmodeling['Age_group'] = dfmodeling['Age'].apply(recode_value_Age)
dfmodeling['Prehospital_GCS']=dfmodeling['phgcs'].apply(recode_value_GCS_bi)
dfmodeling['Outcome'] = dfmodeling.apply(lambda row: 1 if (row['ICU_days'] > 10) or (row['30d_Mortality'] == 1) else 0, axis=1)

dfmodeling['WB_vs_Plasma'] = dfmodeling.apply(lambda row: 1 if ((row['LTOWB_24h']==1) and (row['Plasma_24h']==0))  \
                                              else 0 if ((row['LTOWB_24h']==0) and (row['Plasma_24h']==1)) else 2, axis=1)

###### defning LTOWB vs plasma for admission time
dfmodeling['WB_vs_Plasma_q'] = dfmodeling.apply(lambda row: 1 if ((row['LTOWB_0h']>0) and (row['Plasma_0h']==0))  \
                                              else 0 if ((row['LTOWB_0h']==0) and (row['Plasma_0h']>0)) else 2, axis=1)


#########################################
##################################################################################### modeling

#################################### predicting plasma (or alternatively LTOWB or resolution prediction)

ID='ID'

bis=['Penetrating','Biological_sex','Cryo_24h','Platelet_24h','LTOWB_24h'] ####or alternatively 'Plasma_24h'
quants=['ISS','Head_AIS','Age','phgcs','Shock','phrr','Fluids','Total_blood']
sets=quants+bis

####### for resolving status prediction, use sets=['Penetrating','Biological_sex','ISS','Head_AIS',
	###'Age','phgcs','Shock','phrr',#'ED_CPR','Plasma_share','Fluids','Total_blood']


treat='Plasma_24h' ### alternatively,'LTOWB_24h', or 'Outcome'

dfmodel=dfmodeling.loc[:,['ID','Severity','Age_group','Shock_bi','phhr','phsbp','TBI','Head_AIS_bi',
                          'Total_blood_0h','Fluids_0h','Platelet_0h','Cryo_0h','ED_CPR',
                          'Plasma_0h','LTOWB_0h','PRBC_0h','Prehospital_GCS',
                          'Outcome','30d_Mortality',treat]+sets]

dfmodel.loc[:,quants]=(dfmodel.loc[:,quants]-dfmodel.loc[:,quants].mean())/dfmodel.loc[:,quants].std()

################


dftrain=dfmodel.loc[dfmodel['ID'].isin(df0['ID'])==False,:]
choice=LogisticRegression(solver='lbfgs',penalty='l2',random_state=42)


xmodel=dftrain[sets]
ymodel=dftrain[treat]


model=choice
model.fit(xmodel,ymodel)
ypredicttrain=model.predict(xmodel)
yprobatrain=np.array([p[1] for p in model.predict_proba(xmodel)])



rocutrain=roc_auc_score(ymodel,yprobatrain)

print(rocutrain)

####################################################################### ROC curve train
##############################

fpr,tpr,thresh=roc_curve(ymodel,yprobatrain)

point=np.where((tpr-fpr)==np.max(tpr-fpr))[0][0]
aucf=round(auc(fpr,tpr),2)
plt.figure(figsize=(6,6))
ax=plt.axes()
ax.set_facecolor('whitesmoke')
ax.edgecolor="Black"
ax.linewidth=2
plt.grid(color='grey',linestyle='-',linewidth=0.3)
plt.scatter(x=fpr[point],y=tpr[point],s=100,color='Black',zorder=3)
plt.plot(fpr,tpr,label=f'Logistic regression AuROC: {aucf}',color='Crimson',linewidth=3,zorder=2)
plt.plot([0,1],[0,1],color='navy',linestyle='--',linewidth=2,zorder=1,label="Null model")
plt.title('Roc curve on the training data',fontsize=20)
plt.xlabel('\n1 - Specificity',fontsize=16)
plt.ylabel('Sensativity\n',fontsize=16)
plt.legend(loc='lower right')
plt.show()
print('Best cut-off point:',thresh[point],'\nSensativity:',
      round(100*tpr[point],1),',Specificity:',round(100*(1-fpr[point]),1),
     '\nGeometric mean of diagnostic Odds:',round(np.sqrt((tpr[point]*(1-fpr[point]))/((1-tpr[point])*fpr[point])),1))
print()

############################################################################## calib curve train
prob_true_train, prob_pred_train= calibration_curve(dftrain[treat], yprobatrain, n_bins=8)

plt.figure(figsize=(6,6))
ax=plt.axes()
ax.set_facecolor('whitesmoke')
ax.edgecolor="Black"
ax.linewidth=2
plt.grid(color='grey',linestyle='-',linewidth=0.3)

plt.plot(prob_pred_train, prob_true_train,color='Green',linewidth=3,label="Logistic regression model")
plt.plot([.0,1.0], [.0,1.0],color='navy',linestyle='--',linewidth=2,zorder=1,label="Perfect calibration")
plt.title('Calibration curve on the training set',fontsize=20)
plt.ylabel("Fraction of Positives",fontsize=16)
plt.xlabel("Average Prediction",fontsize=16)
plt.legend()

plt.show()


print('\n')

############################### testing on omic patients

dftest=dfmodel.loc[(dfmodel['ID'].isin(df0['ID'])==True),:]

x=dftest[sets]
y=dftest[treat]

ypredict=model.predict(x)
yproba=np.array([p[1] for p in model.predict_proba(x)])
    
############################################################## weight management


overlaplist=[p[0] for p in model.predict_proba(x)]
overlapw=np.array(overlaplist)

overlapwsum=np.array(overlaplist).sum()
squarer = lambda t: t ** 2
overlapwsquares = np.array([squarer(xi) for xi in overlaplist])
overlapwsumsquares=overlapwsquares.sum()

effective_n=(overlapwsum**2)/overlapwsumsquares
weight_frame=pd.DataFrame({'ID':dftest['ID'].values,'Overlap_weights':overlapw})

rocu=roc_auc_score(y,yproba)
print(rocu)
####################################################################### ROC curve test
##############################

fpr,tpr,thresh=roc_curve(y,yproba)

point=np.where((tpr-fpr)==np.max(tpr-fpr))[0][0]
aucf=round(auc(fpr,tpr),2)
plt.figure(figsize=(6,6))
ax=plt.axes()
ax.set_facecolor('whitesmoke')
ax.edgecolor="Black"
ax.linewidth=2
plt.grid(color='grey',linestyle='-',linewidth=0.3)
plt.scatter(x=fpr[point],y=tpr[point],s=100,color='Black',zorder=3)
plt.plot(fpr,tpr,label=f'Logistic regression AuROC: {aucf}',color='Crimson',linewidth=3,zorder=2)
plt.plot([0,1],[0,1],color='navy',linestyle='--',linewidth=2,zorder=1,label="Null model")
plt.title('Roc curve on the test set',fontsize=20)
plt.xlabel('\n1 - Specificity',fontsize=16)
plt.ylabel('Sensativity\n',fontsize=16)
plt.legend(loc='lower right')
plt.show()
print('Best cut-off point:',thresh[point],'\nSensativity:',
      round(100*tpr[point],1),',Specificity:',round(100*(1-fpr[point]),1),
     '\nGeometric mean of diagnostic Odds:',round(np.sqrt((tpr[point]*(1-fpr[point]))/((1-tpr[point])*fpr[point])),1))
print()

######################################################################## Cal curve test

prob_true_test, prob_pred_test= calibration_curve(dftest[treat], yproba, n_bins=8,strategy='uniform')

plt.figure(figsize=(6,6))
ax=plt.axes()
ax.set_facecolor('whitesmoke')
ax.edgecolor="Black"
ax.linewidth=2
plt.grid(color='grey',linestyle='-',linewidth=0.3)

plt.plot(prob_pred_test, prob_true_test,color='Green',linewidth=3,label="Propensity_based Model")
plt.plot([.0,1.0], [.0,1.0],color='navy',linestyle='--',linewidth=2,zorder=1, label="Perfect calibration ")
plt.title('Calibration curve on the test set',fontsize=20)
plt.ylabel("Fraction of Positives",fontsize=16)
plt.xlabel("Average Prediction",fontsize=16)
plt.legend()

plt.show()

    
for i in range(len(sets)):
    OR= np.exp(model.coef_)[0][i]
    print(sets[i],OR)    

######################## creating a merged dataset and finding effective sample size

print('effective n:',effective_n)
merged_df0=df0.merge(weight_frame, on='ID',how='left'))
