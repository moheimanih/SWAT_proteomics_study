################ Multivariable regression and SEM for blood products and clinical outcomes ##################
############################################################################################################

################################################ creating new variables in SWAT cohort data

dfmodeling['Plasma_share_0h'] =dfmodeling['Plasma_0h']/(dfmodeling['Total_blood_0h']+0.1)
dfmodeling['LTOWB_share_0h'] =dfmodeling['LTOWB_0h']/(dfmodeling['Total_blood_0h']+0.1)
dfmodeling['Platelet_share_0h'] =dfmodeling['Platelet_0h']/(dfmodeling['Total_blood_0h']+0.1)
dfmodeling['Cryo_share_0h'] =dfmodeling['Cryo_0h']/(dfmodeling['Total_blood_0h']+0.1)
dfmodeling['PRBC_share_0h'] =dfmodeling['PRBC_0h']/(dfmodeling['Total_blood_0h']+0.1)

dfmodeling['Plasma_recieved_by_admission']=sp.stats.mstats.winsorize(dfmodeling['Plasma_0h'].values,(0,0.01))
dfmodeling['PRBC_recieved_by_admission']=sp.stats.mstats.winsorize(dfmodeling['PRBC_0h'].values,(0,0.01))
dfmodeling['Platelet_recieved_by_admission']=sp.stats.mstats.winsorize(dfmodeling['Platelet_0h'].values,(0,0.01))
dfmodeling['Cryo_recieved_by_admission']=sp.stats.mstats.winsorize(dfmodeling['Cryo_0h'].values,(0,0.01))
dfmodeling['LTOWB_recieved_by_admission']=sp.stats.mstats.winsorize(dfmodeling['LTOWB_0h'].values,(0,0.01))
dfmodeling['Fluids_recieved_by_admission']=sp.stats.mstats.winsorize(dfmodeling['Fluids_0h'].values,(0,0.01))

dfmodeling['Other_components_by_0h']=dfmodeling['PRBC_recieved_by_admission']+ \
dfmodeling['Platelet_recieved_by_admission']+dfmodeling['Cryo_recieved_by_admission']

dfmodeling['Blood_0_24h']=dfmodeling['Total_blood_24h']-dfmodeling['Total_blood_0h']
dfmodeling['Totals_0_24h']=sp.stats.mstats.winsorize(dfmodeling['Blood_0_24h'].values,(0,0.01))

######################################################### Multivariate regression for first-day transfusion or outcomes

dfshare=dfmodeling.loc[((dfmodeling['Plasma_4h']>=1)&(dfmodeling['LTOWB_4h']>=1)),:].reset_index()
####### use only if filtering patients is desired

constants=['Age','ISS','phgcs','phrr','phsbp','phhr','Shock','Fluids_0h','Total_blood_0h','Total_blood_24h',
           'Plasma_0h','Platelet_0h','Cryo_0h','PRBC_0h','LTOWB_0h','Other_components_by_0h',
           'Plasma_recieved_by_admission','PRBC_recieved_by_admission','Platelet_recieved_by_admission',
           'Cryo_recieved_by_admission','LTOWB_recieved_by_admission','Fluids_recieved_by_admission',
           'Plasma_share_0h','Plasma_share_received_by_0h','Totals_0_24h']

dfshare[constants]=(dfshare[constants]-dfshare[constants].mean())/dfshare[constants].std()
constants=constants+['Center']

X2='Plasma_recieved_by_admission'
Z2='LTOWB_recieved_by_admission'
Y2='Totals_0_24h'

dfshare['Penetrating_X_Plasma']=dfshare['Penetrating']*dfshare[X2]
dfshare['Shock_X_Plasma']=dfshare['phhr']*dfshare[X2]
dfshare['TBI_X_Plasma']=dfshare['TBI']*dfshare[X2]
dfshare['Penetrating_X_TBI_X_Plasma']=dfshare['Penetrating']*dfshare['TBI']*dfshare[X2]


###### if 'outcome is selected as y variable, add the following parameter: ,family=sm.families.Binomial()
share_model=sm.GLM.from_formula("Totals_0_24h ~ \
Plasma_recieved_by_admission+LTOWB_recieved_by_admission+Other_components_by_0h+ \
Age+TBI+ISS+phgcs+phrr+phsbp+phhr+Fluids_0h+Penetrating+Biological_sex",dfshare) 

print(share_model.fit(cov_type='HC3').summary())

share=dfshare[X2]
rest=dfshare[Y2]
########################################## Charting

sns.set(rc={'figure.figsize':(90,60)})
with sns.axes_style({'axes.facecolor':'white','axes.edgecolor':'black',
                     'axes.grid':True,'grid.color':'Grey','grid.linestyle':'--'}):

    plot=sns.lmplot(y=Y2,x=X2,hue='Penetrating',data=dfshare,order=1,col='Shock_bi', ### remove col or hue if desired
                height=8,palette="Set1",markers=['o','x'])
    plot.set_axis_labels(f'\n{X2}',f'{Y2}\n',fontsize=20)
    plt.show(plot)

print('\n')

################## Partial correlation between injury groups and first-day transfusions

dfsharent2=dfshare.loc[(dfshare['TBI']==0) & (dfshare['Penetrating']==1),:].reset_index(drop=True)
print()
print('\nPenet and NBI')
print(pg.partial_corr(data=dfsharent2, x='Plasma_0h', y=Y2,
                      covar=['Age','ISS','phgcs','phhr','LTOWB_0h','Other_products_0h',
                             'Biological_sex','phsbp','phrr'
                             ,'Fluids_0h'],method='spearman').round(3))

dfsharen2=dfshare.loc[(dfshare['TBI']==0) & (dfshare['Penetrating']==0),:].reset_index(drop=True)
print()
print('No Penet and NBI')
print(pg.partial_corr(data=dfsharen2, x='Plasma_0h', y=Y2,
                      covar=['Age','ISS','phgcs','phhr','LTOWB_0h','Other_products_0h',
                             'Biological_sex','phsbp','phrr'
                             ,'Fluids_0h'],method='spearman').round(3))

dfsharen2=dfshare.loc[(dfshare['TBI']==1) & (dfshare['Penetrating']==0),:].reset_index(drop=True)
print()
print('No Penet and TBI')
print(pg.partial_corr(data=dfsharen2, x='Plasma_0h', y=Y2,
                      covar=['Age','ISS','phgcs','phhr','LTOWB_0h','Other_products_0h',
                             'Biological_sex','phsbp','phrr'
                             ,'Fluids_0h'],method='spearman').round(3))

##############################
########################################################## SEM model for first-day transfusions

X='Plasma_share_0h'
Z='LTOWB_share_0h'

dfmodeling['Total_units_0h_24h']=sp.stats.mstats.winsorize(dfmodeling['Blood_0_24h'].values,(0.01,0.01))
Y='Total_units_0h_24h'

dfshare['Penetrating_X_Plasma']=dfshare['Penetrating']*dfshare[X]
dfshare['Shock_X_Plasma']=dfshare['phhr']*dfshare[X]
dfshare['TBI_X_Plasma']=dfshare['TBI']*dfshare[X]
dfshare['Penetrating_X_TBI_X_Plasma']=dfshare['Penetrating']*dfshare['TBI']*dfshare[X]


mod = """ Total_units_0h_24h ~ \
TBI_X_Plasma+Plasma_share_0h+ \
Total_blood_0h+LTOWB_share_0h+Cryo_share_0h+ \
Age+TBI+ISS+phgcs+phrr+phsbp+phhr+Fluids_0h+Penetrating+Biological_sex
      """
semmodel = semopy.ModelEffects(mod)
semmodel.fit(dfshare,group='Center')
print(semmodel.inspect())
