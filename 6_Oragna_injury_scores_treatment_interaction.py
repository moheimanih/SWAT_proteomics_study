##################### Correlation heatmap of modules, organ-specific PCA and PT correlation ###########################
######################################################################################################################

############################# designing a new frame for downstream analysis from community detected frames

importantes_modulo=['Module1','Module3','Module7','Module8','Module10']
## tier one and two modules as identified through leiden community detection

dfmodular=df6.loc[df6['Time point']==0,['ID','platelet factor 4','fibrinogen','Biological_sex','Age','Shock','phrr','phhr',
                                        'phsbp','phgcs','ISS','Head_AIS','TBI','Chest_AIS',
                                        'Penetrating','Plasma','LTOWB','Platelet','Cryo','PRBC',
                                        'Plasma_units','LTOWB_units','Prob_resolving','Total_24h_units',
                                        'Platelet_units','Cryo_units','PRBC_units','Fluids_0h',
                                        'Total_blood_0h','Total_blood_24h',
                                        'INR','PT','Outcome']+importantes_modulo].reset_index(drop=True)

modular4=df6.loc[df6['Time point']==4,['ID','Plasma_units','LTOWB_units',
                                       'Platelet_units','Cryo_units','Total_24h_units',
                                       'PRBC_units','Fluids_4h','Total_blood_4h']].reset_index(drop=True)
                                       
modular24=df6.loc[df6['Time point']==24,['ID','Plasma_units','LTOWB_units',
                                       'Platelet_units','Cryo_units','Total_24h_units',
                                       'PRBC_units','Fluids_0h','Total_blood_24h']].reset_index(drop=True)
                  
newnames=['ID','PF4','FIB','Biological_sex','Age','Shock','phrr','phhr',
          'phsbp','phgcs','ISS','Head_AIS','TBI','Chest_AIS',
                                        'Penetrating','Plasma','LTOWB','Platelet','Cryo','PRBC',
                                        'Plasma_units','LTOWB_units','Prob_resolving','Total_24h_units',
                                        'Platelet_units','Cryo_units','PRBC_units','Fluids_0h',
                                        'Total_blood_0h','Total_blood_24h','INR','PT',
                                        'Non_resolution','FLAME1','FLAME2','CLOT','ALPHA','SURF']
dfmodular.columns=newnames

############################################################### PCA for brain and heart scores
ID='ID'
time=0
sample_set='SWAT'

#number of pcas to be shown
pcas=[0,1,2]

#recoding variable
def recoder(df,name,newname,val1,rec1,val2,rec2):
    newcolumn=[]
    for i in enumerate(df.index):
        row=i[1]
        value=df.loc[row,name]
        if value==val1:
            newcolumn.append(rec1)
        else:
            newcolumn.append(rec2)
    df[newname]=newcolumn    

df1=df0.loc[df0['Sample']==sample_set,:]

for groups in ['Brain','Heart']:
    
    if groups=='Brain':
        f=brain_pros

        df1.loc[:,f]=np.log2(df1.loc[:,f])
        df2=df1.loc[df1['Time point']==time,:].reset_index()
        recoder(df2,variable,'OutcomeN',df1[variable].unique()[0],0,df1[variable].unique()[1],1)
        df3=df2.loc[:,f]
        df3=(df3-df3.mean())/df3.std()
        df3=df3.dropna(axis=1)
        df3=df3.dropna(axis=0)
        df4=df2.iloc[:,:first]
    
    elif groups=='Heart':
        f=cardiac_pros

        df1.loc[:,f]=np.log2(df1.loc[:,f])
        df2=df1.loc[df1['Time point']==time,:].reset_index()
        recoder(df2,variable,'OutcomeN',df1[variable].unique()[0],0,df1[variable].unique()[1],1)

        df3=df2.loc[:,f]
        df3=(df3-df3.mean())/df3.std()
        df3=df3.dropna(axis=1)
        df3=df3.dropna(axis=0)
        df4=df2.iloc[:,:first]
        
    #################################### PCA

    pca=PCA()
    pca.fit(df3.iloc[:,:-1])
    pca_data=pca.transform(df3.iloc[:,:-1])
    per_var=np.round(pca.explained_variance_ratio_*100,decimals=1)
    var=per_var[:10]
    print(var)
    labels=[f'{groups}_PC' + str(x) for x in range (1, len(df3.columns))]

    #visualizing the explained variances:
    sns.set(rc={'figure.figsize':(6,6)})
    ax=plt.axes()
    ax.set_facecolor('whitesmoke')
    ax.edgecolor="Black"
    plt.grid(color='grey',linestyle='-',linewidth=0.3)
    colormap=['maroon' if i%2==0 else 'sandybrown' for i in range(len(var))]
    plt.barh(labels[:10],width=var,height=0.8,color=colormap)
    plt.gca().invert_yaxis()
    plt.xlabel('Percentage of explained variance')
    plt.ylabel('Principal component')
    plt.title('Scree plot',fontsize=20)
    plt.show()

    #Printing highest loadings:
    for comp in pcas:
        loading=pd.Series(pca.components_[comp],index=df3.columns[:-1])
        sortedd=loading.abs().sort_values(ascending=False)
        topsortedd=sortedd[0:10].index.values
        print(time,f'hour {groups}_PC{comp+1}:\n',loading[topsortedd],'\n')

    pca_df=pd.DataFrame(pca_data)

    chartdf=pca_df.iloc[:,0:10]
    chartdf.columns=labels[:10]

    for column in chartdf.columns:
        df4[column]=chartdf[column]
        dfmodular[column]=chartdf[column]


############################### Assessing the first three PCs for organ specificity

print(dfmodular.loc[df6['Time point']==0,['Chest_AIS','Head_AIS','ISS','phgcs','TBI','Heart_PC1',
                                          'Brain_PC1','Brain_PC2','Brain_PC3'
                                   ]].corr(method='spearman'))

print(dfmodular.loc[df6['Time point']==0,['Chest_AIS','Head_AIS','ISS','phgcs','TBI','Brain_PC2',
                                          'Heart_PC1','Heart_PC2','Heart_PC3',
                                   ]].corr(method='spearman'))

######################################## Building data frames for correlation analysis

dfmodular['WB_vs_Plasma'] = dfmodular.apply(lambda row: 1 if ((row['LTOWB']==1) and (row['Plasma']==0))  \
                                              else 0 if ((row['LTOWB']==0) and (row['Plasma']==1)) else 2, axis=1)

dfmodular['Plasma_share'] = dfmodular['Plasma_units']/(dfmodular['Total_blood_0h']+0.1)
dfmodular['LTOWB_share'] = dfmodular['LTOWB_units']/(dfmodular['Total_blood_0h']+0.1)

dfmodularmapdf=dfmodular.loc[:,['LTOWB_share','FLAME2','FLAME1','SURF','ALPHA','CLOT','Plasma_share']]

print(dfmodularmapdf.corr(method='spearman'))
dfmodularmapdf.corr(method='spearman').style.background_gradient(cmap ='coolwarm')
splot=sns.clustermap(dfmodularmapdf.corr(method='spearman'),
                    center=0.0,linewidth=0,figsize=(5,5),dendrogram_ratio=0.05,cbar_pos=None,
                    cmap=sns.diverging_palette(265,-366, n=1001),row_cluster=False,col_cluster=False)

plt.title(f'Correlation between modules and treatments\n',fontsize=14)

plt.show()


for x in ['FLAME1','FLAME2','LTOWB_share','CLOT','ALPHA','SURF']:
    print(sp.stats.spearmanr(dfmodular['Plasma_share'],dfmodular[x]))
    print()


######################################################### partial correlation of admission PT and plasma (alternatively LTOWB)

B_median=np.quantile(dfmodular['Brain_PC2'],0.5))
H_median=np.quantile(dfmodular['Heart_PC1'],0.5))


print('TBI')
dfmodulart=dfmodular.loc[
                         (dfmodular['Brain_PC2']>B_median),:].reset_index(drop=True)
print(pg.partial_corr(data=dfmodulart, x='Plasma_share', y='INR',#'Total_units_0h_24h',PT'
                      covar=['Total_blood_0h','Age','phhr','phgcs',#'Cryo_share','Platelet_share',
                             'ISS','Fluids_0h','Penetrating','phsbp','phrr','TBI','Heart_PC1',
                                                     'Biological_sex'],method='spearman').round(4))

print('Heart')
dfmodulart=dfmodular.loc[
                         (dfmodular['Heart_PC1']>H_median),:].reset_index(drop=True)
print(pg.partial_corr(data=dfmodulart, x='Plasma_share', y='INR',#'Total_units_0h_24h',PT'
                      covar=['Total_blood_0h','Age','phhr','phgcs',#'Cryo_share','Platelet_share',
                             'ISS','Fluids_0h','Penetrating','phsbp','phrr','TBI','Brain_PC2',
                                                     'Biological_sex'],method='spearman').round(4))

#################################### graphing PT and blood product correlation in high and low brain injury groups

interventions = ['LTOWB share', 'Plasma share']
classes = ['Low', 'High']

###################### inputing means and intervals caclulated from above
means = {
    'LTOWB share': [-0.11, -0.0237],
    'Plasma share': [-0.115, -0.2623]
}

ci_low = {
    'LTOWB share': [-0.36, -0.28],
    'Plasma share': [-0.365,-0.49]
}

ci_high = {
    'LTOWB share': [0.15,  0.24],
    'Plasma share': [0.145, -0.0]
}

errors_lower = {
    'LTOWB share': [means['LTOWB share'][i] - ci_low['LTOWB share'][i] for i in range(len(classes))],
    'Plasma share': [means['Plasma share'][i] - ci_low['Plasma share'][i] for i in range(len(classes))]
}

errors_upper = {
    'LTOWB share': [ci_high['LTOWB share'][i] - means['LTOWB share'][i] for i in range(len(classes))],
    'Plasma share': [ci_high['Plasma share'][i] - means['Plasma share'][i] for i in range(len(classes))]
}

bar_width = 0.2
index = np.arange(len(classes))

colors_intervention_1 = [ '#ff7f0e']  # Blue and orange for Intervention 1
colors_intervention_2 = [ '#ff7f0e']  # Similar color scheme for Intervention 2
background_color = '0.97'  # Light gray background

fig, ax = plt.subplots()
fig.patch.set_facecolor(background_color)

ax.set_facecolor(background_color)
ax.errorbar(index , means['Plasma share'], yerr=[errors_lower['Plasma share'], errors_upper['Plasma share']],
            fmt='o', color='green',capsize=5, linestyle='None',elinewidth=7, markersize=18,)

ax.errorbar(index+ bar_width, means['LTOWB share'], yerr=[errors_lower['LTOWB share'], errors_upper['LTOWB share']],
            fmt='o', color='Crimson',capsize=5, linestyle='None',elinewidth=7, markersize=18,)

ax.set_xlabel('\nBrain Injury Score', fontsize=27)
ax.set_ylabel('Partial Spearman correlation\n', fontsize=27)
#ax.set_title('Association of early treatment shares with 0h PT in omic subset\n\n', fontsize=30)
ax.set_xticks(index + bar_width / 5)
ax.set_xticklabels(classes, fontsize=25)

ax.grid(color='grey',linestyle='--',linewidth=0.3)
plt.tight_layout()
plt.show()
