############# codes for adjusted and unadjusted heatmaps and PCA in th omic cohort ############
###############################################################################################


################################################# Unadjusted heatmaps for select protein groups

treat='Plasma'## alternatively LTOWB
dfna=df0.loc[ (df0['Sample']=='Donor') | (df0['Sample']=='SWAT'),:]

dfna.loc[ (dfna['Sample']=='Donor'),treat]=\
np.array(['Donor' for i in range(dfna.loc[ dfna['Sample']=='Donor'].shape[0])])

variables=['Time point',treat]

####################### selected protein groups
f0=coagulation
f1=complement
f2=plateletn
f3=endothel
f4=neutrophil
f5=inflammation
f6=mediator_s
f7=[]

f=list(set(f0+f1+f2+f3+f4+f5+f6))


colors=[]
for protein in f:
    if protein in f0:
        colors.append('Crimson')
    elif protein in f1:
        colors.append('yellowgreen')
    elif protein in f2:
        colors.append('khaki')
    elif protein in f3:
        colors.append('aqua')
    elif protein in f4:
        colors.append('pink')
    elif protein in f5:
        colors.append('grey')
    elif protein in f6:
        colors.append('Green')
    else:
        colors.append('grey')


heatframe=pd.DataFrame({})
dicte={}
if len (variables)==4:
    for val1 in dfna[variables[0]].unique():
        for val2 in dfna[variables[1]].unique():
            for val3 in dfna[variables[2]].unique():
                for val4 in dfna[variables[3]].unique():
                    dicte[f'{val1} {val2} {val3} {val4}']=[]
elif len (variables)==3:
    for val1 in dfna[variables[0]].unique():
        for val2 in dfna[variables[1]].unique():
            for val3 in dfna[variables[2]].unique():
                dicte[f'{val1} {val2} {val3}']=[]
elif len (variables)==2:
    for val1 in dfna[variables[0]].unique():
        for val2 in dfna[variables[1]].unique():
            dicte[f'{val1} {val2}']=[]
elif len (variables)==1:
    for val1 in dfna[variables[0]].unique():
        dicte[f'{val1}']=[]

#calculating z_scores for each groups for all proteins

dfna.loc[:,f]=np.log2(dfna.loc[:,f])


for column in f:
    if dfna[column].isna().sum()>0:
        f.remove(column)
        print(f'{column} was removed')
        
for column in f:
    dfg=dfna.groupby(variables)[str(column)]
    dfz=dfg.mean().apply(lambda r:(r-dfna[str(column)].mean())/dfna[str(column)].std())
    if len (variables)==4:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i][0]}'+' '+f'{dfz.index[i][1]}'+' '+f'{dfz.index[i][2]}'+
                  ' '+f'{dfz.index[i][3]}'].append(dfz.iloc[i])
    if len (variables)==3:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i][0]}'+' '+f'{dfz.index[i][1]}'+' '+f'{dfz.index[i][2]}'].append(dfz.iloc[i])
    elif len (variables)==2:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i][0]}'+' '+f'{dfz.index[i][1]}'].append(dfz.iloc[i])
    elif len(variables)==1:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i]}'].append(dfz.iloc[i])

#Making a new dataframe for plotting

for key in dicte:
    if dicte[key]==[]:
        pass
    else:
        heatframe[key]=dicte[key]



#ploting the clustered heatmap on the dataframe that includes the z-socres for all proteins
plot=sns.clustermap(heatframe,center=0.0,linewidth=0,vmin=-1.3,vmax=1.3,yticklabels=pnames ## list of protein abbreviations
                    ,figsize=(2.1*np.log2(len(f))**1.1,1.8*np.log2(len(f))**1.35)
                    ,metric='euclidean',row_colors=colors
                    ,dendrogram_ratio=0.05,col_cluster=False,cbar_pos=None
                    ,cmap=sns.diverging_palette(258,-372, n=1001))
plt.yticks(fontsize=26,rotation=0)
plt.xticks(fontsize=10,rotation=90)
plt.title(f'Unadjusted levels in {treat}-treated vs. untreated/donors\n',fontsize=30)
plt.show()

################################################# Confounder-adjusted heatmaps for select protein groups

variables=['Time point',treat]
olw='Overlap_weights' ############### from primary overlap modeling

f0=coagulation
f1=complement
f2=platelet
f3=endothelial
f4=neutrophil
f5=inflammation
f6=mediator_s
f7=[]

colors=[]
for protein in f:
    if protein in f0:
        colors.append('Crimson')
    elif protein in f1:
        colors.append('yellowgreen')
    elif protein in f2:
        colors.append('khaki')
    elif protein in f3:
        colors.append('aqua')
    elif protein in f4:
        colors.append('pink')
    elif protein in f5:
        colors.append('grey')
    elif protein in f6:
        colors.append('Green')
    else:
        colors.append('blue')

df1=merged_df0.loc[merged_df0['Sample']=='SWAT',:]
df1.loc[:,f]=np.log2(df1.loc[:,f])

for column in f:
    if merged_df0[column].isna().sum()>0:
        f.remove(column)
        print(f'{column} was removed')


heatframe=pd.DataFrame({})
dicte={}
if len (variables)==4:
    for val1 in df1[variables[0]].unique():
        for val2 in df1[variables[1]].unique():
            for val3 in df1[variables[2]].unique():
                for val4 in df1[variables[3]].unique():
                    dicte[f'{val1} {val2} {val3} {val4}']=[]
elif len (variables)==3:
    for val1 in df1[variables[0]].unique():
        for val2 in df1[variables[1]].unique():
            for val3 in df1[variables[2]].unique():
                dicte[f'{val1} {val2} {val3}']=[]
elif len (variables)==2:
    for val1 in df1[variables[0]].unique():
        for val2 in df1[variables[1]].unique():
            dicte[f'{val1} {val2}']=[]
elif len (variables)==1:
    for val1 in df1[variables[0]].unique():
        dicte[f'{val1}']=[]


for column in f:
    
    weighted_average = np.average(df1[str(column)], weights=df1[olw])
    weighted_stds=weighted_std(df1[str(column)],df1[olw])

    dfg=df1.groupby(variables)
    dfm=dfg.apply(lambda x: np.average(x[str(column)], weights=x[olw]))
    dfz=dfm.apply(lambda r:(r-weighted_average)/weighted_stds)
    
    
    if len (variables)==4:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i][0]}'+' '+f'{dfz.index[i][1]}'+' '+f'{dfz.index[i][2]}'+
                  ' '+f'{dfz.index[i][3]}'].append(dfz.iloc[i])
    if len (variables)==3:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i][0]}'+' '+f'{dfz.index[i][1]}'+' '+f'{dfz.index[i][2]}'].append(dfz.iloc[i])
    elif len (variables)==2:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i][0]}'+' '+f'{dfz.index[i][1]}'].append(dfz.iloc[i])
    elif len(variables)==1:
        for i in range(0,len(dfz.index)):
            dicte[f'{dfz.index[i]}'].append(dfz.iloc[i])

#Making a new dataframe for plotting

for key in dicte:
    if dicte[key]==[]:
        pass
    else:
        heatframe[key]=dicte[key]

#ploting the clustered heatmap on the dataframe that includes the z-socres for all proteins
plot=sns.clustermap(heatframe,center=0.0,linewidth=0,vmin=-1.3,vmax=1.3,yticklabels=pnames ## protein names
                    ,figsize=(2.1*np.log2(len(f))**1.1,1.8*np.log2(len(f))**1.35)
                    ,metric='euclidean',row_colors=colors
                    ,dendrogram_ratio=0.05,col_cluster=False,cbar_pos=None
                    ,cmap=sns.diverging_palette(258,-372, n=1001))
plt.yticks(fontsize=26,rotation=0)
plt.xticks(fontsize=10,rotation=90)
plt.title(f'Adjusted levels in {treat}-treated and untreated\n',fontsize=30)
plt.show()



####################################
######################################################333###33############ PCA

pcas=[0,1]
variable='Plasma'## or any other covariate

f=list(set(coagulation+platelet+complement+inflammation+endothelial+neutrophil+mediator_s))

#Designing the dataframes for downstream analysis

df1=df0.loc[(df0['ID']!=''),:] ### fill '' with potential extreme outlier IDs
df1.loc[:,f]=np.log2(df1.loc[:,f])
df1.loc[:,f]=(df1.loc[:,f]-df1.loc[:,f].mean())/df1.loc[:,f].std()
df2=df1.loc[df1['Sample']!='Donor',:].reset_index()
df3=df2.loc[:,f]
df3=df3.dropna(axis=1)
df3=df3.dropna(axis=0)
df4=df2.loc[:,['Sample','Time point','LTOWB','Plasma','TBI','Penetrating']+new_cats]
          
####PCA

pca=PCA()
pca.fit(df3.iloc[:,:-1])
pca_data=pca.transform(df3.iloc[:,:-1])
per_var=np.round(pca.explained_variance_ratio_*100,decimals=1)
var=per_var[:10]
labels=['PC' + str(x) for x in range (1, len(df3.columns))]

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

#ploting the patients based on the two first PCs:
pca_df=pd.DataFrame(pca_data)

chartdf=pca_df.iloc[:,0:10]
chartdf.columns=labels[:10]

for column in chartdf.columns:
    df4[column]=chartdf[column]

if len(pcas)==2:
    sns.set(rc={'figure.figsize':(8,8)})
    ax=plt.axes()
    ax.set_facecolor('whitesmoke')
    colorlist=[colormap[c] for c in df4[variable].values]
    labell=[l for l in df4[variable].values]
    plt.scatter(df4.loc[:,'PC1'],df4.loc[:,'PC2'],s=100,color=colorlist)
    plt.title(f'PCA using the selected markers\n(Dots present individual patient samples)', size=23)
    plt.xlabel('\nPC1',fontsize=22)
    plt.ylabel('PC2\n',fontsize=22)
    plt.legend()
    plt.grid(color='grey',linestyle='-',linewidth=0.3)
    plt.show()

