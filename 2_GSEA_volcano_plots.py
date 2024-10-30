########################
############################################################ GSEA

############## getting a data frame with gene names in columns and defining the total gene list
address='SWAT_proteins_clinical_genes.csv'
dfG=pd.read_csv(address)
first=list(dfG.columns).index('CRYBB2')
last=list(dfG.columns).index('IRF6')

allofthem=list(dfG.columns[first:last+1])

for column in allofthem:
    if 'fc_mouse' in column:
        allofthem.remove(column)

dfG['WB_vs_Plasma'] = dfG.apply(lambda row: 1 if ((row['LTOWB']==1) and (row['Plasma']==0))  \
                                              else 0 if ((row['LTOWB']==0) and (row['Plasma']==1)) else 2, axis=1)

############

variable='WB_vs_Plasma' ### alternatively sample, time point, etc.

######## filtering for severely injured patients with a low probability of resolution (avoid if full sample desired)

df1=dfG.loc[((dfG['ISS']>24))&((dfG['Time point']==0)&(dfG['Prob_resolving']<0.9)&
                               (dfG['WB_vs_Plasma']!=2)),:]#
  
    
df1.iloc[:,first:last+1]=np.log2(df1.iloc[:,first:last+1])

df1.iloc[:,first:last+1]= \
(df1.iloc[:,first:last+1]-np.mean(df1.iloc[:,first:last+1]))/np.std(df1.iloc[:,first:last+1])

dfgsea=df1.iloc[:,first:last+1].T
phenotypes=df1[variable].values


base=['Reactome_2022','GO_Biological_Process_2023']
enr = gp.gsea(data=dfgsea,
           gene_sets=base,
           cls=phenotypes,
           permutation_type='phenotype', 
           permutation_num=2000, 
           outdir=None, 
           no_plot=True,
           method='signal_to_noise',
           processes=4,
           format='png')

################################################ plotting GSEA

dictform=[]
x=[]
y=[]
countpos=0
countneg=0
numselect=15


threshold=0.05 ## or alternatively 0.1
p_value_choice='FWER p-val' ## alternatively NOM p-val

#marker size correction factor
cf=5

#making arrays of columns:
names=np.array(enr.res2d.Term)
nes=np.array(enr.res2d.NES)
fdr= np.array(enr.res2d[p_value_choice])

#combining them into a dictionary:
for i in range (len(names)):
    dictform.append([names[i],nes[i],fdr[i]])
    if nes[i]<0:
        countneg=countneg+1
        dictform[i].append(0)
    else:
        countpos=countpos+1
        dictform[i].append(1)

#sorting the dictionary based on the NES values
def sorter(lister):
    lister.sort(key=lambda x:x[1],reverse=False)
    return(lister)

sorter(dictform)
paths=dictform
paths=dictform[:min(numselect,countneg)]+dictform[-min(numselect,countpos):]#

#Making the axis:
for path in paths:
    x.append(path[1])
    y.append(path[0])
    

#################### drawing

colormap=['Orangered' if (j[3]==1 and j[2]<threshold) else 'pink' if (j[3]==1 and j[2]>=threshold)
          else 'Blue' if (j[3]==0 and j[2]<threshold) else 'dodgerblue' for j in paths]
sizer=min(numselect,countneg)+min(numselect,countpos)
plt.figure(figsize=(sizer/6.6,sizer/3.5))

ax=plt.axes()
ax.set_facecolor('whitesmoke')

s=[]
for percent in enr.res2d['Gene %'].tolist():
    s.append(np.round(cf*float(percent[:3]),1))
S=np.array(s)
for y_val, width, size, colors in zip(y, x, S, colormap):
        plt.scatter(width, y_val, color=colors, s=size,edgecolor='black')
ax.axvline(0, color='grey', linestyle='--', linewidth=1.5)

plt.xlabel('\nNES',fontsize=18)
plt.xlim([-2.4,2.4])
plt.xticks(fontsize=14)
plt.yticks(rotation=0,fontsize=10)

plt.title('GSEA\n',fontsize=20)
for spine in ax.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(2)
ax.grid(color='grey',linestyle='-',linewidth=0.3)

sizes = [np.quantile(S,0.1), np.quantile(S,0.5),
         np.quantile(S,0.9)]
markers = [plt.Line2D([0,0], [0,0], color='whitesmoke', marker='o',
                      markersize=np.log2(si), markerfacecolor='Grey', markeredgecolor='black') for si in sizes]

color_legend_labels = [f'NES >0, Adj-p <{threshold}','NES >0','NES <0',f'NES <0, Adj-p <{threshold}'] 
colors = ['Orangered','Pink','dodgerblue','Blue']  
patches = [mpatches.Patch(color=color, label=label) for color, label in zip(colors, color_legend_labels)]
legend2 = plt.legend(handles=patches, labels=color_legend_labels, title='',
                     loc='upper left',bbox_to_anchor=(0.48, 0.18),fontsize=10)
ax.add_artist(legend2)
plt.show()


################################################### 
################################################################### volcano plots

variable='Plasma'# or 'LTOWB
ID='ID'
condition=((df1['ISS']>24)&(df1['Prob_resolving']<0.9)) ## condition=(df1['ISS']<75) if all samples are to be included

time=0
f=list(set(allofthem))

df1=dfG.loc[(dfG['Time point']==time),:]
df2=df1.loc[condition,:].reset_index()
df2=df2.dropna(axis=1)
df2=df2.dropna(axis=0)

df2.loc[:,f]=np.log2(df2.loc[:,f])
df2.loc[:,f]=(df2.loc[:,f]-df2.loc[:,f].mean())/df2.loc[:,f].std() #### only choose if standardization is desired

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
    
recoder(df2,variable,'OutcomeN',df1[variable].unique()[0],1,'filler',0)

##################

dictform=[]
pvalues=[]
for column in f:
    columninfo=[]
    columninfo.append(column)
    pvalue=sp.stats.mannwhitneyu(df2.loc[df2['OutcomeN']==1,column],
                            df2.loc[df2['OutcomeN']==0,column])
    pvalues.append(pvalue[1])
    dictform.append(columninfo)

fdr=statsmodels.stats.multitest.multipletests(pvalues,method='fdr_bh')

counter=0
for i in fdr[1]:
    p=0-np.log10(i)
    dictform[counter].append(p)        
    counter=counter+1

counter=0
for column in f:
    fol=np.mean(df2.loc[df2['OutcomeN']==1,column])-np.mean(df2.loc[df2['OutcomeN']==0,column])
    dictform[counter].append(fol)
    counter=counter+1

def sorter(lister):
    lister.sort(key=lambda x:x[2])
    return(lister)

sorter(dictform)

def belong (x,c1,c2):
    if x in df1.columns[c1:(c2+1)]:
        return True
    else:
        return False


rank=1         
for icon in dictform:
    if (icon[1]>=10) or (np.abs(icon[2])>=2)): ########### change the thresholds based on the plot
        
        print(rank,icon[0])          
        icon.append(rank)
        icon.append(1)
        icon.append(500)

        rank=rank+1    
    else:
        icon.append('')
        icon.append(0)

logps=[]
logfolds=[]
ranks=[]
categories=[]
for member in dictform:
    logps.append(member[1])
    logfolds.append(member[2])
    ranks.append(member[3])
    categories.append(member[4])
y=np.array(logps)
x=np.array(logfolds)
rankname=np.array(ranks)
rankname=np.array(['' for x in ranks])
category=np.array(categories)

##########volcano plot

plt.figure(figsize=(26,16))
ax=plt.axes()
ax.set_facecolor('white')
category=np.array(categories)
colormap=np.array(['Grey','Red'])
plot=plt.scatter(x,y,s=230,c=colormap[category])

for i in range(len(dictform)):
    plt.annotate(rankname[i], (logfolds[i],logps[i] + 0.03),color='black',fontsize=20)

plt.title(f'{variable} differences\n',fontsize=40)
plt.xlabel('\nLog-2 fold difference',fontsize=35)
plt.ylabel('-Log10 (P-value)\n',fontsize=35)
plt.xlim([-2.7,2.7])
plt.ylim([0,42])
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
plt.grid(color='Grey',linestyle='--',linewidth=0.3)
plt.rcParams['axes.edgecolor']="Black"
plt.rcParams['axes.linewidth']=4


######## change horizontal and vertical lines based on plot

plt.axvline(2,color="black",linestyle='--')
plt.axvline(-2,color="black",linestyle='--')
plt.axhline(0.01,color="black",linestyle='--')

plt.show()

