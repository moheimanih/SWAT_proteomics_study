################### Louvain detection and affinity propagation for sensitivity analysis ####################
############################################################################################################


############################################## Louvain community detection

##### Inputs
f=list(set(coagulation+complement+platelet+neutrophil+endothelial+inflammation+mediator_s)) # list of selected proteins

######################### adjust treatment and covariates accordingly
treat='Plasma'
treatment='Plasma'
if treat=='LTOWB':
    mode='Plasma'
else:
    mode='LTOWB'
confounder_q=['ISS','Head_AIS','Age','Shock','phgcs','phrr','Fluids',
              'Cryo_units','Platelet_units','PRBC_units','LTOWB_units']#'Total_blood'
confounder_bi=['Penetrating','Biological_sex']#,'Cryo','Platelet','PRBC',mode]
moderator=confounder_q+confounder_bi

powered=2
splits=10
sample_set='SWAT'

reps=500
split=10

time1=0
time2=24
time3=4

colormap={0:'lightcoral',1:'khaki',2:'yellowgreen',3:'orange',16:'blue',5:'pink',18:'yellow',7:'darkgoldenrod',8:'Green'
          ,9:'aqua',17:'Black',11:'Grey',12:'forestgreen',13:'paleturquoise',14:'navy',15:'darkmagenta',
         4:'indigo',10:'salmon',6:'firebrick',19:'dodgerblue'}

######################################## defining functions

def sorter(lister):
    lister.sort(key=lambda x:x[1],reverse=False)
    return(lister)

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

print()

######################################### data frames

df1=merged_df0.loc[(merged_df0['Severity']!='Milds'),:]
df1=df1.loc[df1['PAID']!='',:] ## fill '' with outlier IDs

for fs in f:
    if fs not in df1.columns:
        f.remove(fs)
        print(f'{fs} was removed from list')

df1.loc[:,f]=np.log2(df1.loc[:,f])
df2=df1.loc[:,['ID',treat,'Time point']+confounder_q+confounder_bi+f]

############## robust standardizaton of protein values
qns2=np.apply_along_axis(robustbase.Qn,0,np.array(df2.loc[:,f]))
dictqn2=dict(zip(f,qns2))

df2.loc[:,f]=(df2.loc[:,f]-df2.loc[:,f].median())
for column in df2.loc[:,f].columns:
    df2[column]=df2[column]/dictqn2[column]

df2=df2.dropna(axis=1)
df2=df2.dropna(axis=0)
df3=df2.loc[df2['Time point']==time1,f]
df4=df2.loc[df2['Time point']==time2,f]
  
###################################################### finding differences and puting it in a frame
ind=np.intersect1d(np.unique(np.array(df2.loc[df2['Time point']==time2,'ID'])),
                   np.unique(np.array(df2.loc[df2['Time point']==time1,'ID'])))
dfcs=pd.DataFrame(columns=f,index=ind)
for i,ind in enumerate(dfcs.index):
    for j in f:
        level24=np.array(df2.loc[(df2['ID']==ind) & (df2['Time point']==time2),j])[0]

        level0=np.array(df2.loc[(df2['ID']==ind) & (df2['Time point']==time1),j])[0]
        dfcs.at[ind,j]=level24-level0

dfcs=dfcs.dropna(axis=1)    
dfcs=dfcs.dropna(axis=0).reset_index()   
dfcs=dfcs.drop(columns='index')
dfcs=dfcs.astype(float)

########################################## correlation network draw and community detection

corr3=df3.corr(method='spearman')
corr4=df4.corr(method='spearman')
corrcs=dfcs.corr(method='spearman')

corr=((((corr3+corr4)**powered)/2)+(corrcs**powered))/2

g = nx.Graph()
    
# add the nodes to the graph, using the variable names as labels
g.add_nodes_from(corr.columns)

# add the edges to the graph, using the correlations as weights

dump=[]
for i, row in corr.iterrows():
    for j in row.index:
        if (j,i) in dump or i==j:          
            continue
        else:
            weight=row[j]
            dump.append((i,j))
            g.add_edge(i, j, weight=weight**powered)
        
g.remove_edges_from(nx.selfloop_edges(g))

Colors=[]
shapes={}

##### getting louvain communities
partition=community.best_partition(g,random_state=int(round(100*np.random.random_sample(),0)))
for cluster_id, nodes in partition.items():
    print(f"{cluster_id}: Cluster {nodes}")
    Colors.append(nodes)
    if nodes not in shapes:
        shapes[nodes]=[cluster_id]
    else:
        shapes[nodes].append(cluster_id)
         
colorlist=[colormap[c] for c in Colors]
pos = nx.spring_layout(g, k=20*1/np.sqrt(len(g.nodes())), iterations=300)
nodesize=[60000*(((corr[v]**1).sum())**0) for v in g]

# draw the graph
sns.set(rc={'figure.figsize':(70,70)})
nx.draw(g, with_labels=False,node_color=colorlist,node_shape='o',node_size=nodesize,pos=pos,
        font_size=50,font_color='Black',width=0.2,alpha=0.99)
plt.show()
print()


#############################################
############################################################## Affinity propagation clustering and graphing each module 

detailed_modules=True
types='Agg' # alternatively 'Diff'
show_average_module_levels=True

time1=0
time2=24
time3=4

olw='Overlap_weights' # as calculated from adjustment modeling

reps=500
split=10

f=list(set(coagulation+complement+platelet+endothelial+neutrophil+inflammation+mediator_s)) # list of selected proteins

confounder_q=['ISS','Age','phrr','Shock','phgcs','Head_AIS',
              'Fluids','Cryo_units','Platelet_units','PRBC_units','LTOWB_units']

confounder_bi=['Penetrating','Biological_sex']
treatment='Plasma' ## or 'LTOWB
moderator=confounder_bi+confounder_q

everything=['PAID','Time point','Age_group','Severity',treatment]+confounder_q+confounder_bi


df00=merged_df0.loc[merged_df0['ISS']!=0,:]
df00.loc[:,f]=np.log2(df00.loc[:,f])
df00.loc[:,f]=(df00.loc[:,f]-df00.loc[:,f].mean())/df00.loc[:,f].std()
df1=df00.loc[(df00['Time point']==time1)|(df00['Time point']==time2),:]

param_grid = {'damping': [0.5, 0.6, 0.7, 0.8, 0.9]}
best_score = -1
best_params = None

########################################################## AF modeling 

if types=='Agg':
    df=df1.loc[:,f+everything].reset_index(drop=True)

elif types=='Diff':
    dfk=df1.loc[df1['Time point']==time1,everything].reset_index(drop=True)

    df=df1.loc[df1['Time point']==time2,f].reset_index(drop=True) \
    -df1.loc[df1['Time point']==time1,f].reset_index(drop=True)
    df.loc[:,everything]=dfk.loc[:,everything]

for params in ParameterGrid(param_grid):
    modelafp_test = AffinityPropagation(**params, max_iter=2000, convergence_iter=30,affinity='euclidean')
    modelafp_test.fit(df[f].T)
    modelafp_test_labels = modelafp_test.labels_

    # Note: silhouette_score expects the original data and the labels as parameters
    score = silhouette_score(df[f].T, modelafp_test_labels)

    if score > best_score:
        best_score = score
        best_params = params

print("Best score:", best_score)
print("Best params:", best_params)

afp=AffinityPropagation(damping=best_params['damping'], max_iter=2000, convergence_iter=30,
                        copy=True, preference=None,affinity='euclidean',
                        verbose=False,random_state=int(round(100*np.random.random_sample(),0)))

modelafp=afp.fit(df[f].T)    


examples=modelafp.cluster_centers_indices_#+len(everything)
print(len(examples),examples)

n_clusters=len(modelafp.cluster_centers_indices_)

    
labelsafp=['Cluster ' + str(x) for x in range (0,n_clusters)]
cluster_pros=[[] for x in range (0,n_clusters)]

for clust in range (0,n_clusters):
    for fs in range(len(f)):
        if modelafp.labels_[fs]==clust:
            cluster_pros[clust].append(f[fs])
            
key_value_pairs = zip(labelsafp,cluster_pros)

# convert the list of key-value pairs to a dictionary
afp_clusters= dict(key_value_pairs)
cluster_df=df.iloc[:,examples]#+len(everything)

############################################## getting module level averages into df00 columns
for labeln in range(len(labelsafp)):
    
    label=labelsafp[labeln]
    df00[label]=0
    weights=0
    
    for memb in cluster_pros[labeln]:
        df00[label]=df00[label]+df00[memb]
        weights=weights+1
        
    df00[label]=df00[label]/weights
    df00[label]=(df00[label]-df00[label].mean())/df00[label].std()

print()

####################### Densmap

Colors=[]
for pro in f:
    for j in cluster_pros:
        if pro in j:
            Colors.append(cluster_pros.index(j))
colormap={0:'lightcoral',1:'khaki',2:'yellowgreen',3:'orange',4:'blue',5:'pink',6:'yellow',7:'darkgoldenrod',8:'Green'
          ,9:'aqua',10:'Black',11:'Grey',12:'forestgreen',13:'paleturquoise',14:'navy',15:'darkmagenta',
         16:'indigo',17:'salmon',18:'firebrick',19:'dodgerblue',20:'dimgrey',21:'seashell',22:'darkcyan',23:'magenta'}
colorlist=[colormap[c] for c in Colors]


print("DENS_MAP")

embedding=umap.UMAP(random_state=int(round(100*np.random.random_sample(),0))
                    ,n_components=2,densmap=True,dens_lambda=0.6,n_neighbors=15,
                   learning_rate=0.6, init='spectral',
                    spread=1.0,min_dist=0.2).fit_transform(df[f].T)

sns.set(rc={'figure.figsize':(9,9)})
ax=plt.axes()
ax.set_facecolor('whitesmoke')
plt.scatter(embedding[:, 0],embedding[:, 1],s=100,c=colorlist)#,label=colormap)
plt.title('DensMAP projection of selected proteins\n', fontsize=24)
plt.xlabel('\nD1',fontsize=25)
plt.ylabel('D2\n',fontsize=25)
plt.grid(color='grey',linestyle='-',linewidth=0.3)
plt.legend()
plt.show()  


################################# module trajectory and adjustment with Affinity propagation

#################################### unadjusted plot
if detailed_modules:
   
    def plotinfo (name):
        plt.title(name,fontsize=25)
        plt.xlabel('\nTime point (hours post-admission)',fontsize=20)
        plt.ylabel('Standardized Module score',fontsize=22)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=14)
        ax = fig.gca()
        handles, labeled = ax.get_legend_handles_labels()
        label_mapping = {0: "No", 1: "Yes"}
        labeled = [label_mapping.get(int(float(labela)), labela) for labela in labeled]
        legend_properties = {'size': '15', 'family': 'Arial'}
        plt.legend(handles, labeled,title=f'Received {treat}',   
                   frameon=True, 
                   edgecolor='black', 
                   facecolor='0.99', 
                   framealpha=1.0, 
                   fancybox=True, 
                   shadow=True, 
                   title_fontsize='16', 
                   prop=legend_properties)

        plt.grid(color='0.10',linestyle='--',linewidth=0.3)
        plot.set_facecolor('0.98') 

    ############################################## Choose what to show and analyze, module or median cluster member

    print()
    if show_average_module_levels:   
        choice_of_columns=labelsafp
        print('We going with modules')
    else:
        choice_of_columns=cluster_df.columns
        print('We going with examplar proteins')


    ##############################################    
    cluster_counter=0

    for label in choice_of_columns:

        print('\n\n')
        print(label)
        print()

        if show_average_module_levels:
            def group_printer(name):
                print(cluster_pros[cluster_counter])
        else:
            def group_printer(name):
                for group in cluster_pros:
                    if name in group:
                        print('Members:',group)
        group_printer(label)
        cluster_counter=cluster_counter+1
      
        #creating a fake 14 h time point for plot design only
        df_with_14=df00.groupby(['Time point',treat])[label].median()

        mean_0=(df_with_14[4][0]+df_with_14[24][0])/2
        mean_1=(df_with_14[4][1]+df_with_14[24][1])/2
        new_rows = []


        for k,row in df00.iterrows():
            if df00.at[k,'Time point']==24:
                new_row = row.copy()
                new_row['Time point'] = int(14)

                if df00.at[k,treat]==0:
                    new_row[label] = float(mean_0)

                elif df00.at[k,treat]==1:
                    new_row[label] = float(mean_1)

                new_rows.append(new_row)

        duplicated_rows = pd.DataFrame(new_rows)
        df00_with_14 = pd.concat([df00, duplicated_rows], ignore_index=True)

        fig =plt.figure(figsize=(10,6))
        plot=sns.pointplot(y=df00_with_14.loc[df00_with_14['Sample']=='Donor',label],
                           x=df00_with_14.loc[df00_with_14['Sample']=='Donor','Time point'],
                           order=[-1],#linestyles=['dashed','dotted','dashdot','solid'],markers=['D','o','^','v'],
                           capsize=0.12,errwidth=4.4,estimator=np.median,ci=95,n_boot=10000,dodge=True,color='forestgreen')

        plot=sns.pointplot(y=df00_with_14.loc[df00_with_14['Sample']!='Donor',label],
                           x=df00_with_14.loc[df00_with_14['Sample']!='Donor','Time point'],
                           order=[-1,0,4,14,24],hue=df00_with_14.loc[df00_with_14['Sample']!='Donor',treat],
                           linestyles=['dashed','dotted','dashdot','solid'],markers=['D','o','^','v'],#size=2000,'s'
                           capsize=0.12,errwidth=4.4,estimator=np.median,ci=95,n_boot=10000,dodge=True)

        sns.color_palette('rocket')
        plotinfo(f'{label}\n')
        plt.show()
        print()


        ########################################################### adjusted plot
        def weighted_average_and_se(group):

            data_col = group[label]
            weights_col = group[olw]

            ################################################ new sample size
            overlaplistS=group[olw].to_list()
            overlapwS=np.array(overlaplistS)

            overlapwsumS=np.array(overlaplistS).sum()
            squarer = lambda t: t ** 2
            overlapwsquaresS= np.array([squarer(xi) for xi in overlaplistS])
            overlapwsumsquaresS=overlapwsquaresS.sum()

            effective_nS=(overlapwsumS**2)/overlapwsumsquaresS
        
            weighted_avg = np.average(data_col, weights=weights_col)
            variance = np.average((data_col-weighted_avg)**2, weights=weights_col)
            weighted_se = np.sqrt(variance/effective_nS)

            return pd.Series({'Weighted_Average': weighted_avg, 'Standard_Error': weighted_se})

        df_agg = df00.groupby(['Time point', treat]).apply(weighted_average_and_se).reset_index()       
        sns.set_style("whitegrid")
        unique_treatments = df_agg[treat].unique()
        df_agg = df00.groupby(['Time point', treat]).apply(weighted_average_and_se).reset_index()
        sns.set_style("whitegrid")
        palette = sns.color_palette("pastel", n_colors=len(unique_treatments))
        palette=['forestGreen','Crimson']
        fig = plt.figure(figsize=(8, 8))
        ax = fig.gca()
        dodge_amount =0.19 # adjust this value as needed

        for idx, treatment_value in enumerate(unique_treatments):
            subset = df_agg[df_agg[treat] == treatment_value]
            if treatment_value == 0:  # replace 'SpecificTreatment1' with the actual name
                x_values = subset['Time point'] - dodge_amount
            elif treatment_value == 1:  # replace 'SpecificTreatment2' with the actual name
                x_values = subset['Time point'] + dodge_amount
            else:
                x_values = subset['Time point']
            if treatment_value == 0:  # replace 'TreatmentToBeDashed' with the actual name
                line_style = '--'
            else:
                line_style = '-'

            plt.errorbar(x_values, subset['Weighted_Average'], 
                         yerr=subset['Standard_Error'], label=treatment_value,
                         fmt='-o', capsize=10, color=palette[idx],
                         elinewidth=5.5, capthick=4, linewidth=3.5, markersize=10, linestyle=line_style)

        sns.despine()

        handles, labeled = ax.get_legend_handles_labels()
        label_mapping = {0: "No", 1: "Yes"}
        labeled = [label_mapping.get(int(float(labela)), labela) for labela in labeled]
        legend_properties = {'size': '17', 'family': 'Arial'}

        plt.legend(handles, labeled,title=f'Received {treat}',   
                   frameon=True, 
                   edgecolor='black', 
                   facecolor='0.99', 
                   framealpha=1.0, 
                   fancybox=True, 
                   shadow=True, 
                   title_fontsize='16', 
                   prop=legend_properties)
        plt.xlabel('\nTime point (hours post-admission)', fontsize=20)
        plt.ylabel('Standardized module level\n', fontsize=17)
        plt.title('\n\n\n', fontsize=17)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=15)
        plt.tight_layout()
        plt.show()

        ###################################### Causal modeling
      
        for timer in [time1,time3,time2]:

            print('\n')
            print(label,timer)
            target=[label]
          
            ################################################ Designing the dataframes for downstream analysis
            df10=df00.loc[df00['Time point']==timer,:].reset_index()
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

            df10.loc[:,confounder_q]=(df10.loc[:,confounder_q]-df10.loc[:,confounder_q].mean())/df10.loc[:,confounder_q].std()
            df20=df10.loc[df10['ISS']!=0,:]
            df20=df20.dropna(axis=1)
            df30=df20.loc[:,confounder_bi+confounder_q]
            df40=df20.loc[:,moderator]
            df50=df20.loc[:,target]

            W=np.array(df30)
            X=np.array(df40)
            Y=np.array(df50)
            T=np.array(df20[treatment])

            ####################################################################### DoubleML
            print('EconML with heterogeneity')

            est=econdml.CausalForestDML(model_y=RandomForestRegressor(n_estimators=reps,min_samples_split=split),
                                    model_t=RandomForestClassifier(n_estimators=reps,min_samples_split=split),
                                    cv=split,random_state=42,n_estimators=reps,discrete_treatment=True)

            est.fit(Y, T, X=X, W=W)
            print(est.ate__inference().summary_frame())
            print()


else:
    def plotinfo (name):
        plt.title(name,fontsize=25)
        plt.xlabel('\nTime point (hours post-admission)',fontsize=20)
        plt.ylabel('Standardized Module score',fontsize=22)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=14)
        ax = fig.gca()
        handles, labeled = ax.get_legend_handles_labels()
        label_mapping = {0: "No", 1: "Yes"}
        labeled = [label_mapping.get(int(float(labela)), labela) for labela in labeled]
        legend_properties = {'size': '15', 'family': 'Arial'}
        plt.legend(handles, labeled,title=f'Received {treat}',   

                   frameon=True, 
                   edgecolor='black', 
                   facecolor='0.99', 
                   framealpha=1.0, 
                   fancybox=True, 
                   shadow=True, 
                   title_fontsize='16', 
                   prop=legend_properties)

        plt.grid(color='0.10',linestyle='--',linewidth=0.3)
        plot.set_facecolor('0.98') 

    ############################################## Choose what to show and analyze, module or median cluster member

    print()
    if show_average_module_levels:   
        choice_of_columns=labelsafp
        print('We going with modules')
    else:
        choice_of_columns=cluster_df.columns
        print('We going with examplar proteins')

    ##############################################    
    cluster_counter=0

    for label in choice_of_columns:

        print('\n\n\n')
        print(label)
        print()

        if show_average_module_levels:
            def group_printer(name):
                print(cluster_pros[cluster_counter])
        else:
            def group_printer(name):
                for group in cluster_pros:
                    if name in group:
                        print('Members:',group)
        group_printer(label)
        cluster_counter=cluster_counter+1
       
