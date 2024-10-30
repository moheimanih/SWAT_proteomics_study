############## Effect modification analysis and propensity matching for treatment effect estimation ################
####################################################################################################################

############################################# Causal forest for treatment effect on PT

####################### model specification
reps=500
split=10
Y_prediction=RandomForestRegressor(n_estimators=reps,min_samples_split=split,random_state=42)
T_prediction=RandomForestClassifier(n_estimators=reps,min_samples_split=split,random_state=42)

clf=RandomForestClassifier(n_estimators=reps,max_depth=split,min_samples_split=split,random_state=42)
reg=lgb.LGBMRegressor(n_estimators=reps,max_depth=split,min_samples_split=split,random_state=42,verbose= -100)

############################### data selection

def recode_value_blood(value):
    if value==0:
        return 0
    else:
        return 1

dfmodeling['LTOWB_bi']=dfmodeling['LTOWB_0h'].apply(recode_value_blood)
dfmodeling['Plasma_bi']=dfmodeling['Plasma_0h'].apply(recode_value_blood)
dfmodeling['PRBC_bi']=dfmodeling['PRBC_0h'].apply(recode_value_blood)
dfmodeling['Platelet_bi']=dfmodeling['Platelet_0h'].apply(recode_value_blood)
dfmodeling['Cryo_bi']=dfmodeling['Cryo_0h'].apply(recode_value_blood)

treat='Plasma_bi'
target='PT_0h'
confounder_bi=['Penetrating','Biological_sex','Platelet_bi','Cryo_bi','PRBC_bi','LTOWB_bi']#
confounder_q=['ISS','Age','phgcs','Fluids_0h','Shock','phrr','Head_AIS','Total_blood_0h']       
moderator=confounder_q+confounder_bi

###################################### Double ML with causal forrest

dfclin=dfmodeling.loc[(dfmodeling[treat]!=2),:].reset_index()
dfclin.loc[:,target]=(dfclin.loc[:,target] \
                                    -dfclin.loc[:,target].mean())/dfclin.loc[:,target].std()
dfclin.loc[:,confounder_q]=(dfclin.loc[:,confounder_q] \
                                    -dfclin.loc[:,confounder_q].mean())/dfclin.loc[:,confounder_q].std()


print('EconML with heterogeneity for binary treatments')
df30=dfclin.loc[:,confounder_bi+confounder_q]
df40=dfclin.loc[:,moderator]
df50=dfclin.loc[:,target]

W=np.array(df30)
X=np.array(df40)
Y=np.array(df50)

T=np.array(dfclin[treat])

est=econdml.CausalForestDML(model_y=reg,model_t=clf,#reg
                            cv=split,random_state=42,n_estimators=reps,discrete_treatment=True)

est.fit(Y, T, X=X, W=W)
print(est.ate__inference().summary_frame())

#################################### visualizing effect modification

means={}
conf_intervals={}
for risk in ['Head_AIS_bi','Prehospital_GCS',#'Platelet_bi','Cryo_bi', add as desired
             'Shock_bi','Penetrating','Severity','Age_group','Biological_sex']:
        mo=[x for x in moderator if x!=risk]
        X1=dfclin.loc[dfclin[risk]==1,moderator]
        X0=dfclin.loc[dfclin[risk]==0,moderator]


        summaries1=est.effect_inference(X1).population_summary(alpha=0.05,value=1,decimals=2, tol=0.0001)#
        risk_means1=summaries1.mean_point
        risk_vals1=(summaries1.conf_int_mean()[0],summaries1.conf_int_mean()[1])
        means[f'{risk} on']=risk_means1
        conf_intervals[f'{risk} on']=risk_vals1

        summaries0=est.effect_inference(X0).population_summary(alpha=0.05,value=1,decimals=2, tol=0.0001)# 
        risk_means0=summaries0.mean_point
        risk_vals0=(summaries0.conf_int_mean()[0],summaries0.conf_int_mean()[1])
        means[f'{risk} off']=risk_means0
        conf_intervals[f'{risk} off']=risk_vals0

        print(f'{risk}')
        print('Pval: ',np.round(6*sp.stats.mannwhitneyu(est.effect(X1),est.effect(X0))[1],3))
        print('CohenD: ',np.round(cohend(est.effect(X1),est.effect(X0)),3))
        print()

colors = ['Crimson','green'] * (len(means) // 2)
positions = [0]
narrow_spacing = 1
wide_spacing = 2


for i in range(1, len(means)):
    if i % 2 != 0:
        positions.append(positions[-1] + narrow_spacing)
    else:
        positions.append(positions[-1] + wide_spacing)

sns.set(rc={'figure.figsize':(4,6)})
fig, ax = plt.subplots()

for i, ((key, mean), color) in enumerate(zip(means.items(), colors)):
    lower, upper = conf_intervals[key]
    error = [[mean - lower], [upper - mean]]

    ax.errorbar(mean, positions[i], xerr=error, fmt='o', color=color, ecolor=color,
                label=key if i < 2 else "", capsize=5, capthick=2, linestyle='')

ax.set_facecolor('0.98')
plt.yticks(positions, means.keys())
plt.xlabel(f'\n Standardized plasma effect on {target}',fontsize=16)
plt.ylabel('Potential moderators\n',fontsize=14)
plt.title(f'Assessing effect moodification at 0h\n',fontsize=17)
ax.grid(color='grey',linestyle='--',linewidth=0.3)
plt.show()
print()


###################################
####################################################### propensity matching for PT/INR (Alternatively for platelet count)

def cohend(d1: pd.Series, d2: pd.Series) -> float:
    n1, n2 = len(d1), len(d2)
    s1, s2 = np.var(d1, ddof=1), np.var(d2, ddof=1)
    s = np.sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    u1, u2 = np.mean(d1), np.mean(d2)
    return (u1 - u2) / s

N_treat='WB_vs_Plasma_q' ## exclusive WB vs plasma at admission time
dfclin2=dfmodeling.loc[(dfmodeling[N_treat]!=2),:].reset_index()

indices1 = []
indices2 = []
indices7 = []
indices = []

### misses were indexed before imputation
for value in PT0_misses: 
    indices1.extend(dfclin2.index[dfclin2['index'] == value].tolist())
for value in INR0_misses:
    indices2.extend(dfclin2.index[dfclin2['index'] == value].tolist())
        
for value in Platelets0_misses:
    indices7.extend(dfclin2.index[dfclin2['index'] == value].tolist())

indices=list(set(indices1+indices2)) # alternatively indices7 for platelets

df_droped=dfclin2.drop(index=indices).reset_index()
print(dfclin2.shape[0],df_droped.shape[0],len(PT0_misses),len(INR0_misses),len(Platelets0_misses),len(indices))

############################################ matching

target='PT_0h'
cons=['Head_AIS_bi','Penetrating','Biological_sex','ISS','Age','phgcs','Fluids_0h','Shock','phrr',
      'Platelet_0h','Cryo_0h','LTOWB_plasma_PRBC']#LTOWB_plasma_PRBC signifies sum of products
df_matched=df_droped.loc[df_droped[N_treat]!=2,cons+['ID',N_treat,target]]

sns.set(rc={'figure.figsize':(10,8)}, font_scale = 1.3)
psm = PsmPy(df_matched,
            treatment=N_treat, indx='ID', exclude = [target])
psm.logistic_ps(balance = True)
psm.knn_matched(matcher='propensity_logit', replacement=False, caliper=0.15)
psm.plot_match(Title='Side by side matched controls',
               Ylabel='Number of patients', Xlabel= 'Propensity logit',
               names = ['Componants', 'LTOWB'], save=False)
plt.show()

psm.effect_size_plot(save=False)
plt.show()

matchedIDs=psm.df_matched['ID'].values
print('Selected for each group:',(len(matchedIDs))/2)
filtered_df=dfclin2[dfclin2['ID'].isin(matchedIDs)]

############################################################### Comparing outcomes

PT_p=sp.stats.mannwhitneyu(filtered_df.loc[filtered_df[N_treat]==1,'PT_0h'],
                            filtered_df.loc[filtered_df[N_treat]==0,'PT_0h'])[1]
PT_cohen=round(cohend(filtered_df.loc[filtered_df[N_treat]==1,'PT_0h'],
                      filtered_df.loc[filtered_df[N_treat]==0,'PT_0h']),2)

INR_p=sp.stats.mannwhitneyu(filtered_df.loc[filtered_df[N_treat]==1,'INR_0h'],
                            filtered_df.loc[filtered_df[N_treat]==0,'INR_0h'])[1]
INR_cohen=round(cohend(filtered_df.loc[filtered_df[N_treat]==1,'INR_0h'],
                       filtered_df.loc[filtered_df[N_treat]==0,'INR_0h']),2)

Platelets_p=sp.stats.mannwhitneyu(filtered_df.loc[filtered_df[N_treat]==1,'Platelets_0h'],
                            filtered_df.loc[filtered_df[N_treat]==0,'Platelets_0h'])[1]
Platelets_cohen=round(cohend(filtered_df.loc[filtered_df[N_treat]==1,'Platelets_0h'],
                       filtered_df.loc[filtered_df[N_treat]==0,'Platelets_0h']),2)

print('PT\n','pval',round(PT_p,4),'Cohens d',PT_cohen)     
print('\n\nINR\n','pval',round(INR_p,4),'Cohens d',INR_cohen)        
print('\n\nPlatelet count\n','pval',round(Platelets_p,3),'Cohens d',Platelets_cohen)
print(sp.stats.fisher_exact(pd.crosstab(filtered_df[N_treat],filtered_df['30d_Mortality'])))

############################ errors
for target in ['PT_0h','INR_0h','Platelets_0h']:
   
    print('\n')

    plot.set_facecolor('0.99')
    grouped = filtered_df.groupby(N_treat)[target].agg(['mean', sp.stats.sem]).reset_index()
    grouped.rename(columns={'mean': 'Mean', 'sem': 'SE'}, inplace=True)

    category_numeric = {cat: i for i, cat in enumerate(grouped[N_treat])}
    grouped['Category_numeric'] = grouped[N_treat].map(category_numeric)

    colors = ['green','Crimson']
    labels = ['Plasma', 'LTOWB']

    sns.set(rc={'figure.figsize':(4,6)})
    fig, ax = plt.subplots()
    for i, (index, row) in enumerate(grouped.iterrows()):
        ax.errorbar(row['Category_numeric'], row['Mean'], yerr=row['SE'], label=labels[i],
                    elinewidth=7, markersize=18, fmt='o', capsize=5, linestyle=None, color=colors[i])
    ax.set_facecolor('0.98')

    ax.set_xticks(list(category_numeric.values()))
    ax.set_xticklabels(list(category_numeric.keys()),fontsize=20)
    plt.xticks(ticks=[0, 1], labels=['Component', 'LTOWB'])
    plt.yticks(fontsize=14)
    ax.set_title(f'{target}\n',fontsize=25)
    ax.grid(color='grey',linestyle='--',linewidth=0.3)
    plt.show()

############################### matched summary table

k_filtered=filtered_df.groupby(N_treat)[['Plasma_0h','Platelet_0h','PRBC_0h','LTOWB_0h','Cryo_0h','Total_blood_0h',
                                    'Fluids_0h','Age','ISS','Head_AIS','Face_AIS','Chest_AIS','Abdomen_AIS',
                                    'Extremity_AIS','External_AIS','30d_Mortality','Biological_sex',
                                    'Penetrating','TBI','phgcs','Shock','phsbp','phrr','phhr','Outcome',
                                    'Hospital_days','ICU_days',
                                    'Ventilator_days','INR_0h','PT_0h','Platelets_0h']].mean() ## alternatively quantile(0.5)


################################ refutation

for target in['PT_0h','INR_0h','Platelets_0h']:
    print('\n\n')
    print(target)
   
    dfwhy0=df_droped.loc[df_droped[N_treat]!=2,[target,N_treat]+cons]

    model=dowhy.CausalModel(
            data = dfwhy0,
            treatment=N_treat,
            outcome=target,
            common_causes=cons,
            effect_modifiers=cons
            #graph=causal_graph.replace("\n", " "),
            #instruments=data["instrument_names"]
            )
    #model.view_model()

    identified_estimand = model.identify_effect(proceed_when_unidentifiable=True) 
    random_refuter_match=model.refute_estimate(identified_estimand,
                                             causal_estimate_match,method_name="random_common_cause") 
    print('Match random refuter: ',random_refuter_match)


